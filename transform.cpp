#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstring>

void printUsage(const char *path) {
    std::cerr << "\nProgram: Transform genotypes/dosages to dominance encoding\n";
    std::cerr << "Usage: " << path << " --input <input.vcf.gz> [options]\n";
    std::cerr << "\nRequired Options:\n";
    std::cerr << "  --input/-i <file>          Input VCF/BCF file (supports .vcf, .vcf.gz, .bcf)\n";
    std::cerr << "\nMain Options:\n";
    std::cerr << "  --mode/-m <mode>           Processing mode (default: dominance)\n";
    std::cerr << "                             Supported modes:\n";
    std::cerr << "                               - dominance: Dominance deviation encoding\n";
    std::cerr << "                               - recessive: Set heterozygotes to 0, homozygotes unchanged\n";
    std::cerr << "  --scale-dosages            Enable dosage scaling to [0,2] range\n";
    std::cerr << "  --scale-dosages-factor <f> Apply scaling factor to dosages (default: 1.0)\n";
    std::cerr << "\nAdditional Options:\n";
    std::cerr << "  --set-variant-id           Set variant IDs to chr:pos:ref:alt format\n";
    std::cerr << "  --all-info                 Include additional info in output (frequencies, min/max)\n";
    std::cerr << "  --gene-map/-g <file>       Optional gene mapping file for gene-specific scaling\n";
    std::cerr << "                             Format: tab-separated with headers 'variant' and 'gene'\n";
    std::cerr << "\nInput Requirements:\n";
    std::cerr << "  - Input file must contain either GT (genotype) or DS (dosage) format fields\n";
    std::cerr << "  - DS values should be in range [0,2] where:\n";
    std::cerr << "    0 = homozygous reference (aa)\n";
    std::cerr << "    1 = heterozygous (Aa)\n";
    std::cerr << "    2 = homozygous alternate (AA)\n";
    std::cerr << "\nExample Usage:\n";
    std::cerr << "  Basic usage:\n";
    std::cerr << "    " << path << " --input sample.vcf.gz --scale-dosages\n";
    std::cerr << "\n  With gene mapping:\n";
    std::cerr << "    " << path << " --input sample.vcf.gz --scale-dosages --gene-map genes.txt\n";
    std::cerr << "\n  Full output with variant IDs:\n";
    std::cerr << "    " << path << " --input sample.vcf.gz --scale-dosages --all-info --set-variant-id\n";
    std::cerr << "\n  Using recessive mode:\n";
    std::cerr << "    " << path << " --input sample.vcf.gz --mode recessive\n";
}

std::vector<std::string> sortChromosomes(const std::set<std::string> &contigs) {
    std::vector<std::string> chromosomes(contigs.begin(), contigs.end());
    std::sort(chromosomes.begin(), chromosomes.end(), [](const std::string &a, const std::string &b) {
        std::string a_num = a.substr(0, 3) == "chr" ? a.substr(3) : a;
        std::string b_num = b.substr(0, 3) == "chr" ? b.substr(3) : b;

        if (isdigit(a_num[0]) && isdigit(b_num[0])) {
            return std::stoi(a_num) < std::stoi(b_num);
        }
        if (isdigit(a_num[0]) && !isdigit(b_num[0])) {
            return true;
        }
        if (!isdigit(a_num[0]) && isdigit(b_num[0])) {
            return false;
        }
        return a < b;
    });
    return chromosomes;
}

bool parseArguments(int argc, char *argv[], std::string &pathInput, std::string &mode, double &scalingFactor, bool &scaleDosage, bool &setVariantId, bool &allInfo, std::string &geneMapPath) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return false;
        } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
            pathInput = argv[++i];
        } else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
            mode = argv[++i];
        } else if (arg == "--scale-dosages-factor" && i + 1 < argc) {
            scalingFactor = std::stod(argv[++i]);
        } else if (arg == "--scale-dosages") {
            scaleDosage = true;
        } else if (arg == "--set-variant-id") {
            setVariantId = true;
        } else if (arg == "--all-info") {
            allInfo = true;
        } else if ((arg == "--gene-map" || arg == "-g") && i + 1 < argc) {
            geneMapPath = argv[++i];
        } else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return false;
        }
    }

    if (pathInput.empty()) {
        std::cerr << "Error! --input must be provided." << std::endl;
        printUsage(argv[0]);
        return false;
    }
    return true;
}

std::map<std::string, std::string> readGeneMap(const std::string &geneMapPath) {
    std::map<std::string, std::string> geneMap;
    std::ifstream infile(geneMapPath);
    if (!infile) {
        std::cerr << "Error: Cannot open gene map file for reading: " << geneMapPath << std::endl;
        exit(1);
    }

    std::string line, variant, gene;
    std::getline(infile, line); // skip header
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> variant >> gene)) { break; }
        geneMap[variant] = gene;
    }

    infile.close();
    return geneMap;
}


bool hasFormat(const bcf_hdr_t *hdr, const char *format) {
    return bcf_hdr_id2int(hdr, BCF_DT_ID, format) >= 0;
}

void processSampleGenotypes(int *gt_ptr, float ds_value, bool has_gt, bool has_ds, 
                          int &aa_count, int &Aa_count, int &AA_count) {
    if (has_gt && (!bcf_gt_is_missing(gt_ptr[0]) && !bcf_gt_is_missing(gt_ptr[1]))) {
        // Use GT if available
        int allele1 = bcf_gt_allele(gt_ptr[0]);
        int allele2 = bcf_gt_allele(gt_ptr[1]);
        if (allele1 == 0 && allele2 == 0) aa_count++;
        else if (allele1 != allele2) Aa_count++;
        else if (allele1 > 0 && allele2 > 0) AA_count++;
    } else if (has_ds && !std::isnan(ds_value)) {
        // Round DS value to nearest integer
        int rounded_ds = std::round(ds_value);
        if (rounded_ds < 0 || rounded_ds > 2) {
            std::cerr << "Warning: DS value " << ds_value << " out of range [0,2]. Skipping.\n";
            return;
        }
        if (rounded_ds == 0) aa_count++;
        else if (rounded_ds == 1) Aa_count++;
        else if (rounded_ds == 2) AA_count++;
    }
}

// Modify calculateGlobalAndGeneDosages function
void calculateGlobalAndGeneDosages(htsFile *fp, bcf_hdr_t *hdr, double &globalMinDomDosage, 
    double &globalMaxDomDosage, int n_samples, std::set<std::string> &chromosomes, 
    const std::map<std::string, std::string> &geneMap, 
    std::map<std::string, std::pair<double, double>> &geneDosages) {
    
    bcf1_t *rec = bcf_init();
    int *gt_arr = NULL, ngt_arr = 0;
    float *ds_arr = NULL;
    int nds_arr = 0;
    
    bool has_gt = hasFormat(hdr, "GT");
    bool has_ds = hasFormat(hdr, "DS");

    while (bcf_read(fp, hdr, rec) == 0) {
        chromosomes.insert(bcf_hdr_id2name(hdr, rec->rid));
        bcf_unpack(rec, BCF_UN_STR);
        
        int ngt = has_gt ? bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) : 0;
        int nds = has_ds ? bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds_arr) : 0;
        
        int aa_count = 0, Aa_count = 0, AA_count = 0;
        
        for (int i = 0; i < n_samples; ++i) {
            int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
            float ds_value = has_ds ? ds_arr[i] : 0.0f;
            
            processSampleGenotypes(gt_ptr, ds_value, has_gt, has_ds, aa_count, Aa_count, AA_count);
        }

        if (AA_count > 0) {
            double r = static_cast<double>(aa_count) / n_samples;
            double h = static_cast<double>(Aa_count) / n_samples;
            double a = static_cast<double>(AA_count) / n_samples;

            double dom_dosage_aa = -h * a;
            double dom_dosage_Aa = 2 * a * r;
            double dom_dosage_AA = -h * r;

            globalMinDomDosage = std::min({globalMinDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
            globalMaxDomDosage = std::max({globalMaxDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

            std::string variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" + std::to_string(rec->pos + 1) + ":" + rec->d.allele[0] + ":" + rec->d.allele[1];
            if (geneMap.find(variantId) != geneMap.end()) {
                std::string gene = geneMap.at(variantId);
                double &geneMinDosage = geneDosages[gene].first;
                double &geneMaxDosage = geneDosages[gene].second;
                geneMinDosage = std::min({geneMinDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
                geneMaxDosage = std::max({geneMaxDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
            }
        }
    }

    free(gt_arr);
    free(ds_arr);
    bcf_destroy(rec);
}

void printHeader(const bcf_hdr_t *hdr, const std::vector<std::string> &sortedContigs, const std::string &mode, double globalMinDomDosage, double globalMaxDomDosage, bool allInfo, bool scaleDosage, bool geneMapSpecified) {
    int n_samples = bcf_hdr_nsamples(hdr);

    std::cout << "##fileformat=VCFv4.2\n";
    std::cout << "##EncodingMode=" << mode << "\n";
    // output chromosomes
    for (const auto& chr : sortedContigs) {
        std::cout << "##contig=<ID=" << chr << ">\n";
    }
    std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n";
    std::cout << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
    
    if (mode == "dominance") {
        std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the alternate allele based on dominance model\">\n";
    } else if (mode == "recessive") {
        std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the alternate allele with recessive encoding (het=0, hom=2)\">\n";
    }
    
    if (allInfo) {
        std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of bi-allelic references (aa)\">\n";
        std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of heterozygotes (Aa)\">\n";
        std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of bi-allelic alternates (AA)\">\n";
    }
    if (scaleDosage) {
        if (geneMapSpecified) {
            std::cout << "##INFO=<ID=localMinDomDosage,Number=1,Type=Float,Description=\"Local minimum dominance dosage\">\n";
            std::cout << "##INFO=<ID=localMaxDomDosage,Number=1,Type=Float,Description=\"Local maximum dominance dosage\">\n";
            std::cout << "##INFO=<ID=group,Number=1,Type=String,Description=\"Group (gene) used for local scaling\">\n";
        } else {
            std::cout << "##INFO=<ID=globalMinDomDosage,Number=1,Type=Float,Description=\"Global minimum dominance dosage\">\n";
            std::cout << "##INFO=<ID=globalMaxDomDosage,Number=1,Type=Float,Description=\"Global maximum dominance dosage\">\n";
        }
    }
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 0; i < n_samples; ++i) {
        std::cout << "\t" << bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
    }
    std::cout << "\n";
}


void processVcfFile(htsFile *fp, bcf_hdr_t *hdr, const std::string &mode, double globalMinDomDosage, double globalMaxDomDosage, bool allInfo, bool scaleDosage, bool setVariantId, double scalingFactor, const std::map<std::string, std::string> &geneMap, const std::map<std::string, std::pair<double, double>> &geneDosages) {
    bcf1_t *rec = bcf_init();
    int *gt_arr = NULL, ngt_arr = 0;
    float *ds_arr = NULL;
    int nds_arr = 0;
    int n_samples = bcf_hdr_nsamples(hdr);
    const double epsilon = 1e-8;
    int countVariantsWithoutHomAlt = 0;
    const int limit = 5;
    int discardedVariantsCount = 0;
    const int warningLimit = 5;

    bool has_gt = hasFormat(hdr, "GT");
    bool has_ds = hasFormat(hdr, "DS");

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        
        int ngt = has_gt ? bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) : 0;
        int nds = has_ds ? bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds_arr) : 0;
        
        if (!has_gt && !has_ds) continue;

        int aa_count = 0, Aa_count = 0, AA_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
            float ds_value = has_ds ? ds_arr[i] : 0.0f;
            
            processSampleGenotypes(gt_ptr, ds_value, has_gt, has_ds, aa_count, Aa_count, AA_count);
        }

        // For dominance mode, we need at least one homozygous alternate
        if ((mode == "dominance" && AA_count == 0) || (mode != "dominance" && mode != "recessive")) {
            if (mode == "dominance" && AA_count == 0) {
                countVariantsWithoutHomAlt++;
                if (countVariantsWithoutHomAlt <= limit) {
                    std::cerr << "variant '" << bcf_hdr_id2name(hdr, rec->rid) << ":" << (rec->pos + 1) << ":"
                             << rec->d.allele[0] << ":" << (rec->n_allele > 1 ? rec->d.allele[1] : ".")  
                             << "' has no homozygous alternate alleles. Skipping..\n";
                }
            }
            continue;
        }

        double r = static_cast<double>(aa_count) / n_samples;
        double h = static_cast<double>(Aa_count) / n_samples;
        double a = static_cast<double>(AA_count) / n_samples;

        double dom_dosage_aa = -h * a;
        double dom_dosage_Aa = 2 * a * r;
        double dom_dosage_AA = -h * r;

        std::string variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" + 
                              std::to_string(rec->pos + 1) + ":" + 
                              rec->d.allele[0] + ":" + rec->d.allele[1];
                              
        if (!geneMap.empty() && geneMap.find(variantId) == geneMap.end()) {
            if (discardedVariantsCount < warningLimit) {
                std::cerr << "Warning: Variant '" << variantId << "' not found in gene map. Discarding.\n";
            }
            discardedVariantsCount++;
            continue;
        }

        double localMinDomDosage = globalMinDomDosage;
        double localMaxDomDosage = globalMaxDomDosage;
        std::string group;

        if (!geneMap.empty() && geneMap.find(variantId) != geneMap.end()) {
            std::string gene = geneMap.at(variantId);
            localMinDomDosage = geneDosages.at(gene).first;
            localMaxDomDosage = geneDosages.at(gene).second;
            group = gene;
        }

        if (setVariantId) {
            variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" + 
                       std::to_string(rec->pos + 1) + ":" + 
                       rec->d.allele[0] + ":" + 
                       (rec->n_allele > 1 ? rec->d.allele[1] : ".");
        }

        // Output variant information
        std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" 
                 << (rec->pos + 1) << "\t" 
                 << variantId << "\t"
                 << rec->d.allele[0] << "\t";
        
        if (rec->n_allele > 1) {
            std::cout << rec->d.allele[1];
        } else {
            std::cout << ".";
        }

        std::cout << "\t.\t.\t";
        std::cout << "AC=" << Aa_count + 2 * AA_count << ";AN=" << 2 * n_samples;
        
        if (scaleDosage) {
            if (!geneMap.empty()) {
                std::cout << ";localMinDomDosage=" << localMinDomDosage
                         << ";localMaxDomDosage=" << localMaxDomDosage
                         << ";group=" << group;
            } else {
                std::cout << ";globalMinDomDosage=" << globalMinDomDosage
                         << ";globalMaxDomDosage=" << globalMaxDomDosage;
            }
        }

        if (allInfo) {
            std::cout << ";r=" << r
                     << ";h=" << h
                     << ";a=" << a;
        }

        std::cout << "\tDS";

        // Output dosage values for each sample
        for (int i = 0; i < n_samples; ++i) {
            double dosage = 0.0;
            int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
            float ds_value = has_ds ? ds_arr[i] : 0.0f;
            
            if ((has_gt && !bcf_gt_is_missing(gt_ptr[0]) && !bcf_gt_is_missing(gt_ptr[1])) || 
                (has_ds && !std::isnan(ds_value))) {
                
                int geno;
                if (has_gt && !bcf_gt_is_missing(gt_ptr[0]) && !bcf_gt_is_missing(gt_ptr[1])) {
                    int allele1 = bcf_gt_allele(gt_ptr[0]);
                    int allele2 = bcf_gt_allele(gt_ptr[1]);
                    if (allele1 == 0 && allele2 == 0) geno = 0;
                    else if (allele1 != allele2) geno = 1;
                    else if (allele1 > 0 && allele2 > 0) geno = 2;
                    else continue;
                } else {
                    geno = std::round(ds_value);
                    if (geno < 0 || geno > 2) continue;
                }

                if (mode == "dominance") {
                    // Dominance mode encoding
                    if (geno == 0) dosage = -h * a;
                    else if (geno == 1) dosage = 2 * a * r;
                    else if (geno == 2) dosage = -h * r;

                    if (scaleDosage) {
                        dosage = 2 * ((dosage - localMinDomDosage) / (localMaxDomDosage - localMinDomDosage));
                        if (dosage < 0) {
                            dosage = 0;
                        } else if (dosage > 2 + epsilon) {
                            std::cerr << "Error: Dosage value " << dosage << " exceeds 2 + epsilon (" 
                                     << (2 + epsilon) << ")" << std::endl;
                            exit(1);
                        }
                    }
                } else if (mode == "recessive") {
                    // Recessive mode encoding - set heterozygotes to 0, homozygotes unchanged
                    if (geno == 0) dosage = 0.0;
                    else if (geno == 1) dosage = 0.0; // Set heterozygotes to 0
                    else if (geno == 2) dosage = 2.0;
                }

                if (scalingFactor != 0) {
                    dosage *= scalingFactor;
                }
            }

            std::cout << "\t" << dosage;
        }

        std::cout << std::endl;
    }

    if (discardedVariantsCount > warningLimit) {
        std::cerr << "Note: Total number of variants discarded due to not being present in the gene map: " 
                  << discardedVariantsCount << std::endl;
    }

    if (mode == "dominance") {
        std::cerr << "Note: Total discarded variants without homozygous alternate alleles: " 
                  << countVariantsWithoutHomAlt << std::endl;
    }

    free(gt_arr);
    free(ds_arr);
    bcf_destroy(rec);
}

int main(int argc, char *argv[]) {
    std::string pathInput;
    std::string mode = "dominance";
    double scalingFactor = 1.0;
    bool scaleDosage = false;
    bool setVariantId = false;
    bool allInfo = false;
    std::string geneMapPath;

    if (!parseArguments(argc, argv, pathInput, mode, scalingFactor, scaleDosage, setVariantId, allInfo, geneMapPath)) {
        return 1;
    }
    
    // Validate mode
    if (mode != "dominance" && mode != "recessive") {
        std::cerr << "Error: Invalid mode '" << mode << "'. Only 'dominance' or 'recessive' modes are supported." << std::endl;
        return 1;
    }

    htsFile *fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: Cannot read header from VCF/BCF file." << std::endl;
        bcf_close(fp);
        return 1;
    }

    int n_samples = bcf_hdr_nsamples(hdr);

    double globalMinDomDosage = std::numeric_limits<double>::max();
    double globalMaxDomDosage = std::numeric_limits<double>::lowest();
    std::set<std::string> chromosomes;
    std::map<std::string, std::string> geneMap;
    std::map<std::string, std::pair<double, double>> geneDosages;

    if (!geneMapPath.empty()) {
        geneMap = readGeneMap(geneMapPath);
        for (const auto &pair : geneMap) {
            geneDosages[pair.second] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
        }
    }

    calculateGlobalAndGeneDosages(fp, hdr, globalMinDomDosage, globalMaxDomDosage, n_samples, chromosomes, geneMap, geneDosages);

    std::vector<std::string> sortedContigs = sortChromosomes(chromosomes);
    bcf_close(fp);

    // Reopen file for actual processing
    fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    hdr = bcf_hdr_read(fp);

    printHeader(hdr, sortedContigs, mode, globalMinDomDosage, globalMaxDomDosage, allInfo, scaleDosage, !geneMap.empty());

    processVcfFile(fp, hdr, mode, globalMinDomDosage, globalMaxDomDosage, allInfo, scaleDosage, setVariantId, scalingFactor, geneMap, geneDosages);

    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return 0;
}

