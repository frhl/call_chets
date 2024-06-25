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

void printUsage(const char *path) {
    std::cerr << "\nProgram: transform genotypes to dominance encoding\n";
    std::cerr << "Usage: " << path << " --input <input.vcf.gz> --mode dominance [--output <output.vcf.gz>] [--format vcf|bcf|vcf.gz] [--scale-dosages] [--scale-dosages-factor <factor>] [--all-info]\n";
    std::cerr << "Options:\n";
    std::cerr << "  --input/-i                : Input VCF/BCF file.\n";
    std::cerr << "  --mode/-m                 : Specify mode for processing genotypes. Supports 'dominance'.\n";
    std::cerr << "  --scale-dosages           : Enable dosage scaling. Scales dosages to be between 0 and 2. \n";
    std::cerr << "  --scale-dosages-factor    : Apply a scaling factor to the dosages. Default is 1.0.\n";
    std::cerr << "  --set-variant-id          : Set the variant ID to chr:pos:ref:alt.\n";
    std::cerr << "  --all-info                : Include all calculated information (e.g., allele frequencies, min/max dosage) in the INFO field.\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << path << " --input example.vcf.gz --scale-dosages | gzip > output.vcf.gz \n";
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

bool parseArguments(int argc, char *argv[], std::string &pathInput, std::string &mode, double &scalingFactor, bool &scaleDosage, bool &setVariantId, bool &allInfo) {
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

void calculateGlobalDosages(htsFile *fp, bcf_hdr_t *hdr, double &globalMinDomDosage, double &globalMaxDomDosage, int n_samples, std::set<std::string> &chromosomes) {
    bcf1_t *rec = bcf_init();
    int *gt_arr = NULL, ngt_arr = 0;

    while (bcf_read(fp, hdr, rec) == 0) {
        chromosomes.insert(bcf_hdr_id2name(hdr, rec->rid));
        bcf_unpack(rec, BCF_UN_STR);
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue;

        int aa_count = 0, Aa_count = 0, AA_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            int *gt_ptr = gt_arr + i * 2;
            if (bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1])) continue;

            int allele1 = bcf_gt_allele(gt_ptr[0]);
            int allele2 = bcf_gt_allele(gt_ptr[1]);
            if (allele1 == 0 && allele2 == 0) aa_count++;
            else if (allele1 != allele2) Aa_count++;
            else if (allele1 > 0 && allele2 > 0) AA_count++;
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
        }
    }

    free(gt_arr);
    bcf_destroy(rec);
}

void printHeader(const bcf_hdr_t *hdr, const std::vector<std::string> &sortedContigs, const std::string &mode, double globalMinDomDosage, double globalMaxDomDosage, bool allInfo, bool scaleDosage) {
    int n_samples = bcf_hdr_nsamples(hdr);

    std::cout << "##fileformat=VCFv4.2\n";
    // output chromosomes
    for (const auto& chr : sortedContigs) {
        std::cout << "##contig=<ID=" << chr << ">\n";
    }
    std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n";
    std::cout << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
    std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the alternate allele based on dominance model\">\n";
    if (allInfo) {
        std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of bi-allelic references (aa)\">\n";
        std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of heterozygotes (Aa)\">\n";
        std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of bi-allelic alternates (AA)\">\n";
    }
    if (scaleDosage) {
        std::cout << "##INFO=<ID=globalMinDomDosage,Number=1,Type=Float,Description=\"Global minimum dominance dosage\">\n";
        std::cout << "##INFO=<ID=globalMaxDomDosage,Number=1,Type=Float,Description=\"Global maximum dominance dosage\">\n";
    }
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 0; i < n_samples; ++i) {
        std::cout << "\t" << bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
    }
    std::cout << "\n";
}

void processVcfFile(htsFile *fp, bcf_hdr_t *hdr, const std::string &mode, double globalMinDomDosage, double globalMaxDomDosage, bool allInfo, bool scaleDosage, bool setVariantId, double scalingFactor) {
    bcf1_t *rec = bcf_init();
    int *gt_arr = NULL, ngt_arr = 0;
    int n_samples = bcf_hdr_nsamples(hdr);
    const double epsilon = 1e-8;
    int countVariantsWithoutHomAlt = 0;
    const int limit = 5;

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue;

        int aa_count = 0, Aa_count = 0, AA_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            int *gt_ptr = gt_arr + i * 2;
            if (bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1])) continue;

            int allele1 = bcf_gt_allele(gt_ptr[0]);
            int allele2 = bcf_gt_allele(gt_ptr[1]);
            if (allele1 == 0 && allele2 == 0) aa_count++;
            else if (allele1 != allele2) Aa_count++;
            else if (allele1 > 0 && allele2 > 0) AA_count++;
        }

        if (AA_count > 0) {
            double r = static_cast<double>(aa_count) / n_samples;
            double h = static_cast<double>(Aa_count) / n_samples;
            double a = static_cast<double>(AA_count) / n_samples;

            double dom_dosage_aa = -h * a;
            double dom_dosage_Aa = 2 * a * r;
            double dom_dosage_AA = -h * r;

            std::string variantId = ".";
            if (setVariantId) {
                variantId = bcf_hdr_id2name(hdr, rec->rid) + std::string(":") + std::to_string(rec->pos + 1) + std::string(":") + rec->d.allele[0] + std::string(":");
                if (rec->n_allele > 1) {
                    variantId += rec->d.allele[1];
                } else {
                    variantId += ".";
                }
            }

            std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << variantId << "\t"
                      << rec->d.allele[0] << "\t";
            if (rec->n_allele > 1) {
                std::cout << rec->d.allele[1];
            } else {
                std::cout << ".";
            }

            std::cout << "\t.\t.\t";
            std::cout << "AC=" << Aa_count + 2 * AA_count << ";AN=" << 2 * n_samples;
            if (scaleDosage) {
                std::cout << ";globalMinDomDosage=" << globalMinDomDosage
                          << ";globalMaxDomDosage=" << globalMaxDomDosage;
            }

            if (allInfo) {
                std::cout << ";r=" << r
                          << ";h=" << h
                          << ";a=" << a;
            }

            std::cout << "\tDS";

            for (int i = 0; i < n_samples; ++i) {
                int *gt_ptr = gt_arr + i * 2;
                double dosage = 0.0;

                if (!bcf_gt_is_missing(gt_ptr[0]) && !bcf_gt_is_missing(gt_ptr[1])) {
                    int allele1 = bcf_gt_allele(gt_ptr[0]);
                    int allele2 = bcf_gt_allele(gt_ptr[1]);

                    if (allele1 == 0 && allele2 == 0) {
                        dosage = -h * a;
                    } else if (allele1 != allele2) {
                        dosage = 2 * a * r;
                    } else if (allele1 > 0 && allele2 > 0) {
                        dosage = -h * r;
                    }

                    if (scaleDosage) {
                        dosage = 2 * ((dosage - globalMinDomDosage) / (globalMaxDomDosage - globalMinDomDosage));
                        if (dosage < 0) {
                            dosage = 0;
                        } else if (dosage > 2 + epsilon) {
                            std::cerr << "Error: Dosage value " << dosage << " exceeds 2 + epsilon (" << (2 + epsilon) << ")" << std::endl;
                            exit(1);
                        }
                    }

                    if (scalingFactor != 0) {
                        dosage *= scalingFactor;
                    }
                }

                std::cout << "\t" << dosage;
            }

            std::cout << std::endl;
        } else {
            countVariantsWithoutHomAlt++;
            if (countVariantsWithoutHomAlt <= limit) {
                std::cerr << "variant '" << bcf_hdr_id2name(hdr, rec->rid) << ":" << (rec->pos + 1) << ":"
                          << rec->d.allele[0] << ":" << (rec->n_allele > 1 ? rec->d.allele[1] : ".") << "'";
                std::cerr << " has no homozygous alternate alleles. Skipping..\n";
            }
        }
    }
    std::cerr << "Note: Total discarded variants without homozygous alternate alleles: " << countVariantsWithoutHomAlt << std::endl;

    free(gt_arr);
    bcf_destroy(rec);
}

int main(int argc, char *argv[]) {
    std::string pathInput;
    std::string mode = "dominance";
    double scalingFactor = 1.0;
    bool scaleDosage = false;
    bool setVariantId = false;
    bool allInfo = false;

    if (!parseArguments(argc, argv, pathInput, mode, scalingFactor, scaleDosage, setVariantId, allInfo)) {
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
    calculateGlobalDosages(fp, hdr, globalMinDomDosage, globalMaxDomDosage, n_samples, chromosomes);

    std::vector<std::string> sortedContigs = sortChromosomes(chromosomes);
    bcf_close(fp);

    // Reopen file for actual processing
    fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    hdr = bcf_hdr_read(fp);

    printHeader(hdr, sortedContigs, mode, globalMinDomDosage, globalMaxDomDosage, allInfo, scaleDosage);

    processVcfFile(fp, hdr, mode, globalMinDomDosage, globalMaxDomDosage, allInfo, scaleDosage, setVariantId, scalingFactor);

    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return 0;
}

