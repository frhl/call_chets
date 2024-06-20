#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <zlib.h>
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
    std::cerr << "  --all-info                : Include all calculated information (e.g., allele frequencies, min/max dosage) in the INFO field.\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << path << " --input example.vcf.gz --scale-dosages | gzip > output.vcf.gz \n";
}

int main(int argc, char *argv[]) {
    std::string pathInput;
    std::string mode = "dominance";
    double scalingFactor = 1.0;
    bool scaleDosage = false;
    bool allInfo = false;
    int countVariantsWithoutHomAlt = 0;
    const int limit = 5;
    const double epsilon = 1e-8;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        }
        else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
            pathInput = argv[++i];
        }
        else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
            mode = argv[++i];
        }
        else if (arg == "--scale-dosages-factor" && i + 1 < argc) {
            scalingFactor = std::stod(argv[++i]);
        }
        else if (arg == "--scale-dosages") {
            scaleDosage = true;
        }
        else if (arg == "--all-info") {
            allInfo = true;
        }
        else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (pathInput.empty()) {
        std::cerr << "Error! --input must be provided." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Open VCF/BCF file for first pass
    htsFile *fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: Cannot read header from VCF/BCF file." << std::endl;
        bcf_close(fp);
        return 1;
    }

    bcf1_t *rec = bcf_init();
    if (!rec) {
        std::cerr << "Error: Cannot initialize BCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return 1;
    }

    int *gt_arr = NULL, ngt_arr = 0; // Genotype array
    int n_samples = bcf_hdr_nsamples(hdr); // Number of samples in the VCF file

    // Variables for global min and max dominance dosages
    double globalMinDomDosage = std::numeric_limits<double>::max();
    double globalMaxDomDosage = std::numeric_limits<double>::lowest();

    // First pass: Calculate global min and max dominance dosages
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue; // Skip if no genotype information

        // Calculate AC, AN, aa_count, Aa_count, AA_count
        int aa_count = 0, Aa_count = 0, AA_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            // Point to the current sample
            int *gt_ptr = gt_arr + i * 2;

            // Check for missing genotype data
            if (bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1])) continue;

            int allele1 = bcf_gt_allele(gt_ptr[0]);
            int allele2 = bcf_gt_allele(gt_ptr[1]);

            // Count hom alt, het, hom ref
            if (allele1 == 0 && allele2 == 0) aa_count++;
            else if (allele1 != allele2) Aa_count++;
            else if (allele1 > 0 && allele2 > 0) AA_count++;
        }

        // Only proceed with encoding if at least one hom_alt is present
        if (AA_count > 0) {
            // Calculate allele frequencies based on genotype counts
            double r = static_cast<double>(aa_count) / n_samples;
            double h = static_cast<double>(Aa_count) / n_samples;
            double a = static_cast<double>(AA_count) / n_samples;

            // Calculate dominance dosages
            double dom_dosage_aa = -h * a;
            double dom_dosage_Aa = 2 * a * r;
            double dom_dosage_AA = -h * r;

            // Update global min and max dominance dosage values
            globalMinDomDosage = std::min({globalMinDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
            globalMaxDomDosage = std::max({globalMaxDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
        }
    }

    // Close the file from the first pass
    bcf_close(fp);

    // Open VCF/BCF file for second pass
    fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    // Read the VCF header
    hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: Cannot read header from VCF/BCF file." << std::endl;
        bcf_close(fp);
        return 1;
    }

    rec = bcf_init();
    if (!rec) {
        std::cerr << "Error: Cannot initialize BCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return 1;
    }

    gt_arr = NULL;
    ngt_arr = 0; // Genotype array
    n_samples = bcf_hdr_nsamples(hdr); // Number of samples in the VCF file

    // Print VCF header
    std::cout << "##fileformat=VCFv4.2\n";
    std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n";
    std::cout << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
    std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the alternate allele based on dominance model\">\n";
    if (allInfo == true) {
        std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of bi-allelic references (aa)\">\n";
        std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of heterozygotes (Aa)\">\n";
        std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of bi-allelic alternates (AA)\">\n";
    }
    if (scaleDosage) {
        std::cout << "##INFO=<ID=globalMinDomDosage,Number=1,Type=Float,Description=\"Global minimum dominance dosage\">\n";
        std::cout << "##INFO=<ID=globalMaxDomDosage,Number=1,Type=Float,Description=\"Global maximum dominance dosage\">\n";
    }
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 0; i < n_samples; i++) {
        std::cout << "\t" << bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
    }
    std::cout << "\n";

    // Second pass: Scale and print dosages
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue; // Skip if no genotype information

        // Calculate AC, AN, aa_count, Aa_count, AA_count
        int aa_count = 0, Aa_count = 0, AA_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            // Point to the current sample
            int *gt_ptr = gt_arr + i * 2;

            // Check for missing genotype data
            if (bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1])) continue;

            int allele1 = bcf_gt_allele(gt_ptr[0]);
            int allele2 = bcf_gt_allele(gt_ptr[1]);

            // Count hom alt, het, hom ref
            if (allele1 == 0 && allele2 == 0) aa_count++;
            else if (allele1 != allele2) Aa_count++;
            else if (allele1 > 0 && allele2 > 0) AA_count++;
        }

        // Only proceed with encoding if at least one hom_alt is present
        if (AA_count > 0) {
            // Calculate allele frequencies based on genotype counts
            double r = static_cast<double>(aa_count) / n_samples;
            double h = static_cast<double>(Aa_count) / n_samples;
            double a = static_cast<double>(AA_count) / n_samples;

            // Calculate dominance dosages
            double dom_dosage_aa = -h * a;
            double dom_dosage_Aa = 2 * a * r;
            double dom_dosage_AA = -h * r;

            // Output VCF line for this variant
            std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t.\t"
                      << rec->d.allele[0] << "\t"; // REF allele
            if (rec->n_allele > 1) {
                std::cout << rec->d.allele[1]; // ALT allele
            } else {
                std::cout << ".";
            }

            std::cout << "\t.\t.\t"; // QUAL and FILTER columns
            std::cout << "AC=" << Aa_count + 2 * AA_count << ";AN=" << 2 * n_samples; // INFO column
            if (scaleDosage) {
                std::cout << ";globalMinDomDosage=" << globalMinDomDosage
                          << ";globalMaxDomDosage=" << globalMaxDomDosage;
            }

            if ((mode == "dominance") & (allInfo == true)) {
                std::cout << ";r=" << r
                          << ";h=" << h
                          << ";a=" << a;
            }

            std::cout << "\tDS"; // FORMAT column

            // Output transformed genotypes based on dominance model
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

                    // Scale dosage to be between 0 and 2 using global min and max values
                    if (scaleDosage) {
                        dosage = 2 * ((dosage - globalMinDomDosage) / (globalMaxDomDosage - globalMinDomDosage));
                        // Ensure dosage is within [0, 2] and check for errors
                        if (dosage < 0) {
                            dosage = 0;
                        } else if (dosage > 2 + epsilon) {
                            std::cerr << "Error: Dosage value " << dosage << " exceeds 2 + epsilon (" << (2 + epsilon) << ")" << std::endl;
                            exit(1);
                        }
                    }

                    // Apply scaling factor
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
                          << rec->d.allele[0] << ":" << (rec->n_allele > 1 ? rec->d.allele[1] : ".") << "";
                std::cerr << "' has no homozygous alternate alleles. Skipping..\n";
            }
        }
    }
    std::cerr << "Note: Total discarded variants without homozygous alternate alleles: " << countVariantsWithoutHomAlt << std::endl;

    // Cleanup
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return 0;
}





