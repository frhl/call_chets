#include <iostream>
#include <fstream>
#include <string>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <cmath>

int missing_count = 0;
bool verbose=false;

void printUsage(const char *path) {
    std::cerr << "\nProgram: VCF PP Filter\n\n";
    std::cerr << "Usage: " << path << " --input <VCF/BCF file> --output <output file> --pp-threshold <threshold>\n\n";

    std::cerr << "Description:\n";
    std::cerr << "  This program filters variants in a VCF/BCF file based on the PP (posterior probability) field.\n";
    std::cerr << "  If the PP value for a genotype is below the specified threshold, that genotype is set to missing.\n\n";

    std::cerr << "Arguments:\n";
    std::cerr << "  --input/-i <VCF/BCF file>   : Required. Path to the VCF/BCF file to be processed.\n";
    std::cerr << "  --output/-o <output file>   : Required. Path to the output file where the filtered VCF/BCF data will be written.\n";
    std::cerr << "  --pp-threshold/-p <float>   : Required. Threshold value for the PP field. Genotypes with PP below this value will be set to missing.\n\n";
    std::cerr << "  --verbose/-v                 : Optional. Enables verbose mode, printing progress every 100 variants processed.\n\n";

    std::cerr << "Notes:\n";
    std::cerr << "  Ensure that the VCF/BCF file is properly formatted.\n";
    std::cerr << "  The program reads the PP field from the FORMAT column of the VCF/BCF file.\n";
    std::cerr << "  The output file format (VCF or BCF) is determined by the extension of the output file path.\n";
}

void set_genotype_to_missing(bcf_fmt_t *fmt_pp, int i) {
    float missing_value = bcf_float_missing;
    ((float*)fmt_pp->p)[i] = missing_value;
    missing_count++;
}

int main(int argc, char *argv[]) {
    std::string inputFilePath;
    std::string outputFilePath;
    float ppThreshold = 0.0;
    int variant_count = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
            inputFilePath = argv[++i];
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            outputFilePath = argv[++i];
        } else if ((arg == "--pp-threshold" || arg == "-p") && i + 1 < argc) {
            ppThreshold = std::stof(argv[++i]);
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = true;
        } else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (inputFilePath.empty() || outputFilePath.empty() || ppThreshold <= 0) {
        std::cerr << "Error! Missing required arguments." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    htsFile *inFile = bcf_open(inputFilePath.c_str(), "r");
    if (!inFile) {
        std::cerr << "Could not open input file: " << inputFilePath << std::endl;
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(inFile);

    htsFile *outFile;
    if (outputFilePath.size() >= 4 && outputFilePath.substr(outputFilePath.size() - 4) == ".bcf") {
        outFile = bcf_open(outputFilePath.c_str(), "wb");  // BCF format (binary)
    } else if (outputFilePath.size() >= 7 && outputFilePath.substr(outputFilePath.size() - 7) == ".vcf.gz") {
        outFile = bcf_open(outputFilePath.c_str(), "wz");  // Compressed VCF format
    } else {
        outFile = bcf_open(outputFilePath.c_str(), "w");   // Non-compressed VCF format
    }


    if (!outFile) {
        std::cerr << "Could not open output file: " << outputFilePath << std::endl;
        return 1;
    }

    bcf_hdr_write(outFile, hdr);

   bcf1_t *rec = bcf_init();

    while (bcf_read(inFile, hdr, rec) == 0) {
        variant_count++;
        bcf_unpack(rec, BCF_UN_ALL);

        // Get the PP field for each sample
        bcf_fmt_t *fmt_pp = bcf_get_fmt(hdr, rec, "PP");
        bcf_fmt_t *fmt_gt = bcf_get_fmt(hdr, rec, "GT");
        if (fmt_pp) {
            int num_samples = bcf_hdr_nsamples(hdr);
            for (int i = 0; i < num_samples; ++i) {
                std::string genotype = std::to_string(bcf_gt_allele(fmt_gt->p[i*2])) + "|"
                                   + std::to_string(bcf_gt_allele(fmt_gt->p[i*2+1])); 
            
                // Accessing the PP value for each sample
                if (genotype != "0|0" && genotype != "1|1") {
                    float* pp_array = (float*)(fmt_pp->p);
                    float pp_value = pp_array[i];

                    // Check PP value and set genotype to missing if below the threshold
                    if (pp_value != bcf_float_missing && pp_value < ppThreshold) {
                        set_genotype_to_missing(fmt_pp, i);
                        missing_count++;
                    }
	        }
            }
            if (verbose && variant_count % 100 == 0) {
                std::cerr << "Processed " << variant_count << " variants. "
                          << "Heterozygotes set to missing so far: " << missing_count << std::endl;
            }
        }

        // Write the processed record to output file
        bcf_write(outFile, hdr, rec);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(inFile);
    bcf_close(outFile);

    std::cerr << "Total genotypes set to missing: " << missing_count << std::endl;

    return 0;
}

