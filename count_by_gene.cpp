#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <zlib.h> // For gzipped file operations
#include <limits>

// Define a structure to hold genotype counts for a gene
struct GenotypeCounts {
    int homRefCount = 0;
    int hetCount = 0;
    int homAltCount = 0;

    void addCounts(int homRef, int het, int homAlt) {
        homRefCount += homRef;
        hetCount += het;
        homAltCount += homAlt;
    }
};

void printUsage(const char* path) {
    std::cerr << "\nProgram: VCF Gene Counter\n";
    std::cerr << "Usage: " << path << " --input <input.vcf.gz> --map <variant_to_gene.map.gz> [--max-af <max allele frequency>] [--output <output.txt>]\n";
    std::cerr << "Options:\n";
    std::cerr << "  --input/-i                : Input VCF file.\n";
    std::cerr << "  --map/-m                  : Gzipped variant to gene mapping file.\n";
    std::cerr << "  --max-af                  : Maximum allele frequency threshold (default includes all).\n";
}

int main(int argc, char* argv[]) {
    std::string pathInput, outputPath = "gene_counts.txt", mappingFilePath;
    float maxAf = std::numeric_limits<float>::max();
    int totalSamples = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
            pathInput = argv[++i];
        } else if ((arg == "--map" || arg == "-m") && i + 1 < argc) {
            mappingFilePath = argv[++i];
        } else if (arg == "--max-af" && i + 1 < argc) {
            maxAf = std::stof(argv[++i]);
        } else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (pathInput.empty() || mappingFilePath.empty()) {
        std::cerr << "Error! --input and --map must be provided." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Load mapping file
    gzFile mappingFile = gzopen(mappingFilePath.c_str(), "r");
    if (!mappingFile) {
        std::cerr << "Error: Cannot open mapping file for reading: " << mappingFilePath << std::endl;
        return 1;
    }

    char buf[1024];
    std::string variant, gene;
    std::map<std::string, std::vector<std::string>> variantToGene;
    bool isFirstLineMappingFile = true;
    while (gzgets(mappingFile, buf, sizeof(buf))) {
        std::stringstream ss(buf);
        ss >> variant >> gene;

        if(ss.fail() || variant.empty() || gene.empty()) {
            if(!isFirstLineMappingFile) {
                std::cerr << "Error: Failed to extract two columns from line (variant gene): " << buf << ". Please, fix this line in --gene-map and retry." << std::endl;
                gzclose(mappingFile);
                return 1;
            }
        } else {
            variantToGene[variant].push_back(gene);
        }
        isFirstLineMappingFile = false;
    }
    gzclose(mappingFile);

    // Open the VCF file
    htsFile *vcfFile = bcf_open(pathInput.c_str(), "r");
    if (!vcfFile) {
        std::cerr << "Error: Cannot open VCF file for reading: " << pathInput << std::endl;
        return 1;
    }

    bcf_hdr_t *header = bcf_hdr_read(vcfFile);
    if (!header) {
        std::cerr << "Error: Cannot read header from VCF file." << std::endl;
        bcf_close(vcfFile);
        return 1;
    }

    bcf1_t *record = bcf_init();
    if (!record) {
        std::cerr << "Error: Cannot initialize BCF record." << std::endl;
        bcf_hdr_destroy(header);
        bcf_close(vcfFile);
        return 1;
    }

    totalSamples = bcf_hdr_nsamples(header);
    std::map<std::string, GenotypeCounts> geneCounts;
    int ngt, *gt_arr = nullptr, ngt_arr = 0;

    // Iterate over all records (variants) in the VCF file
    while (bcf_read(vcfFile, header, record) == 0) {
        bcf_unpack(record, BCF_UN_ALL);
        std::string chrom = bcf_hdr_id2name(header, record->rid);
        int pos = record->pos + 1;
        std::string ref = record->d.allele[0];
        
        // If there are alternate alleles, just consider the first one for simplicity
        std::string alt = record->n_allele > 1 ? record->d.allele[1] : ".";
        std::string variantKey = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;

        // Check if the current variant is in our mapping
        if (variantToGene.find(variantKey) != variantToGene.end()) {
            // Get genotype information
            ngt = bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
            if (ngt <= 0) continue; // Skip if no genotype information

            int n_samples = bcf_hdr_nsamples(header);
            int homRefCount = 0, hetCount = 0, homAltCount = 0;

            for (int i = 0; i < n_samples; ++i) {
                // Check each allele of the genotype (diploid assumption: two alleles per genotype)
                if (gt_arr[i * 2] == bcf_gt_missing || gt_arr[i * 2 + 1] == bcf_gt_missing) continue; // Skip missing genotypes
                
                int allele1 = bcf_gt_allele(gt_arr[i * 2]);
                int allele2 = bcf_gt_allele(gt_arr[i * 2 + 1]);

                if (allele1 == 0 && allele2 == 0) homRefCount++;
                else if (allele1 != allele2) hetCount++;
                else homAltCount++; // Assuming only one alt allele for simplicity
            }

	    // Calculate allele count (AC) and allele frequency (AF)
            int alleleCount = 2 * homAltCount + hetCount;
            float alleleFrequency = static_cast<float>(alleleCount) / (2 * totalSamples);

            // Skip variant if allele frequency exceeds maxAf
            if (alleleFrequency > maxAf) continue;

            // Add counts to all genes this variant is associated with
            for (const auto& gene : variantToGene[variantKey]) {
                geneCounts[gene].addCounts(homRefCount, hetCount, homAltCount);
            }
        }
    }

    // Free resources
    if (gt_arr) free(gt_arr);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(vcfFile);

    std::cout << "region\taa\tAa\tAA\tsamples\n";
    for (const auto& pair : geneCounts) {
	std::cout << pair.first << "\t" 
                << pair.second.homRefCount << "\t" 
                << pair.second.hetCount << "\t" 
                << pair.second.homAltCount << "\t"
		<< totalSamples << "\n";
    }

    return 0;
}
       


