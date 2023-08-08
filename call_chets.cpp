#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <limits>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <Genotype File> <Mapping File>" << std::endl;
        return 1;
    }

    std::ifstream genotypeFile(argv[1]);
    std::ifstream mappingFile(argv[2]);

    if (!genotypeFile.is_open() || !mappingFile.is_open()) {
        std::cerr << "Error opening files." << std::endl;
        return 1;
    }

    // Map variant to gene
    std::map<std::string, std::string> variantToGene;
    std::string variant, gene;

    while (mappingFile >> variant >> gene) {
        variantToGene[variant] = gene;
    }

    // Process genotype file
    struct GeneData {
        std::vector<std::string> variants;
        int config10 = 0;
        int config01 = 0;
        int config11 = 0;
    };

    std::map<std::string, std::map<std::string, GeneData>> sampleGeneVariants;

    std::string line;
    // Skip the header line of genotype file
    std::getline(genotypeFile, line);

    std::string sample, variantIndex, genotype;
    while (genotypeFile >> sample >> variantIndex >> variant >> genotype) {
        // Populate the data structure with variant-genotype
        sampleGeneVariants[sample][variantToGene[variant]].variants.push_back(variant + "-" + genotype);
        if (genotype == "1|0") {
            sampleGeneVariants[sample][variantToGene[variant]].config10++;
        } else if (genotype == "0|1") {
            sampleGeneVariants[sample][variantToGene[variant]].config01++;
        } else if (genotype == "1|1") {
            sampleGeneVariants[sample][variantToGene[variant]].config11++;
        }
    }

    // Print the output to stdout
    for (const auto &samplePair : sampleGeneVariants) {
        for (const auto &genePair : samplePair.second) {
            std::cout << samplePair.first << "\t" << genePair.first;
            
            bool isFirst = true;
            for (const auto &var : genePair.second.variants) {
                if (isFirst) {
                    std::cout << "\t" << var;
                    isFirst = false;
                } else {
                    std::cout << ";" << var;
                }
            }
            std::cout << "\t" << genePair.second.config10 << ":" << genePair.second.config01 << ":" << genePair.second.config11 << std::endl;
        }
    }

    return 0;
}

