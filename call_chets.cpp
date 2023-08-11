#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>


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

    // Mapping variant to gene
    std::map<std::string, std::string> variantToGene;
    std::string variant, gene;

    while (mappingFile >> variant >> gene) {
        variantToGene[variant] = gene;
    }

    // Process genotype file
    std::map<std::string, std::map<std::string, std::vector<std::string>>> sampleGeneVariants;
    std::string line;
    // Skip the header of genotype file
    std::getline(genotypeFile, line);

    std::string sample, variantIndex, genotype;

    while (genotypeFile >> sample >> variantIndex >> variant >> genotype) {
        // Update our data structure with this variant-genotype
        sampleGeneVariants[sample][variantToGene[variant]].push_back(variant + "-" + genotype);
    }


    // Print the output
    for (const auto &samplePair : sampleGeneVariants) {
        for (const auto &genePair : samplePair.second) {
	    std::stringstream ss(variant);
            std::string chromosome;
            std::getline(ss, chromosome, ':');  // Split at the first ':'
            std::cout << samplePair.first << "\t" << chromosome << "\t" << genePair.first << "\t";
 
            // Counting configurations
            int count10 = 0;
            int count01 = 0;
            int count11 = 0;
            for (const auto &var : genePair.second) {
                size_t pos = var.find('-');
                std::string genotype = var.substr(pos + 1);
                if (genotype == "1|0") count10++;
                else if (genotype == "0|1") count01++;
                else if (genotype == "1|1") count11++;
            }
            
            //std::cout << count10 << ":" << count01 << ":" << count11 << "\t";
	    std::string callValue;
	    int dosage = 0;

            // Determine the "call" value
            if ((count10 > 0 && count01 == 0 && count11 == 0) || (count10 == 0 && count01 > 0 && count11 == 0)) {
                callValue = "het";
	        dosage = 1;
            } else if (count10 == 0 && count01 == 0 && count11 > 0) {
                callValue = "hom" ;
		dosage = 2;
	    } else if (count10 >= 0 && count01 >= 0 && count11 > 0) {
	    	callValue = "hom+het";
		dosage = 2;
	    } else if (count10 > 0 && count01 > 0 && count11 == 0) {
            	callValue="chet";
		dosage = 2;
	    } else {
            	callValue = "na";
		dosage = 0;
	    }

	    std::cout << callValue << "\t" << dosage << "\t"; // Print call and dosage

            for (size_t i = 0; i < genePair.second.size(); ++i) {
                std::cout << genePair.second[i];
                if (i != genePair.second.size() - 1) {
                    std::cout << ";";
                }
            }
            std::cout << std::endl;
        }
    }

    return 0;
}

