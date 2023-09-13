#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <zlib.h>
#include <set>

std::string getVersion() {
    std::ifstream versionFile(".version");
    if (!versionFile.is_open()) {
        std::cerr << "Warning: Unable to open .version file." << std::endl;
    }
    std::string version;
    std::getline(versionFile, version);
    versionFile.close();

    return version;
}



void printUsage(const char* path) {
	std::string version = getVersion();
	std::cerr << "\nProgram: chet tools v" << version << "\n" << std::endl;
	std::cerr << "Usage: " << path << " --geno <Genotype File> --map <Mapping File>\n" << std::endl;
	std::cerr << "\nDescription:" << std::endl;
	std::cerr << "  Call co-occuring variants using a variant-to-gene mapping file" << std::endl;
	std::cerr << "  alongside a file containing alternate genotypes." << std::endl;
	std::cerr << "\nInput:" << std::endl; 
	std::cerr << "  --geno/-g : .gz file containing alternate genotype data in the" << std::endl;
	std::cerr << "              form of sample id, variant id and genotype that is" << std::endl;
	std::cerr << "              seperated by whitespace. For example a line may look." << std::endl;
	std::cerr << "              like the following: 'Sample1 chr21:12314:C:T 1|0'." << std::endl;
	std::cerr << "  --map/-m  : File mapping variants to genes. This is expected" << std::endl;
	std::cerr << "              to contain at least two columns with a header"<< std::endl;
	std::cerr << "              in the format: variant, gene, info (optional)." << std::endl;
	std::cerr << "Output Format:" << std::endl;
	std::cerr << "  Sample Chromosome Gene Call Dosage Variant-Genotype..." << std::endl;
	std::cerr << "\nNotes:" << std::endl;
	std::cerr << "  See README for how to generate --geno file from a VCF/BCF." << std::endl; 

}


int main(int argc, char* argv[]) {

   std::string pathGeno;
   std::string pathMap;

   for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
       } else if ((arg == "--geno" || arg == "-g" ) && i + 1 < argc) {
            pathGeno = argv[++i];
       } else if ((arg == "--map" || arg == "-m" ) && i + 1 < argc) {
            pathMap = argv[++i];
       } else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
       }
    }

    // check if needed flags 
    if (pathGeno.empty() || pathMap.empty()) {
        std::cerr << "Error! Both --geno and --map must be provided." << std::endl;
        printUsage(argv[0]);
        return 1;
    } 

    // expect .gz file no matter what
    gzFile genotypeFile = gzopen(pathGeno.c_str(), "rb");
    gzFile mappingFile = gzopen(pathMap.c_str(), "rb");

    if (!genotypeFile) {
        std::cerr << "Error: Cannot open --geno file for reading: " << pathGeno << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    if (!mappingFile) {
        std::cerr << "Error: Cannot open --mapping file for reading: " << pathMap << std::endl;
        printUsage(argv[0]);
        return 1;
    }


    // Mapping variant to gene
    std::map<std::string, std::vector<std::pair<std::string, std::string>>> variantToGene;
    std::string variant, gene, thirdColumn;
    std::set<std::string> uniqueVariantsKept;
    std::set<std::string> uniqueVariantsDiscarded;
    std::set<std::string> multiGeneVariants;
    char buf[1024];

    //std::cerr << "* Processing mapping and genotype file..\n";
    while (gzgets(mappingFile, buf, sizeof(buf))) {  // Use gzgets to read a line
        std::stringstream ss(buf);
        ss >> variant >> gene;
        std::pair<std::string, std::string> geneInfo = ss >> thirdColumn ? std::make_pair(gene, thirdColumn) : std::make_pair(gene, "");
        variantToGene[variant].push_back(geneInfo);	
    }

    gzclose(mappingFile);

    // Process genotype file
    std::map<std::string, std::map<std::string, std::vector<std::string>>> sampleGeneVariants;
    std::string line;
    std::string sample,  genotype;
    gzgets(genotypeFile, buf, sizeof(buf));

    int discardedVariants = 0;
    int keptVariants = 0; 

    //std::cerr << "* Processing genotype file.. This can take a while\n";
    while (gzgets(genotypeFile, buf, sizeof(buf))) {
        std::stringstream ss(buf);
        ss >> sample >> variant >> genotype;

        // ensure that the input genotype is always alt
        if (genotype != "1|0" && genotype != "0|1" && genotype != "1|1") {
               std::cerr << "Error: Unexpected genotype value (" << genotype << ") in sample " << sample << " for variant " << variant << "." << std::endl;
               return 1;
	}

	// continue to find corresponding gene
        if (variantToGene.find(variant) != variantToGene.end()) {
		if (variantToGene[variant].size() > 1) {
			multiGeneVariants.insert(variant);
		}
		for (const auto& geneInfoPair : variantToGene[variant]) {
			std::string geneValue = geneInfoPair.first;
			std::string modifiedVariant = variant + "-" + genotype;
			if (!geneInfoPair.second.empty()) {
			    modifiedVariant += "-" + geneInfoPair.second;
			}
			sampleGeneVariants[sample][geneValue].push_back(modifiedVariant);
			uniqueVariantsKept.insert(variant);
		}
    	} else {
                uniqueVariantsDiscarded.insert(variant);
        }
     }

    gzclose(genotypeFile);

    // Check if there are no matching variants
    if (uniqueVariantsKept.empty()) {
        std::cerr << "Error: No matching variants found between the --map and --geno file!" << std::endl;
        return 1;
    }

    std::cerr << "Variants in mapping file (kept): " << uniqueVariantsKept.size() << std::endl;
    std::cerr << "Variants not in mapping file (Discarded): " << uniqueVariantsDiscarded.size() << std::endl;
    std::cerr << "Variants mapping to more than one gene: " << multiGeneVariants.size() << std::endl;
    std::cerr << "* Generating sample-gene-variant file.." << std::endl;

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
                std::string genotype = var.substr(pos + 1, 3);
                if (genotype == "1|0") count10++;
                else if (genotype == "0|1") count01++;
                else if (genotype == "1|1") count11++;
            }
            
            //std::cout << count10 << ":" << count01 << ":" << count11 << "\t";
	    std::string callValue;
	    int dosage = 0;

            // Determine the "call" value
            if ((count10 == 1 && count01 == 0 && count11 == 0) || (count10 == 0 && count01 == 1 && count11 == 0)) {
                callValue = "het";
	        dosage = 1;
            } else if (count10 == 0 && count01 > 1 && count11 == 0) {
                callValue = "cis";
                dosage = 1;
	    } else if (count10 > 1 && count01 == 0 && count11 == 0) {
                callValue = "cis";
                dosage = 1;
	    } else if (count10 == 0 && count01 == 0 && count11 > 0) {
                callValue = "hom" ;
		dosage = 2;
	    } else if (count10 >= 0 && count01 >= 0 && count11 > 0) {
	    	callValue = "hom"; // previous "hom+het"
		dosage = 2;
	    } else if (count10 > 0 && count01 > 0 && count11 == 0) {
            	callValue="chet";
		dosage = 2;
	    } else {
            	callValue = "na";
		dosage = 0;
	    }
	    std::cout << callValue << "\t" << dosage << "\t";
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



