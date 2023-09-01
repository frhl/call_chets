#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <zlib.h>
#include <set>

int main(int argc, char* argv[]) {
    if (argc != 3) {
	std::cerr << "\nProgram: chet tools v0.0.1 (encode_vcf)\n" << std::endl;
	std::cerr << "Usage: " << argv[0] << " <Genotype File> <Mapping File>\n" << std::endl;
	std::cerr << "\nDescription:" << std::endl;
	std::cerr << "  Maps variants to genes using a self-generated mapping file and" << std::endl;
	std::cerr << "  takes in the output from 'get_non_ref_sites'. Will identify co-" << std::endl;
	std::cerr << "  occurring variants in cis and trans." << std::endl;
	std::cerr << "\nInput:" << std::endl; 
	std::cerr << "  <Genotype File> : .gz file containing sample genotype data," << std::endl;
	std::cerr << "                    typically output from 'get_non_ref_sites'.\n" << std::endl;
	std::cerr << "  <Mapping File>  : File mapping variants to genes. Expected" << std::endl;
	std::cerr << "                    to contain at least two columns with a header"<< std::endl;
	std::cerr << "                    in the format: variant, gene, info (optional).\n" << std::endl;
	std::cerr << "Output Format:" << std::endl;
	std::cerr << "  Sample Chromosome Gene Call Dosage Variant-Genotype..." << std::endl;
	std::cerr << "\nNotes:" << std::endl;
	std::cerr << "  Call values can be 'het', 'cis', 'hom', 'hom+het', 'chet." << std::endl; 
	return 1;
    }

    // expect .gz file no matter what
    gzFile genotypeFile = gzopen(argv[1], "rb");
    gzFile mappingFile = gzopen(argv[2], "rb");

    // Mapping variant to gene
    //std::map<std::string, std::string> variantToGene;
    std::map<std::string, std::pair<std::string, std::string>> variantToGene; 
    std::string variant, gene, thirdColumn;
    std::set<std::string> uniqueVariantsKept;
    std::set<std::string> uniqueVariantsDiscarded;
    char buf[1024];

    //std::cerr << "* Processing mapping and genotype file..\n";
    while (gzgets(mappingFile, buf, sizeof(buf))) {  // Use gzgets to read a line
        std::stringstream ss(buf);
        ss >> variant >> gene;
        if (ss >> thirdColumn) {  // Try to read the third column
                variantToGene[variant] = {gene, thirdColumn};
        } else {
                variantToGene[variant] = {gene, ""};  // No third column
        }
    }

    gzclose(mappingFile);

    // Process genotype file
    std::map<std::string, std::map<std::string, std::vector<std::string>>> sampleGeneVariants;
    std::string line;
    std::string sample, variantIndex, genotype;
    gzgets(genotypeFile, buf, sizeof(buf));

    int discardedVariants = 0;
    int keptVariants = 0; 

    //std::cerr << "* Processing genotype file.. This can take a while\n";
    while (gzgets(genotypeFile, buf, sizeof(buf))) {
        std::stringstream ss(buf);
        ss >> sample >> variantIndex >> variant >> genotype;
        if (variantToGene.find(variant) != variantToGene.end()) {
                std::string geneValue = variantToGene[variant].first;
                std::string modifiedVariant = variant + "-" + genotype;
                if (!variantToGene[variant].second.empty()) {
                        modifiedVariant += "-" + variantToGene[variant].second;
                }
                sampleGeneVariants[sample][geneValue].push_back(modifiedVariant);
                uniqueVariantsKept.insert(variant);
        } else {
                uniqueVariantsDiscarded.insert(variant);
        }
     }

    gzclose(genotypeFile);

    std::cerr << "Variants in mapping file (kept): " << uniqueVariantsKept.size() << std::endl;
    std::cerr << "Variants not in mapping file (Discarded): " << uniqueVariantsDiscarded.size() << std::endl;
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
	    	callValue = "hom+het";
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

