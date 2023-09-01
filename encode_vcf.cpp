#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "\nchet tools v0.0.1\n" << std::endl;
        std::cerr << "\nUsage: " << argv[0] << " <input> <samples> <additive|recessive>\n" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << " <input> \tOutput file from 'get_non_ref_sites'." << std::endl;
        std::cerr << " <samples> \tFile path to list of samples seperated by lines. No header required." << std::endl;
        std::cerr << " <mode> \tEither 'additive' to encode dosages as 0,1 and 2. Heterozygotes and cis" << std::endl;
        std::cerr << "        \tvariants with be kept. Use 'recessive' to encode dosages as 0 and 2 and" << std::endl;
        std::cerr << "        \tThus only keep compound heterozygotes and homozygotes downstream.\n" << std::endl;
	return 1;
    }

    std::ifstream longFile(argv[1]);
    std::ifstream sampleFile(argv[2]);
    std::string mode(argv[3]);

    if (mode != "additive" && mode != "recessive") {
        std::cerr << "Error: Invalid dosage encoding mode provided. Expecting 'additive' or 'recessive'." << std::endl;
        return 1;
    }


    // Collect all samples into a set
    std::set<std::string> samples;
    std::string sampleLine;
    while (std::getline(sampleFile, sampleLine)) {
        // Remove trailing and leading whitespaces
        sampleLine.erase(0, sampleLine.find_first_not_of(" \t\n\r"));
        sampleLine.erase(sampleLine.find_last_not_of(" \t\n\r") + 1);

        // If line is not empty, add the sample to our set
        if (!sampleLine.empty()) {
            samples.insert(sampleLine);
        }
    }

    std::map<std::string, std::map<std::string, int>> geneSampleDosage;

    // Map the chrom
    std::map<std::string, std::string> geneToChromosome;

    // Process long format file
    std::string line, sample, chromosome, gene, configuration, variantInfo;
    std::set<std::string> contigs;
    int dosage;
    while (std::getline(longFile, line)) {
        std::istringstream iss(line);
        iss >> sample >> chromosome >> gene >> configuration >> dosage >> variantInfo;

	// only 2 allowed in recessive encoding
	if (mode == "recessive" && dosage == 1) {
            dosage = 0;
        }
        geneToChromosome[gene] = chromosome;
        geneSampleDosage[gene][sample] = dosage;
	contigs.insert(chromosome);
    }

    // Print the output header
    std::cout << "##fileformat=VCFv4.2\n";
    std::cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    for (const auto& chr : contigs) {
        std::cout << "##contig=<ID=" << chr << ">\n";
    }
    std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sampleName : samples) {
        std::cout << "\t" << sampleName;
    }
    std::cout << std::endl;
    int rowIndex = 1;
    // Print the output data
    for (const auto &genePair : geneSampleDosage) {
        std::cout << geneToChromosome[genePair.first] << "\t" << rowIndex << "\t" << genePair.first << "\tA\tB\t.\t.\t.\tDS";  // Repeated columns filled with dots, except for ID which gets the gene name

        for (const auto& sample : samples) {
            std::cout << "\t";

            if (genePair.second.find(sample) != genePair.second.end()) {
                std::cout << genePair.second.at(sample);
            } else {
                std::cout << "0";  // Default to 0 if sample-gene combo doesn't exist
            }
        }
        std::cout << std::endl;
        rowIndex++;
    }

    return 0;
}

