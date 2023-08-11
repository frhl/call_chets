#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <Long Format File> <Sample List File> <Output File>" << std::endl;
        return 1;
    }

    std::ifstream longFile(argv[1]);
    std::ifstream sampleFile(argv[2]);
    std::ofstream outputFile(argv[3]);

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

    // Process long format file
    std::string line, sample, gene, configuration, variantInfo;
    int dosage;
    while (std::getline(longFile, line)) {
        std::istringstream iss(line);
        iss >> sample >> gene >> configuration >> dosage >> variantInfo;

        geneSampleDosage[gene][sample] = dosage;
    }

    // Print the output header
    outputFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sampleName : samples) {
        outputFile << "\t" << sampleName;
    }
    outputFile << std::endl;

    // Print the output data
    for (const auto &genePair : geneSampleDosage) {
        outputFile << ".\t.\t" << genePair.first << "\t.\t.\t.\t.\t.\t.";  // Repeated columns filled with dots, except for ID which gets the gene name

        for (const auto& sample : samples) {
            outputFile << "\t";

            if (genePair.second.find(sample) != genePair.second.end()) {
                outputFile << genePair.second.at(sample);
            } else {
                outputFile << "0";  // Default to 0 if sample-gene combo doesn't exist
            }
        }
        outputFile << std::endl;
    }

    return 0;
}

