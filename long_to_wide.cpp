#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <Input File> <Output File>" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    std::ofstream outputFile(argv[2]);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    // Data structure to store dosages
    std::map<std::string, std::map<std::string, int>> geneSampleDosage;
    std::set<std::string> allSamples;
    std::set<std::string> allGenes;

    // Read and parse the input file
    std::string line, sample, gene, config, call, dosage_str;
    int dosage;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        iss >> sample >> gene >> config >> dosage_str >> call;
        dosage = std::stoi(dosage_str);

        allSamples.insert(sample);
        allGenes.insert(gene);

        geneSampleDosage[gene][sample] = dosage;
    }

    // Initialize with zeros for every sample-gene combination not already present
    for (const auto &gene : allGenes) {
        for (const auto &sample : allSamples) {
            if (geneSampleDosage[gene].find(sample) == geneSampleDosage[gene].end()) {
                geneSampleDosage[gene][sample] = 0;
            }
        }
    }

    // Write the matrix to the output file
    outputFile << "Gene";
    for (const auto &sample : allSamples) {
        outputFile << "\t" << sample;
    }
    outputFile << std::endl;

    for (const auto &genePair : geneSampleDosage) {
        outputFile << genePair.first;
        for (const auto &sample : allSamples) {
            outputFile << "\t" << genePair.second.at(sample);
        }
        outputFile << std::endl;
    }

    return 0;
}

