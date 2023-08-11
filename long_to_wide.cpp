#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <Long Format File> <Sample List File>" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    std::ifstream sampleFile(argv[2]);

    if (!inputFile.is_open() || !sampleFile.is_open()) {
        std::cerr << "Error opening the files." << std::endl;
        return 1;
    }

    // This map will store our wide-format data. It's a map of gene to a map of sample to dosage.
    std::map<std::string, std::map<std::string, int>> data;

    // Read the samples from the sample list file and store them in a set.
    std::set<std::string> samples;
    std::string sampleId;
    std::string sampleLine;
    while (std::getline(sampleFile, sampleLine)) {
       // Remove trailing and leading whitespaces
       sampleLine.erase(0, sampleLine.find_first_not_of(" \t\n\r"));
       sampleLine.erase(sampleLine.find_last_not_of(" \t\n\r") + 1);

       // If line is not empty, add the sample to our list
       if (!sampleLine.empty()) {
          samples.insert(sampleLine);
       }
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string sample, gene, call, variantDetails;
        int dosage;

        iss >> sample >> gene >> call >> dosage >> variantDetails;

        samples.insert(sample); // This ensures any new samples from the main table are also included.
        data[gene][sample] = dosage;
    }

    // Print header
    std::cout << "Gene";
    for (const auto& sample : samples) {
        std::cout << "\t" << sample;
    }
    std::cout << std::endl;

    // Print data
    for (const auto& genePair : data) {
       std::cout << genePair.first;
        for (const auto& sample : samples) {
           auto it = genePair.second.find(sample);
           std::cout << "\t";
            if (it != genePair.second.end()) {
                std::cout << it->second;
            } else {
                std::cout << "0";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}


