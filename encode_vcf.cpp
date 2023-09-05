#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <zlib.h>
#include <vector>
#include <sstream>
#include <string>

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
	std::cerr << "\n\nUsage: " << path << " --input <input> --samples <samples> --mode [<additive|recessive>]" << std::endl;
	std::cerr << "\nDescription:" << std::endl;
	std::cerr << "  Converts 'call_chets' output to VCF for downstream analysis."<< std::endl;
	std::cerr << "  Results are streamed to standard output." << std::endl;
	std::cerr << "\nOptions:" << std::endl;
	std::cerr << "  --input/-i   : Output from 'call_chets'.\n";
	std::cerr << "  --samples/-s : List of samples. One per line. No header.\n";
	std::cerr << "  --mode/-m    : Use 'additive' for dosages of 0, 1, and 2. This keeps" << std::endl;
	std::cerr << "                 heterozygotes and cis variants. Use 'recessive' for"  << std::endl;
	std::cerr << "                 dosages of 0 and 2, targeting compound heterozygotes"  << std::endl;
	std::cerr << "                 and homozygotes. (Default='additive')"  << std::endl;
	std::cerr << "\nExample:" << std::endl;
	std::cerr << "  ./encode_vcf called_chets.txt.gz samples.txt additive | bgzip > out.vcf.gz\n\n";
}

int main(int argc, char* argv[]) {
    
    std::string pathInput;
    std::string pathSamples;
    std::string mode = "additive";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
	    return 0;
        } else if ((arg == "--input" || arg == "-i" ) && i + 1 < argc) {
            pathInput = argv[++i];
        } else if ((arg == "--samples" || arg == "-s" ) && i + 1 < argc) {
            pathSamples = argv[++i];
        } else if ((arg == "--mode" || arg == "-m" ) && i + 1 < argc) {
            mode = argv[++i];
       } else {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (pathInput.empty() || pathSamples.empty()) {
        std::cerr << "Error! Both --input and --samples must be provided." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    gzFile longFile = gzopen(pathInput.c_str(), "rb");
    std::ifstream sampleFile(pathSamples);

    if (!longFile) {
        std::cerr << "Error: Cannot open gzipped <input> file for reading: " << pathInput << std::endl;
        return 1;
    }

    if (!sampleFile) {
        std::cerr << "Error: Cannot open <samples> file for reading: " << pathSamples << std::endl;
        return 1;
    }

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

        // Split the line using white spaces
        std::istringstream iss(sampleLine);
        std::string word;
        std::vector<std::string> wordsInLine;

        while (iss >> word) {
            wordsInLine.push_back(word);
        }

        // Check if line contains more than one word or is empty
        if (wordsInLine.size() != 1 || sampleLine.empty()) {
            std::cerr << "Error: Invalid samples file format. Expected one sample per line but detected more!" << std::endl;
            printUsage(argv[0]);
            sampleFile.close(); // Close the sample file when done
            return 1;
         }

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
    char buffer[4096];
    while (gzgets(longFile, buffer, sizeof(buffer))) {
        std::string line(buffer);
        std::istringstream iss(line);
        iss >> sample >> chromosome >> gene >> configuration >> dosage >> variantInfo;
	if (mode == "recessive" && dosage == 1) {
            dosage = 0;
        }
        geneToChromosome[gene] = chromosome;
        geneSampleDosage[gene][sample] = dosage;
	contigs.insert(chromosome);
    }

    gzclose(longFile);

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

