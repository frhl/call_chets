#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <zlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>


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
	std::cerr << "  --min-ac     : Filters to genes with sum of DS >= argument" << std::endl;
	std::cerr << "  --max-ac     : Filters to genes with sum of DS < argument" << std::endl;
	std::cerr << "\nExample:" << std::endl;
	std::cerr << "  ./encode_vcf called_chets.txt.gz samples.txt additive | bgzip > out.vcf.gz\n\n";
}


std::vector<std::string> sortChromosomes(const std::set<std::string>& contigs) {
    std::vector<std::string> chromosomes(contigs.begin(), contigs.end());
    std::sort(chromosomes.begin(), chromosomes.end(), [](const std::string& a, const std::string& b) {
        // Extract the part after 'chr' prefix
        std::string a_num = a.substr(0, 3) == "chr" ? a.substr(3) : a;
        std::string b_num = b.substr(0, 3) == "chr" ? b.substr(3) : b;

        // If both are numeric chromosomes
        if (isdigit(a_num[0]) && isdigit(b_num[0])) {
            return std::stoi(a_num) < std::stoi(b_num);
        }
        // If a is numeric but b is not (e.g., a = chr2, b = chrX)
        if (isdigit(a_num[0]) && !isdigit(b_num[0])) {
            return true;
        }
        // If b is numeric but a is not
        if (!isdigit(a_num[0]) && isdigit(b_num[0])) {
            return false;
        }
        // If neither are numeric (e.g., chrX, chrY, etc.), then just use lexicographical order
        return a < b;
    });
    return chromosomes;
}


int main(int argc, char* argv[]) {
    
    std::string pathInput;
    std::string pathSamples;
    std::string mode = "additive";
    int minAC = 0;
    int maxAC = INT_MAX;

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
        } else if (arg == "--min-ac") {
	    minAC = std::stoi(argv[++i]);
            if (minAC < 0) {
                 std::cerr << "Error! Minimum Allele Count (AC) cannot be negative." << std::endl;
                 return 1;
            }
        } else if (arg == "--max-ac") {
	    maxAC = std::stoi(argv[++i]);
	    if (maxAC < 0) {
                 std::cerr << "Error! Maximum Allele Count (AC) cannot be negative." << std::endl;
                 return 1;
             }
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
    std::map<std::string, std::string> geneToChromosome;
    std::map<std::string, int> geneAC;
    std::map<std::string, int> geneBI;
    std::map<std::string, int> geneChet;
    std::map<std::string, int> geneHom;
    std::map<std::string, int> geneHet;
    std::map<std::string, int> geneCis;

    // Process long format file
    std::string line, sample, chromosome, gene, configuration, variantInfo;
    std::set<std::string> contigs;
    int dosage;
    char buffer[4096];
    while (gzgets(longFile, buffer, sizeof(buffer))) {
        std::string line(buffer);
        std::istringstream iss(line);
        iss >> sample >> chromosome >> gene >> configuration >> dosage >> variantInfo;
	
	// count instances of mono-allelic
	if (mode == "recessive" && dosage == 1) {
            dosage = 0;
        }

	// count individual sites
	if (configuration == "chet" && dosage == 2) {
	    geneChet[gene] += 1;
	} else if (configuration == "hom" && dosage == 2) {
	    geneHom[gene] += 1;
	} else if (configuration == "cis" && dosage == 1) {
	    geneCis[gene] += 1;
	} else if (configuration == "het" && dosage == 1) {
	    geneHet[gene] += 1;
	}
	
	// count instances of bi-allelic
	if (dosage == 2) {
            geneBI[gene] += 1;
	} 
        geneToChromosome[gene] = chromosome;
        geneSampleDosage[gene][sample] = dosage;
        geneAC[gene] += dosage;
	contigs.insert(chromosome);
    }

    gzclose(longFile);

    // sort chromosomes
    std::vector<std::string> sortedContigs = sortChromosomes(contigs);

    // Print the output header
    std::cout << "##fileformat=VCFv4.2\n";
    std::cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    for (const auto& chr : sortedContigs) {
        std::cout << "##contig=<ID=" << chr << ">\n";
    }
    std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
    std::cout << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n";
    std::cout << "##INFO=<ID=BI,Number=1,Type=Integer,Description=\"Bi-allelic Count\">\n";
    std::cout << "##INFO=<ID=CHET,Number=1,Type=Integer,Description=\"Compound heterozygous pseudo Count\">\n";
    std::cout << "##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Homozygous Count\">\n";
    std::cout << "##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Heterozygous Count\">\n";
    std::cout << "##INFO=<ID=CIS,Number=1,Type=Integer,Description=\"Cis pseudo Count \">\n";
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sampleName : samples) {
        std::cout << "\t" << sampleName;
    }

    std::cout << std::endl;
    // Print the output data
    for (const auto& chr : sortedContigs) {
        int rowIndex = 0;
	for (const auto &genePair : geneSampleDosage) {
            if (geneToChromosome[genePair.first] == chr) {  // Only process genes on the current chromosome
		rowIndex++;
		int currentAC = geneAC[genePair.first];
		int currentBI = geneBI[genePair.first];
		int currentChet = geneChet[genePair.first];
		int currentHom = geneHom[genePair.first];
		int currentHet = geneHet[genePair.first];
		int currentCis = geneCis[genePair.first];
		// only output variants that fit our criteria
		if (currentAC >= minAC && currentAC < maxAC) {

		    // get left side of VCF body
		    std::cout << geneToChromosome[genePair.first] 
			  << "\t" << rowIndex 
			  << "\t" << genePair.first 
			  << "\tA\tB\t.\t.\t"
			  << "AC=" << currentAC 
			  << ";BI=" << currentBI 
			  << ";CHET=" << currentChet 
			  << ";HOM=" << currentHom 
			  << ";CIS=" << currentCis 
			  << ";HET=" << currentHet 
			  << "\tDS";  
		    
		    // get right side of body VCF
		    for (const auto& sample : samples) {
			std::cout << "\t";

			if (genePair.second.find(sample) != genePair.second.end()) {
			    std::cout << genePair.second.at(sample);
			} else {
			    std::cout << "0";  // Default to 0 if sample-gene combo doesn't exist
			}
		    }
		    std::cout << std::endl;
		}
	    }
       }
    }

    return 0;
}
