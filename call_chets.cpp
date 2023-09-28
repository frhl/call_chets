#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <zlib.h>
#include <set>

std::string getVersion()
{
    std::ifstream versionFile(".version");
    if (!versionFile.is_open())
    {
        std::cerr << "Warning: Unable to open .version file." << std::endl;
    }
    std::string version;
    std::getline(versionFile, version);
    versionFile.close();

    return version;
}

void printUsage(const char *path)
{
    std::string version = getVersion();
    std::cerr << "\nProgram: chet tools v" << version << "\n"
              << std::endl;
    std::cerr << "Usage: " << path << " --geno <Genotype File> --map <Mapping File>\n"
              << std::endl;
    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  Call co-occuring variants using a variant-to-gene mapping file" << std::endl;
    std::cerr << "  alongside a file containing alternate genotypes." << std::endl;
    std::cerr << "\nInput:" << std::endl;
    std::cerr << "  --geno/-g : .gz file containing alternate genotype data in the" << std::endl;
    std::cerr << "              form of sample id, variant id and genotype that is" << std::endl;
    std::cerr << "              seperated by whitespace. For example a line may look." << std::endl;
    std::cerr << "              like the following: 'Sample1 chr21:12314:C:T 1|0'." << std::endl;
    std::cerr << "  --map/-m  : File mapping variants to genes. This is expected" << std::endl;
    std::cerr << "              to contain at least two columns with a header" << std::endl;
    std::cerr << "              in the format: variant, gene, info (optional)." << std::endl;
    std::cerr << "Output Format:" << std::endl;
    std::cerr << "  Sample Chromosome Gene Call Dosage Variant-Genotype..." << std::endl;
    std::cerr << "\nNotes:" << std::endl;
    std::cerr << "  See README for how to generate --geno file from a VCF/BCF." << std::endl;
}

int main(int argc, char *argv[])
{

    std::string pathGeno;
    std::string pathMap;
    std::string pathInfoMap;
    std::string pathDosageMap;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if ((arg == "--geno" || arg == "-g") && i + 1 < argc)
        {
            pathGeno = argv[++i];
        }
        else if ((arg == "--gene-map" || arg == "-m") && i + 1 < argc)
        {
            pathMap = argv[++i];
        }
        else if ((arg == "--info-map" || arg == "-i") && i + 1 < argc)
        {
            pathInfoMap = argv[++i];
        }

        else if ((arg == "--dosage-map" || arg == "-p") && i + 1 < argc)
        {
            pathDosageMap = argv[++i];
        }
        else
        {
            std::cerr << "Error! Unknown or incomplete argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    // check if needed flags
    if (pathGeno.empty() || pathMap.empty())
    {
        std::cerr << "Error! Both --geno and --map must be provided." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // expect .gz file no matter what
    gzFile genotypeFile = gzopen(pathGeno.c_str(), "rb");
    gzFile mappingFile = gzopen(pathMap.c_str(), "rb");

    if (!genotypeFile)
    {
        std::cerr << "Error: Cannot open --geno file for reading: " << pathGeno << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    if (!mappingFile)
    {
        std::cerr << "Error: Cannot open --mapping file for reading: " << pathMap << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Mapping variant to gene
    std::map<std::string, std::vector<std::pair<std::string, std::string>>> variantToGene;
    std::map<std::string, std::map<std::string, double>> variantGeneToPathogenicity;
    std::string variant, gene, thirdColumn;
    char buf[1024];

    // Read the mapping file
    while (gzgets(mappingFile, buf, sizeof(buf)))
    {
        std::stringstream ss(buf);
        ss >> variant >> gene;
        std::pair<std::string, std::string> geneInfo = ss >> thirdColumn ? std::make_pair(gene, thirdColumn) : std::make_pair(gene, "");
        variantToGene[variant].push_back(geneInfo);
    }
    gzclose(mappingFile);

    // read the info-map file if provided
    std::map<std::pair<std::string, std::string>, std::string> infoMap;
    if (!pathInfoMap.empty())
    {
        gzFile infoMapFile = gzopen(pathInfoMap.c_str(), "rb");
        if (!infoMapFile)
        {
            std::cerr << "Warning: Cannot open --info-map file for reading: " << pathInfoMap << ". Continuing without info data." << std::endl;
        }
        else
        {
            char infoBuf[1024];
            while (gzgets(infoMapFile, infoBuf, sizeof(infoBuf)))
            {
                std::stringstream ss(infoBuf);
                std::string variant, gene, info;
                ss >> variant >> gene >> info;
                infoMap[std::make_pair(variant, gene)] = info;
            }
            gzclose(infoMapFile);
        }
    }

    // read the dosage-map file if provided
    std::map<std::pair<std::string, std::string>, float> variantGeneDosage;
    if (!pathDosageMap.empty())
    {
        gzFile dosageMapFile = gzopen(pathDosageMap.c_str(), "rb");
        if (!dosageMapFile)
        {
            std::cerr << "Warning: Cannot open --dosage-map file for reading: " << pathDosageMap << ". Continuing without dosage data." << std::endl;
        }
        else
        {
            char pathoBuf[1024];
            while (gzgets(dosageMapFile, pathoBuf, sizeof(pathoBuf)))
            {
                std::stringstream ss(pathoBuf);
                std::string variant, gene;
                float newDosage;
                ss >> variant >> gene >> newDosage;
                variantGeneDosage[std::make_pair(variant, gene)] = newDosage;
            }
            gzclose(dosageMapFile);
        }
    }

    // Process genotype file
    std::map<std::string, std::map<std::string, std::vector<std::string>>> sampleGeneVariants;
    std::map<std::string, std::map<std::string, float>> sampleVariantDosage;
    std::map<std::string, std::map<std::string, std::string>> sampleVariantInfo;
    std::map<std::string, std::map<std::string, std::map<int, std::vector<std::string>>>> sampleGeneHaplotypeVariant;

    std::set<std::string> uniqueVariantsKept;
    std::set<std::string> uniqueVariantsDiscarded;
    std::set<std::string> multiGeneVariants;
    
    bool isFirstLine = true;
    while (gzgets(genotypeFile, buf, sizeof(buf)))
    {
        std::stringstream ss(buf);
        std::string sample, variant, genotype;
        ss >> sample >> variant >> genotype;

        // if header skip line
        if ((genotype != "1|0" && genotype != "0|1" && genotype != "1|1") && (isFirstLine))
        {
            continue;
        }
        isFirstLine = false;

        // Warn user that skipping is happening
        if (genotype != "1|0" && genotype != "0|1" && genotype != "1|1")
        {
            std::cerr << "Warning: Skipping unexpected genotype value (" << genotype << ") in sample " << sample << " for variant " << variant << "." << std::endl;
            continue;
        }

        if (variantToGene.find(variant) != variantToGene.end())
        {
            if (variantToGene[variant].size() > 1)
            {
                multiGeneVariants.insert(variant);
            }

            for (const auto &geneInfoPair : variantToGene[variant])
            {
                std::string gene = geneInfoPair.first;
                uniqueVariantsKept.insert(variant);

                // save haplotype
                if (genotype == "1|0") {
                    sampleGeneHaplotypeVariant[sample][gene][1].push_back(variant); // Haplotype 1
                } else if (genotype == "0|1") {
                    sampleGeneHaplotypeVariant[sample][gene][2].push_back(variant); // Haplotype 2
                } else if (genotype == "1|1") {
                    sampleGeneHaplotypeVariant[sample][gene][1].push_back(variant); // Haplotype 1
                    sampleGeneHaplotypeVariant[sample][gene][2].push_back(variant); // Haplotype 2
                }
            }
        }
        else
        {
            uniqueVariantsDiscarded.insert(variant);
        }
    }


    // Check if there are no matching variants
    if (uniqueVariantsKept.empty())
    {
        std::cerr << "Error: No matching variants found between the --map and --geno file!" << std::endl;
        return 1;
    }

    std::cerr << "Variants in mapping file (kept): " << uniqueVariantsKept.size() << std::endl;
    std::cerr << "Variants not in mapping file (Discarded): " << uniqueVariantsDiscarded.size() << std::endl;
    std::cerr << "Variants mapping to more than one gene: " << multiGeneVariants.size() << std::endl;
    std::cerr << "* Generating sample-gene-variant file.." << std::endl;

    std::cerr << "* Generating output.." << std::endl;

    for (const auto &samplePair : sampleGeneHaplotypeVariant)
    {
        const std::string &sample = samplePair.first;
        for (const auto &genePair : samplePair.second)
        {
            const std::string &gene = genePair.first;
            const auto &haplotypeVariantMap = genePair.second;

            std::set<std::string> haplotype1Variants = haplotypeVariantMap.count(1) ? std::set<std::string>(haplotypeVariantMap.at(1).begin(), haplotypeVariantMap.at(1).end()) : std::set<std::string>();
            std::set<std::string> haplotype2Variants = haplotypeVariantMap.count(2) ? std::set<std::string>(haplotypeVariantMap.at(2).begin(), haplotypeVariantMap.at(2).end()) : std::set<std::string>();

            std::string callValue;
            int dosage = 0;

            if (!haplotype1Variants.empty() && haplotype2Variants.empty())
            {
                // Variants only on haplotype 1
                callValue = (haplotype1Variants.size() == 1) ? "het" : "cis";
                dosage = 1;
            }
            else if (!haplotype2Variants.empty() && haplotype1Variants.empty())
            {
                // Variants only on haplotype 2
                callValue = (haplotype2Variants.size() == 1) ? "het" : "cis";
                dosage = 1;
            }
            else if (!haplotype1Variants.empty() && !haplotype2Variants.empty())
            {
                // Variants on both haplotypes
                std::set<std::string> intersection;
                std::set_intersection(haplotype1Variants.begin(), haplotype1Variants.end(), haplotype2Variants.begin(), haplotype2Variants.end(), std::inserter(intersection, intersection.begin()));

                if (!intersection.empty())
                {
                    // There is at least one variant in common between the two haplotypes
                    callValue = "hom";
                    dosage = 2;
                }
                else
                {
                    // Different variants on each haplotype
                    callValue = "chet";
                    dosage = 2;
                }
            }
            else
            {
                // Handle the case of no variants, or another unexpected case
                callValue = "na";
                dosage = 0;
            }

            // Print output here
            std::cout << sample << "\t" << gene << "\t" << callValue << "\t" << dosage;
            // Add additional output fields as necessary
            std::cout << std::endl;
        }
    }

    return 0;
}

