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
    std::cerr << "Usage: " << path << " --geno <Genotype File> --map <Mapping File> [--info-map <Info Mapping File>] [--dosage-map <Dosage Mapping File>]\n"
              << std::endl;

    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  The program calls co-occurring variants using a variant-to-gene mapping file alongside a file containing alternate genotypes." << std::endl;
    std::cerr << "  It calculates the call value and dosage based on the haplotype information of each variant." << std::endl;

    std::cerr << "\nInput:" << std::endl;
    std::cerr << "  --geno/-g <Genotype File>   : Required. A gzipped file containing alternate genotype data." << std::endl;
    std::cerr << "                                The file should have sample id, variant id, and genotype separated by whitespace." << std::endl;
    std::cerr << "                                Each line might look like the following: 'Sample1 chr21:12314:C:T 1|0'." << std::endl;
    std::cerr << "  --gene-map/-m <Map File>    : Required. File mapping variants to genes. It should contain at least" << std::endl;
    std::cerr << "                                two columns with a header in the format: variant, gene, info (optional)." << std::endl;
    std::cerr << "  --info-map/-i               : Optional. File mapping variants to their corresponding information." << std::endl;
    std::cerr << "                                The file should contain variant, gene, and info columns, mapping each variant " << std::endl;
    std::cerr << "                                to infomation provided. The file must contain variant, gene, and info columns." << std::endl;
    std::cerr << "  --score-map/-p              : Optional. File mapping variants to their score information." << std::endl;
    std::cerr << "                                The file should contain variant, gene, and score columns." << std::endl;
    std::cerr << "  --haplotype-collapse/-hc     : Optional. Specifies the rule for combining variant dosages within a haplotype." << std::endl;
    std::cerr << "                                Options are 'product', 'additive', 'min' or 'max'. Please, note that" << std::endl;
    std::cerr << "                                when 'product' is chosen the results is P(affected) by determining  " << std::endl;
    std::cerr << "                                P(affected) = 1-((1-P1_score)*(1-P2_score) .. * (1-Pn_score)) for n variants" << std::endl;
    std::cerr << "                                across that haplotype. Default value 'product'. " << std::endl;
    std::cerr << "  --gene-collapse/-gc         : Optional. Specifies the rule for combining resulting variant scores across" << std::endl;
    std::cerr << "                                haplotypes. The options are 'additive', 'product','min' or 'max'. The" << std::endl;
    std::cerr << "                                default is 'product'." << std::endl;
    std::cerr << "  --show-haplotype-score/-shs : Optional. Prints the haplotype-specific scores (when '-c' is specified)." << std::endl;
    std::cerr << "  --show-variants/-s          : Optional. Print variants involved in encoding as an extra column." << std::endl;
    std::cerr << "  --debug                     : Optional. Print out more information during run" << std::endl;

    std::cerr << "\nOutput Format:" << std::endl;
    std::cerr << "  The output will be printed to the console with the following columns:" << std::endl;
    std::cerr << "  Sample, Gene, Call, Dosage, Score (Will be zero unless --score-map is defined), info (with --info-map)" << std::endl;

    std::cerr << "\nNotes:" << std::endl;
    std::cerr << "  Ensure that the --geno file is appropriately formatted (optionally gzipped)." << std::endl;
    std::cerr << "  The program will attempt to match variants between the --gene-map and --geno file, discarding those that do not match." << std::endl;
}

int main(int argc, char *argv[])
{

    std::string pathGeno;
    std::string pathMap;
    std::string pathInfoMap;
    std::string pathScoreMap;
    std::string haplotypeCollapseRule = "product";
    std::string geneCollapseRule = "product";
    bool haplotypeCollapseRuleSet = false;
    bool geneCollapseRuleSet = false;
    bool showVariants = false;
    bool showHaplotypeScore = false;
    bool verbose = false;

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
        else if ((arg == "--score-map" || arg == "-p") && i + 1 < argc)
        {
            pathScoreMap = argv[++i];
        }
        else if ((arg == "--haplotype-collapse" || arg == "-hc") && i + 1 < argc)
        {
            haplotypeCollapseRule = argv[++i];
            haplotypeCollapseRuleSet = true;
            if (haplotypeCollapseRule != "product" && haplotypeCollapseRule != "additive" && haplotypeCollapseRule != "max" && haplotypeCollapseRule != "min")
            {
                std::cerr << "Error! Invalid haplotype-collapse rule: " << haplotypeCollapseRule << std::endl;
                printUsage(argv[0]);
                return 1;
            }
            
        }
        else if ((arg == "--gene-collapse" || arg == "-gc") && i + 1 < argc)
        {
            geneCollapseRule = argv[++i];
            geneCollapseRuleSet = true;
            if (geneCollapseRule != "product" && geneCollapseRule != "additive" && geneCollapseRule != "max" && geneCollapseRule != "min")
            {
                std::cerr << "Error! Invalid gene-collapse rule: " << geneCollapseRule << std::endl;
                printUsage(argv[0]);
                return 1;
            }
            
        }
        else if (arg == "--show-haplotype-score" || arg == "-shs")
        {
            showHaplotypeScore = true;
        }
        else if (arg == "--show-variants" || arg == "-sv")
        {
            showVariants = true;
        }
        else if (arg == "--debug")
        {
            verbose = true;
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

    // print some warnings in case score map is not specified but other arguments that depend on it are specifed
    if (pathScoreMap.empty() && showHaplotypeScore){
        std::cerr << "Warning: ignored --show-haplotype-score argument. This can only be used when --score-map is also specified!" << std::endl;
    }
    if (pathScoreMap.empty() && geneCollapseRuleSet){
        std::cerr << "Warning: ignored --gene-collapse argument. This can only be used when --score-map is also specified!" << std::endl;
    }
    if (pathScoreMap.empty() && haplotypeCollapseRuleSet){
        std::cerr << "Warning: ignored --haplotype-collapse argument. This can only be used when --score-map is also specified!" << std::endl;
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

    // read the score-map file if provided
    std::map<std::pair<std::string, std::string>, float> variantGeneScore;
    if (!pathScoreMap.empty())
    {
        gzFile scoreMapFile = gzopen(pathScoreMap.c_str(), "rb");
        if (!scoreMapFile)
        {
            std::cerr << "Warning: Cannot open --score-map file for reading: " << pathScoreMap << ". Continuing without score data." << std::endl;
        }
        else
        {
            char pathoBuf[1024];
            while (gzgets(scoreMapFile, pathoBuf, sizeof(pathoBuf)))
            {
                std::stringstream ss(pathoBuf);
                std::string variant, gene;
                float newScore;
                ss >> variant >> gene >> newScore;
                variantGeneScore[std::make_pair(variant, gene)] = newScore;
            }
            gzclose(scoreMapFile);
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
                if (genotype == "1|0")
                {
                    sampleGeneHaplotypeVariant[sample][gene][1].push_back(variant); // Haplotype 1
                }
                else if (genotype == "0|1")
                {
                    sampleGeneHaplotypeVariant[sample][gene][2].push_back(variant); // Haplotype 2
                }
                else if (genotype == "1|1")
                {
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

    // only print stuff if user wants it
    if (verbose)
    {
        std::cerr << "Variants in mapping file (kept): " << uniqueVariantsKept.size() << std::endl;
        std::cerr << "Variants not in mapping file (Discarded): " << uniqueVariantsDiscarded.size() << std::endl;
        std::cerr << "Variants mapping to more than one gene: " << multiGeneVariants.size() << std::endl;
    }

    // Print results
    for (const auto &samplePair : sampleGeneHaplotypeVariant)
    {
        const std::string &sample = samplePair.first;
        for (const auto &genePair : samplePair.second)
        {
            const std::string &gene = genePair.first;
            const auto &haplotypeVariantMap = genePair.second;
            std::string chromosome;

            std::set<std::string> haplotype1Variants = haplotypeVariantMap.count(1) ? std::set<std::string>(haplotypeVariantMap.at(1).begin(), haplotypeVariantMap.at(1).end()) : std::set<std::string>();
            std::set<std::string> haplotype2Variants = haplotypeVariantMap.count(2) ? std::set<std::string>(haplotypeVariantMap.at(2).begin(), haplotypeVariantMap.at(2).end()) : std::set<std::string>();

            // extract chromosome from string
            std::set<std::string> mergedVariants;
            std::set_union(haplotype1Variants.begin(), haplotype1Variants.end(), haplotype2Variants.begin(), haplotype2Variants.end(), std::inserter(mergedVariants, mergedVariants.begin()));
            const std::string& variant = *mergedVariants.begin(); // taking the first variant
            std::stringstream ss(variant);
            std::getline(ss, chromosome, ':'); // extracting chromosome part

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

            // print the basic information
            std::cout << sample << "\t" << chromosome << "\t" << gene << "\t" << callValue << "\t" << dosage;

            int totalMapped = 0;
            int totalVariants = 0;
            float geneScore = 0.0f;
            float haplotype1Score = 0.0f;
            float haplotype2Score = 0.0f;

            // calculate probability of knockout for each haplotype
            if (!pathScoreMap.empty())
            {

                // Only for product we the default
                if (haplotypeCollapseRule == "product")
                {
                    haplotype1Score = 1.0f;
                    haplotype2Score = 1.0f;
                }
                else if (haplotypeCollapseRule == "min")
                {
                    // variants in cis on haplotype 1
                    if (!haplotype1Variants.empty() && haplotype2Variants.empty())
                    {
                        haplotype1Score = 1.0f;
                        haplotype2Score = 0.0f;

                        // variants in cis on haplotype 2
                    }
                    else if (!haplotype2Variants.empty() && haplotype1Variants.empty())
                    {
                        haplotype1Score = 0.0f;
                        haplotype2Score = 1.0f;

                        // variants in chets
                    }
                    else
                    {
                        haplotype1Score = 1.0f;
                        haplotype2Score = 1.0f;
                    }
                }
                else
                {
                    haplotype1Score = 0.0f;
                    haplotype2Score = 0.0f;
                }

                // keep track of how many variants are mapped here
                int counthaplotype1ScoreMapped = 0;
                int counthaplotype2ScoreMapped = 0;
                int totalHaplotype1Variants = haplotype1Variants.size();
                int totalHaplotype2Variants = haplotype2Variants.size();

                // Run for each haplotype separately (haplotype 1)
                for (const auto &variant : haplotype1Variants)
                {
                    auto it = variantGeneScore.find(std::make_pair(variant, gene));
                    float mappedScore = (it != variantGeneScore.end()) ? it->second : 1.0f; // default to 1.0 if not found
                    if (mappedScore != 1.0f)
                        counthaplotype1ScoreMapped++;
                    if (haplotypeCollapseRule == "product")
                    {
                        haplotype1Score *= (1 - mappedScore);
                    }
                    else if (haplotypeCollapseRule == "max")
                    {
                        haplotype1Score = std::max(haplotype1Score, mappedScore);
                    }
                    else if (haplotypeCollapseRule == "min")
                    {
                        haplotype1Score = std::min(haplotype1Score, mappedScore);
                    }
                    else if (haplotypeCollapseRule == "additive")
                    {
                        haplotype1Score += mappedScore;
                    }
                }

                // Run for each haplotype separately (haplotype 2)
                for (const auto &variant : haplotype2Variants)
                {
                    auto it = variantGeneScore.find(std::make_pair(variant, gene));
                    float mappedScore = (it != variantGeneScore.end()) ? it->second : 1.0f; // default to 1.0 if not found
                    if (mappedScore != 1.0f)
                        counthaplotype2ScoreMapped++;
                    if (haplotypeCollapseRule == "product")
                    {
                        haplotype2Score *= (1 - mappedScore);
                    }
                    else if (haplotypeCollapseRule == "max")
                    {
                        haplotype2Score = std::max(haplotype2Score, mappedScore);
                    }
                    else if (haplotypeCollapseRule == "min")
                    {
                        haplotype2Score = std::min(haplotype2Score, mappedScore);
                    }
                    else if (haplotypeCollapseRule == "additive")
                    {
                        haplotype2Score += mappedScore;
                    }
                }

                // 1 minus the haplotype score is the probability of knockout for each haplotype
                if (haplotypeCollapseRule == "product")
                {
                    haplotype2Score = 1 - haplotype2Score;
                    haplotype1Score = 1 - haplotype1Score;
                }

                // count total variants mapped by gene-sample
                totalMapped = counthaplotype1ScoreMapped + counthaplotype2ScoreMapped;
                totalVariants = totalHaplotype1Variants + totalHaplotype2Variants;

                // compute final score that depending on gene collapse rules
                if (geneCollapseRule == "product")
                {
                    geneScore = haplotype1Score * haplotype2Score;
                }
                else if (geneCollapseRule == "max")
                {
                    geneScore = std::max(haplotype1Score, haplotype2Score);
                }
                else if (geneCollapseRule == "min")
                {
                    geneScore = std::min(haplotype1Score, haplotype2Score);
                }
                else if (geneCollapseRule == "additive")
                {
                    geneScore = haplotype1Score + haplotype2Score;
                }

                // Print results here
                std::cout << "\tg=" << geneCollapseRule << "\t" << geneScore;
                if (showHaplotypeScore)
                {
                    std::cout << "\th=" << haplotypeCollapseRule << "\t" << haplotype1Score << "\t" << haplotype2Score;
                }
            }

            // allow user to show variants even if info file is not specified
            if (showVariants && pathInfoMap.empty())
            {

                std::cout << "\t";

                // For Haplotype 1 variants
                if (!haplotype1Variants.empty())
                {
                    for (auto it = haplotype1Variants.begin(); it != haplotype1Variants.end(); ++it)
                    {
                        std::cout << *it;
                        if (std::next(it) != haplotype1Variants.end())
                        {
                            std::cout << ";";
                        }
                    }
                }

                if (!haplotype1Variants.empty() && !haplotype2Variants.empty())
                {
                    std::cout << "|";
                }

                // For Haplotype 2 variants
                if (!haplotype2Variants.empty())
                {
                    for (auto it = haplotype2Variants.begin(); it != haplotype2Variants.end(); ++it)
                    {
                        std::cout << *it;
                        if (std::next(it) != haplotype2Variants.end())
                        {
                            std::cout << ";";
                        }
                    }
                }
            }

            if (!pathInfoMap.empty())
            {
                {

                    std::cout << "\t";

                    // For Haplotype 1 variants
                    if (!haplotype1Variants.empty())
                    {
                        for (auto it = haplotype1Variants.begin(); it != haplotype1Variants.end(); ++it)
                        {
                            std::pair<std::string, std::string> key = std::make_pair(*it, gene);
                            std::string info = (infoMap.find(key) != infoMap.end()) ? infoMap[key] : "NA";
                            std::cout << *it << ":" << info;
                            if (std::next(it) != haplotype1Variants.end())
                            {
                                std::cout << ";";
                            }
                        }
                    }

                    if (!haplotype1Variants.empty() && !haplotype2Variants.empty())
                    {
                        std::cout << "|";
                    }

                    // For Haplotype 2 variants
                    if (!haplotype2Variants.empty())
                    {
                        for (auto it = haplotype2Variants.begin(); it != haplotype2Variants.end(); ++it)
                        {
                            std::pair<std::string, std::string> key = std::make_pair(*it, gene);
                            std::string info = (infoMap.find(key) != infoMap.end()) ? infoMap[key] : "NA";
                            std::cout << *it << ":" << info;
                            if (std::next(it) != haplotype2Variants.end())
                            {
                                std::cout << ";";
                            }
                        }
                    }
                }
            }

            std::cout << std::endl;
        }
    }
    return 0;
}

