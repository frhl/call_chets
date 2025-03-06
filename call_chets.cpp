#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <zlib.h>
#include <set>
#include <regex>

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

#ifndef GIT_DATE
#define GIT_DATE "unknown"
#endif

#ifndef VERSION
#define VERSION "0.3.0"
#endif

std::string getVersion()
{
    // Use preprocessor macros for version and git info
    std::string version = VERSION;
    std::string gitCommit = GIT_COMMIT;
    std::string gitDate = GIT_DATE;
    
    // Try to read from .version file if it exists
    std::ifstream versionFile(".version");
    if (versionFile.is_open())
    {
        std::string fileVersion;
        if (std::getline(versionFile, fileVersion) && !fileVersion.empty())
        {
            version = fileVersion;
        }
        versionFile.close();
    }
    
    // Construct full version string
    std::string versionInfo = version + " / commit = " + gitCommit + " / release = " + gitDate;
    return versionInfo;
}

#include <ctime>

void printUsage(const char *path)
{
    // Get version info
    std::string version = getVersion();
    
    // Get current date and time
    std::time_t now = std::time(nullptr);
    char timestr[100];
    std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S", std::localtime(&now));
    
    // Print header with version and run date
    std::cerr << "\n[CALL_CHETS] compound heterozygosity and cis variant caller"
              << "\n  * Version       : " << version
              << "\n  * Run date      : " << timestr
              << "\n" << std::endl;
              
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
    std::cerr << "  --haplotype-collapse/-hc    : Optional. Specifies the rule for combining variant dosages within a single" << std::endl;
    std::cerr << "                              : haplotype. This is only relevant when +1 variants are present on the haplotype." << std::endl;
    std::cerr << "                                Options are 'product', 'additive', 'min' or 'max'. Please, note that" << std::endl;
    std::cerr << "                                when 'product' is chosen the results is P(affected) by determining  " << std::endl;
    std::cerr << "                                P(affected) = 1-((1-P1_score)*(1-P2_score) .. * (1-Pn_score)) for n variants" << std::endl;
    std::cerr << "                                across that haplotype. Default value 'product'. " << std::endl;
    std::cerr << "  --gene-collapse/-gc         : Optional. Specifies the rule for combining resulting variant scores across" << std::endl;
    std::cerr << "                                haplotypes. The options are 'additive', 'product','min' or 'max'. The" << std::endl;
    std::cerr << "                                default is 'product'. This argument is only relevant when --score-map is defined." << std::endl;
    std::cerr << "  --show-haplotype-scores/-shs: Optional. Prints the haplotype-specific scores (when '-c' is specified)." << std::endl;
    std::cerr << "  --show-variants/-s          : Optional. Print variants involved in encoding as an extra column." << std::endl;
    std::cerr << "  --debug                     : Optional. Print out more information during run" << std::endl;

    std::cerr << "\nOutput Format:" << std::endl;
    std::cerr << "  The output will be printed to the console with the following columns:" << std::endl;
    std::cerr << "  Sample, Gene, Call, Dosage, Score (Will be zero unless --score-map is defined), info (with --info-map)" << std::endl;

    std::cerr << "\nNotes:" << std::endl;
    std::cerr << "  Ensure that the --geno file is appropriately formatted (optionally gzipped)." << std::endl;
    std::cerr << "  The program will attempt to match variants between the --gene-map and --geno file, discarding those that do not match." << std::endl;
}

// Function to validate variant format (chr:pos:ref:alt)
bool isValidVariantFormat(const std::string& variant) {
    std::regex pattern("^(chr)?[0-9XYM]{1,2}:[0-9]+:[ACGT]+:[ACGT]+$");
    return std::regex_match(variant, pattern);
}

// Function to validate score value (between 0 and 1)
bool isValidScore(float score) {
    return score >= 0.0f && score <= 1.0f;
}

// Function to check if a file is empty
bool isFileEmpty(gzFile file) {
    char testBuf[2];
    int bytesRead = gzread(file, testBuf, 1);
    gzrewind(file);  // Reset file position
    return bytesRead <= 0;
}

// Function to count columns in a string
size_t countColumns(const std::string& line) {
    std::stringstream ss(line);
    std::string temp;
    size_t count = 0;
    while (ss >> temp) {
        count++;
    }
    return count;
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
        else if (arg == "--show-haplotype-scores" || arg == "-shs")
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
        gzclose(genotypeFile);
        printUsage(argv[0]);
        return 1;
    }

    // Check if files are empty
    if (isFileEmpty(genotypeFile)) {
        std::cerr << "Error: Genotype file is empty: " << pathGeno << std::endl;
        gzclose(genotypeFile);
        gzclose(mappingFile);
        return 1;
    }

    if (isFileEmpty(mappingFile)) {
        std::cerr << "Error: Mapping file is empty: " << pathMap << std::endl;
        gzclose(genotypeFile);
        gzclose(mappingFile);
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
    std::map<std::string, std::vector<std::string>> variantToGene;
    std::map<std::string, std::map<std::string, double>> variantGeneToPathogenicity;
    std::string variant, gene;
    char buf[4096]; // Increased buffer size for potentially larger lines

    bool isFirstLineMappingFile = true;
    int mappingLineCount = 0;
    int validMappingLines = 0;
    int invalidFormatVariants = 0;

    while (gzgets(mappingFile, buf, sizeof(buf)))
    {
        mappingLineCount++;
        std::string line(buf);
        
        // Remove trailing newline characters
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Check column count
        size_t columnCount = countColumns(line);
        if (columnCount < 2 && !isFirstLineMappingFile) {
            std::cerr << "Error: Line " << mappingLineCount << " in mapping file has fewer than required 2 columns: '" << line << "'" << std::endl;
            gzclose(mappingFile);
            gzclose(genotypeFile);
            return 1;
        }

        std::stringstream ss(line);
        ss >> variant >> gene;

        // Check whether data extraction was successful.
        if(ss.fail() || variant.empty() || gene.empty())
        {
            if(!isFirstLineMappingFile)
            {
                std::cerr << "Error: Failed to extract two columns from line " << mappingLineCount 
                          << " (variant gene): '" << line << "'. Please, fix this line in --gene-map and retry." << std::endl;
                gzclose(mappingFile);
                gzclose(genotypeFile);
                return 1;
            }
        }
        else
        {
            // Validate variant format (only for non-header lines)
            if (!isFirstLineMappingFile && !isValidVariantFormat(variant)) {
                invalidFormatVariants++;
                if (verbose) {
                    std::cerr << "Warning: Line " << mappingLineCount << " - Invalid variant format: " << variant 
                              << ". Expected format: chr:pos:ref:alt" << std::endl;
                }
            }
            
            variantToGene[variant].push_back(gene);
            validMappingLines++;
        }
        isFirstLineMappingFile = false;
    }
    gzclose(mappingFile);

    // Check if any valid mapping lines were found after header
    if (validMappingLines == 0) {
        std::cerr << "Error: No valid mapping data found in: " << pathMap << std::endl;
        gzclose(genotypeFile);
        return 1;
    }

    if (invalidFormatVariants > 0) {
        std::cerr << "Warning: " << invalidFormatVariants << " variants in mapping file have invalid format." << std::endl;
    }

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
            if (isFileEmpty(infoMapFile)) {
                std::cerr << "Warning: Info mapping file is empty: " << pathInfoMap << ". Continuing without info data." << std::endl;
                gzclose(infoMapFile);
            }
            else {
                char infoBuf[4096];
                bool isFirstLine = true;
                int infoLineCount = 0;
                int validInfoLines = 0;

                while (gzgets(infoMapFile, infoBuf, sizeof(infoBuf)))
                {
                    infoLineCount++;
                    std::string line(infoBuf);
                    
                    // Remove trailing newline characters
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
                    
                    // Skip empty lines
                    if (line.empty()) {
                        continue;
                    }

                    // Check column count
                    size_t columnCount = countColumns(line);
                    if (columnCount < 3 && !isFirstLine) {
                        std::cerr << "Error: Line " << infoLineCount << " in info-map file has fewer than required 3 columns: '" 
                                << line << "'" << std::endl;
                        gzclose(infoMapFile);
                        gzclose(genotypeFile);
                        return 1;
                    }

                    std::stringstream ss(line);
                    std::string variant, gene, info;
                    ss >> variant >> gene >> info;

                    // Check whether data extraction was successful.
                    if(ss.fail() || ss.bad())
                    {
                        if(!isFirstLine)
                        {
                            std::cerr << "Error: Failed to extract variant, gene, and info from line " << infoLineCount 
                                    << ": '" << line << "'. Please, fix this line in --info-map and retry." << std::endl;
                            gzclose(infoMapFile);
                            gzclose(genotypeFile);
                            return 1;
                        }
                    }
                    else
                    {
                        // Validate variant format (only for non-header lines)
                        if (!isFirstLine && !isValidVariantFormat(variant)) {
                            if (verbose) {
                                std::cerr << "Warning: Line " << infoLineCount << " in info-map - Invalid variant format: " 
                                        << variant << ". Expected format: chr:pos:ref:alt" << std::endl;
                            }
                        }
                        
                        infoMap[std::make_pair(variant, gene)] = info;
                        validInfoLines++;
                    }
                    isFirstLine = false;
                }

                // Check if any valid info lines were found after header
                if (validInfoLines == 0) {
                    std::cerr << "Warning: No valid info mapping data found in: " << pathInfoMap 
                            << ". Continuing without info data." << std::endl;
                }

                gzclose(infoMapFile);
            }
        }
    }

    // read the score-map file if provided
    std::map<std::pair<std::string, std::string>, float> variantGeneScore;
    if (!pathScoreMap.empty())
    {
        gzFile scoreMapFile = gzopen(pathScoreMap.c_str(), "rb");
        if (!scoreMapFile)
        {
            std::cerr << "Warning: Cannot open --score-map file for reading: "
                    << pathScoreMap << ". Continuing without score data." << std::endl;
        }
        else
        {
            if (isFileEmpty(scoreMapFile)) {
                std::cerr << "Warning: Score mapping file is empty: " << pathScoreMap << ". Continuing without score data." << std::endl;
                gzclose(scoreMapFile);
            }
            else {
                char pathoBuf[4096];
                bool isFirstLine = true;
                int scoreLineCount = 0;
                int validScoreLines = 0;
                int invalidScoreValues = 0;

                while (gzgets(scoreMapFile, pathoBuf, sizeof(pathoBuf)))
                {
                    scoreLineCount++;
                    std::string line(pathoBuf);
                    
                    // Remove trailing newline characters
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
                    
                    // Skip empty lines
                    if (line.empty()) {
                        continue;
                    }

                    // Check column count
                    size_t columnCount = countColumns(line);
                    if (columnCount < 3 && !isFirstLine) {
                        std::cerr << "Error: Line " << scoreLineCount << " in score-map file has fewer than required 3 columns: '" 
                                << line << "'" << std::endl;
                        gzclose(scoreMapFile);
                        gzclose(genotypeFile);
                        return 1;
                    }

                    std::stringstream ss(line);
                    std::string variant, gene;
                    float newScore;
                    ss >> variant >> gene >> newScore;

                    // Check whether data extraction was successful.
                    if(ss.fail() || ss.bad())
                    {
                        if(!isFirstLine)
                        {
                            std::cerr << "Error: Failed to extract variant, gene, and score from line " << scoreLineCount 
                                    << ": '" << line << "'. Please, fix this line in --score-map and retry." << std::endl;
                            gzclose(scoreMapFile);
                            gzclose(genotypeFile);
                            return 1;
                        }
                    }
                    else
                    {
                        // Validate variant format (only for non-header lines)
                        if (!isFirstLine && !isValidVariantFormat(variant)) {
                            if (verbose) {
                                std::cerr << "Warning: Line " << scoreLineCount << " in score-map - Invalid variant format: " 
                                        << variant << ". Expected format: chr:pos:ref:alt" << std::endl;
                            }
                        }
                        
                        // Validate score range
                        if (!isFirstLine && !isValidScore(newScore)) {
                            invalidScoreValues++;
                            if (verbose) {
                                std::cerr << "Warning: Line " << scoreLineCount << " - Score value " << newScore 
                                        << " is outside valid range [0-1]" << std::endl;
                            }
                        }
                        
                        variantGeneScore[std::make_pair(variant, gene)] = newScore;
                        validScoreLines++;
                    }
                    isFirstLine = false;
                }

                // Check if any valid score lines were found after header
                if (validScoreLines == 0) {
                    std::cerr << "Warning: No valid score mapping data found in: " << pathScoreMap 
                            << ". Continuing without score data." << std::endl;
                }

                if (invalidScoreValues > 0) {
                    std::cerr << "Warning: " << invalidScoreValues << " scores in score-map file are outside the valid range [0-1]." << std::endl;
                }

                gzclose(scoreMapFile);
            }
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
    std::set<std::string> uniqueSamples;

    bool isFirstLine = true;
    int genoLineCount = 0;
    int validGenoLines = 0;
    int skippedGenoLines = 0;
    int invalidFormatGenoVariants = 0;
    int invalidGenotypeFormat = 0;

    while (gzgets(genotypeFile, buf, sizeof(buf)))
    {
        genoLineCount++;
        std::string line(buf);
        
        // Remove trailing newline characters
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Check column count
        size_t columnCount = countColumns(line);
        if (columnCount < 3) {
            std::cerr << "Error: Line " << genoLineCount << " in genotype file has fewer than required 3 columns: '" 
                      << line << "'" << std::endl;
            skippedGenoLines++;
            continue;
        }

        std::stringstream ss(line);
        std::string sample, variant, genotype;
        ss >> sample >> variant >> genotype;

        // if header skip line
        if ((genotype != "1|0" && genotype != "0|1" && genotype != "1|1") && (isFirstLine))
        {
            isFirstLine = false;
            continue;
        }

        uniqueSamples.insert(sample);

        // Validate variant format
        if (!isValidVariantFormat(variant)) {
            invalidFormatGenoVariants++;
            if (verbose && invalidFormatGenoVariants <= 10) { // Only show first 10 warnings to avoid flooding
                std::cerr << "Warning: Line " << genoLineCount << " - Invalid variant format: " << variant 
                          << ". Expected format: chr:pos:ref:alt" << std::endl;
            }
        }

        // Warn user that skipping is happening
        if (genotype != "1|0" && genotype != "0|1" && genotype != "1|1")
        {
            invalidGenotypeFormat++;
            if (invalidGenotypeFormat <= 10) { // Limit number of warnings
                std::cerr << "Warning: Line " << genoLineCount << " - Skipping unexpected genotype value (" << genotype 
                          << ") in sample " << sample << " for variant " << variant 
                          << ". Expected values: 1|0, 0|1, or 1|1." << std::endl;
            }
            else if (invalidGenotypeFormat == 11) {
                std::cerr << "Warning: Additional invalid genotype formats found. Suppressing further warnings." << std::endl;
            }
            skippedGenoLines++;
            continue;
        }

        validGenoLines++;

        if (variantToGene.find(variant) != variantToGene.end())
        {
            if (variantToGene[variant].size() > 1)
            {
                multiGeneVariants.insert(variant);
            }

            for (const auto &geneInfoPair : variantToGene[variant])
            {
                std::string gene = geneInfoPair;
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
    gzclose(genotypeFile);

    // Check for processing statistics
    if (skippedGenoLines > 0) {
        std::cerr << "Warning: Skipped " << skippedGenoLines << " lines in genotype file due to format issues." << std::endl;
    }

    if (validGenoLines == 0) {
        std::cerr << "Error: No valid genotype data found in file: " << pathGeno << std::endl;
        return 1;
    }

    // Check if there are no matching variants
    if (uniqueVariantsKept.empty())
    {
        std::cerr << "Error: No matching variants found between the --map and --geno file!" << std::endl;
        return 1;
    }

    // Report on sample count
    if (verbose) {
        std::cerr << "Found " << uniqueSamples.size() << " unique samples in genotype file." << std::endl;
    }

    // only print stuff if user wants it
    if (verbose)
    {
        std::cerr << "Variants in mapping file (kept): " << uniqueVariantsKept.size() << std::endl;
        std::cerr << "Variants not in mapping file (Discarded): " << uniqueVariantsDiscarded.size() << std::endl;
        std::cerr << "Variants mapping to more than one gene: " << multiGeneVariants.size() << std::endl;
    }

    // Print results
    int resultsCount = 0;
    for (const auto &samplePair : sampleGeneHaplotypeVariant)
    {
        const std::string &sample = samplePair.first;
        for (const auto &genePair : samplePair.second)
        {
            resultsCount++;
            const std::string &gene = genePair.first;
            const auto &haplotypeVariantMap = genePair.second;
            std::string chromosome;

            std::set<std::string> haplotype1Variants = haplotypeVariantMap.count(1) ? std::set<std::string>(haplotypeVariantMap.at(1).begin(), haplotypeVariantMap.at(1).end()) : std::set<std::string>();
            std::set<std::string> haplotype2Variants = haplotypeVariantMap.count(2) ? std::set<std::string>(haplotypeVariantMap.at(2).begin(), haplotypeVariantMap.at(2).end()) : std::set<std::string>();

            // Safety check for empty variant sets
            if (haplotype1Variants.empty() && haplotype2Variants.empty()) {
                if (verbose) {
                    std::cerr << "Warning: Empty variant sets for sample " << sample << " and gene " << gene << ". Skipping." << std::endl;
                }
                continue;
            }

            // extract chromosome from string
            std::set<std::string> mergedVariants;
            std::set_union(haplotype1Variants.begin(), haplotype1Variants.end(), haplotype2Variants.begin(), haplotype2Variants.end(), std::inserter(mergedVariants, mergedVariants.begin()));
            
            if (mergedVariants.empty()) {
                std::cerr << "Error: No valid variants found for sample " << sample << " and gene " << gene << std::endl;
                continue;
            }
            
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

    // Final summary
    if (resultsCount == 0) {
        std::cerr << "Warning: No results were generated. Check your input files and parameters." << std::endl;
    }
    else if (verbose) {
        std::cerr << "Successfully generated " << resultsCount << " result entries." << std::endl;
    }

    return 0;
}
