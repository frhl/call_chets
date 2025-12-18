#include "logging.hpp"
#include "version.hpp"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

// Struct for encoding preview examples
struct EncodingExample {
  std::string gene;
  float r, h, a;
  float raw_aa, raw_Aa, raw_AA;
  float scaled_aa, scaled_Aa, scaled_AA;
  int ac;
};

// Print encoding preview to stderr
void printEncodingPreview(const std::vector<EncodingExample> &examples,
                          bool hasScaling, bool hasMore) {
  if (examples.empty())
    return;

  std::cerr << "\n  * Encoding Preview (first " << examples.size()
            << " genes):" << std::endl;

  if (hasScaling) {
    std::cerr
        << "      Gene        | r      | h      | a      | Raw: aa   | Raw: Aa "
           "  | Raw: AA   | Scaled: aa | Scaled: Aa | Scaled: AA | AC"
        << std::endl;
    std::cerr << "      "
                 "------------|--------|--------|--------|-----------|---------"
                 "--|-----------|------------|------------|------------|----"
              << std::endl;

    for (const auto &ex : examples) {
      std::cerr << "      " << std::left << std::setw(12)
                << ex.gene.substr(0, 12) << "| " << std::fixed
                << std::setprecision(4) << std::setw(6) << ex.r << " | "
                << std::setw(6) << ex.h << " | " << std::setw(6) << ex.a
                << " | " << std::setprecision(6) << std::setw(9) << ex.raw_aa
                << " | " << std::setw(9) << ex.raw_Aa << " | " << std::setw(9)
                << ex.raw_AA << " | " << std::setw(10) << ex.scaled_aa << " | "
                << std::setw(10) << ex.scaled_Aa << " | " << std::setw(10)
                << ex.scaled_AA << " | " << ex.ac << std::endl;
    }
  } else {
    std::cerr << "      Gene        | r      | h      | a      | aa (hom ref) "
                 "| Aa (het)  | AA (hom alt) | AC"
              << std::endl;
    std::cerr << "      "
                 "------------|--------|--------|--------|--------------|------"
                 "-----|--------------|----"
              << std::endl;

    for (const auto &ex : examples) {
      std::cerr << "      " << std::left << std::setw(12)
                << ex.gene.substr(0, 12) << "| " << std::fixed
                << std::setprecision(4) << std::setw(6) << ex.r << " | "
                << std::setw(6) << ex.h << " | " << std::setw(6) << ex.a
                << " | " << std::setprecision(6) << std::setw(12) << ex.raw_aa
                << " | " << std::setw(9) << ex.raw_Aa << " | " << std::setw(12)
                << ex.raw_AA << " | " << ex.ac << std::endl;
    }
  }

  if (hasMore) {
    std::cerr << "      (truncated...)" << std::endl;
  }
  std::cerr << std::endl;
}

void printUsage(const char *path) {
  // Get version info
  std::string version = getFullVersion();

  // Get current date and time
  std::time_t now = std::time(nullptr);
  char timestr[100];
  std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S",
                std::localtime(&now));

  // Print header with version and run date
  std::cerr << "\n[ENCODE_VCF] vcf encoder for CALL_CHETS output"
            << "\n  * Version       : " << version
            << "\n  * Run date      : " << timestr << "\n";

  std::cerr << "\nUsage: " << path << " --input <input> --samples <samples>"
            << " --mode [<additive|recessive|dominance|001|012|010|011>]\n";
  std::cerr << "\nDescription:";
  std::cerr
      << "\n  Converts 'call_chets' output to VCF for downstream analysis.";
  std::cerr << "\n  Results are streamed to standard output.\n";
  std::cerr << "\nOptions:";
  std::cerr << "\n  --input/-i   : Output from 'call_chets'.";
  std::cerr << "\n  --samples/-s : List of samples. One per line. No header.";
  std::cerr << "\n  --mode/-m    : Specify genotype/dosage encoding. Use "
               "'additive' or '012'";
  std::cerr << "\n                 for dosages of 0, 1, and 2. Use 'recessive' "
               "or '001' for";
  std::cerr << "\n                 dosages of 0 and 2. Use 'dominance' to "
               "encode orthogonal";
  std::cerr << "\n                 contribution for dominance effects. '010' "
               "and '011' represent";
  std::cerr << "\n                 custom modes setting bi-allelics to zero or "
               "one respectively.";
  std::cerr
      << "\n  --min-ac     : Filters to genes with sum of DS >= argument.";
  std::cerr
      << "\n  --max-ac     : Filters to genes with sum of DS < argument.\n";
  std::cerr << "\n  --all-info   : Populate INFO column with all details "
               "(relevant when mode='dominance').\n";
  std::cerr << "\n  --global-dom-dosage : Use global min/max dominance dosage "
               "for scaling. \n";
  std::cerr << "\n  --no-dosage-scaling : Disables any scaling of dosages from "
               "0 to 2. \n";
  std::cerr << "\n  --force-chr-out-name <chr_name> : Forces the output "
               "chromosome to be the specified";
  std::cerr
      << "\n                 name regardless of the gene's actual chromosome.";
  std::cerr << "\nExample:";
  std::cerr << "\n  ./encode_vcf called_chets.txt.gz samples.txt additive | "
               "bgzip > out.vcf.gz\n\n";
}

std::vector<std::string> sortChromosomes(const std::set<std::string> &contigs) {
  std::vector<std::string> chromosomes(contigs.begin(), contigs.end());
  std::sort(chromosomes.begin(), chromosomes.end(),
            [](const std::string &a, const std::string &b) {
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
              // If neither are numeric (e.g., chrX, chrY, etc.), then just use
              // lexicographical order
              return a < b;
            });
  return chromosomes;
}

int main(int argc, char *argv[]) {
  if (argc > 0) {
    std::string programName = argv[0];
    if (programName.find("encode_vcf") != std::string::npos) {
      std::cerr << "Warning: 'encode_vcf' is deprecated and will be removed in "
                   "a future release. Please use 'make_pseudo_vcf' instead."
                << std::endl;
    }
  }

  std::string pathInput;
  std::string pathSamples;
  std::string mode = "additive";
  std::string suffix = "";
  int minAC = 0;
  int maxAC = INT_MAX;
  int rowIndex = 0;
  std::string forcedChromosomeName;
  float scalingFactor = 1.0;
  bool scaleDosage = true;
  bool allInfo = false;
  bool globalDomDosage = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h" || arg == "-?" || arg == "-help") {
      printUsage(argv[0]);
      return 0;
    } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
      pathInput = argv[++i];
    } else if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
      pathSamples = argv[++i];
    } else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
      mode = argv[++i];
    } else if (arg == "--global-dom-dosage") {
      globalDomDosage = true;
    } else if (arg == "--force-chr-out-name" && i + 1 < argc) {
      forcedChromosomeName = argv[++i];
    } else if (arg == "--min-ac" && i + 1 < argc) {
      minAC = std::stoi(argv[++i]);
      if (minAC < 0) {
        std::cerr << "Error! Minimum Allele Count (AC) cannot be negative."
                  << std::endl;
        return 1;
      }
    } else if (arg == "--max-ac" && i + 1 < argc) {
      maxAC = std::stoi(argv[++i]);
      if (maxAC < 0) {
        std::cerr << "Error! Maximum Allele Count (AC) cannot be negative."
                  << std::endl;
        return 1;
      }
    } else if (arg == "--scaling-factor" && i + 1 < argc) {
      scalingFactor = std::stof(argv[++i]);
    } else if (arg == "--suffix" && i + 1 < argc) {
      suffix = argv[++i];
    } else if (arg == "--no-dosage-scaling") {
      scaleDosage = false; // Disable dosage scaling which is only relevant when
                           // mode='dominance'
    } else if (arg == "--all-info") {
      allInfo = true;
    } else {
      std::cerr << "Error! Unknown argument: " << arg << std::endl;
      printUsage(argv[0]);
      return 1;
    }
  }

  if (pathInput.empty() || pathSamples.empty()) {
    std::cerr << "Error! --input and --samples are required." << std::endl;
    printUsage(argv[0]);
    return 1;
  }

  std::map<std::string, std::string> files;
  files["Input"] = pathInput;
  files["Samples"] = pathSamples;

  std::map<std::string, std::string> params;
  params["Mode"] = mode;
  params["Scale Dosage"] = scaleDosage ? "Yes" : "No";
  if (minAC > 0)
    params["Min AC"] = std::to_string(minAC);
  if (maxAC < INT_MAX)
    params["Max AC"] = std::to_string(maxAC);
  if (!forcedChromosomeName.empty())
    params["Force Chr"] = forcedChromosomeName;
  if (globalDomDosage)
    params["Global Dom"] = "Yes";

  call_chets::printHeader("ENCODE_VCF", "vcf encoder for CALL_CHETS output",
                          files, params);

  gzFile longFile = gzopen(pathInput.c_str(), "rb");
  std::ifstream sampleFile(pathSamples);

  if (!longFile) {
    std::cerr << "Error: Cannot open gzipped <input> file for reading: "
              << pathInput << std::endl;
    return 1;
  }

  if (!sampleFile) {
    std::cerr << "Error: Cannot open <samples> file for reading: "
              << pathSamples << std::endl;
    gzclose(longFile);
    return 1;
  }

  // normalise mode input
  if (mode == "recessive")
    mode = "001";
  else if (mode == "additive")
    mode = "012";

  if (mode != "001" && mode != "012" && mode != "010" && mode != "011" &&
      mode != "dominance") {
    std::cerr << "Error: Invalid dosage encoding mode provided. Only "
                 "'012|additive, '001|recessive', '010', '011' or 'dominance'."
              << std::endl;
    printUsage(argv[0]);
    gzclose(longFile);
    sampleFile.close();
    return 1;
  }

  // Collect all samples into a set
  std::set<std::string> samples;
  std::string sampleLine;
  int sampleCount = 0;
  while (std::getline(sampleFile, sampleLine)) {
    // Remove trailing and leading whitespaces
    sampleLine.erase(0, sampleLine.find_first_not_of(" \t\n\r"));
    sampleLine.erase(sampleLine.find_last_not_of(" \t\n\r") + 1);

    // Skip empty lines
    if (sampleLine.empty()) {
      continue;
    }

    // Split the line using white spaces
    std::istringstream iss(sampleLine);
    std::string word;
    std::vector<std::string> wordsInLine;

    while (iss >> word) {
      wordsInLine.push_back(word);
    }

    // Check if line contains more than one word
    if (wordsInLine.size() != 1) {
      std::cerr << "Error: Invalid samples file format. Expected one sample "
                   "per line but detected more on line "
                << (sampleCount + 1) << "!" << std::endl;
      printUsage(argv[0]);
      sampleFile.close();
      gzclose(longFile);
      return 1;
    }

    // Add the sample to our set
    samples.insert(sampleLine);
    sampleCount++;
  }

  // Check if sample file is empty
  if (samples.empty()) {
    std::cerr << "Error: No valid samples found in file: " << pathSamples
              << std::endl;
    gzclose(longFile);
    sampleFile.close();
    return 1;
  }

  // set total alleles
  int totalSamples = samples.size();
  int geneAN = totalSamples * 2;

  std::map<std::string, std::map<std::string, int>> geneSampleDosage;
  std::map<std::string, std::string> geneToChromosome;
  std::map<std::string, int> geneAC;
  std::map<std::string, int> geneBI;
  std::map<std::string, int> geneChet;
  std::map<std::string, int> geneHom;
  std::map<std::string, int> geneHet;
  std::map<std::string, int> geneCis;

  // Process long format file
  std::string line, sample, chromosome, gene, configuration;
  std::set<std::string> contigs;
  char buffer[4096];
  float dosage;
  std::string variantInfo;
  int lineCount = 0;
  int matchingSampleCount = 0;
  bool isFirstLine = true;

  while (gzgets(longFile, buffer, sizeof(buffer))) {
    lineCount++;
    std::string line(buffer);
    std::istringstream iss(line);

    // Check for empty lines or comment lines
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // Count expected columns
    std::vector<std::string> columns;
    std::string field;
    while (iss >> field) {
      columns.push_back(field);
    }

    // Validate number of columns
    if (columns.size() < 6) {
      std::cerr << "Error: Line " << lineCount
                << " has fewer than the expected 6 columns. Found "
                << columns.size() << " columns." << std::endl;
      gzclose(longFile);
      return 1;
    }

    // Reset stringstream to parse columns individually
    iss.clear();
    iss.str(line);

    // Parse columns
    iss >> sample >> chromosome >> gene >> configuration >> dosage >>
        variantInfo;

    // Check if configuration is one of the expected values
    if (configuration != "chet" && configuration != "het" &&
        configuration != "cis" && configuration != "hom") {
      std::cerr << "Error: 4th column '" << configuration
                << "' is not one of the expected values (chet, het, cis, hom) "
                   "in line "
                << lineCount << ": " << line << std::endl;
      gzclose(longFile);
      return 1;
    }

    // Validate dosage value
    if (dosage < 0.0f || dosage > 2.0f) {
      std::cerr << "Error: Dosage value '" << dosage
                << "' is outside the expected range [0-2] in line " << lineCount
                << ": " << line << std::endl;
      gzclose(longFile);
      return 1;
    }

    // Skip processing this line if the sample is not in the set
    if (samples.find(sample) == samples.end())
      continue;

    // If we reach here, we have at least one matching sample
    matchingSampleCount++;

    // for recessive mode, set hets to zero
    if (mode == "001" && dosage == 1.0f) {
      dosage = 0;
    }

    // for het mode, set bi-allelic as zero
    if (mode == "010" && dosage == 2.0f) {
      dosage = 0;
    }

    // for '011' mode, set bi-allelic to 1
    if (mode == "011" && dosage == 2.0f) {
      dosage = 1;
    }

    // count indiviual sites
    if (configuration == "chet") {
      geneChet[gene] += 1;
      geneAC[gene] += 2;
    } else if (configuration == "hom") {
      geneHom[gene] += 1;
      geneAC[gene] += 2;
    } else if (configuration == "cis") {
      geneCis[gene] += 1;
      geneAC[gene] += 1;
    } else if (configuration == "het") {
      geneHet[gene] += 1;
      geneAC[gene] += 1;
    }

    // count instances of bi-allelic
    if (configuration == "chet" || configuration == "hom") {
      geneBI[gene] += 1;
    }

    // count instances of bi-allelic
    geneToChromosome[gene] = chromosome;
    geneSampleDosage[gene][sample] = dosage;
    contigs.insert(chromosome);
  }

  // Check if the input file was empty
  if (lineCount == 0) {
    std::cerr << "Error: Input file is empty: " << pathInput << std::endl;
    gzclose(longFile);
    return 1;
  }

  // Check if there was no overlap between samples in input file and sample list
  if (matchingSampleCount == 0) {
    std::cerr << "Error: No overlap found between samples in input file and "
                 "provided sample list."
              << std::endl;
    std::cerr << "       Please check that sample names match between files."
              << std::endl;
    gzclose(longFile);
    return 1;
  }

  gzclose(longFile);

  // compute global minimum or maximum dominance scaling
  // which will be applied when the global-dom-dosage arugment is on
  float globalMinDomDosage = std::numeric_limits<float>::max();
  float globalMaxDomDosage = std::numeric_limits<float>::min();

  // Collect first 5 genes for encoding preview
  std::vector<EncodingExample> encodingExamples;
  const int MAX_EXAMPLES = 5;
  int totalGenesProcessed = 0;

  if (globalDomDosage && mode == "dominance") {

    for (const auto &gene : geneAC) {
      int currentAN = samples.size() * 2;
      float aa_count = static_cast<float>(
          (currentAN / 2) -
          (geneCis[gene.first] + geneHet[gene.first] + geneBI[gene.first]));
      float Aa_count =
          static_cast<float>(geneHet[gene.first] + geneCis[gene.first]);
      float AA_count = static_cast<float>(geneBI[gene.first]);

      float r = aa_count / (currentAN / 2);
      float h = Aa_count / (currentAN / 2);
      float a = AA_count / (currentAN / 2);

      float dom_dosage_aa = -h * a;
      float dom_dosage_Aa = 2 * a * r;
      float dom_dosage_AA = -h * r;

      globalMinDomDosage = std::min(
          {globalMinDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
      globalMaxDomDosage = std::max(
          {globalMaxDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

      totalGenesProcessed++;
    }

    // Second pass to collect examples with scaled values (now that we have
    // global min/max)
    for (const auto &gene : geneAC) {
      if (encodingExamples.size() >= MAX_EXAMPLES)
        break;

      // Skip genes with no homozygotes (they won't be output anyway)
      if (geneBI[gene.first] == 0)
        continue;

      int currentAN = samples.size() * 2;
      float aa_count = static_cast<float>(
          (currentAN / 2) -
          (geneCis[gene.first] + geneHet[gene.first] + geneBI[gene.first]));
      float Aa_count =
          static_cast<float>(geneHet[gene.first] + geneCis[gene.first]);
      float AA_count = static_cast<float>(geneBI[gene.first]);

      float r = aa_count / (currentAN / 2);
      float h = Aa_count / (currentAN / 2);
      float a = AA_count / (currentAN / 2);

      float dom_dosage_aa = -h * a;
      float dom_dosage_Aa = 2 * a * r;
      float dom_dosage_AA = -h * r;

      EncodingExample ex;
      ex.gene = gene.first;
      ex.r = r;
      ex.h = h;
      ex.a = a;
      ex.raw_aa = dom_dosage_aa;
      ex.raw_Aa = dom_dosage_Aa;
      ex.raw_AA = dom_dosage_AA;
      ex.scaled_aa = 2 * ((dom_dosage_aa - globalMinDomDosage) /
                          (globalMaxDomDosage - globalMinDomDosage));
      ex.scaled_Aa = 2 * ((dom_dosage_Aa - globalMinDomDosage) /
                          (globalMaxDomDosage - globalMinDomDosage));
      ex.scaled_AA = 2 * ((dom_dosage_AA - globalMinDomDosage) /
                          (globalMaxDomDosage - globalMinDomDosage));
      ex.ac = gene.second;
      encodingExamples.push_back(ex);
    }

    // Print global scaling info
    std::cerr << "\n  * Global Dominance Scaling:" << std::endl;
    std::cerr << "      + Min Dosage: " << std::fixed << std::setprecision(6)
              << globalMinDomDosage << std::endl;
    std::cerr << "      + Max Dosage: " << globalMaxDomDosage << std::endl;
    std::cerr << "      + Target Range: [0.0, 2.0]" << std::endl;

    // Print encoding preview
    printEncodingPreview(encodingExamples, true,
                         totalGenesProcessed > MAX_EXAMPLES);
  }

  // For non-global dominance mode or other modes, collect examples separately
  if (!globalDomDosage && mode == "dominance") {
    for (const auto &gene : geneAC) {
      // Skip genes with no homozygotes (they won't be output anyway)
      if (geneBI[gene.first] == 0)
        continue;

      // Skip if not meeting AC filter criteria
      if (gene.second < minAC || gene.second >= maxAC)
        continue;

      if (encodingExamples.size() >= MAX_EXAMPLES) {
        totalGenesProcessed++;
        continue;
      }

      int currentAN = samples.size() * 2;
      float aa_count = static_cast<float>(
          (currentAN / 2) -
          (geneCis[gene.first] + geneHet[gene.first] + geneBI[gene.first]));
      float Aa_count =
          static_cast<float>(geneHet[gene.first] + geneCis[gene.first]);
      float AA_count = static_cast<float>(geneBI[gene.first]);

      float r = aa_count / (currentAN / 2);
      float h = Aa_count / (currentAN / 2);
      float a = AA_count / (currentAN / 2);

      float dom_dosage_aa = -h * a;
      float dom_dosage_Aa = 2 * a * r;
      float dom_dosage_AA = -h * r;

      // Per-variant min/max
      float minDom = std::min({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
      float maxDom = std::max({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

      EncodingExample ex;
      ex.gene = gene.first;
      ex.r = r;
      ex.h = h;
      ex.a = a;
      ex.raw_aa = dom_dosage_aa;
      ex.raw_Aa = dom_dosage_Aa;
      ex.raw_AA = dom_dosage_AA;
      ex.scaled_aa = 2 * ((dom_dosage_aa - minDom) / (maxDom - minDom));
      ex.scaled_Aa = 2 * ((dom_dosage_Aa - minDom) / (maxDom - minDom));
      ex.scaled_AA = 2 * ((dom_dosage_AA - minDom) / (maxDom - minDom));
      ex.ac = gene.second;
      encodingExamples.push_back(ex);
      totalGenesProcessed++;
    }

    std::cerr << "\n  * Per-Gene Dominance Scaling (each gene scaled "
                 "independently to [0, 2])"
              << std::endl;
    printEncodingPreview(encodingExamples, true,
                         totalGenesProcessed > MAX_EXAMPLES);
  }

  // sort chromosomes
  std::vector<std::string> sortedContigs = sortChromosomes(contigs);

  // Print the output header
  std::cout << "##fileformat=VCFv4.2\n";
  std::cout << "##EncodingMode=" << mode << "\n"; // Add this lin
  std::cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";

  // Check if forced chromosome name is provided and add it to the header
  if (!forcedChromosomeName.empty()) {
    std::cout << "##contig=<ID=" << forcedChromosomeName << ">\n";
  } else {
    // If no forced chromosome name, add all known contigs
    for (const auto &chr : sortedContigs) {
      std::cout << "##contig=<ID=" << chr << ">\n";
    }
  }

  std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
  std::cout
      << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n";
  std::cout
      << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Number\">\n";
  std::cout << "##INFO=<ID=BI,Number=1,Type=Integer,Description=\"Bi-allelic "
               "Count\">\n";
  std::cout << "##INFO=<ID=CHET,Number=1,Type=Integer,Description=\"Compound "
               "heterozygous pseudo Count\">\n";
  std::cout << "##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Homozygous "
               "Count\">\n";
  std::cout << "##INFO=<ID=HET,Number=1,Type=Integer,Description="
               "\"Heterozygous Count\">\n";
  std::cout << "##INFO=<ID=CIS,Number=1,Type=Integer,Description=\"Cis pseudo "
               "Count\">\n";
  // add dominance relevant information
  if (mode == "dominance") {
    std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic references (aa)\">\n";
    std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of "
                 "heterozygotes (Aa)\">\n";
    std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic alternates (AA)\">\n";
    if (allInfo == true) {
      std::cout << "##INFO=<ID=minDosage,Number=1,Type=Float,Description="
                   "\"Minimum scaled dosage value for dominance mode\">\n";
      std::cout << "##INFO=<ID=maxDosage,Number=1,Type=Float,Description="
                   "\"Maximum scaled dosage value for dominance mode\">\n";
      std::cout << "##INFO=<ID=DS0,Number=1,Type=Float,Description=\"Scaled "
                   "dosage for genotype aa in dominance mode\">\n";
      std::cout << "##INFO=<ID=DS1,Number=1,Type=Float,Description=\"Scaled "
                   "dosage for genotype Aa in dominance mode\">\n";
      std::cout << "##INFO=<ID=DS2,Number=1,Type=Float,Description=\"Scaled "
                   "dosage for genotype AA in dominance mode\">\n";
    }
  }
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (const auto &sampleName : samples) {
    std::cout << "\t" << sampleName;
  }

  std::cout << std::endl;
  // Print the output data
  for (const auto &chr : sortedContigs) {
    if (forcedChromosomeName.empty()) {
      rowIndex = 0;
    }

    for (const auto &genePair : geneSampleDosage) {
      if (geneToChromosome[genePair.first] ==
          chr) { // Only process genes on the current chromosome
        rowIndex++;
        int currentAN = geneAN;
        int currentAC = geneAC[genePair.first];
        int currentBI = geneBI[genePair.first];
        int currentChet = geneChet[genePair.first];
        int currentHom = geneHom[genePair.first];
        int currentHet = geneHet[genePair.first];
        int currentCis = geneCis[genePair.first];

        // Skip processing this gene/pseudo-variant if
        // in dominance mode and no homozygotes present
        if (mode == "dominance" && currentBI == 0) {
          continue;
        }

        // derive reference alleles
        float aa_count = static_cast<float>(
            (currentAN / 2) - (currentCis + currentHet + currentBI));
        float Aa_count = static_cast<float>(currentHet + currentCis);
        float AA_count = static_cast<float>(currentBI);

        // rename accordinly
        float r = static_cast<float>(aa_count / (currentAN / 2));
        float h = static_cast<float>(Aa_count / (currentAN / 2));
        float a = static_cast<float>(AA_count / (currentAN / 2));

        // Pre-compute the three possible dosage configurations
        // when performing dominance deviation
        float dom_dosage_aa = -h * a;    // Dosage for genotype aa
        float dom_dosage_Aa = 2 * a * r; // Dosage for genotype Aa
        float dom_dosage_AA = -h * r;    // Dosage for genotype AA

        // Find min and max dosage values for scaling
        float minDomDosage =
            std::min({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
        float maxDomDosage =
            std::max({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

        // overwrite dominance dosage
        if (globalDomDosage && mode == "dominance") {
          minDomDosage = globalMinDomDosage;
          maxDomDosage = globalMaxDomDosage;
        }

        // only output variants that fit our criteria
        if (currentAC >= minAC && currentAC < maxAC) {

          // get left side of VCF body and force chromosome name if the
          // --force-chr-name is given
          std::cout << (forcedChromosomeName.empty()
                            ? geneToChromosome[genePair.first]
                            : forcedChromosomeName)
                    << "\t" << rowIndex << "\t" << genePair.first + suffix
                    << "\tA\tB\t.\t.\t"
                    << "AC=" << currentAC << ";AN=" << currentAN
                    << ";BI=" << currentBI << ";CHET=" << currentChet
                    << ";HOM=" << currentHom << ";CIS=" << currentCis;

          // be verbose at times'
          if ((mode == "dominance") & (allInfo == true)) {
            std::cout << ";r=" << r << ";h=" << h << ";a=" << a
                      << ";minDosage=" << minDomDosage
                      << ";maxDosage=" << maxDomDosage << ";DS0="
                      << 2 * (((-h * a) - minDomDosage) /
                              (maxDomDosage - minDomDosage))
                      << ";DS1="
                      << 2 * (((2 * a * r) - minDomDosage) /
                              (maxDomDosage - minDomDosage))
                      << ";DS2="
                      << 2 * (((-h * r) - minDomDosage) /
                              (maxDomDosage - minDomDosage));
          }

          std::cout << "\tDS"; //

          // get right side of body VCF
          for (const auto &sample : samples) {
            std::cout << "\t";
            float dosage = 0;
            if (genePair.second.find(sample) != genePair.second.end()) {
              dosage = genePair.second.at(sample);
            }

            if (mode == "dominance") {
              if (dosage == 0.0) {
                dosage = -h * a;
              } else if (dosage == 1.0) {
                dosage = 2 * a * r;
              } else if (dosage == 2.0) {
                dosage = -h * r;
              }
              if (scaleDosage) {
                dosage = 2 * ((dosage - minDomDosage) /
                              (maxDomDosage - minDomDosage));
              }

              if (scalingFactor != 1.0) {
                dosage *= scalingFactor;
              }
            }
            std::cout << dosage;
          }
          std::cout << std::endl;
        }
      }
    }
  }

  return 0;
}
