#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <zlib.h>

#include "ChetCaller.hpp"
#include "logging.hpp"
#include "version.hpp"

// Function to print usage
void printUsage(const char *path) {
  std::string version = getFullVersion();

  std::time_t now = std::time(nullptr);
  char timestr[100];
  std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S",
                std::localtime(&now));

  std::cerr << "\n[CALL_CHETS] compound heterozygosity and cis variant caller"
            << "\n  * Version       : " << version
            << "\n  * Run date      : " << timestr << "\n"
            << std::endl;

  std::cerr
      << "Usage: " << path
      << " --geno <Genotype File> --map <Mapping File> [options]\n\n"
      << "Required Arguments:\n"
      << "  --geno/-g <file>           Path to phased/unphased genotype file "
         "(gzipped)\n"
      << "  --gene-map/-m <file>       Path to variant-to-gene mapping file\n\n"
      << "Optional Arguments:\n"
      << "  --info-map/-i <file>       Path to variant info mapping file (AF, "
         "AC, etc.)\n"
      << "  --score-map/-p <file>      Path to variant score mapping file\n"
      << "  --unphased                 Run in unphased mode (treats all "
         "variants as unphased)\n"
      << "  --show-variants/-sv        Show detailed variant information in "
         "output\n"
      << "  --verbose/-v               Enable verbose logging\n\n"
      << "Score/Collapse Options (requires --score-map):\n"
      << "  --haplotype-collapse/-hc   Rule for collapsing haplotype scores "
         "[product|min|max|additive] (default: product)\n"
      << "  --gene-collapse/-gc        Rule for collapsing gene scores "
         "[product|min|max|additive] (default: product)\n"
      << "  --show-haplotype-scores/-shs Show calculated scores for each "
         "haplotype\n"
      << std::endl;
}

// Main function
int main(int argc, char *argv[]) {
  if (argc > 0) {
    std::string programName = argv[0];
    if (programName.find("call_chets") != std::string::npos) {
      std::cerr << "Warning: 'call_chets' is deprecated and will be removed in "
                   "a future release. Please use 'interpret_phase' instead."
                << std::endl;
    }
  }

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
  bool unphasedMode = false;

  // Argument parsing
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      // Re-implement full printUsage or keep the simplified one
      // Ideally we keep the logic but for now let's just minimal
      printUsage(argv[0]);
      return 0;
    } else if ((arg == "--geno" || arg == "-g") && i + 1 < argc)
      pathGeno = argv[++i];
    else if ((arg == "--gene-map" || arg == "-m") && i + 1 < argc)
      pathMap = argv[++i];
    else if ((arg == "--info-map" || arg == "-i") && i + 1 < argc)
      pathInfoMap = argv[++i];
    else if ((arg == "--score-map" || arg == "-p") && i + 1 < argc)
      pathScoreMap = argv[++i];
    else if ((arg == "--haplotype-collapse" || arg == "-hc") && i + 1 < argc) {
      haplotypeCollapseRule = argv[++i];
      haplotypeCollapseRuleSet = true;
    } else if ((arg == "--gene-collapse" || arg == "-gc") && i + 1 < argc) {
      geneCollapseRule = argv[++i];
      geneCollapseRuleSet = true;
    } else if (arg == "--show-haplotype-scores" || arg == "-shs")
      showHaplotypeScore = true;
    else if (arg == "--show-variants" || arg == "-sv")
      showVariants = true;
    else if (arg == "--unphased")
      unphasedMode = true;
    else if (arg == "--debug")
      verbose = true;
    else {
      std::cerr << "Error! Unknown or incomplete argument: " << arg
                << std::endl;
      // printUsage(argv[0]);
      return 1;
    }
  }

  if (pathGeno.empty() || pathMap.empty()) {
    std::cerr << "Error! Both --geno and --map must be provided." << std::endl;
    return 1;
  }

  // Initialize ChetCaller
  call_chets::ChetCaller caller;
  caller.setVerbose(verbose);
  caller.setUnphasedMode(unphasedMode);
  caller.setShowVariants(showVariants);
  caller.setShowHaplotypeScores(showHaplotypeScore);

  if (haplotypeCollapseRuleSet)
    caller.setHaplotypeCollapseRule(haplotypeCollapseRule);
  if (geneCollapseRuleSet)
    caller.setGeneCollapseRule(geneCollapseRule);

  // Prepare logging
  std::map<std::string, std::string> files;
  files["Genotypes"] = pathGeno;
  files["Gene Map"] = pathMap;
  if (!pathInfoMap.empty())
    files["Info Map"] = pathInfoMap;
  if (!pathScoreMap.empty())
    files["Score Map"] = pathScoreMap;

  std::map<std::string, std::string> params;
  params["Mode"] = unphasedMode ? "Unphased" : "Phased";
  params["Verbose"] = verbose ? "Yes" : "No";
  params["Show Variants"] = showVariants ? "Yes" : "No";

  if (!unphasedMode && !pathScoreMap.empty()) {
    params["Haplo Collapse"] = haplotypeCollapseRule;
    params["Gene Collapse"] = geneCollapseRule;
    params["Show Scores"] = showHaplotypeScore ? "Yes" : "No";
  }

  call_chets::printHeader("CALL_CHETS",
                          "compound heterozygosity and cis variant caller",
                          files, params);

  std::cerr << "Reading genotype data..." << std::endl;

  // Validate and set collapse rules
  if (haplotypeCollapseRuleSet) {
    if (haplotypeCollapseRule != "product" &&
        haplotypeCollapseRule != "additive" && haplotypeCollapseRule != "max" &&
        haplotypeCollapseRule != "min") {
      std::cerr << "Error! Invalid haplotype-collapse rule: "
                << haplotypeCollapseRule << std::endl;
      return 1;
    }
    caller.setHaplotypeCollapseRule(haplotypeCollapseRule);
  }

  if (geneCollapseRuleSet) {
    if (geneCollapseRule != "product" && geneCollapseRule != "additive" &&
        geneCollapseRule != "max" && geneCollapseRule != "min") {
      std::cerr << "Error! Invalid gene-collapse rule: " << geneCollapseRule
                << std::endl;
      return 1;
    }
    caller.setGeneCollapseRule(geneCollapseRule);
  }

  caller.setShowHaplotypeScores(showHaplotypeScore);
  caller.setShowVariants(showVariants);

  // Warnings for ignored args
  if (pathScoreMap.empty()) {
    if (showHaplotypeScore)
      std::cerr << "Warning: ignored --show-haplotype-score argument. This can "
                   "only be used when --score-map is also specified!"
                << std::endl;
    if (geneCollapseRuleSet)
      std::cerr << "Warning: ignored --gene-collapse argument. This can only "
                   "be used when --score-map is also specified!"
                << std::endl;
    if (haplotypeCollapseRuleSet)
      std::cerr << "Warning: ignored --haplotype-collapse argument. This can "
                   "only be used when --score-map is also specified!"
                << std::endl;
  }

  if (unphasedMode) {
    if (!pathScoreMap.empty())
      std::cerr << "Warning: --score-map is ignored in --unphased mode."
                << std::endl;
    if (showHaplotypeScore)
      std::cerr
          << "Warning: --show-haplotype-scores is ignored in --unphased mode."
          << std::endl;
    if (haplotypeCollapseRuleSet || geneCollapseRuleSet)
      std::cerr << "Warning: collapse rules are ignored in --unphased mode."
                << std::endl;
  }

  // Load data
  if (!caller.loadGeneMap(pathMap))
    return 1;

  if (!pathInfoMap.empty()) {
    if (!caller.loadInfoMap(pathInfoMap))
      return 1;
  }

  if (!pathScoreMap.empty()) {
    if (!caller.loadScoreMap(pathScoreMap))
      return 1;
  }

  // Process
  if (!caller.processGenotypes(pathGeno))
    return 1;

  // Print run stats
  caller.printStats();

  // Output
  caller.printResults();

  return 0;
}