#ifndef CHET_CALLER_HPP
#define CHET_CALLER_HPP

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace call_chets {

/**
 * @brief Structure to hold the result of a single chet/cis call.
 */
struct ChetResult {
  std::string sample;
  std::string chromosome;
  std::string gene;
  std::string call; // "het", "hom", "cis", "chet", "na"
  int dosage;

  float geneScore = 0.0f;
  float haplotype1Score = 0.0f;
  float haplotype2Score = 0.0f;

  // Variants involved
  std::vector<std::string> haplotype1Variants;
  std::vector<std::string> haplotype2Variants;
  std::vector<std::string> unphasedVariants;
};

/**
 * @brief core class for calling compound heterozygotes and cis variants.
 *
 * The ChetCaller class handles the processing of genotype, mapping, and scoring
 * data to identify and qualify compound heterozygous and cis variants variants.
 * It supports various configuration options for collapsing rules (haplotype and
 * gene level), phased/unphased data, and output formatting.
 */
class ChetCaller {
public:
  /**
   * @brief Construct a new Chet Caller object with default settings.
   */
  ChetCaller();

  /**
   * @brief Destroy the Chet Caller object.
   */
  ~ChetCaller();

  // section: Configuration

  /**
   * @brief Set the rule for compressing dosages within a haplotype.
   *
   * @param rule Rule name: "product", "additive", "min", "max". Default:
   * "product".
   */
  void setHaplotypeCollapseRule(const std::string &rule);

  /**
   * @brief Set the rule for combining scores across haplotypes (at gene level).
   *
   * @param rule Rule name: "product", "additive", "min", "max". Default:
   * "product".
   */
  void setGeneCollapseRule(const std::string &rule);

  /**
   * @brief Enable or disable printing of individual haplotype scores.
   *
   * @param show If true, haplotype scores are included in the output.
   */
  void setShowHaplotypeScores(bool show);

  /**
   * @brief Enable or disable printing of the specific variants involved in the
   * call.
   *
   * @param show If true, variant IDs are included in the output.
   */
  void setShowVariants(bool show);

  /**
   * @brief Set the mode to treat input data as unphased.
   *
   * @param unphased If true, phase information is ignored, and only 'het'/'hom'
   * are called.
   */
  void setUnphasedMode(bool unphased);

  /**
   * @brief Enable verbose logging to stderr.
   *
   * @param verbose If true, progress and warnings are printed to stderr.
   */
  void setVerbose(bool verbose);

  // section: Data Loading

  /**
   * @brief Load variant-to-gene mapping from a file.
   *
   * @param path Path to the mapping file (can be gzipped).
   * @return true If loading was successful.
   * @return false If loading failed or file format was invalid.
   */
  bool loadGeneMap(const std::string &path);

  /**
   * @brief Load variant info annotations from a file (optional).
   *
   * @param path Path to the info map file (can be gzipped).
   * @return true If loading was successful (or file missing but optional).
   * @return false If file exists but format is invalid.
   */
  bool loadInfoMap(const std::string &path);

  /**
   * @brief Load variant scores from a file (optional).
   *
   * @param path Path to the score map file (can be gzipped).
   * @return true If loading was successful.
   * @return false If file exists but format is invalid.
   */
  bool loadScoreMap(const std::string &path);

  // section: Processing

  /**
   * @brief Process a genotype file and compute calls.
   *
   * reads the genotype file, matches variants to genes using the loaded map,
   * stores relevant data in internal structures. Call `printResults` after
   * this.
   *
   * @param path Path to the genotype file (can be gzipped).
   * @return true If processing was successful.
   * @return false If file could not be opened or had critical errors.
   */
  bool processGenotypes(const std::string &path);

  // section: Output

  /**
   * @brief Retrieve the calculated results.
   *
   * @return std::vector<ChetResult> A list of all calculated results.
   */
  std::vector<ChetResult> getResults() const;

  /**
   * @brief Calculate final scores and print results to stdout.
   *
   * Iterates through the processed data, applies the collapse rules, and prints
   * the resulting calls/dosages/scores to standard output in tab-separated
   * format.
   */
  void printResults() const;

private:
  // internal state
  std::string defaultHaplotypeCollapseRule;
  std::string defaultGeneCollapseRule;
  bool showHaplotypeScore;
  bool showVariants;
  bool unphasedMode;
  bool verbose;

  // Data storage
  std::map<std::string, std::vector<std::string>> variantToGene;
  std::map<std::pair<std::string, std::string>, std::string> infoMap;
  std::map<std::pair<std::string, std::string>, float> variantGeneScore;

  // Result storage: sample -> gene -> haplotype -> variants
  std::map<std::string,
           std::map<std::string, std::map<int, std::vector<std::string>>>>
      sampleGeneHaplotypeVariant;

  // Stats
  std::set<std::string> uniqueSamples;
  std::set<std::string> uniqueVariantsKept;
  std::set<std::string> uniqueVariantsDiscarded;
  std::set<std::string> multiGeneVariants;
  int invalidFormatVariants;

  // Helper methods
  static bool isValidVariantFormat(const std::string &variant);
  static bool isValidScore(float score);
  static size_t countColumns(const std::string &line);
};

} // namespace call_chets

#endif // CHET_CALLER_HPP
