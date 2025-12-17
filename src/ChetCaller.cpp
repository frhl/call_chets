#include "ChetCaller.hpp"
#include <algorithm>
#include <iostream>
#include <regex>
#include <sstream>
#include <zlib.h>

namespace call_chets {

ChetCaller::ChetCaller()
    : defaultHaplotypeCollapseRule("product"),
      defaultGeneCollapseRule("product"), showHaplotypeScore(false),
      showVariants(false), unphasedMode(false), verbose(false),
      invalidFormatVariants(0) {}

ChetCaller::~ChetCaller() {}

void ChetCaller::setHaplotypeCollapseRule(const std::string &rule) {
  defaultHaplotypeCollapseRule = rule;
}

void ChetCaller::setGeneCollapseRule(const std::string &rule) {
  defaultGeneCollapseRule = rule;
}

void ChetCaller::setShowHaplotypeScores(bool show) {
  showHaplotypeScore = show;
}

void ChetCaller::setShowVariants(bool show) { showVariants = show; }

void ChetCaller::setUnphasedMode(bool unphased) { unphasedMode = unphased; }

void ChetCaller::setVerbose(bool v) { verbose = v; }

bool ChetCaller::isValidVariantFormat(const std::string &variant) {
  static const std::regex pattern(
      "^(chr)?[0-9XYM]{1,2}:[0-9]+:[ACGT]+:[ACGT]+$");
  return std::regex_match(variant, pattern);
}

bool ChetCaller::isValidScore(float score) {
  return score >= 0.0f && score <= 1.0f;
}

size_t ChetCaller::countColumns(const std::string &line) {
  std::stringstream ss(line);
  std::string temp;
  size_t count = 0;
  while (ss >> temp) {
    count++;
  }
  return count;
}

bool ChetCaller::loadGeneMap(const std::string &path) {
  gzFile mappingFile = gzopen(path.c_str(), "rb");
  if (!mappingFile) {
    std::cerr << "Error: Cannot open --mapping file for reading: " << path
              << std::endl;
    return false;
  }

  char buf[4096];
  bool isFirstLineMappingFile = true;
  int mappingLineCount = 0;
  int validMappingLines = 0;
  invalidFormatVariants = 0;

  std::string variant, gene;

  while (gzgets(mappingFile, buf, sizeof(buf))) {
    mappingLineCount++;
    std::string line(buf);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    if (line.empty())
      continue;

    size_t columnCount = countColumns(line);
    if (columnCount < 2 && !isFirstLineMappingFile) {
      std::cerr << "Error: Line " << mappingLineCount
                << " in mapping file has fewer than required 2 columns: '"
                << line << "'" << std::endl;
      gzclose(mappingFile);
      return false;
    }

    std::stringstream ss(line);
    ss >> variant >> gene;

    if (ss.fail() || variant.empty() || gene.empty()) {
      if (!isFirstLineMappingFile) {
        std::cerr << "Error: Failed to extract two columns from line "
                  << mappingLineCount << " (variant gene): '" << line
                  << "'. Please, fix this line in --gene-map and retry."
                  << std::endl;
        gzclose(mappingFile);
        return false;
      }
    } else {
      if (!isFirstLineMappingFile) {
        if (!isValidVariantFormat(variant)) {
          invalidFormatVariants++;
          if (verbose) {
            std::cerr << "Warning: Line " << mappingLineCount
                      << " - Invalid variant format: " << variant
                      << ". Expected format: chr:pos:ref:alt" << std::endl;
          }
        }
        variantToGene[variant].push_back(gene);
        validMappingLines++;
      }
    }
    isFirstLineMappingFile = false;
  }
  gzclose(mappingFile);

  if (validMappingLines == 0) {
    std::cerr << "Error: No valid mapping data found in: " << path << std::endl;
    return false;
  }

  if (invalidFormatVariants > 0) {
    std::cerr << "Warning: " << invalidFormatVariants
              << " variants in mapping file have invalid format." << std::endl;
  }

  return true;
}

bool ChetCaller::loadInfoMap(const std::string &path) {
  gzFile infoMapFile = gzopen(path.c_str(), "rb");
  if (!infoMapFile) {
    std::cerr << "Warning: Cannot open --info-map file for reading: " << path
              << ". Continuing without info data." << std::endl;
    return true; // Not a fatal error
  }

  // Check empty handled by simple read check loop or externally, skipping
  // explicit check for brevity as usage pattern suggests simple read loop is
  // fine usually. But let's check basic read.

  char infoBuf[4096];
  bool isFirstLine = true;
  int infoLineCount = 0;
  int validInfoLines = 0;

  while (gzgets(infoMapFile, infoBuf, sizeof(infoBuf))) {
    infoLineCount++;
    std::string line(infoBuf);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    if (line.empty())
      continue;

    size_t columnCount = countColumns(line);
    if (columnCount < 3 && !isFirstLine) {
      std::cerr << "Error: Line " << infoLineCount
                << " in info-map file has fewer than required 3 columns: '"
                << line << "'" << std::endl;
      gzclose(infoMapFile);
      return false;
    }

    std::stringstream ss(line);
    std::string variant, gene, info;
    ss >> variant >> gene >> info;

    if (ss.fail() || ss.bad()) {
      if (!isFirstLine) {
        std::cerr
            << "Error: Failed to extract variant, gene, and info from line "
            << infoLineCount << ": '" << line
            << "'. Please, fix this line in --info-map and retry." << std::endl;
        gzclose(infoMapFile);
        return false;
      }
    } else {
      if (!isFirstLine) {
        if (!isValidVariantFormat(variant)) {
          if (verbose) {
            std::cerr << "Warning: Line " << infoLineCount
                      << " in info-map - Invalid variant format: " << variant
                      << ". Expected format: chr:pos:ref:alt" << std::endl;
          }
        }
        infoMap[std::make_pair(variant, gene)] = info;
        validInfoLines++;
      }
    }
    isFirstLine = false;
  }

  if (validInfoLines == 0) {
    std::cerr << "Warning: No valid info mapping data found in: " << path
              << ". Continuing without info data." << std::endl;
  }

  gzclose(infoMapFile);
  return true;
}

bool ChetCaller::loadScoreMap(const std::string &path) {
  gzFile scoreMapFile = gzopen(path.c_str(), "rb");
  if (!scoreMapFile) {
    std::cerr << "Warning: Cannot open --score-map file for reading: " << path
              << ". Continuing without score data." << std::endl;
    return true;
  }

  char pathoBuf[4096];
  bool isFirstLine = true;
  int scoreLineCount = 0;
  int validScoreLines = 0;
  int invalidScoreValues = 0;

  while (gzgets(scoreMapFile, pathoBuf, sizeof(pathoBuf))) {
    scoreLineCount++;
    std::string line(pathoBuf);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    if (line.empty())
      continue;

    size_t columnCount = countColumns(line);
    if (columnCount < 3 && !isFirstLine) {
      std::cerr << "Error: Line " << scoreLineCount
                << " in score-map file has fewer than required 3 columns: '"
                << line << "'" << std::endl;
      gzclose(scoreMapFile);
      return false;
    }

    std::stringstream ss(line);
    std::string variant, gene;
    float newScore;
    ss >> variant >> gene >> newScore;

    if (ss.fail() || ss.bad()) {
      if (!isFirstLine) {
        std::cerr
            << "Error: Failed to extract variant, gene, and score from line "
            << scoreLineCount << ": '" << line
            << "'. Please, fix this line in --score-map and retry."
            << std::endl;
        gzclose(scoreMapFile);
        return false;
      }
    } else {
      if (!isFirstLine) {
        if (!isValidVariantFormat(variant)) {
          if (verbose) {
            std::cerr << "Warning: Line " << scoreLineCount
                      << " in score-map - Invalid variant format: " << variant
                      << ". Expected format: chr:pos:ref:alt" << std::endl;
          }
        }

        if (!isValidScore(newScore)) {
          invalidScoreValues++;
          if (verbose) {
            std::cerr << "Warning: Line " << scoreLineCount << " - Score value "
                      << newScore << " is outside valid range [0-1]"
                      << std::endl;
          }
        }

        variantGeneScore[std::make_pair(variant, gene)] = newScore;
        validScoreLines++;
      }
    }
    isFirstLine = false;
  }

  if (validScoreLines == 0) {
    std::cerr << "Warning: No valid score mapping data found in: " << path
              << ". Continuing without score data." << std::endl;
  }

  if (invalidScoreValues > 0) {
    std::cerr << "Warning: " << invalidScoreValues
              << " scores in score-map file are outside the valid range [0-1]."
              << std::endl;
  }

  gzclose(scoreMapFile);
  return true;
}

bool ChetCaller::processGenotypes(const std::string &path) {
  gzFile genotypeFile = gzopen(path.c_str(), "rb");
  if (!genotypeFile) {
    std::cerr << "Error: Cannot open --geno file for reading: " << path
              << std::endl;
    return false;
  }

  char buf[4096];
  bool isFirstLine = true;
  int genoLineCount = 0;
  int validGenoLines = 0;
  int skippedGenoLines = 0;
  int invalidFormatGenoVariants = 0;
  int invalidGenotypeFormat = 0;

  uniqueSamples.clear();
  uniqueVariantsKept.clear();
  uniqueVariantsDiscarded.clear();
  multiGeneVariants.clear();

  while (gzgets(genotypeFile, buf, sizeof(buf))) {
    genoLineCount++;
    std::string line(buf);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    if (line.empty())
      continue;

    size_t columnCount = countColumns(line);
    if (columnCount < 3) {
      std::cerr << "Error: Line " << genoLineCount
                << " in genotype file has fewer than required 3 columns: '"
                << line << "'" << std::endl;
      skippedGenoLines++;
      continue;
    }

    std::stringstream ss(line);
    std::string sample, variant, genotype;
    ss >> sample >> variant >> genotype;

    if (isFirstLine) {
      bool isPhasedGenotype = (genotype == "1|0" || genotype == "0|1" ||
                               genotype == "1|1" || genotype == "0|0");
      bool isUnphasedGenotype = (genotype == "1/0" || genotype == "0/1" ||
                                 genotype == "1/1" || genotype == "0/0");

      if (!isPhasedGenotype && !isUnphasedGenotype) {
        isFirstLine = false;
        continue;
      }
      isFirstLine = false;
    }

    uniqueSamples.insert(sample);

    if (!isValidVariantFormat(variant)) {
      invalidFormatGenoVariants++;
      if (verbose && invalidFormatGenoVariants <= 10) {
        std::cerr << "Warning: Line " << genoLineCount
                  << " - Invalid variant format: " << variant
                  << ". Expected format: chr:pos:ref:alt" << std::endl;
      }
    }

    bool isUnphasedGenotype = (genotype.find('/') != std::string::npos);
    if (!unphasedMode && isUnphasedGenotype) {
      invalidGenotypeFormat++;
      if (invalidGenotypeFormat <= 10) {
        std::cerr << "Warning: Line " << genoLineCount
                  << " - Unphased genotype '" << genotype
                  << "' detected. Use --unphased flag to process unphased "
                     "data. Skipping."
                  << std::endl;
      } else if (invalidGenotypeFormat == 11) {
        std::cerr << "Warning: Additional unphased genotypes found. "
                     "Suppressing further warnings."
                  << std::endl;
      }
      skippedGenoLines++;
      continue;
    }

    bool validGenotype = false;
    std::string normalizedGenotype = genotype;

    if (unphasedMode) {
      if (genotype == "0/1" || genotype == "1/0" || genotype == "0|1" ||
          genotype == "1|0") {
        validGenotype = true;
        normalizedGenotype = "0/1";
      } else if (genotype == "1/1" || genotype == "1|1") {
        validGenotype = true;
        normalizedGenotype = "1/1";
      }
    } else {
      validGenotype =
          (genotype == "1|0" || genotype == "0|1" || genotype == "1|1");
      normalizedGenotype = genotype;
    }

    if (!validGenotype) {
      invalidGenotypeFormat++;
      if (invalidGenotypeFormat <= 10) {
        std::string expectedFormat =
            unphasedMode ? "0/1, 1/0, 0|1, 1|0, or 1/1" : "1|0, 0|1, or 1|1";
        std::cerr << "Warning: Line " << genoLineCount
                  << " - Skipping unexpected genotype value (" << genotype
                  << ") in sample " << sample << " for variant " << variant
                  << ". Expected values: " << expectedFormat << "."
                  << std::endl;
      } else if (invalidGenotypeFormat == 11) {
        std::cerr << "Warning: Additional invalid genotype formats found. "
                     "Suppressing further warnings."
                  << std::endl;
      }
      skippedGenoLines++;
      continue;
    }

    validGenoLines++;

    if (variantToGene.find(variant) != variantToGene.end()) {
      if (variantToGene[variant].size() > 1) {
        multiGeneVariants.insert(variant);
      }

      for (const auto &geneInfoPair : variantToGene[variant]) {
        std::string gene = geneInfoPair;
        uniqueVariantsKept.insert(variant);

        if (unphasedMode) {
          if (normalizedGenotype == "0/1") {
            sampleGeneHaplotypeVariant[sample][gene][0].push_back(
                variant); // Mark as het
          } else if (normalizedGenotype == "1/1") {
            sampleGeneHaplotypeVariant[sample][gene][1].push_back(
                variant); // Mark as hom
          }
        } else {
          if (normalizedGenotype == "1|0") {
            sampleGeneHaplotypeVariant[sample][gene][1].push_back(
                variant); // H1
          } else if (normalizedGenotype == "0|1") {
            sampleGeneHaplotypeVariant[sample][gene][2].push_back(
                variant); // H2
          } else if (normalizedGenotype == "1|1") {
            sampleGeneHaplotypeVariant[sample][gene][1].push_back(variant);
            sampleGeneHaplotypeVariant[sample][gene][2].push_back(variant);
          }
        }
      }
    } else {
      uniqueVariantsDiscarded.insert(variant);
    }
  }
  gzclose(genotypeFile);

  if (skippedGenoLines > 0) {
    std::cerr << "Warning: Skipped " << skippedGenoLines
              << " lines in genotype file due to format issues." << std::endl;
  }

  if (validGenoLines == 0) {
    std::cerr << "Error: No valid genotype data found in file: " << path
              << std::endl;
    return false;
  }

  if (uniqueVariantsKept.empty()) {
    std::cerr << "Error: No matching variants found between the --map and "
                 "--geno file!"
              << std::endl;
    return false;
  }

  if (verbose) {
    std::cerr << "Found " << uniqueSamples.size()
              << " unique samples in genotype file." << std::endl;
    std::cerr << "Variants in mapping file (kept): "
              << uniqueVariantsKept.size() << std::endl;
    std::cerr << "Variants not in mapping file (Discarded): "
              << uniqueVariantsDiscarded.size() << std::endl;
    std::cerr << "Variants mapping to more than one gene: "
              << multiGeneVariants.size() << std::endl;
  }

  return true;
}

std::vector<ChetResult> ChetCaller::getResults() const {
  std::vector<ChetResult> results;

  for (const auto &samplePair : sampleGeneHaplotypeVariant) {
    const std::string &sample = samplePair.first;
    for (const auto &genePair : samplePair.second) {
      const std::string &gene = genePair.first;
      const auto &haplotypeVariantMap = genePair.second;
      std::string chromosome;

      std::string callValue;
      int dosage = 0;
      std::set<std::string> haplotype1Variants, haplotype2Variants, allVariants;
      std::vector<std::string> unphasedVariantsList, h1VariantsList,
          h2VariantsList;

      if (unphasedMode) {
        std::set<std::string> hetVariants =
            haplotypeVariantMap.count(0)
                ? std::set<std::string>(haplotypeVariantMap.at(0).begin(),
                                        haplotypeVariantMap.at(0).end())
                : std::set<std::string>();
        std::set<std::string> homVariants =
            haplotypeVariantMap.count(1)
                ? std::set<std::string>(haplotypeVariantMap.at(1).begin(),
                                        haplotypeVariantMap.at(1).end())
                : std::set<std::string>();

        allVariants.insert(hetVariants.begin(), hetVariants.end());
        allVariants.insert(homVariants.begin(), homVariants.end());
        unphasedVariantsList.assign(allVariants.begin(), allVariants.end());

        if (allVariants.empty()) {
          if (verbose)
            std::cerr << "Warning: Empty variant sets for sample " << sample
                      << " and gene " << gene << ". Skipping." << std::endl;
          continue;
        }

        const std::string &variant = *allVariants.begin();
        std::stringstream ss(variant);
        std::getline(ss, chromosome, ':');

        if (!homVariants.empty()) {
          callValue = "hom";
          dosage = 2;
        } else if (!hetVariants.empty()) {
          callValue = "het";
          dosage = 1;
        } else {
          callValue = "na";
          dosage = 0;
        }

      } else { // Phased mode
        haplotype1Variants =
            haplotypeVariantMap.count(1)
                ? std::set<std::string>(haplotypeVariantMap.at(1).begin(),
                                        haplotypeVariantMap.at(1).end())
                : std::set<std::string>();
        haplotype2Variants =
            haplotypeVariantMap.count(2)
                ? std::set<std::string>(haplotypeVariantMap.at(2).begin(),
                                        haplotypeVariantMap.at(2).end())
                : std::set<std::string>();

        h1VariantsList.assign(haplotype1Variants.begin(),
                              haplotype1Variants.end());
        h2VariantsList.assign(haplotype2Variants.begin(),
                              haplotype2Variants.end());

        if (haplotype1Variants.empty() && haplotype2Variants.empty()) {
          if (verbose)
            std::cerr << "Warning: Empty variant sets for sample " << sample
                      << " and gene " << gene << ". Skipping." << std::endl;
          continue;
        }

        std::set<std::string> mergedVariants;
        std::set_union(haplotype1Variants.begin(), haplotype1Variants.end(),
                       haplotype2Variants.begin(), haplotype2Variants.end(),
                       std::inserter(mergedVariants, mergedVariants.begin()));

        if (mergedVariants.empty()) {
          std::cerr << "Error: No valid variants found for sample " << sample
                    << " and gene " << gene << std::endl;
          continue;
        }

        const std::string &variant = *mergedVariants.begin();
        std::stringstream ss(variant);
        std::getline(ss, chromosome, ':');

        if (!haplotype1Variants.empty() && haplotype2Variants.empty()) {
          callValue = (haplotype1Variants.size() == 1) ? "het" : "cis";
          dosage = 1;
        } else if (!haplotype2Variants.empty() && haplotype1Variants.empty()) {
          callValue = (haplotype2Variants.size() == 1) ? "het" : "cis";
          dosage = 1;
        } else if (!haplotype1Variants.empty() && !haplotype2Variants.empty()) {
          std::set<std::string> intersection;
          std::set_intersection(
              haplotype1Variants.begin(), haplotype1Variants.end(),
              haplotype2Variants.begin(), haplotype2Variants.end(),
              std::inserter(intersection, intersection.begin()));
          if (!intersection.empty()) {
            callValue = "hom";
            dosage = 2;
          } else {
            callValue = "chet";
            dosage = 2;
          }
        } else {
          callValue = "na";
          dosage = 0;
        }
      }

      ChetResult result;
      result.sample = sample;
      result.chromosome = chromosome;
      result.gene = gene;
      result.call = callValue;
      result.dosage = dosage;
      result.unphasedVariants = unphasedVariantsList;
      result.haplotype1Variants = h1VariantsList;
      result.haplotype2Variants = h2VariantsList;

      // Score calculation
      float geneScore = 0.0f;
      float haplotype1Score = 0.0f;
      float haplotype2Score = 0.0f;

      bool hasScores = !variantGeneScore.empty();

      if (hasScores && !unphasedMode) {
        if (defaultHaplotypeCollapseRule == "product") {
          haplotype1Score = 1.0f;
          haplotype2Score = 1.0f;
        } else if (defaultHaplotypeCollapseRule == "min") {
          if (!haplotype1Variants.empty() && haplotype2Variants.empty()) {
            haplotype1Score = 1.0f;
            haplotype2Score = 0.0f;
          } else if (!haplotype2Variants.empty() &&
                     haplotype1Variants.empty()) {
            haplotype1Score = 0.0f;
            haplotype2Score = 1.0f;
          } else {
            haplotype1Score = 1.0f;
            haplotype2Score = 1.0f;
          }
        } else {
          haplotype1Score = 0.0f;
          haplotype2Score = 0.0f;
        }

        auto processHaplotype = [&](const std::set<std::string> &vars,
                                    float &hScore) {
          for (const auto &variant : vars) {
            auto it = variantGeneScore.find(std::make_pair(variant, gene));
            float mappedScore =
                (it != variantGeneScore.end()) ? it->second : 1.0f;
            if (defaultHaplotypeCollapseRule == "product")
              hScore *= (1 - mappedScore);
            else if (defaultHaplotypeCollapseRule == "max")
              hScore = std::max(hScore, mappedScore);
            else if (defaultHaplotypeCollapseRule == "min")
              hScore = std::min(hScore, mappedScore);
            else if (defaultHaplotypeCollapseRule == "additive")
              hScore += mappedScore;
          }
        };

        processHaplotype(haplotype1Variants, haplotype1Score);
        processHaplotype(haplotype2Variants, haplotype2Score);

        if (defaultHaplotypeCollapseRule == "product") {
          haplotype2Score = 1 - haplotype2Score;
          haplotype1Score = 1 - haplotype1Score;
        }

        if (defaultGeneCollapseRule == "product")
          geneScore = haplotype1Score * haplotype2Score;
        else if (defaultGeneCollapseRule == "max")
          geneScore = std::max(haplotype1Score, haplotype2Score);
        else if (defaultGeneCollapseRule == "min")
          geneScore = std::min(haplotype1Score, haplotype2Score);
        else if (defaultGeneCollapseRule == "additive")
          geneScore = haplotype1Score + haplotype2Score;
      }

      result.geneScore = geneScore;
      result.haplotype1Score = haplotype1Score;
      result.haplotype2Score = haplotype2Score;

      results.push_back(result);
    }
  }
  return results;
}

void ChetCaller::printResults() const {
  std::vector<ChetResult> results = getResults();

  if (results.empty()) {
    std::cerr << "Warning: No results were generated. Check your input files "
                 "and parameters."
              << std::endl;
    return;
  }

  bool hasScores = !variantGeneScore.empty();

  for (const auto &res : results) {
    std::cout << res.sample << "\t" << res.chromosome << "\t" << res.gene
              << "\t" << res.call << "\t" << res.dosage;

    if (hasScores && !unphasedMode) {
      std::cout << "\tg=" << defaultGeneCollapseRule << "\t" << res.geneScore;
      if (showHaplotypeScore) {
        std::cout << "\th=" << defaultHaplotypeCollapseRule << "\t"
                  << res.haplotype1Score << "\t" << res.haplotype2Score;
      }
    }

    // Printing variants
    if (showVariants || !infoMap.empty()) {
      std::cout << "\t";

      auto printVars = [&](const std::vector<std::string> &vars, bool useInfo) {
        for (size_t i = 0; i < vars.size(); ++i) {
          const auto &variant = vars[i];
          std::cout << variant;
          if (useInfo) {
            std::pair<std::string, std::string> key =
                std::make_pair(variant, res.gene);
            auto infoIt = infoMap.find(key);
            std::string info =
                (infoIt != infoMap.end()) ? infoIt->second : "NA";
            std::cout << ":" << info;
          }
          if (i < vars.size() - 1)
            std::cout << ";";
        }
      };

      bool useInfo = !infoMap.empty();

      if (unphasedMode) {
        printVars(res.unphasedVariants, useInfo);
      } else {
        if (!res.haplotype1Variants.empty()) {
          printVars(res.haplotype1Variants, useInfo);
        }

        if (!res.haplotype1Variants.empty() &&
            !res.haplotype2Variants.empty()) {
          std::cout << "|";
        }

        if (!res.haplotype2Variants.empty()) {
          printVars(res.haplotype2Variants, useInfo);
        }
      }
    }
    std::cout << std::endl;
  }

  if (verbose) {
    std::cerr << "Successfully generated " << results.size()
              << " result entries." << std::endl;
  }
}

} // namespace call_chets
