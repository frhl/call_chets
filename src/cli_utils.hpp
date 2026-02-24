#ifndef CLI_UTILS_HPP
#define CLI_UTILS_HPP

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace call_chets {

/**
 * @brief Sort chromosome names in natural order.
 *
 * Sorts chromosome names such that numeric chromosomes come first (1, 2, ..., 22),
 * followed by non-numeric chromosomes (X, Y, MT) in lexicographical order.
 * Handles both "chr" prefixed and non-prefixed chromosome names.
 *
 * @param contigs Set of chromosome/contig names to sort.
 * @return std::vector<std::string> Sorted vector of chromosome names.
 */
inline std::vector<std::string>
sortChromosomes(const std::set<std::string> &contigs) {
  std::vector<std::string> chromosomes(contigs.begin(), contigs.end());
  std::sort(chromosomes.begin(), chromosomes.end(),
            [](const std::string &a, const std::string &b) {
              // Extract the part after 'chr' prefix
              std::string a_num = a.substr(0, 3) == "chr" ? a.substr(3) : a;
              std::string b_num = b.substr(0, 3) == "chr" ? b.substr(3) : b;

              // If both are numeric chromosomes
              if (!a_num.empty() && !b_num.empty() && std::isdigit(a_num[0]) &&
                  std::isdigit(b_num[0])) {
                return std::stoi(a_num) < std::stoi(b_num);
              }
              // If a is numeric but b is not (e.g., a = chr2, b = chrX)
              if (!a_num.empty() && std::isdigit(a_num[0]) &&
                  (b_num.empty() || !std::isdigit(b_num[0]))) {
                return true;
              }
              // If b is numeric but a is not
              if (!b_num.empty() && std::isdigit(b_num[0]) &&
                  (a_num.empty() || !std::isdigit(a_num[0]))) {
                return false;
              }
              // If neither are numeric (e.g., chrX, chrY, etc.), then just use
              // lexicographical order
              return a < b;
            });
  return chromosomes;
}

/**
 * @brief Read a two-column mapping file (variant/gene to group).
 *
 * Reads a tab-separated file with a header line. The first column is the key
 * (variant or gene ID), the second column is the group/value.
 *
 * @param groupMapPath Path to the mapping file.
 * @param errorMessage Output parameter for error message if reading fails.
 * @return std::map<std::string, std::string> Map from key to group.
 *         Returns empty map on error (check errorMessage).
 */
inline std::map<std::string, std::string>
readGroupMap(const std::string &groupMapPath, std::string &errorMessage) {
  std::map<std::string, std::string> groupMap;
  std::ifstream infile(groupMapPath);
  if (!infile) {
    errorMessage =
        "Error: Cannot open group map file for reading: " + groupMapPath;
    return groupMap;
  }

  std::string line, key, group;
  std::getline(infile, line); // skip header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (!(iss >> key >> group)) {
      break;
    }
    groupMap[key] = group;
  }

  infile.close();
  errorMessage.clear();
  return groupMap;
}

/**
 * @brief Read a two-column mapping file (variant/gene to group).
 *
 * Overload that throws std::runtime_error on failure.
 *
 * @param groupMapPath Path to the mapping file.
 * @return std::map<std::string, std::string> Map from key to group.
 * @throws std::runtime_error if file cannot be opened.
 */
inline std::map<std::string, std::string>
readGroupMapOrThrow(const std::string &groupMapPath) {
  std::string errorMessage;
  auto result = readGroupMap(groupMapPath, errorMessage);
  if (!errorMessage.empty()) {
    throw std::runtime_error(errorMessage);
  }
  return result;
}

} // namespace call_chets

#endif // CLI_UTILS_HPP
