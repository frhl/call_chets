#include "logging.hpp"
#include "version.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

void printUsage(const char *path) {
  // Get current date and time
  std::time_t now = std::time(nullptr);
  char timestr[100];
  std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S",
                std::localtime(&now));

  std::cerr << "\n[ENCODE_VCF_BY_GROUP] VCF encoder for group-based analysis"
            << "\n  * Version       : " << getFullVersion()
            << "\n  * Run date      : " << timestr << "\n";

  std::cerr << "\nUsage: " << path << " --input <input> --samples <samples>"
            << " --mode [<additive|recessive|dominance|001|012|010|011>]"
            << " [--group-map <group-map>]\n";
  std::cerr << "\nDescription:";
  std::cerr
      << "\n  Converts 'call_chets' output to VCF for downstream analysis.";
  std::cerr << "\n  Results are streamed to standard output.\n";
  std::cerr << "\nOptions:";
  std::cerr << "\n  --input/-i       : Output from 'call_chets'.";
  std::cerr
      << "\n  --samples/-s     : List of samples. One per line. No header.";
  std::cerr << "\n  --mode/-m        : Specify genotype/dosage encoding. Use "
               "'additive' or '012'";
  std::cerr << "\n                     for dosages of 0, 1, and 2. Use "
               "'recessive' or '001' for";
  std::cerr << "\n                     dosages of 0 and 2. Use 'dominance' to "
               "encode orthogonal";
  std::cerr << "\n                     contribution for dominance effects. "
               "'010' and '011' represent";
  std::cerr << "\n                     custom modes setting bi-allelics to "
               "zero or one respectively.";
  std::cerr
      << "\n  --min-ac         : Filters to genes with sum of DS >= argument.";
  std::cerr
      << "\n  --max-ac         : Filters to genes with sum of DS < argument.";
  std::cerr << "\n  --scaling-factor : Apply a scaling factor to the dosages. "
               "Default is 1.0.";
  std::cerr
      << "\n  --suffix         : Suffix to append to gene names in the output.";
  std::cerr << "\n  --no-dosage-scaling : Disable dosage scaling (only "
               "relevant when mode='dominance').";
  std::cerr << "\n  --all-info       : Populate INFO column with all details "
               "(relevant when mode='dominance').";
  std::cerr << "\n  --group-map      : File mapping genes to groups. Two "
               "columns: gene, group.\n";
  std::cerr << "\nExample:";
  std::cerr << "\n  ./encode_vcf called_chets.txt.gz samples.txt additive | "
               "bgzip > out.vcf.gz\n\n";
}

std::vector<std::string> sortChromosomes(const std::set<std::string> &contigs) {
  std::vector<std::string> chromosomes(contigs.begin(), contigs.end());
  std::sort(chromosomes.begin(), chromosomes.end(),
            [](const std::string &a, const std::string &b) {
              std::string a_num = a.substr(0, 3) == "chr" ? a.substr(3) : a;
              std::string b_num = b.substr(0, 3) == "chr" ? b.substr(3) : b;
              if (isdigit(a_num[0]) && isdigit(b_num[0])) {
                return std::stoi(a_num) < std::stoi(b_num);
              }
              if (isdigit(a_num[0]) && !isdigit(b_num[0])) {
                return true;
              }
              if (!isdigit(a_num[0]) && isdigit(b_num[0])) {
                return false;
              }
              return a < b;
            });
  return chromosomes;
}

std::map<std::string, std::string>
readGroupMap(const std::string &groupMapPath) {
  std::map<std::string, std::string> groupMap;
  std::ifstream infile(groupMapPath);
  if (!infile) {
    std::cerr << "Error: Cannot open group map file for reading: "
              << groupMapPath << std::endl;
    exit(1);
  }
  std::string line, gene, group;
  std::getline(infile, line); // skip header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (!(iss >> gene >> group)) {
      break;
    }
    groupMap[gene] = group;
  }
  infile.close();
  return groupMap;
}

int main(int argc, char *argv[]) {

  std::string pathInput;
  std::string pathSamples;
  std::string mode = "additive";
  std::string suffix = "";
  std::string groupMapPath;
  int minAC = 0;
  int maxAC = INT_MAX;
  float scalingFactor = 1.0;
  bool scaleDosage = true;
  bool allInfo = false;

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
    } else if ((arg == "--group-map") && i + 1 < argc) {
      groupMapPath = argv[++i];
    } else {
      std::cerr << "Error! Unknown or incomplete argument: " << arg
                << std::endl;
      printUsage(argv[0]);
      return 1;
    }
  }

  if (pathInput.empty() || pathSamples.empty()) {
    std::cerr << "Error! Both --input and --samples must be provided."
              << std::endl;
    printUsage(argv[0]);
    return 1;
  }

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
    return 1;
  }

  // Normalise mode input
  if (mode == "recessive")
    mode = "001";
  else if (mode == "additive")
    mode = "012";

  if (mode != "001" && mode != "012" && mode != "010" && mode != "011" &&
      mode != "dominance") {
    std::cerr
        << "Error: Invalid dosage encoding mode provided. Only '012|additive', "
           "'001|recessive', '010', '011', or 'dominance'."
        << std::endl;
    printUsage(argv[0]);
    return 1;
  }

  // Collect all samples into a set
  std::set<std::string> samples;
  std::string sampleLine;
  while (std::getline(sampleFile, sampleLine)) {
    sampleLine.erase(0, sampleLine.find_first_not_of(" \t\n\r"));
    sampleLine.erase(sampleLine.find_last_not_of(" \t\n\r") + 1);
    if (!sampleLine.empty()) {
      samples.insert(sampleLine);
    }
  }

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

  std::map<std::string, std::string> geneToGroup;
  std::map<std::string, std::pair<float, float>> groupDosages;
  std::set<std::string> contigs;

  if (!groupMapPath.empty()) {
    geneToGroup = readGroupMap(groupMapPath);
    for (const auto &pair : geneToGroup) {
      groupDosages[pair.second] = {std::numeric_limits<float>::max(),
                                   std::numeric_limits<float>::lowest()};
    }
  }

  std::string line, sample, chromosome, gene, configuration, variantInfo;
  char buffer[4096];
  float dosage;

  while (gzgets(longFile, buffer, sizeof(buffer))) {
    std::string line(buffer);
    std::istringstream iss(line);
    iss >> sample >> chromosome >> gene >> configuration >> dosage >>
        variantInfo;

    if (configuration != "chet" && configuration != "het" &&
        configuration != "cis" && configuration != "hom") {
      std::cerr << "Error: 4th column '" << configuration
                << "' is not one of the expected values (chet, het, cis, hom) "
                   "in line: "
                << line << std::endl;
      gzclose(longFile);
      printUsage(argv[0]);
      return 1;
    }

    if (samples.find(sample) == samples.end())
      continue;

    if (mode == "001" && dosage == 1.0f) {
      dosage = 0;
    }
    if (mode == "010" && dosage == 2.0f) {
      dosage = 0;
    }
    if (mode == "011" && dosage == 2.0f) {
      dosage = 1;
    }

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

    if (configuration == "chet" || configuration == "hom") {
      geneBI[gene] += 1;
    }

    geneToChromosome[gene] = chromosome;
    geneSampleDosage[gene][sample] = dosage;
    contigs.insert(chromosome);

    if (!groupMapPath.empty() && geneToGroup.find(gene) != geneToGroup.end()) {
      float r =
          static_cast<float>((geneAN / 2) -
                             (geneCis[gene] + geneHet[gene] + geneBI[gene])) /
          (geneAN / 2);
      float h =
          static_cast<float>(geneHet[gene] + geneCis[gene]) / (geneAN / 2);
      float a = static_cast<float>(geneBI[gene]) / (geneAN / 2);
      float dom_dosage_aa = -h * a;
      float dom_dosage_Aa = 2 * a * r;
      float dom_dosage_AA = -h * r;
      groupDosages[geneToGroup[gene]].first =
          std::min({groupDosages[geneToGroup[gene]].first, dom_dosage_aa,
                    dom_dosage_Aa, dom_dosage_AA});
      groupDosages[geneToGroup[gene]].second =
          std::max({groupDosages[geneToGroup[gene]].second, dom_dosage_aa,
                    dom_dosage_Aa, dom_dosage_AA});
    }
  }

  gzclose(longFile);

  std::vector<std::string> sortedContigs = sortChromosomes(contigs);

  std::cout << "##fileformat=VCFv4.2\n";
  std::cout << "##EncodingMode=" << mode << "\n";
  std::cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
  for (const auto &chr : sortedContigs) {
    std::cout << "##contig=<ID=" << chr << ">\n";
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
  std::cout << "##INFO=<ID=HET,Number=1,Type=Integer,Description=\""
               "Heterozygous Count\">\n";
  std::cout << "##INFO=<ID=CIS,Number=1,Type=Integer,Description=\"Cis pseudo "
               "Count \">\n";
  if (mode == "dominance") {
    std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic references (aa)\">\n";
    std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of "
                 "heterozygotes (Aa)\">\n";
    std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic alternates (AA)\">\n";
    if (allInfo) {
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

  int discardedGenesCount = 0;
  const int warningLimit = 5;
  for (const auto &chr : sortedContigs) {
    int rowIndex = 0;
    for (const auto &genePair : geneSampleDosage) {
      if (geneToChromosome[genePair.first] == chr) {
        rowIndex++;
        if (!groupMapPath.empty() &&
            geneToGroup.find(genePair.first) == geneToGroup.end()) {
          if (discardedGenesCount < warningLimit) {
            std::cerr << "Warning: Gene '" << genePair.first
                      << "' not found in group map. Discarding.\n";
          }
          discardedGenesCount++;
          continue;
        }

        int currentAN = geneAN;
        int currentAC = geneAC[genePair.first];
        int currentBI = geneBI[genePair.first];
        int currentChet = geneChet[genePair.first];
        int currentHom = geneHom[genePair.first];
        int currentHet = geneHet[genePair.first];
        int currentCis = geneCis[genePair.first];

        if (mode == "dominance" && currentBI == 0) {
          continue;
        }

        float aa_count = static_cast<float>(
            (currentAN / 2) - (currentCis + currentHet + currentBI));
        float Aa_count = static_cast<float>(currentHet + currentCis);
        float AA_count = static_cast<float>(currentBI);

        float r = static_cast<float>(aa_count / (currentAN / 2));
        float h = static_cast<float>(Aa_count / (currentAN / 2));
        float a = static_cast<float>(AA_count / (currentAN / 2));
        float totalFrequency = r + h + a;

        float dom_dosage_aa = -h * a;
        float dom_dosage_Aa = 2 * a * r;
        float dom_dosage_AA = -h * r;

        float minDomDosage =
            std::min({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
        float maxDomDosage =
            std::max({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

        if (currentAC >= minAC && currentAC < maxAC) {
          std::cout << geneToChromosome[genePair.first] << "\t" << rowIndex
                    << "\t" << genePair.first + suffix << "\tA\tB\t.\t.\t"
                    << "AC=" << currentAC << ";AN=" << currentAN
                    << ";BI=" << currentBI << ";CHET=" << currentChet
                    << ";HOM=" << currentHom << ";CIS=" << currentCis;

          if (mode == "dominance") {
            std::cout << ";r=" << r << ";h=" << h << ";a=" << a;
            if (!groupMapPath.empty()) {
              std::string group = geneToGroup[genePair.first];
              minDomDosage = groupDosages[group].first;
              maxDomDosage = groupDosages[group].second;
              std::cout << ";minDosage=" << minDomDosage
                        << ";maxDosage=" << maxDomDosage << ";group=" << group;
            } else if (allInfo) {
              std::cout << ";minDosage=" << minDomDosage
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
          }
          std::cout << "\tDS";
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

  if (discardedGenesCount > warningLimit) {
    std::cerr << "Note: Total number of genes discarded due to not being "
                 "present in the group map: "
              << discardedGenesCount << std::endl;
  }

  return 0;
}
