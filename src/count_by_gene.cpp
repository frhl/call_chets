#include "logging.hpp"
#include "version.hpp"
#include <ctime>
#include <fstream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

struct IndividualGenotype {
  bool hasHomAlt = false;
  bool hasHet = false;
  bool hasHomRef = false;
};

using GeneToIndivMap = std::map<std::string, std::map<int, IndividualGenotype>>;

void printUsage(const char *path) {
  // Get current date and time
  std::time_t now = std::time(nullptr);
  char timestr[100];
  std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S",
                std::localtime(&now));

  std::cerr << "\n[COUNT_BY_GENE] Count genotypes by gene"
            << "\n  * Version       : " << getFullVersion()
            << "\n  * Run date      : " << timestr << "\n";

  std::cerr << "\nUsage: " << path
            << " --input <input.vcf.gz> --map <variant_to_gene.map.gz> "
               "[--max-af <max allele frequency>] [--output <output.txt>]\n";
  std::cerr << "Options:\n";
  std::cerr << "  --input/-i                : Input VCF file.\n";
  std::cerr << "  --map/-m                  : Gzipped variant to gene mapping "
               "file.\n";
  std::cerr << "  --max-af                  : Maximum allele frequency "
               "threshold (default includes all).\n";
  std::cerr << "  --max-maf                 : Maximum minor allele frequency "
               "threshold (default includes all).\n";
}

int main(int argc, char *argv[]) {
  std::string pathInput, outputPath = "gene_counts.txt", mappingFilePath;
  float maxAf = std::numeric_limits<float>::max();
  float maxMaf = std::numeric_limits<float>::max();
  int totalSamples = 0;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printUsage(argv[0]);
      return 0;
    } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
      pathInput = argv[++i];
    } else if ((arg == "--map" || arg == "-m") && i + 1 < argc) {
      mappingFilePath = argv[++i];
    } else if (arg == "--max-af" && i + 1 < argc) {
      maxAf = std::stof(argv[++i]);
    } else if (arg == "--max-maf" && i + 1 < argc) {
      maxMaf = std::stof(argv[++i]);
    } else {
      std::cerr << "Error! Unknown or incomplete argument: " << arg
                << std::endl;
      printUsage(argv[0]);
      return 1;
    }
  }

  if (pathInput.empty() || mappingFilePath.empty()) {
    std::cerr << "Error! --input and --map must be provided." << std::endl;
    printUsage(argv[0]);
    return 1;
  }

  // Load mapping file
  gzFile mappingFile = gzopen(mappingFilePath.c_str(), "r");
  if (!mappingFile) {
    std::cerr << "Error: Cannot open mapping file for reading: "
              << mappingFilePath << std::endl;
    return 1;
  }

  char buf[1024];
  std::string variant, gene;
  std::map<std::string, std::vector<std::string>> variantToGene;
  bool isFirstLineMappingFile = true;
  while (gzgets(mappingFile, buf, sizeof(buf))) {
    std::stringstream ss(buf);
    ss >> variant >> gene;

    if (ss.fail() || variant.empty() || gene.empty()) {
      if (!isFirstLineMappingFile) {
        std::cerr
            << "Error: Failed to extract two columns from line (variant gene): "
            << buf << ". Please, fix this line in --gene-map and retry."
            << std::endl;
        gzclose(mappingFile);
        return 1;
      }
    } else {
      variantToGene[variant].push_back(gene);
    }
    isFirstLineMappingFile = false;
  }
  gzclose(mappingFile);

  // Open the VCF file
  htsFile *vcfFile = bcf_open(pathInput.c_str(), "r");
  if (!vcfFile) {
    std::cerr << "Error: Cannot open VCF file for reading: " << pathInput
              << std::endl;
    return 1;
  }

  bcf_hdr_t *header = bcf_hdr_read(vcfFile);
  if (!header) {
    std::cerr << "Error: Cannot read header from VCF file." << std::endl;
    bcf_close(vcfFile);
    return 1;
  }

  bcf1_t *record = bcf_init();
  if (!record) {
    std::cerr << "Error: Cannot initialize BCF record." << std::endl;
    bcf_hdr_destroy(header);
    bcf_close(vcfFile);
    return 1;
  }

  totalSamples = bcf_hdr_nsamples(header);
  GeneToIndivMap geneToIndividuals;
  int ngt, *gt_arr = nullptr, ngt_arr = 0;

  // iterate over all variants
  while (bcf_read(vcfFile, header, record) == 0) {
    bcf_unpack(record, BCF_UN_ALL);
    std::string chrom = bcf_hdr_id2name(header, record->rid);
    int pos = record->pos + 1;
    std::string ref = record->d.allele[0];
    std::string alt = record->n_allele > 1 ? record->d.allele[1] : ".";
    std::string variantKey =
        chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;

    // Check if the current variant is in mapping
    if (variantToGene.find(variantKey) != variantToGene.end()) {
      ngt = bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
      if (ngt <= 0)
        continue; // Skip if no genotype information

      int n_samples = bcf_hdr_nsamples(header);
      int alleleCount = 0;

      for (int i = 0; i < n_samples; ++i) {
        if (gt_arr[i * 2] == bcf_gt_missing ||
            gt_arr[i * 2 + 1] == bcf_gt_missing)
          continue;

        int allele1 = bcf_gt_allele(gt_arr[i * 2]);
        int allele2 = bcf_gt_allele(gt_arr[i * 2 + 1]);

        // Count allele occurrences
        alleleCount += (allele1 > 0) + (allele2 > 0);
      }

      // Calculate allele frequency
      float alleleFrequency = static_cast<float>(alleleCount) / (2 * n_samples);
      float minorAlleleFrequency =
          std::min(alleleFrequency, 1.0f - alleleFrequency);

      // Skip variant if allele frequency exceeds maxAf
      if (alleleFrequency > maxAf || minorAlleleFrequency > maxMaf)
        continue;

      for (int i = 0; i < totalSamples; ++i) {
        if (gt_arr[i * 2] == bcf_gt_missing ||
            gt_arr[i * 2 + 1] == bcf_gt_missing)
          continue;

        int allele1 = bcf_gt_allele(gt_arr[i * 2]);
        int allele2 = bcf_gt_allele(gt_arr[i * 2 + 1]);
        bool isHomRef = allele1 == 0 && allele2 == 0;
        bool isHet = allele1 != allele2;
        bool isHomAlt = allele1 == allele2 && allele1 != 0;

        for (const auto &gene : variantToGene[variantKey]) {
          auto &indiv = geneToIndividuals[gene][i];
          if (isHomAlt)
            indiv.hasHomAlt = true;
          else if (isHet)
            indiv.hasHet = true;
          else if (isHomRef)
            indiv.hasHomRef = !indiv.hasHomAlt && !indiv.hasHet;
        }
      }
    }
  }

  // Free resources
  if (gt_arr)
    free(gt_arr);
  bcf_destroy(record);
  bcf_hdr_destroy(header);
  bcf_close(vcfFile);

  std::cout << "region\taa\tAa\tAA\tsamples\n";
  // Output preparation and writing
  for (const auto &genePair : geneToIndividuals) {
    const auto &gene = genePair.first;
    const auto &individuals = genePair.second;

    int homRefCount = 0, hetCount = 0, homAltCount = 0;
    for (const auto &indivPair : individuals) {
      if (indivPair.second.hasHomAlt)
        homAltCount++;
      else if (indivPair.second.hasHet)
        hetCount++;
      else if (indivPair.second.hasHomRef)
        homRefCount++;
    }
    std::cout << gene << "\t" << homAltCount << "\t" << hetCount << "\t"
              << homRefCount << "\t" << totalSamples << "\n";
  }

  return 0;
}
