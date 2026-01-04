#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "logging.hpp"
#include "version.hpp"

// Struct for variant encoding preview
struct VariantEncodingExample {
  std::string variantId;
  float r, h, a;
  float raw_aa, raw_Aa, raw_AA;
  float scaled_aa, scaled_Aa, scaled_AA;
  float variantMin, variantMax;
};

// Print variant encoding preview to stderr
void printVariantEncodingPreview(
    const std::vector<VariantEncodingExample> &examples,
    const std::string &scalingMode, bool hasMore) {
  if (examples.empty())
    return;

  std::cerr << "\n  * Encoding Preview (first " << examples.size()
            << " variants):" << std::endl;
  std::cerr << "      Variant          | r      | h      | a      | Raw: aa   "
               "| Raw: Aa   | Raw: AA   | Scaled: aa | Scaled: Aa | Scaled: AA"
            << std::endl;
  std::cerr << "      "
               "-----------------|--------|--------|--------|-----------|------"
               "-----|-----------|------------|------------|------------"
            << std::endl;

  for (const auto &ex : examples) {
    // Truncate variant ID to 17 chars for display
    std::string displayId =
        ex.variantId.length() > 17 ? ex.variantId.substr(0, 17) : ex.variantId;

    std::cerr << "      " << std::left << std::setw(17) << displayId << "| "
              << std::fixed << std::setprecision(4) << std::setw(6) << ex.r
              << " | " << std::setw(6) << ex.h << " | " << std::setw(6) << ex.a
              << " | " << std::setprecision(6) << std::setw(9) << ex.raw_aa
              << " | " << std::setw(9) << ex.raw_Aa << " | " << std::setw(9)
              << ex.raw_AA << " | " << std::setw(10) << ex.scaled_aa << " | "
              << std::setw(10) << ex.scaled_Aa << " | " << ex.scaled_AA
              << std::endl;
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

  std::cerr
      << "\n[ORTHOGONALIZE] Transform genotypes/dosages to dominance encoding"
      << "\n  * Version       : " << version
      << "\n  * Run date      : " << timestr << "\n";

  std::cerr << "\nUsage: " << path << " --input <input.vcf.gz> [options]\n";
  std::cerr << "\nRequired Options:\n";
  std::cerr << "  --input/-i <file>          Input VCF/BCF file (supports "
               ".vcf, .vcf.gz, .bcf)\n";
  std::cerr << "\nMain Options:\n";
  std::cerr
      << "  --mode/-m <mode>           Processing mode (default: dominance)\n";
  std::cerr << "                             Supported modes:\n";
  std::cerr << "                               - dominance: Dominance "
               "deviation encoding (alias: nonadditive)\n";
  std::cerr << "                               - recessive: Set heterozygotes "
               "to 0, homozygotes unchanged\n";
  std::cerr << "\nScaling Options (mutually exclusive):\n";
  std::cerr << "  --scale-per-variant        Scale each variant independently "
               "to [0,2] (single pass)\n";
  std::cerr << "                             Use for: Individual variant "
               "testing (faster)\n";
  std::cerr << "                             Note: Betas/SEs not comparable "
               "across variants\n";
  std::cerr << "  --scale-globally           Scale using global min/max across "
               "all variants (two passes)\n";
  std::cerr << "                             Use for: Group-based testing with "
               "comparable betas/SEs\n";
  std::cerr << "  --scale-by-group <file>    Scale using group-specific "
               "min/max from mapping file (two passes)\n";
  std::cerr << "                             Format: tab-separated with "
               "headers 'variant' and 'gene'\n";
  std::cerr << "                             Use for: Gene-based testing with "
               "within-gene comparability\n";
  std::cerr << "  --scale-factor <f>         Apply additional scaling factor "
               "(default: 1.0)\n";
  std::cerr << "                             Can be combined with any scaling "
               "option above\n";
  std::cerr << "\n  --min-hom-count <n>        Minimum number of homozygous "
               "alternate alleles\n                             required "
               "(default: 1)\n\nAdditional Options:\n";
  std::cerr << "  --set-variant-id           Set variant IDs to "
               "chr:pos:ref:alt format\n";
  std::cerr << "  --all-info                 Include additional info in output "
               "(frequencies, min/max)\n";
  std::cerr << "\nInput Requirements:\n";
  std::cerr << "  - Input file must contain either GT (genotype) or DS "
               "(dosage) format fields\n";
  std::cerr << "  - DS values should be in range [0,2] where:\n";
  std::cerr << "    0 = homozygous reference (aa)\n";
  std::cerr << "    1 = heterozygous (Aa)\n";
  std::cerr << "    2 = homozygous alternate (AA)\n";
  std::cerr << "\nExample Usage:\n";
  std::cerr << "  Per-variant scaling (default, fastest):\n";
  std::cerr << "    " << path << " --input sample.vcf.gz --scale-per-variant\n";
  std::cerr << "\n  Global scaling across all variants:\n";
  std::cerr << "    " << path << " --input sample.vcf.gz --scale-globally\n";
  std::cerr << "\n  Group-based scaling:\n";
  std::cerr << "    " << path
            << " --input sample.vcf.gz --scale-by-group genes.txt\n";
  std::cerr << "\n  Full output with variant IDs:\n";
  std::cerr << "    " << path
            << " --input sample.vcf.gz --scale-per-variant --all-info "
               "--set-variant-id\n";
  std::cerr << "\n  Using recessive mode:\n";
  std::cerr << "    " << path << " --input sample.vcf.gz --mode recessive\n";
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

bool parseArguments(int argc, char *argv[], std::string &pathInput,
                    std::string &mode, double &scalingFactor,
                    bool &scalePerVariant, bool &scaleGlobally,
                    bool &setVariantId, bool &allInfo,
                    std::string &scaleByGroupPath, int &minHomCount) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printUsage(argv[0]);
      exit(0);
    } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
      pathInput = argv[++i];
    } else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
      mode = argv[++i];
      if (mode == "nonadditive")
        mode = "dominance";
    } else if (arg == "--min-hom-count" && i + 1 < argc) {
      minHomCount = std::stoi(argv[++i]);
    } else if (arg == "--scale-factor" && i + 1 < argc) {
      scalingFactor = std::stod(argv[++i]);
    } else if (arg == "--scale-per-variant") {
      scalePerVariant = true;
    } else if (arg == "--scale-globally") {
      scaleGlobally = true;
    } else if (arg == "--set-variant-id") {
      setVariantId = true;
    } else if (arg == "--all-info") {
      allInfo = true;
    } else if (arg == "--scale-by-group" && i + 1 < argc) {
      scaleByGroupPath = argv[++i];
    } else {
      std::cerr << "Error! Unknown or incomplete argument: " << arg
                << std::endl;
      printUsage(argv[0]);
      return false;
    }
  }

  if (pathInput.empty()) {
    std::cerr << "Error: --input is required\n";
    printUsage(argv[0]);
    return 1;
  }

  std::map<std::string, std::string> files;
  files["Input"] = pathInput;
  if (!scaleByGroupPath.empty())
    files["Group File"] = scaleByGroupPath;

  std::map<std::string, std::string> params;
  params["Mode"] = mode;
  params["Min Hom Count"] = std::to_string(minHomCount);
  params["Scale Variant"] = scalePerVariant ? "Yes" : "No";
  params["Scale Global"] = scaleGlobally ? "Yes" : "No";
  params["Scale Factor"] = std::to_string(scalingFactor);
  params["All Info"] = allInfo ? "Yes" : "No";

  call_chets::printHeader("RECODE",
                          "Transform genotypes/dosages to dominance encoding",
                          files, params);

  // Check for mutually exclusive scaling options
  int scalingModes = (scalePerVariant ? 1 : 0) + (scaleGlobally ? 1 : 0) +
                     (!scaleByGroupPath.empty() ? 1 : 0);
  if (scalingModes > 1) {
    std::cerr << "Error! --scale-per-variant, --scale-globally, and "
                 "--scale-by-group are mutually exclusive."
              << std::endl;
    printUsage(argv[0]);
    return false;
  }

  return true;
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

  std::string line, variant, group;
  std::getline(infile, line); // skip header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (!(iss >> variant >> group)) {
      break;
    }
    groupMap[variant] = group;
  }

  infile.close();
  return groupMap;
}

bool hasFormat(const bcf_hdr_t *hdr, const char *format) {
  return bcf_hdr_id2int(hdr, BCF_DT_ID, format) >= 0;
}

void processSampleGenotypes(int *gt_ptr, float ds_value, bool has_gt,
                            bool has_ds, int &aa_count, int &Aa_count,
                            int &AA_count, int &non_missing_count,
                            long &roundedDosageCount) {
  if (has_gt &&
      (!bcf_gt_is_missing(gt_ptr[0]) && !bcf_gt_is_missing(gt_ptr[1]))) {
    // Use GT if available
    int allele1 = bcf_gt_allele(gt_ptr[0]);
    int allele2 = bcf_gt_allele(gt_ptr[1]);
    non_missing_count++;
    if (allele1 == 0 && allele2 == 0)
      aa_count++;
    else if (allele1 != allele2)
      Aa_count++;
    else if (allele1 > 0 && allele2 > 0)
      AA_count++;
  } else if (has_ds && !std::isnan(ds_value)) {
    // Round DS value to nearest integer
    int rounded_ds = std::round(ds_value);
    if (std::abs(ds_value - rounded_ds) > 1e-4) {
      roundedDosageCount++;
    }
    if (rounded_ds < 0 || rounded_ds > 2) {
      std::cerr << "Warning: DS value " << ds_value
                << " out of range [0,2]. Skipping.\n";
      return;
    }
    non_missing_count++;
    if (rounded_ds == 0)
      aa_count++;
    else if (rounded_ds == 1)
      Aa_count++;
    else if (rounded_ds == 2)
      AA_count++;
  }
}

// Calculate global and group-specific dosages
void calculateGlobalAndGroupDosages(
    htsFile *fp, bcf_hdr_t *hdr, double &globalMinDomDosage,
    double &globalMaxDomDosage, int n_samples,
    std::set<std::string> &chromosomes,
    const std::map<std::string, std::string> &groupMap,
    std::map<std::string, std::pair<double, double>> &groupDosages,
    bool &hasMissingValues) {

  bcf1_t *rec = bcf_init();
  int *gt_arr = NULL, ngt_arr = 0;
  float *ds_arr = NULL;
  int nds_arr = 0;
  long roundedDosageCount = 0; // Local counter for this pass
  long processedVariants = 0;

  bool has_gt_format = hasFormat(hdr, "GT");
  bool has_ds_format = hasFormat(hdr, "DS");

  auto pass1_start_time = std::chrono::steady_clock::now();

  while (bcf_read(fp, hdr, rec) == 0) {
    processedVariants++;
    if (processedVariants % 10000 == 0) {
      auto now = std::chrono::steady_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
          now - pass1_start_time);
      int rate = elapsed.count() > 0 ? processedVariants / elapsed.count() : 0;
      std::cerr << "\r    [Pass 1] " << processedVariants
                << " variants processed (" << elapsed.count() << "s elapsed, ~"
                << rate << " var/sec)" << std::flush;
    }
    chromosomes.insert(bcf_hdr_id2name(hdr, rec->rid));
    bcf_unpack(rec, BCF_UN_STR);

    int ngt =
        has_gt_format ? bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) : 0;
    int nds = has_ds_format
                  ? bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds_arr)
                  : 0;

    // Per-variant flags for successful retrieval
    bool has_gt = has_gt_format && ngt >= 0;
    bool has_ds = has_ds_format && nds >= 0;

    int aa_count = 0, Aa_count = 0, AA_count = 0, non_missing_count = 0;

    for (int i = 0; i < n_samples; ++i) {
      int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
      float ds_value = has_ds ? ds_arr[i] : 0.0f;

      processSampleGenotypes(gt_ptr, ds_value, has_gt, has_ds, aa_count,
                             Aa_count, AA_count, non_missing_count,
                             roundedDosageCount);
    }

    if (non_missing_count < n_samples) {
      hasMissingValues = true;
    }

    if (AA_count > 0 && non_missing_count > 0) {
      double r = static_cast<double>(aa_count) / non_missing_count;
      double h = static_cast<double>(Aa_count) / non_missing_count;
      double a = static_cast<double>(AA_count) / non_missing_count;

      double dom_dosage_aa = -h * a;
      double dom_dosage_Aa = 2 * a * r;
      double dom_dosage_AA = -h * r;

      globalMinDomDosage = std::min(
          {globalMinDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
      globalMaxDomDosage = std::max(
          {globalMaxDomDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});

      std::string variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) +
                              ":" + std::to_string(rec->pos + 1) + ":" +
                              rec->d.allele[0] + ":" + rec->d.allele[1];
      if (groupMap.find(variantId) != groupMap.end()) {
        std::string group = groupMap.at(variantId);
        double &groupMinDosage = groupDosages[group].first;
        double &groupMaxDosage = groupDosages[group].second;
        groupMinDosage = std::min(
            {groupMinDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
        groupMaxDosage = std::max(
            {groupMaxDosage, dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
      }
    }
  }
  std::cerr << std::endl;

  free(gt_arr);
  free(ds_arr);
  bcf_destroy(rec);
}

void printHeader(const bcf_hdr_t *hdr,
                 const std::vector<std::string> &sortedContigs,
                 const std::string &mode, double globalMinDomDosage,
                 double globalMaxDomDosage, bool allInfo, bool scalePerVariant,
                 bool scaleGlobally, bool scaleByGroup) {
  int n_samples = bcf_hdr_nsamples(hdr);

  std::cout << "##fileformat=VCFv4.2\n";
  std::cout << "##EncodingMode=" << mode << "\n";
  // output chromosomes (only if provided)
  for (const auto &chr : sortedContigs) {
    std::cout << "##contig=<ID=" << chr << ">\n";
  }
  std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count "
               "in genotypes\">\n";
  std::cout << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number "
               "of alleles in called genotypes\">\n";

  if (mode == "dominance") {
    std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of "
                 "the alternate allele based on dominance model\">\n";
  } else if (mode == "recessive") {
    std::cout
        << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the "
           "alternate allele with recessive encoding (het=0, hom=2)\">\n";
  }

  if (allInfo) {
    std::cout << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic references (aa)\">\n";
    std::cout << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of "
                 "heterozygotes (Aa)\">\n";
    std::cout << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of "
                 "bi-allelic alternates (AA)\">\n";
  }
  if (scalePerVariant) {
    std::cout << "##INFO=<ID=variantMinDomDosage,Number=1,Type=Float,"
                 "Description=\"Per-variant minimum dominance dosage\">\n";
    std::cout << "##INFO=<ID=variantMaxDomDosage,Number=1,Type=Float,"
                 "Description=\"Per-variant maximum dominance dosage\">\n";
  } else if (scaleGlobally) {
    std::cout << "##INFO=<ID=globalMinDomDosage,Number=1,Type=Float,"
                 "Description=\"Global minimum dominance dosage\">\n";
    std::cout << "##INFO=<ID=globalMaxDomDosage,Number=1,Type=Float,"
                 "Description=\"Global maximum dominance dosage\">\n";
  } else if (scaleByGroup) {
    std::cout << "##INFO=<ID=groupMinDomDosage,Number=1,Type=Float,Description="
                 "\"Group-specific minimum dominance dosage\">\n";
    std::cout << "##INFO=<ID=groupMaxDomDosage,Number=1,Type=Float,Description="
                 "\"Group-specific maximum dominance dosage\">\n";
    std::cout << "##INFO=<ID=group,Number=1,Type=String,Description=\"Group "
                 "used for scaling\">\n";
  }
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (int i = 0; i < n_samples; ++i) {
    std::cout << "\t" << bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
  }
  std::cout << "\n";
}

void processVcfFile(
    htsFile *fp, bcf_hdr_t *hdr, const std::string &mode,
    double globalMinDomDosage, double globalMaxDomDosage, bool allInfo,
    bool scalePerVariant, bool scaleGlobally, bool scaleByGroup,
    bool setVariantId, double scalingFactor,
    const std::map<std::string, std::string> &groupMap,
    const std::map<std::string, std::pair<double, double>> &groupDosages,
    std::set<std::string> &chromosomes, int minHomCount) {
  bcf1_t *rec = bcf_init();
  int *gt_arr = NULL, ngt_arr = 0;
  float *ds_arr = NULL;
  int nds_arr = 0;
  int n_samples = bcf_hdr_nsamples(hdr);
  const double epsilon = 1e-8;
  int countVariantsWithoutHomAlt = 0;
  const int limit = 5;
  int discardedVariantsCount = 0;
  const int warningLimit = 5;
  int lowVarianceVariantsCount = 0;
  long totalVariants = 0;
  long keptVariants = 0;
  const double lowVarianceThreshold =
      0.0001; // Warn if any two genotype dosages differ by less than this
  long roundedDosageCount = 0;
  long haploidWarningCount = 0;

  // Variables for output preview
  std::vector<VariantEncodingExample> encodingExamples;
  bool previewPrinted = false;
  const int MAX_EXAMPLES = 5;

  bool has_gt_format = hasFormat(hdr, "GT");
  bool has_ds_format = hasFormat(hdr, "DS");

  // Track processing time (pass2_start is already defined in main)
  auto processing_start_time = std::chrono::steady_clock::now();

  while (bcf_read(fp, hdr, rec) == 0) {
    totalVariants++;
    bcf_unpack(rec, BCF_UN_STR);

    // Collect chromosome names as we process
    chromosomes.insert(bcf_hdr_id2name(hdr, rec->rid));

    int ngt =
        has_gt_format ? bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) : 0;
    int nds = has_ds_format
                  ? bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds_arr)
                  : 0;

    // Per-variant flags for successful retrieval
    bool has_gt = has_gt_format && ngt >= 0;
    bool has_ds = has_ds_format && nds >= 0;

    if (!has_gt && !has_ds)
      continue;

    // Progress logging every 10000 variants
    if (keptVariants > 0 && keptVariants % 10000 == 0) {
      auto now = std::chrono::steady_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
          now - processing_start_time);
      int rate = elapsed.count() > 0 ? keptVariants / elapsed.count() : 0;
      int discarded = discardedVariantsCount + countVariantsWithoutHomAlt;
      std::cerr << "\r    [Processing] " << keptVariants
                << " variants written (" << elapsed.count() << "s elapsed, ~"
                << rate << " var/sec, " << discarded << " discarded)"
                << std::flush;
    }

    int aa_count = 0, Aa_count = 0, AA_count = 0, non_missing_count = 0;
    for (int i = 0; i < n_samples; ++i) {
      int *gt_ptr = (has_gt && gt_arr) ? gt_arr + i * 2 : NULL;
      float ds_value = (has_ds && ds_arr) ? ds_arr[i] : 0.0f;

      processSampleGenotypes(gt_ptr, ds_value, has_gt, has_ds, aa_count,
                             Aa_count, AA_count, non_missing_count,
                             roundedDosageCount);
    }

    // For dominance mode, we need at least minHomCount homozygous alternates
    if ((mode == "dominance" && AA_count < minHomCount) ||
        (mode != "dominance" && mode != "recessive")) {
      if (mode == "dominance" && AA_count < minHomCount) {
        countVariantsWithoutHomAlt++;
        if (countVariantsWithoutHomAlt <= limit) {
          std::cerr << "variant '" << bcf_hdr_id2name(hdr, rec->rid) << ":"
                    << (rec->n_allele > 1 ? rec->d.allele[1] : ".") << "' has "
                    << AA_count
                    << " homozygous alternate alleles (min required: "
                    << minHomCount << "). Skipping..\n";
        }
      }
      continue;
    }

    // Skip variants with no non-missing samples
    if (non_missing_count == 0)
      continue;

    keptVariants++;

    std::string variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                            std::to_string(rec->pos + 1) + ":" +
                            rec->d.allele[0] + ":" + rec->d.allele[1];

    // Only calculate frequencies and dominance dosages if needed
    double r = 0.0, h = 0.0, a = 0.0;
    double dom_dosage_aa = 0.0, dom_dosage_Aa = 0.0, dom_dosage_AA = 0.0;
    double localMinDomDosage = 0.0, localMaxDomDosage = 0.0;
    std::string group;

    if (mode == "dominance" || allInfo) {
      r = static_cast<double>(aa_count) / non_missing_count;
      h = static_cast<double>(Aa_count) / non_missing_count;
      a = static_cast<double>(AA_count) / non_missing_count;
    }

    if (mode == "dominance") {
      dom_dosage_aa = -h * a;
      dom_dosage_Aa = 2 * a * r;
      dom_dosage_AA = -h * r;

      if (scaleByGroup && groupMap.find(variantId) == groupMap.end()) {
        if (discardedVariantsCount < warningLimit) {
          std::cerr << "Warning: Variant '" << variantId
                    << "' not found in group map. Discarding.\n";
        }
        discardedVariantsCount++;
        continue;
      }

      localMinDomDosage = globalMinDomDosage;
      localMaxDomDosage = globalMaxDomDosage;

      if (scalePerVariant) {
        // For per-variant scaling, calculate local min/max on the fly
        localMinDomDosage =
            std::min({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
        localMaxDomDosage =
            std::max({dom_dosage_aa, dom_dosage_Aa, dom_dosage_AA});
      } else if (scaleByGroup && groupMap.find(variantId) != groupMap.end()) {
        // For group-based scaling, use group-specific min/max
        std::string groupName = groupMap.at(variantId);
        localMinDomDosage = groupDosages.at(groupName).first;
        localMaxDomDosage = groupDosages.at(groupName).second;
        group = groupName;
      }
      // else: scaleGlobally uses globalMinDomDosage/globalMaxDomDosage (already
      // set above)
    }

    if (setVariantId) {
      variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                  std::to_string(rec->pos + 1) + ":" + rec->d.allele[0] + ":" +
                  (rec->n_allele > 1 ? rec->d.allele[1] : ".");
    }

    // Collection for encoding preview
    if (mode == "dominance" && !previewPrinted &&
        encodingExamples.size() < MAX_EXAMPLES) {
      VariantEncodingExample ex;
      ex.variantId = variantId;
      ex.r = r;
      ex.h = h;
      ex.a = a;
      ex.raw_aa = -h * a;
      ex.raw_Aa = 2 * a * r;
      ex.raw_AA = -h * r;
      ex.variantMin = localMinDomDosage;
      ex.variantMax = localMaxDomDosage;

      if (scalePerVariant || scaleGlobally || scaleByGroup) {
        ex.scaled_aa = 2 * ((ex.raw_aa - localMinDomDosage) /
                            (localMaxDomDosage - localMinDomDosage));
        ex.scaled_Aa = 2 * ((ex.raw_Aa - localMinDomDosage) /
                            (localMaxDomDosage - localMinDomDosage));
        ex.scaled_AA = 2 * ((ex.raw_AA - localMinDomDosage) /
                            (localMaxDomDosage - localMinDomDosage));
      } else {
        ex.scaled_aa = ex.raw_aa;
        ex.scaled_Aa = ex.raw_Aa;
        ex.scaled_AA = ex.raw_AA;
      }

      if (scalingFactor != 0 && scalingFactor != 1.0) {
        ex.scaled_aa *= scalingFactor;
        ex.scaled_Aa *= scalingFactor;
        ex.scaled_AA *= scalingFactor;
      }
      encodingExamples.push_back(ex);

      if (encodingExamples.size() >= MAX_EXAMPLES) {
        std::string modeStr = "";
        if (scaleGlobally)
          modeStr = "Global";
        else if (scalePerVariant)
          modeStr = "Per-Variant";
        else if (scaleByGroup)
          modeStr = "By-Group";
        else
          modeStr = "None"; // Or Raw

        printVariantEncodingPreview(encodingExamples, modeStr, true);
        previewPrinted = true;
      }
    }

    // Output variant information
    std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1)
              << "\t" << variantId << "\t" << rec->d.allele[0] << "\t";

    if (rec->n_allele > 1) {
      std::cout << rec->d.allele[1];
    } else {
      std::cout << ".";
    }

    std::cout << "\t.\t.\t";
    std::cout << "AC=" << Aa_count + 2 * AA_count
              << ";AN=" << 2 * non_missing_count;

    if (scalePerVariant) {
      std::cout << ";variantMinDomDosage=" << localMinDomDosage
                << ";variantMaxDomDosage=" << localMaxDomDosage;
    } else if (scaleGlobally) {
      std::cout << ";globalMinDomDosage=" << globalMinDomDosage
                << ";globalMaxDomDosage=" << globalMaxDomDosage;
    } else if (scaleByGroup) {
      std::cout << ";groupMinDomDosage=" << localMinDomDosage
                << ";groupMaxDomDosage=" << localMaxDomDosage
                << ";group=" << group;
    }

    if (allInfo) {
      std::cout << ";r=" << r << ";h=" << h << ";a=" << a;
    }

    // Check for low variance if scaling is enabled (only in dominance mode)
    if (mode == "dominance" &&
        (scalePerVariant || scaleGlobally || scaleByGroup)) {
      double scaled_aa = 2 * ((dom_dosage_aa - localMinDomDosage) /
                              (localMaxDomDosage - localMinDomDosage));
      double scaled_Aa = 2 * ((dom_dosage_Aa - localMinDomDosage) /
                              (localMaxDomDosage - localMinDomDosage));
      double scaled_AA = 2 * ((dom_dosage_AA - localMinDomDosage) /
                              (localMaxDomDosage - localMinDomDosage));

      // Clamp values
      if (scaled_aa < 0)
        scaled_aa = 0;
      if (scaled_Aa < 0)
        scaled_Aa = 0;
      if (scaled_AA < 0)
        scaled_AA = 0;

      // Check if any two genotype dosages are too close (pairwise comparison)
      double diff_aa_Aa = std::abs(scaled_aa - scaled_Aa);
      double diff_aa_AA = std::abs(scaled_aa - scaled_AA);
      double diff_Aa_AA = std::abs(scaled_Aa - scaled_AA);

      double minDiff = std::min({diff_aa_Aa, diff_aa_AA, diff_Aa_AA});

      if (minDiff < lowVarianceThreshold) {
        lowVarianceVariantsCount++;
        if (lowVarianceVariantsCount <= warningLimit) {
          std::cerr << "Warning: Variant '" << variantId
                    << "' has genotypes with very similar scaled dosages (min "
                       "difference="
                    << minDiff << "). Scaled dosages: aa=" << scaled_aa
                    << ", Aa=" << scaled_Aa << ", AA=" << scaled_AA << "\n";
        }
      }
    }

    std::cout << "\tDS";

    // Pre-allocate buffer for output line (estimate: ~10 bytes per sample on
    // average)
    std::string output_line;
    output_line.reserve(n_samples * 10);

    // Output dosage values for each sample
    for (int i = 0; i < n_samples; ++i) {
      double dosage = 0.0;
      bool is_missing = false;
      int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
      float ds_value = has_ds ? ds_arr[i] : 0.0f;

      if (has_gt) {
        // Read genotype values once and cache them
        int gt0 = gt_ptr[0];
        int gt1 = gt_ptr[1];
        int allele1_missing = bcf_gt_is_missing(gt0);
        int allele2_missing = bcf_gt_is_missing(gt1);

        if (allele1_missing && allele2_missing) {
          is_missing = true;
        } else if (allele1_missing || allele2_missing) {
          // Haploid (one missing, one present)
          dosage = 0.0;
          haploidWarningCount++;
          if (haploidWarningCount <= 5) {
            std::cerr << "Warning: Haploid genotype detected for sample "
                      << i + 1 << " at " << variantId << ". Setting to 0."
                      << std::endl;
          }
        } else {
          // Diploid (both present) - use cached values
          int allele1 = bcf_gt_allele(gt0);
          int allele2 = bcf_gt_allele(gt1);
          int geno = -1;

          if (allele1 == 0 && allele2 == 0)
            geno = 0;
          else if (allele1 != allele2)
            geno = 1;
          else if (allele1 > 0 && allele2 > 0)
            geno = 2;

          if (geno >= 0) {
            if (mode == "dominance") {
              if (geno == 0)
                dosage = -h * a;
              else if (geno == 1)
                dosage = 2 * a * r;
              else if (geno == 2)
                dosage = -h * r;

              if (scalePerVariant || scaleGlobally || scaleByGroup) {
                dosage = 2 * ((dosage - localMinDomDosage) /
                              (localMaxDomDosage - localMinDomDosage));
                if (dosage < 0) {
                  dosage = 0;
                } else if (dosage > 2 + epsilon) {
                  std::cerr << "Error: Dosage value " << dosage
                            << " exceeds 2 + epsilon (" << (2 + epsilon) << ")"
                            << std::endl;
                  exit(1);
                }
              }
            } else if (mode == "recessive") {
              if (geno == 0)
                dosage = 0.0;
              else if (geno == 1)
                dosage = 0.0;
              else if (geno == 2)
                dosage = 2.0;
            }

            if (scalingFactor != 0)
              dosage *= scalingFactor;
          }
        }
      } else if (has_ds) {
        if (std::isnan(ds_value)) {
          is_missing = true;
        } else {
          int geno = std::round(ds_value);
          if (geno >= 0 && geno <= 2) {
            if (mode == "dominance") {
              if (geno == 0)
                dosage = -h * a;
              else if (geno == 1)
                dosage = 2 * a * r;
              else if (geno == 2)
                dosage = -h * r;

              if (scalePerVariant || scaleGlobally || scaleByGroup) {
                dosage = 2 * ((dosage - localMinDomDosage) /
                              (localMaxDomDosage - localMinDomDosage));
                if (dosage < 0) {
                  dosage = 0;
                } else if (dosage > 2 + epsilon) {
                  std::cerr << "Error: Dosage value " << dosage
                            << " exceeds 2 + epsilon (" << (2 + epsilon) << ")"
                            << std::endl;
                  exit(1);
                }
              }
            } else if (mode == "recessive") {
              if (geno == 0)
                dosage = 0.0;
              else if (geno == 1)
                dosage = 0.0;
              else if (geno == 2)
                dosage = 2.0;
            }
            if (scalingFactor != 0)
              dosage *= scalingFactor;
          }
        }
      } else {
        is_missing = true;
      }

      // Append to buffer instead of immediate output
      output_line += '\t';
      if (is_missing) {
        output_line += '.';
      } else if (mode == "recessive" && scalingFactor == 1.0) {
        // Fast path for recessive mode with no scaling: only outputs "0" or "2"
        output_line += (dosage == 0.0) ? '0' : '2';
      } else {
        // Convert dosage to string and append
        char buf[32];
        snprintf(buf, sizeof(buf), "%.6g", dosage);
        output_line += buf;
      }
    }

    // Write entire line at once
    std::cout << output_line << '\n';
  }

  free(gt_arr);
  free(ds_arr);
  bcf_destroy(rec);

  // If we collected fewer than MAX_EXAMPLES, print them now
  if (mode == "dominance" && !previewPrinted && !encodingExamples.empty()) {
    // ... logic for printing preview ...
    std::string modeStr = "";
    if (scaleGlobally)
      modeStr = "Global";
    else if (scalePerVariant)
      modeStr = "Per-Variant";
    else if (scaleByGroup)
      modeStr = "By-Group";
    else
      modeStr = "None"; // Or Raw

    if (!scaleGlobally && !scalePerVariant && !scaleByGroup)
      modeStr = "No Scaling";

    if (scalingFactor != 1.0 && scalingFactor != 0) {
      modeStr += " (Factor: " + std::to_string(scalingFactor) + ")";
    }

    printVariantEncodingPreview(encodingExamples, modeStr, false);
  }

  // Print runtime statistics
  std::cerr << std::endl; // Clear progress line
  std::cerr << "\nProcessing complete:" << std::endl;
  std::cerr << "  * Variants processed       : " << totalVariants << std::endl;
  std::cerr << "  * Variants kept            : " << keptVariants << std::endl;

  if (countVariantsWithoutHomAlt > 0 || discardedVariantsCount > 0) {
    std::cerr << "  * Variants discarded       : "
              << (countVariantsWithoutHomAlt + discardedVariantsCount)
              << std::endl;
    if (countVariantsWithoutHomAlt > 0)
      std::cerr << "      - No homozygous alt    : "
                << countVariantsWithoutHomAlt << std::endl;
    if (discardedVariantsCount > 0)
      std::cerr << "      - Not in group map     : " << discardedVariantsCount
                << std::endl;
  }

  if (lowVarianceVariantsCount > 0) {
    std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Low variance dosages : " << lowVarianceVariantsCount
              << std::endl;
  }
  if (roundedDosageCount > 0) {
    if (lowVarianceVariantsCount == 0)
      std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Rounded dosages      : " << roundedDosageCount
              << std::endl;
    std::cerr << "        (DS values were rounded to nearest integer)"
              << std::endl;
  }

  if (haploidWarningCount > 0) {
    if (lowVarianceVariantsCount == 0 && roundedDosageCount == 0)
      std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Haploid genotypes    : " << haploidWarningCount
              << std::endl;
    std::cerr << "        (set to 0, warnings truncated after 5)" << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc > 0) {
    std::string programName = argv[0];
    if (programName.find("transform") != std::string::npos) {
      std::cerr << "Warning: 'transform' is deprecated and will be removed in "
                   "a future release. Please use 'orthogonalize' instead."
                << std::endl;
    }
  }

  std::string pathInput;
  std::string mode = "dominance";
  double scalingFactor = 1.0;
  bool scalePerVariant = false;
  bool scaleGlobally = false;
  bool setVariantId = false;
  bool allInfo = false;
  std::string scaleByGroupPath;
  int minHomCount = 1;

  if (!parseArguments(argc, argv, pathInput, mode, scalingFactor,
                      scalePerVariant, scaleGlobally, setVariantId, allInfo,
                      scaleByGroupPath, minHomCount)) {
    return 1;
  }

  // Validate mode
  if (mode != "dominance" && mode != "recessive") {
    std::cerr << "Error: Invalid mode '" << mode
              << "'. Only 'dominance' or 'recessive' modes are supported."
              << std::endl;
    return 1;
  }

  // Validate minHomCount for dominance mode
  if (mode == "dominance" && minHomCount < 1) {
    std::cerr << "Error: --min-hom-count must be at least 1 when using "
                 "dominance mode."
              << std::endl;
    return 1;
  }

  htsFile *fp = bcf_open(pathInput.c_str(), "r");
  if (!fp) {
    std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput
              << std::endl;
    return 1;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  if (!hdr) {
    std::cerr << "Error: Cannot read header from VCF/BCF file." << std::endl;
    bcf_close(fp);
    return 1;
  }

  int n_samples = bcf_hdr_nsamples(hdr);

  // Validate that VCF has samples
  if (n_samples == 0) {
    std::cerr << "Error: VCF file has no samples." << std::endl;
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return 1;
  }

  // Calculate number of passes needed
  bool needTwoPass = scaleGlobally || !scaleByGroupPath.empty();
  int numPasses = needTwoPass ? 2 : 1;

  // Print recoding transformation (append to Parameters section)
  std::cerr << "  * Recoding       : ";
  if (mode == "recessive") {
    double scaledValue = 2.0 * scalingFactor;
    // Format the scaled value nicely
    if (scaledValue == static_cast<int>(scaledValue)) {
      std::cerr << "[0, 1, 2] -> [0, 0, " << static_cast<int>(scaledValue)
                << "]" << std::endl;
    } else {
      std::cerr << "[0, 1, 2] -> [0, 0, " << scaledValue << "]" << std::endl;
    }
  } else if (mode == "dominance") {
    std::cerr << "[0, 1, 2] -> [-h*a, 2*a*r, -h*r]" << std::endl;
  }
  std::cerr << "  * Number of passes: " << numPasses << std::endl;
  std::cerr << "  * Number of samples: " << n_samples << std::endl;
  std::cerr << std::endl; // Blank line after Parameters section

  // Validate that VCF has GT or DS format fields
  bool has_gt = hasFormat(hdr, "GT");
  bool has_ds = hasFormat(hdr, "DS");
  if (!has_gt && !has_ds) {
    std::cerr << "Error: VCF file must contain either GT (genotype) or DS "
                 "(dosage) format fields."
              << std::endl;
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return 1;
  }

  double globalMinDomDosage = std::numeric_limits<double>::max();
  double globalMaxDomDosage = std::numeric_limits<double>::lowest();
  std::set<std::string> chromosomes;
  std::map<std::string, std::string> groupMap;
  std::map<std::string, std::pair<double, double>> groupDosages;
  bool hasMissingValues = false;

  // Validate and load group map if specified
  if (!scaleByGroupPath.empty()) {
    groupMap = readGroupMap(scaleByGroupPath);
    if (groupMap.empty()) {
      std::cerr << "Error: Group map file is empty or contains no valid data: "
                << scaleByGroupPath << std::endl;
      bcf_hdr_destroy(hdr);
      bcf_close(fp);
      return 1;
    }
    for (const auto &pair : groupMap) {
      groupDosages[pair.second] = {std::numeric_limits<double>::max(),
                                   std::numeric_limits<double>::lowest()};
    }
  }

  if (needTwoPass) {
    std::cerr << "  * Parsing specified genomic regions (Pass 1/2)..."
              << std::endl;

    auto pass1_start = std::chrono::steady_clock::now();

    calculateGlobalAndGroupDosages(fp, hdr, globalMinDomDosage,
                                   globalMaxDomDosage, n_samples, chromosomes,
                                   groupMap, groupDosages, hasMissingValues);

    if (scaleGlobally) {
      std::cerr << "  * Global Min Dosage: " << globalMinDomDosage << std::endl;
      std::cerr << "  * Global Max Dosage: " << globalMaxDomDosage << std::endl;
    }
    if (hasMissingValues) {
      std::cerr << "  * Warning: Input contains missing values. These are "
                   "skipped/imputed during calculation."
                << std::endl;
    }

    auto pass1_end = std::chrono::steady_clock::now();
    auto pass1_duration = std::chrono::duration_cast<std::chrono::seconds>(
        pass1_end - pass1_start);
    std::cerr << "  * Pass 1 completed in " << pass1_duration.count()
              << " seconds" << std::endl;

    // Validate that we found at least one variant
    if (chromosomes.empty()) {
      std::cerr << "Error: VCF file contains no variants." << std::endl;
      bcf_hdr_destroy(hdr);
      bcf_close(fp);
      return 1;
    }
  }

  // For dominance mode with global/group scaling, validate we found
  // variants with homozygous alternates
  if (mode == "dominance" && needTwoPass &&
      globalMinDomDosage == std::numeric_limits<double>::max()) {
    std::cerr << "Error: No variants with homozygous alternate alleles found. "
                 "Dominance encoding requires at least one variant with AA "
                 "genotype."
              << std::endl;
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return 1;
  }

  bool scaleByGroup = !scaleByGroupPath.empty();

  // If we did two passes, we need to reopen the file
  if (needTwoPass) {
    std::vector<std::string> sortedContigs = sortChromosomes(chromosomes);
    bcf_close(fp);

    // Reopen file for actual processing
    fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
      std::cerr << "Error: Cannot open VCF/BCF file for reading: " << pathInput
                << std::endl;
      return 1;
    }

    hdr = bcf_hdr_read(fp);

    printHeader(hdr, sortedContigs, mode, globalMinDomDosage,
                globalMaxDomDosage, allInfo, scalePerVariant, scaleGlobally,
                scaleByGroup);

    std::cerr << "Processing variants (two passes)..." << std::endl;
  } else {
    // Single pass mode - print header without contigs (we'll collect them
    // during processing)
    std::vector<std::string> emptyContigs;
    printHeader(hdr, emptyContigs, mode, globalMinDomDosage, globalMaxDomDosage,
                allInfo, scalePerVariant, scaleGlobally, scaleByGroup);

    std::cerr << "Processing variants (single pass)..." << std::endl;
  }

  auto pass2_start = std::chrono::steady_clock::now();

  processVcfFile(fp, hdr, mode, globalMinDomDosage, globalMaxDomDosage, allInfo,
                 scalePerVariant, scaleGlobally, scaleByGroup, setVariantId,
                 scalingFactor, groupMap, groupDosages, chromosomes,
                 minHomCount);

  auto pass2_end = std::chrono::steady_clock::now();
  auto pass2_duration =
      std::chrono::duration_cast<std::chrono::seconds>(pass2_end - pass2_start);
  if (needTwoPass) {
    std::cerr << "  * Pass 2 completed in " << pass2_duration.count()
              << " seconds" << std::endl;
  } else {
    std::cerr << "  * Processing completed in " << pass2_duration.count()
              << " seconds" << std::endl;
  }

  bcf_hdr_destroy(hdr);
  bcf_close(fp);

  if (hasMissingValues) {
    std::cerr << "Note: VCF contained missing genotype values. Frequencies "
                 "were calculated using only non-missing samples."
              << std::endl;
  }

  return 0;
}
