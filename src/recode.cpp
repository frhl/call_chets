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

#include "cli_utils.hpp"
#include "hts_raii.hpp"
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

/**
 * Scaling Modes Overview:
 *
 * This tool supports three mutually exclusive scaling modes for dominance
 * encoding:
 *
 * 1. --scale-per-variant (Single Pass, Fast)
 *    - Each variant is scaled independently to [0, 2] range
 *    - Uses variant-specific min/max dominance dosages
 *    - Betas/SEs from regression are NOT comparable across variants
 *    - Best for: Individual variant testing, speed-critical applications
 *
 * 2. --scale-globally (Two Passes, Comparable)
 *    - All variants scaled using global min/max across entire file
 *    - Pass 1: Scan all variants to find global min/max
 *    - Pass 2: Apply global scaling and output
 *    - Betas/SEs ARE comparable across variants
 *    - Best for: Burden tests, group-based testing requiring comparable effects
 *
 * 3. --scale-by-group <file> (Two Passes, Group-Level Comparability)
 *    - Variants scaled using group-specific (e.g., gene-level) min/max
 *    - Requires mapping file: tab-separated with 'variant' and 'gene' columns
 *    - Betas/SEs comparable within groups, not across groups
 *    - Best for: Gene-based burden tests where within-gene comparability needed
 *
 * The raw dominance deviation formula:
 *   aa (hom ref):  -h * a
 *   Aa (het):       2 * a * r
 *   AA (hom alt):  -h * r
 * Where r = freq(aa), h = freq(Aa), a = freq(AA)
 *
 * Scaling formula: scaled = 2 * ((raw - min) / (max - min))
 */
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
  std::cerr
      << "\n  --min-hom-count <n>        Minimum number of minor homozygous "
         "alleles required\n                             (default: 1). "
         "Applies to min(AA,aa) in dominance mode\n";
  std::cerr << "  --min-het-count <n>        Minimum number of heterozygous "
               "alleles required\n                             (default: "
               "1)\n\nAdditional Options:\n";
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
  std::cerr << "    " << path << " --input sample.vcf.gz --mode recessive\n";

  std::cerr << "\nFrequency & Count Filtering:\n";
  std::cerr
      << "  --min-[aac|mac] <n>        Min Alternate/Minor Allele Count\n";
  std::cerr
      << "  --max-[aac|mac] <n>        Max Alternate/Minor Allele Count\n";
  std::cerr
      << "  --min-[aaf|maf] <f>        Min Alternate/Minor Allele Frequency\n";
  std::cerr
      << "  --max-[aaf|maf] <f>        Max Alternate/Minor Allele Frequency\n";
}

// sortChromosomes is now provided by cli_utils.hpp as call_chets::sortChromosomes

bool parseArguments(int argc, char *argv[], std::string &pathInput,
                    std::string &mode, double &scalingFactor,
                    bool &scalePerVariant, bool &scaleGlobally,
                    bool &setVariantId, bool &allInfo,
                    std::string &scaleByGroupPath, int &minHomCount,
                    int &minHetCount, int &minAAC, int &maxAAC, int &minMAC,
                    int &maxMAC, double &minAAF, double &maxAAF, double &minMAF,
                    double &maxMAF) {
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
    } else if (arg == "--min-het-count" && i + 1 < argc) {
      minHetCount = std::stoi(argv[++i]);
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
    } else if (arg == "--min-aac" && i + 1 < argc) {
      minAAC = std::stoi(argv[++i]);
    } else if (arg == "--max-aac" && i + 1 < argc) {
      maxAAC = std::stoi(argv[++i]);
    } else if (arg == "--min-mac" && i + 1 < argc) {
      minMAC = std::stoi(argv[++i]);
    } else if (arg == "--max-mac" && i + 1 < argc) {
      maxMAC = std::stoi(argv[++i]);
    } else if (arg == "--min-aaf" && i + 1 < argc) {
      minAAF = std::stod(argv[++i]);
    } else if (arg == "--max-aaf" && i + 1 < argc) {
      maxAAF = std::stod(argv[++i]);
    } else if (arg == "--min-maf" && i + 1 < argc) {
      minMAF = std::stod(argv[++i]);
    } else if (arg == "--max-maf" && i + 1 < argc) {
      maxMAF = std::stod(argv[++i]);
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
    return false;
  }

  std::map<std::string, std::string> files;
  files["Input"] = pathInput;
  if (!scaleByGroupPath.empty())
    files["Group File"] = scaleByGroupPath;

  std::map<std::string, std::string> params;
  params["Mode"] = mode;
  params["Min Hom Count"] = std::to_string(minHomCount);
  params["Min Het Count"] = std::to_string(minHetCount);
  if (minAAC >= 0)
    params["Min AAC"] = std::to_string(minAAC);
  if (maxAAC != std::numeric_limits<int>::max())
    params["Max AAC"] = std::to_string(maxAAC);
  if (minMAC >= 0)
    params["Min MAC"] = std::to_string(minMAC);
  if (maxMAC != std::numeric_limits<int>::max())
    params["Max MAC"] = std::to_string(maxMAC);
  if (minAAF >= 0.0)
    params["Min AAF"] = std::to_string(minAAF);
  if (maxAAF != std::numeric_limits<double>::max())
    params["Max AAF"] = std::to_string(maxAAF);
  if (minMAF >= 0.0)
    params["Min MAF"] = std::to_string(minMAF);
  if (maxMAF != std::numeric_limits<double>::max())
    params["Max MAF"] = std::to_string(maxMAF);
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

// readGroupMap is now provided by cli_utils.hpp as call_chets::readGroupMap

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

// RAII Wrappers are now provided by hts_raii.hpp
using call_chets::Bcf1UPtr;
using call_chets::BcfHdrUPtr;
using call_chets::HtsFileUPtr;
using call_chets::openVcf;

// Processing stats structure
struct ProcessingStats {
  long totalVariants = 0;
  long keptVariants = 0;
  int countVariantsWithoutHomRef = 0;
  int countVariantsWithoutHomAlt = 0;
  int countVariantsWithoutHet = 0;
  int discardedVariantsCount = 0;
  int lowVarianceVariantsCount = 0;
  long roundedDosageCount = 0;
  long haploidWarningCount = 0;
  int countVariantsMinAAC = 0;
  int countVariantsMaxAAC = 0;
  int countVariantsMinMAC = 0;
  int countVariantsMaxMAC = 0;
  int countVariantsMinAAF = 0;
  int countVariantsMaxAAF = 0;
  int countVariantsMinMAF = 0;
  int countVariantsMaxMAF = 0;
};

// Filter constants
constexpr int WARNING_LIMIT = 5;
constexpr double LOW_VARIANCE_THRESHOLD = 0.0001;
constexpr double EPSILON = 1e-8;

struct VariantFilterParams {
  int minHomCount;
  int minHetCount;
  int minAAC;
  int maxAAC;
  int minMAC;
  int maxMAC;
  double minAAF;
  double maxAAF;
  double minMAF;
  double maxMAF;
  const std::string &mode;
};

// --- Logic Components ---

// Check if variant passes frequency/count filters
bool passesFrequencyFilters(int aac_count, int ref_count,
                            const VariantFilterParams &params,
                            ProcessingStats &stats) {
  int total_alleles = aac_count + ref_count;
  if (total_alleles == 0)
    return false; // Should not happen for non-missing data

  int mac_count = std::min(aac_count, ref_count);
  double aaf_val = (double)aac_count / total_alleles;
  double maf_val = (double)mac_count / total_alleles;

  if ((params.minAAC >= 0 && aac_count < params.minAAC)) {
    stats.countVariantsMinAAC++;
    return false;
  }
  if (params.maxAAC != std::numeric_limits<int>::max() &&
      aac_count > params.maxAAC) {
    stats.countVariantsMaxAAC++;
    return false;
  }
  if ((params.minMAC >= 0 && mac_count < params.minMAC)) {
    stats.countVariantsMinMAC++;
    return false;
  }
  if (params.maxMAC != std::numeric_limits<int>::max() &&
      mac_count > params.maxMAC) {
    stats.countVariantsMaxMAC++;
    return false;
  }
  if ((params.minAAF >= 0.0 && aaf_val < params.minAAF)) {
    stats.countVariantsMinAAF++;
    return false;
  }
  if (params.maxAAF != std::numeric_limits<double>::max() &&
      aaf_val > params.maxAAF) {
    stats.countVariantsMaxAAF++;
    return false;
  }
  if ((params.minMAF >= 0.0 && maf_val < params.minMAF)) {
    stats.countVariantsMinMAF++;
    return false;
  }
  if (params.maxMAF != std::numeric_limits<double>::max() &&
      maf_val > params.maxMAF) {
    stats.countVariantsMaxMAF++;
    return false;
  }

  return true;
}

// Check genotype counts (hom/het) logic
bool passesGenotypeCountFilters(bcf_hdr_t *hdr, bcf1_t *rec, int aa_count,
                                int Aa_count, int AA_count,
                                const VariantFilterParams &params,
                                ProcessingStats &stats) {
  if (params.mode != "dominance")
    return true;

  if (aa_count < params.minHomCount) {
    stats.countVariantsWithoutHomRef++;
    if (stats.countVariantsWithoutHomRef <= WARNING_LIMIT) {
      std::string varId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                          std::to_string(rec->pos + 1) + ":" +
                          rec->d.allele[0] + ":" +
                          (rec->n_allele > 1 ? rec->d.allele[1] : ".");
      std::cerr << "variant '" << varId << "' has " << aa_count
                << " (aa) homozygous alleles (min required: "
                << params.minHomCount << " for minor hom). Skipping..\n";
    } else if (stats.countVariantsWithoutHomRef == WARNING_LIMIT + 1) {
      std::cerr << "(truncated...)" << std::endl;
    }
    return false;
  }

  if (AA_count < params.minHomCount) {
    stats.countVariantsWithoutHomAlt++;
    if (stats.countVariantsWithoutHomAlt <= WARNING_LIMIT) {
      std::string varId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                          std::to_string(rec->pos + 1) + ":" +
                          rec->d.allele[0] + ":" +
                          (rec->n_allele > 1 ? rec->d.allele[1] : ".");
      std::cerr << "variant '" << varId << "' has " << AA_count
                << " (AA) homozygous alleles (min required: "
                << params.minHomCount << " for minor hom). Skipping..\n";
    } else if (stats.countVariantsWithoutHomAlt == WARNING_LIMIT + 1) {
      std::cerr << "(truncated...)" << std::endl;
    }
    return false;
  }

  if (Aa_count < params.minHetCount) {
    stats.countVariantsWithoutHet++;
    if (stats.countVariantsWithoutHet <= WARNING_LIMIT) {
      std::string varId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                          std::to_string(rec->pos + 1) + ":" +
                          rec->d.allele[0] + ":" +
                          (rec->n_allele > 1 ? rec->d.allele[1] : ".");
      std::cerr << "variant '" << varId << "' has " << Aa_count
                << " heterozygotes (min required: " << params.minHetCount
                << "). Skipping..\n";
    } else if (stats.countVariantsWithoutHet == WARNING_LIMIT + 1) {
      std::cerr << "(truncated...)" << std::endl;
    }
    return false;
  }

  return true;
}

// Calculate Dominance Dosages
struct DominanceDosages {
  double aa, Aa, AA;
};

DominanceDosages calculateRawDominance(double r, double h, double a) {
  return {-h * a, 2 * a * r, -h * r};
}

double scaleDosage(double val, double minVal, double maxVal) {
  if (maxVal == minVal)
    return 0.0;
  return 2.0 * ((val - minVal) / (maxVal - minVal));
}

// Main logic for processing a single sample's dosage
double computeSampleDosage(int geno, const std::string &mode, double h,
                           double a, double r, double localMin, double localMax,
                           bool scalingEnabled, double scalingFactor) {
  double dosage = 0.0;
  if (mode == "dominance") {
    if (geno == 0)
      dosage = -h * a;
    else if (geno == 1)
      dosage = 2 * a * r;
    else if (geno == 2)
      dosage = -h * r;

    if (scalingEnabled) {
      dosage = scaleDosage(dosage, localMin, localMax);
      if (dosage < 0)
        dosage = 0;
      else if (dosage > 2 + EPSILON) {
        throw std::runtime_error("Dosage value " + std::to_string(dosage) +
                                 " exceeds 2 + epsilon");
      }
    }
  } else if (mode == "recessive") {
    if (geno == 0 || geno == 1)
      dosage = 0.0;
    else if (geno == 2)
      dosage = 2.0;
  }

  if (scalingFactor != 0 && scalingFactor != 1.0) {
    dosage *= scalingFactor;
  }
  return dosage;
}

void processVariant(
    bcf_hdr_t *hdr, bcf1_t *rec, const std::string &mode, bool allInfo,
    bool scalePerVariant, bool scaleGlobally, bool scaleByGroup,
    bool setVariantId, double scalingFactor, double globalMinDomDosage,
    double globalMaxDomDosage,
    const std::map<std::string, std::pair<double, double>> &groupDosages,
    const std::map<std::string, std::string> &groupMap,
    const VariantFilterParams &filterParams, ProcessingStats &stats,
    std::vector<VariantEncodingExample> &encodingExamples, bool &previewPrinted,
    int MAX_EXAMPLES) {

  // Get Genotypes/Dosages
  int *gt_arr = NULL, ngt_arr = 0;
  float *ds_arr = NULL;
  int nds_arr = 0;

  bool has_gt_format = hasFormat(hdr, "GT");
  bool has_ds_format = hasFormat(hdr, "DS");

  int ngt = has_gt_format ? bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) : 0;
  int nds = has_ds_format
                ? bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds_arr)
                : 0;

  bool has_gt = has_gt_format && ngt >= 0;
  bool has_ds = has_ds_format && nds >= 0;

  if (!has_gt && !has_ds) {
    free(gt_arr);
    free(ds_arr);
    return;
  }

  int n_samples = bcf_hdr_nsamples(hdr);
  int aa_count = 0, Aa_count = 0, AA_count = 0, non_missing_count = 0;

  // First pass over samples: calculate stats
  for (int i = 0; i < n_samples; ++i) {
    int *gt_ptr = (has_gt && gt_arr) ? gt_arr + i * 2 : NULL;
    float ds_value = (has_ds && ds_arr) ? ds_arr[i] : 0.0f;
    processSampleGenotypes(gt_ptr, ds_value, has_gt, has_ds, aa_count, Aa_count,
                           AA_count, non_missing_count,
                           stats.roundedDosageCount);
  }

  // Check genotype count filters
  if (!passesGenotypeCountFilters(hdr, rec, aa_count, Aa_count, AA_count,
                                  filterParams, stats)) {
    free(gt_arr);
    free(ds_arr);
    return;
  }

  if (non_missing_count == 0) {
    free(gt_arr);
    free(ds_arr);
    return;
  }

  // Calculate stats for frequency filters
  int aac_count = Aa_count + 2 * AA_count;
  int ref_count = 2 * aa_count + Aa_count;

  if (!passesFrequencyFilters(aac_count, ref_count, filterParams, stats)) {
    free(gt_arr);
    free(ds_arr);
    return;
  }

  // Passed all filters!
  stats.keptVariants++;

  std::string variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                          std::to_string(rec->pos + 1) + ":" +
                          rec->d.allele[0] + ":" + rec->d.allele[1];

  double r = 0.0, h = 0.0, a = 0.0;
  if (mode == "dominance" || allInfo) {
    r = static_cast<double>(aa_count) / non_missing_count;
    h = static_cast<double>(Aa_count) / non_missing_count;
    a = static_cast<double>(AA_count) / non_missing_count;
  }

  double localMin = 0.0, localMax = 0.0;
  std::string group;
  DominanceDosages rawDom = {0, 0, 0};

  if (mode == "dominance") {
    rawDom = calculateRawDominance(r, h, a);

    if (scaleByGroup && groupMap.find(variantId) == groupMap.end()) {
      if (stats.discardedVariantsCount < WARNING_LIMIT) {
        std::cerr << "Warning: Variant '" << variantId
                  << "' not found in group map. Discarding.\n";
      }
      stats.discardedVariantsCount++;
      free(gt_arr);
      free(ds_arr);
      return;
    }

    localMin = globalMinDomDosage;
    localMax = globalMaxDomDosage;

    if (scalePerVariant) {
      localMin = std::min({rawDom.aa, rawDom.Aa, rawDom.AA});
      localMax = std::max({rawDom.aa, rawDom.Aa, rawDom.AA});
    } else if (scaleByGroup) {
      group = groupMap.at(variantId);
      localMin = groupDosages.at(group).first;
      localMax = groupDosages.at(group).second;
    }
  }

  if (setVariantId) {
    variantId = std::string(bcf_hdr_id2name(hdr, rec->rid)) + ":" +
                std::to_string(rec->pos + 1) + ":" + rec->d.allele[0] + ":" +
                (rec->n_allele > 1 ? rec->d.allele[1] : ".");
  }

  // --- Preview Logic ---
  if (mode == "dominance" && !previewPrinted &&
      encodingExamples.size() < MAX_EXAMPLES) {
    VariantEncodingExample ex;
    ex.variantId = variantId;
    ex.r = r;
    ex.h = h;
    ex.a = a;
    ex.raw_aa = rawDom.aa;
    ex.raw_Aa = rawDom.Aa;
    ex.raw_AA = rawDom.AA;
    ex.variantMin = localMin;
    ex.variantMax = localMax;

    bool scaling = scalePerVariant || scaleGlobally || scaleByGroup;
    if (scaling) {
      ex.scaled_aa = scaleDosage(rawDom.aa, localMin, localMax);
      ex.scaled_Aa = scaleDosage(rawDom.Aa, localMin, localMax);
      ex.scaled_AA = scaleDosage(rawDom.AA, localMin, localMax);
    } else {
      ex.scaled_aa = rawDom.aa;
      ex.scaled_Aa = rawDom.Aa;
      ex.scaled_AA = rawDom.AA;
    }
    if (scalingFactor != 0 && scalingFactor != 1.0) {
      ex.scaled_aa *= scalingFactor;
      ex.scaled_Aa *= scalingFactor;
      ex.scaled_AA *= scalingFactor;
    }
    encodingExamples.push_back(ex);
    // (Check for printing preview done in caller or periodically)
  }

  // --- Output ---
  std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t"
            << variantId << "\t" << rec->d.allele[0] << "\t";
  if (rec->n_allele > 1)
    std::cout << rec->d.allele[1];
  else
    std::cout << ".";
  std::cout << "\t.\t.\tAC=" << aac_count << ";AN=" << 2 * non_missing_count;

  if (scalePerVariant)
    std::cout << ";variantMinDomDosage=" << localMin
              << ";variantMaxDomDosage=" << localMax;
  else if (scaleGlobally)
    std::cout << ";globalMinDomDosage=" << localMin
              << ";globalMaxDomDosage=" << localMax;
  else if (scaleByGroup)
    std::cout << ";groupMinDomDosage=" << localMin
              << ";groupMaxDomDosage=" << localMax << ";group=" << group;

  if (allInfo)
    std::cout << ";r=" << r << ";h=" << h << ";a=" << a;

  // Check low variance
  bool scaling = scalePerVariant || scaleGlobally || scaleByGroup;
  if (mode == "dominance" && scaling) {
    double s_aa = scaleDosage(rawDom.aa, localMin, localMax);
    double s_Aa = scaleDosage(rawDom.Aa, localMin, localMax);
    double s_AA = scaleDosage(rawDom.AA, localMin, localMax);
    if (s_aa < 0)
      s_aa = 0;
    if (s_Aa < 0)
      s_Aa = 0;
    if (s_AA < 0)
      s_AA = 0;

    double diff1 = std::abs(s_aa - s_Aa);
    double diff2 = std::abs(s_aa - s_AA);
    double diff3 = std::abs(s_Aa - s_AA);
    double minDiff = std::min({diff1, diff2, diff3});
    if (minDiff < LOW_VARIANCE_THRESHOLD) {
      stats.lowVarianceVariantsCount++;
      if (stats.lowVarianceVariantsCount <= WARNING_LIMIT) {
        std::cerr
            << "Warning: Variant '" << variantId
            << "' has genotypes with very similar scaled dosages (min diff="
            << minDiff << ").\n";
      }
    }
  }

  std::cout << "\tDS";

  // Build output line
  std::string output_line;
  output_line.reserve(n_samples * 10);

  for (int i = 0; i < n_samples; ++i) {
    double dosage = 0.0;
    bool is_missing = false;
    int *gt_ptr = has_gt ? gt_arr + i * 2 : NULL;
    float ds_value = has_ds ? ds_arr[i] : 0.0f;

    if (has_gt) {
      int allele1_missing = bcf_gt_is_missing(gt_ptr[0]);
      int allele2_missing = bcf_gt_is_missing(gt_ptr[1]);
      if (allele1_missing && allele2_missing) {
        is_missing = true;
      } else if (allele1_missing || allele2_missing) {
        // Haploid
        dosage = 0.0;
        stats.haploidWarningCount++;
        if (stats.haploidWarningCount <= 5) {
          std::cerr << "Warning: Haploid genotype sample " << i + 1 << " at "
                    << variantId << ". Setting to 0.\n";
        }
      } else {
        int g1 = bcf_gt_allele(gt_ptr[0]);
        int g2 = bcf_gt_allele(gt_ptr[1]);
        int geno = -1;
        if (g1 == 0 && g2 == 0)
          geno = 0;
        else if (g1 != g2)
          geno = 1;
        else if (g1 > 0 && g2 > 0)
          geno = 2;

        if (geno >= 0) {
          try {
            dosage = computeSampleDosage(geno, mode, h, a, r, localMin,
                                         localMax, scaling, scalingFactor);
          } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            free(gt_arr);
            free(ds_arr);
            throw;
          }
        }
      }
    } else if (has_ds) {
      if (std::isnan(ds_value)) {
        is_missing = true;
      } else {
        int geno = std::round(ds_value);
        if (geno >= 0 && geno <= 2) {
          try {
            dosage = computeSampleDosage(geno, mode, h, a, r, localMin,
                                         localMax, scaling, scalingFactor);
          } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            free(gt_arr);
            free(ds_arr);
            throw;
          }
        }
      }
    } else {
      is_missing = true;
    }

    output_line += '\t';
    if (is_missing) {
      output_line += '.';
    } else if (mode == "recessive" && scalingFactor == 1.0) {
      output_line += (dosage == 0.0) ? '0' : '2';
    } else {
      char buf[32];
      snprintf(buf, sizeof(buf), "%.6g", dosage);
      output_line += buf;
    }
  }

  std::cout << output_line << '\n';

  free(gt_arr);
  free(ds_arr);
}

void processVcfFile(
    htsFile *fp, bcf_hdr_t *hdr, const std::string &mode,
    double globalMinDomDosage, double globalMaxDomDosage, bool allInfo,
    bool scalePerVariant, bool scaleGlobally, bool scaleByGroup,
    bool setVariantId, double scalingFactor,
    const std::map<std::string, std::string> &groupMap,
    const std::map<std::string, std::pair<double, double>> &groupDosages,
    std::set<std::string> &chromosomes, int minHomCount, int minHetCount,
    int minAAC, int maxAAC, int minMAC, int maxMAC, double minAAF,
    double maxAAF, double minMAF, double maxMAF) {

  Bcf1UPtr rec(bcf_init());
  ProcessingStats stats;
  VariantFilterParams filterParams{minHomCount, minHetCount, minAAC, maxAAC,
                                   minMAC,      maxMAC,      minAAF, maxAAF,
                                   minMAF,      maxMAF,      mode};

  // Preview Logic
  std::vector<VariantEncodingExample> encodingExamples;
  bool previewPrinted = false;
  const int MAX_EXAMPLES = 5;

  auto processing_start_time = std::chrono::steady_clock::now();

  while (bcf_read(fp, hdr, rec.get()) == 0) {
    stats.totalVariants++;
    bcf_unpack(rec.get(), BCF_UN_STR);
    chromosomes.insert(bcf_hdr_id2name(hdr, rec->rid));

    // Progress logging
    if (stats.keptVariants > 0 && stats.keptVariants % 10000 == 0) {
      auto now = std::chrono::steady_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
          now - processing_start_time);
      int rate = elapsed.count() > 0 ? stats.keptVariants / elapsed.count() : 0;
      int discarded =
          stats.discardedVariantsCount + stats.countVariantsWithoutHomRef +
          stats.countVariantsWithoutHomAlt + stats.countVariantsWithoutHet;
      std::cerr << "\r    [Processing] " << stats.keptVariants
                << " variants written (" << elapsed.count() << "s elapsed, ~"
                << rate << " var/sec, " << discarded << " discarded)"
                << std::flush;
    }

    processVariant(hdr, rec.get(), mode, allInfo, scalePerVariant,
                   scaleGlobally, scaleByGroup, setVariantId, scalingFactor,
                   globalMinDomDosage, globalMaxDomDosage, groupDosages,
                   groupMap, filterParams, stats, encodingExamples,
                   previewPrinted, MAX_EXAMPLES);

    if (mode == "dominance" && !previewPrinted &&
        encodingExamples.size() >= MAX_EXAMPLES) {
      std::string modeStr = "";
      if (scaleGlobally)
        modeStr = "Global";
      else if (scalePerVariant)
        modeStr = "Per-Variant";
      else if (scaleByGroup)
        modeStr = "By-Group";
      else
        modeStr = "None";

      printVariantEncodingPreview(encodingExamples, modeStr, true);
      previewPrinted = true;
    }
  }

  // Final preview check
  if (mode == "dominance" && !previewPrinted && !encodingExamples.empty()) {
    std::string modeStr = "";
    if (scaleGlobally)
      modeStr = "Global";
    else if (scalePerVariant)
      modeStr = "Per-Variant";
    else if (scaleByGroup)
      modeStr = "By-Group";
    else
      modeStr = "None";

    if (!scaleGlobally && !scalePerVariant && !scaleByGroup)
      modeStr = "No Scaling";
    if (scalingFactor != 1.0 && scalingFactor != 0) {
      modeStr += " (Factor: " + std::to_string(scalingFactor) + ")";
    }
    printVariantEncodingPreview(encodingExamples, modeStr, false);
  }

  // Print Summary
  std::cerr << std::endl;
  std::cerr << "\nProcessing complete:" << std::endl;
  std::cerr << "  * Variants processed       : " << stats.totalVariants
            << std::endl;
  std::cerr << "  * Variants kept            : " << stats.keptVariants
            << std::endl;

  int totalDiscarded = stats.totalVariants - stats.keptVariants;

  if (totalDiscarded > 0) {
    std::cerr << "  * Variants discarded       : " << totalDiscarded
              << std::endl;
    if (stats.countVariantsWithoutHomRef > 0)
      std::cerr << "      - Hom ref count < " << std::left << std::setw(6)
                << minHomCount << ": " << stats.countVariantsWithoutHomRef
                << std::endl;
    if (stats.countVariantsWithoutHomAlt > 0)
      std::cerr << "      - Hom alt count < " << std::left << std::setw(6)
                << minHomCount << ": " << stats.countVariantsWithoutHomAlt
                << std::endl;
    if (stats.countVariantsWithoutHet > 0)
      std::cerr << "      - Het count < " << std::left << std::setw(10)
                << minHetCount << ": " << stats.countVariantsWithoutHet
                << std::endl;

    if (stats.countVariantsMinAAC > 0)
      std::cerr << "      - Min AAC < " << std::left << std::setw(10) << minAAC
                << ": " << stats.countVariantsMinAAC << std::endl;
    if (stats.countVariantsMaxAAC > 0)
      std::cerr << "      - Max AAC > " << std::left << std::setw(10) << maxAAC
                << ": " << stats.countVariantsMaxAAC << std::endl;
    if (stats.countVariantsMinMAC > 0)
      std::cerr << "      - Min MAC < " << std::left << std::setw(10) << minMAC
                << ": " << stats.countVariantsMinMAC << std::endl;
    if (stats.countVariantsMaxMAC > 0)
      std::cerr << "      - Max MAC > " << std::left << std::setw(10) << maxMAC
                << ": " << stats.countVariantsMaxMAC << std::endl;

    if (stats.countVariantsMinAAF > 0)
      std::cerr << "      - Min AAF < " << std::left << std::setw(10) << minAAF
                << ": " << stats.countVariantsMinAAF << std::endl;
    if (stats.countVariantsMaxAAF > 0)
      std::cerr << "      - Max AAF > " << std::left << std::setw(10) << maxAAF
                << ": " << stats.countVariantsMaxAAF << std::endl;
    if (stats.countVariantsMinMAF > 0)
      std::cerr << "      - Min MAF < " << std::left << std::setw(10) << minMAF
                << ": " << stats.countVariantsMinMAF << std::endl;
    if (stats.countVariantsMaxMAF > 0)
      std::cerr << "      - Max MAF > " << std::left << std::setw(10) << maxMAF
                << ": " << stats.countVariantsMaxMAF << std::endl;

    if (stats.discardedVariantsCount > 0)
      std::cerr << "      - Not in group map    : "
                << stats.discardedVariantsCount << std::endl;
  }

  if (stats.lowVarianceVariantsCount > 0) {
    std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Low variance dosages : "
              << stats.lowVarianceVariantsCount << std::endl;
  }
  if (stats.roundedDosageCount > 0) {
    if (stats.lowVarianceVariantsCount == 0)
      std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Rounded dosages      : " << stats.roundedDosageCount
              << std::endl;
    std::cerr << "        (DS values were rounded to nearest integer)"
              << std::endl;
  }

  if (stats.haploidWarningCount > 0) {
    if (stats.lowVarianceVariantsCount == 0 && stats.roundedDosageCount == 0)
      std::cerr << "  * Warnings:" << std::endl;
    std::cerr << "      - Haploid genotypes    : " << stats.haploidWarningCount
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
  int minHetCount = 1;
  int minAAC = -1, maxAAC = std::numeric_limits<int>::max();
  int minMAC = -1, maxMAC = std::numeric_limits<int>::max();
  double minAAF = -1.0, maxAAF = std::numeric_limits<double>::max();
  double minMAF = -1.0, maxMAF = std::numeric_limits<double>::max();

  if (!parseArguments(argc, argv, pathInput, mode, scalingFactor,
                      scalePerVariant, scaleGlobally, setVariantId, allInfo,
                      scaleByGroupPath, minHomCount, minHetCount, minAAC,
                      maxAAC, minMAC, maxMAC, minAAF, maxAAF, minMAF, maxMAF)) {
    return 1;
  }

  // Validate mode
  if (mode != "dominance" && mode != "recessive") {
    std::cerr << "Error: Invalid mode '" << mode
              << "'. Only 'dominance' or 'recessive' modes are supported."
              << std::endl;
    return 1;
  }

  // Validate minHomCount and minHetCount for dominance mode
  if (mode == "dominance" && minHomCount < 1) {
    std::cerr << "Error: --min-hom-count must be at least 1 when using "
                 "dominance mode."
              << std::endl;
    return 1;
  }

  if (minHetCount < 1) {
    std::cerr << "Error: --min-het-count must be at least 1." << std::endl;
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
    std::string errorMsg;
    groupMap = call_chets::readGroupMap(scaleByGroupPath, errorMsg);
    if (!errorMsg.empty()) {
      std::cerr << errorMsg << std::endl;
      bcf_hdr_destroy(hdr);
      bcf_close(fp);
      return 1;
    }
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
    std::vector<std::string> sortedContigs = call_chets::sortChromosomes(chromosomes);
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
                 minHomCount, minHetCount, minAAC, maxAAC, minMAC, maxMAC,
                 minAAF, maxAAF, minMAF, maxMAF);

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
