#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <iostream>
#include <map>
#include <string>

#include "hts_raii.hpp"
#include "logging.hpp"
#include "version.hpp"

static void printUsage(const char *path) {
  std::cerr << "\nProgram: VCF PP Filter\n\n";
  std::cerr << "Usage: " << path
            << " --input <VCF/BCF file> --output <output file> --pp-threshold "
               "<threshold>\n\n";

  std::cerr << "Description:\n";
  std::cerr << "  This program filters variants in a VCF/BCF file based on the "
               "PP (posterior probability) field.\n";
  std::cerr << "  If the PP value for a genotype is below the specified "
               "threshold, that genotype is set to missing.\n\n";

  std::cerr << "Arguments:\n";
  std::cerr << "  --input/-i <VCF/BCF file>   : Required. Path to the VCF/BCF "
               "file to be processed.\n";
  std::cerr << "  --output/-o <output file>   : Required. Path to the output "
               "file where the filtered VCF/BCF data will be written.\n";
  std::cerr << "  --pp-threshold/-p <float>   : Required. Threshold value for "
               "the PP field. Genotypes with PP below this value will be set "
               "to missing.\n\n";
  std::cerr << "  --verbose/-v                 : Optional. Enables verbose "
               "mode, printing progress every 100 variants processed.\n\n";

  std::cerr << "Notes:\n";
  std::cerr << "  Ensure that the VCF/BCF file is properly formatted.\n";
  std::cerr << "  The program reads the PP field from the FORMAT column of the "
               "VCF/BCF file.\n";
  std::cerr << "  The output file format (VCF or BCF) is determined by the "
               "extension of the output file path.\n";
}

int main(int argc, char *argv[]) {
  if (argc > 0) {
    std::string programName = argv[0];
    if (programName.find("filter_vcf_by_pp") != std::string::npos) {
      std::cerr
          << "Warning: 'filter_vcf_by_pp' is deprecated and will be removed in "
             "a future release. Please use 'filter_pp' instead."
          << std::endl;
    }
  }

  std::string inputPath;
  std::string outputFilePath;
  float ppThreshold = 0.0;
  bool verbose = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printUsage(argv[0]);
      return 0;
    } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
      inputPath = argv[++i];
    } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
      outputFilePath = argv[++i];
    } else if ((arg == "--pp-threshold" || arg == "-p") && i + 1 < argc) {
      ppThreshold = std::stof(argv[++i]);
    } else if (arg == "--verbose" || arg == "-v") {
      verbose = true;
    } else {
      std::cerr << "Error! Unknown or incomplete argument: " << arg
                << std::endl;
      printUsage(argv[0]);
      return 1;
    }
  }

  if (inputPath.empty() || outputFilePath.empty() || ppThreshold <= 0) {
    std::cerr
        << "Error: --input, --output, and --pp-threshold (>0) are required\n";
    printUsage(argv[0]);
    return 1;
  }

  std::map<std::string, std::string> files;
  files["Input"] = inputPath;
  files["Output"] = outputFilePath;

  std::map<std::string, std::string> params;
  params["PP Threshold"] = std::to_string(ppThreshold);
  params["Verbose"] = verbose ? "Yes" : "No";

  call_chets::printHeader("FILTER_PP", "Filter VCF by PP threshold", files,
                          params);

  // Use RAII for input file and header
  call_chets::HtsFileUPtr inFile(bcf_open(inputPath.c_str(), "r"));
  if (!inFile) {
    std::cerr << "Could not open input file: " << inputPath << std::endl;
    return 1;
  }

  call_chets::BcfHdrUPtr hdr(bcf_hdr_read(inFile.get()));
  if (!hdr) {
    std::cerr << "Error: Cannot read header from input file: " << inputPath
              << std::endl;
    return 1;
  }

  // Determine output format from extension
  const char *outMode = "w";
  if (outputFilePath.size() >= 4 &&
      outputFilePath.substr(outputFilePath.size() - 4) == ".bcf") {
    outMode = "wb";
  } else if (outputFilePath.size() >= 7 &&
             outputFilePath.substr(outputFilePath.size() - 7) == ".vcf.gz") {
    outMode = "wz";
  }

  call_chets::HtsFileUPtr outFile(bcf_open(outputFilePath.c_str(), outMode));
  if (!outFile) {
    std::cerr << "Could not open output file: " << outputFilePath << std::endl;
    return 1;
  }

  if (bcf_hdr_write(outFile.get(), hdr.get()) < 0) {
    std::cerr << "Error: Failed to write VCF header" << std::endl;
    return 1;
  }

  call_chets::Bcf1UPtr rec(bcf_init());
  int missing_count = 0;
  int variant_count = 0;

  while (bcf_read(inFile.get(), hdr.get(), rec.get()) == 0) {
    variant_count++;
    bcf_unpack(rec.get(), BCF_UN_ALL);

    // Get the PP field for each sample
    bcf_fmt_t *fmt_pp = bcf_get_fmt(hdr.get(), rec.get(), "PP");
    if (fmt_pp) {
      int num_samples = bcf_hdr_nsamples(hdr.get());

      // Use bcf_get_genotypes for safe genotype access
      int *gt_arr = nullptr;
      int ngt_arr = 0;
      int ngt = bcf_get_genotypes(hdr.get(), rec.get(), &gt_arr, &ngt_arr);

      if (ngt > 0 && gt_arr) {
        for (int i = 0; i < num_samples; ++i) {
          int allele0 = bcf_gt_is_missing(gt_arr[i * 2]) ? -1 : bcf_gt_allele(gt_arr[i * 2]);
          int allele1 = bcf_gt_is_missing(gt_arr[i * 2 + 1]) ? -1 : bcf_gt_allele(gt_arr[i * 2 + 1]);

          // Skip hom-ref (0|0) and hom-alt (1|1) — only filter heterozygotes
          if (allele0 == allele1) {
            continue;
          }

          // Accessing the PP value for each sample
          float *pp_array = reinterpret_cast<float *>(fmt_pp->p);
          float pp_value = pp_array[i];

          // Check PP value and set genotype to missing if below the threshold
          if (pp_value != bcf_float_missing && pp_value < ppThreshold) {
            float missing_value = bcf_float_missing;
            reinterpret_cast<float *>(fmt_pp->p)[i] = missing_value;
            missing_count++;
          }
        }
        free(gt_arr);
      }

      if (verbose && variant_count % 100 == 0) {
        std::cerr << "Processed " << variant_count << " variants. "
                  << "Heterozygotes set to missing so far: " << missing_count
                  << std::endl;
      }
    }

    // Write the processed record to output file
    bcf_write(outFile.get(), hdr.get(), rec.get());
  }

  std::cerr << "Total genotypes set to missing: " << missing_count << std::endl;

  return 0;
}
