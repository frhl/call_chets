#ifndef HTS_RAII_HPP
#define HTS_RAII_HPP

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <memory>
#include <stdexcept>
#include <string>

namespace call_chets {

/**
 * @brief RAII wrapper deleters and smart pointer typedefs for HTSlib resources.
 *
 * These wrappers ensure proper cleanup of HTSlib resources when they go out
 * of scope, preventing memory leaks even in the presence of exceptions.
 */

/**
 * @brief Deleter for htsFile* (VCF/BCF file handle).
 */
struct HtsFileDeleter {
  void operator()(htsFile *fp) const {
    if (fp)
      bcf_close(fp);
  }
};

/**
 * @brief RAII smart pointer for htsFile*.
 */
using HtsFileUPtr = std::unique_ptr<htsFile, HtsFileDeleter>;

/**
 * @brief Deleter for bcf_hdr_t* (VCF/BCF header).
 */
struct BcfHdrDeleter {
  void operator()(bcf_hdr_t *hdr) const {
    if (hdr)
      bcf_hdr_destroy(hdr);
  }
};

/**
 * @brief RAII smart pointer for bcf_hdr_t*.
 */
using BcfHdrUPtr = std::unique_ptr<bcf_hdr_t, BcfHdrDeleter>;

/**
 * @brief Deleter for bcf1_t* (VCF/BCF record).
 */
struct Bcf1Deleter {
  void operator()(bcf1_t *rec) const {
    if (rec)
      bcf_destroy(rec);
  }
};

/**
 * @brief RAII smart pointer for bcf1_t*.
 */
using Bcf1UPtr = std::unique_ptr<bcf1_t, Bcf1Deleter>;

/**
 * @brief Safely open a VCF/BCF file with RAII.
 *
 * @param path Path to the VCF/BCF file.
 * @param mode File mode ("r" for read, "w" for write, etc.).
 * @return HtsFileUPtr RAII-wrapped file handle.
 * @throws std::runtime_error if the file cannot be opened.
 */
inline HtsFileUPtr openVcf(const std::string &path, const char *mode = "r") {
  HtsFileUPtr fp(bcf_open(path.c_str(), mode));
  if (!fp) {
    throw std::runtime_error("Cannot open VCF/BCF file: " + path);
  }
  return fp;
}

/**
 * @brief Safely read VCF/BCF header with RAII.
 *
 * @param fp RAII-wrapped file handle.
 * @return BcfHdrUPtr RAII-wrapped header.
 * @throws std::runtime_error if the header cannot be read.
 */
inline BcfHdrUPtr readVcfHeader(htsFile *fp) {
  BcfHdrUPtr hdr(bcf_hdr_read(fp));
  if (!hdr) {
    throw std::runtime_error("Cannot read header from VCF/BCF file.");
  }
  return hdr;
}

/**
 * @brief Create a new BCF record with RAII.
 *
 * @return Bcf1UPtr RAII-wrapped BCF record.
 * @throws std::runtime_error if the record cannot be created.
 */
inline Bcf1UPtr createBcfRecord() {
  Bcf1UPtr rec(bcf_init());
  if (!rec) {
    throw std::runtime_error("Cannot initialize BCF record.");
  }
  return rec;
}

} // namespace call_chets

#endif // HTS_RAII_HPP
