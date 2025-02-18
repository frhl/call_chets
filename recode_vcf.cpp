#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

void printUsage(const char *path) {
    std::cerr << "Usage: " << path << " --input <input.vcf.gz> [--mode dominance|recessive] [--scale-dosages] [--scale-dosages-factor <factor>]\n";
}

bool parseArguments(int argc, char *argv[], std::string &pathInput, std::string &mode, 
                   double &scalingFactor, bool &scaleDosage) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]); 
            return false;
        } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
            pathInput = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            mode = argv[++i];
            if (mode != "dominance" && mode != "recessive") {
                std::cerr << "Error: mode must be either 'dominance' or 'recessive'\n";
                return false;
            }
        } else if (arg == "--scale-dosages-factor" && i + 1 < argc) {
            scalingFactor = std::stod(argv[++i]);
        } else if (arg == "--scale-dosages") {
            scaleDosage = true;
        }
    }
    
    if (pathInput.empty()) {
        std::cerr << "Error: --input required\n";
        return false;
    }
    
    if (mode == "recessive" && scaleDosage) {
        std::cerr << "Error: --scale-dosages cannot be used with recessive mode\n";
        return false;
    }
    
    return true;
}

void printHeader(const bcf_hdr_t *hdr, bool scaleDosage, const std::string &mode) {
    std::cout << "##fileformat=VCFv4.2\n";
    int ncontigs = hdr->n[BCF_DT_CTG];
    for (int i = 0; i < ncontigs; ++i) {
        std::cout << "##contig=<ID=" << bcf_hdr_int2id(hdr, BCF_DT_CTG, i) << ">\n";
    }

    if (mode == "dominance") {
        std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                  << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles\">\n"
                  << "##INFO=<ID=r,Number=1,Type=Float,Description=\"Frequency of bi-allelic references (aa)\">\n"
                  << "##INFO=<ID=h,Number=1,Type=Float,Description=\"Frequency of heterozygotes (Aa)\">\n"
                  << "##INFO=<ID=a,Number=1,Type=Float,Description=\"Frequency of bi-allelic alternates (AA)\">\n";
        
        if (scaleDosage) {
            std::cout << "##INFO=<ID=localMinDomDosage,Number=1,Type=Float,Description=\"Local minimum dominance dosage\">\n"
                      << "##INFO=<ID=localMaxDomDosage,Number=1,Type=Float,Description=\"Local maximum dominance dosage\">\n";
        }
        
        std::cout << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dominance scaled dosage\">\n";
    } else {
        std::cout << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                  << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles\">\n"
                  << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Recessive encoded dosage (0=ref or het, 1=hom alt)\">\n";
    }
    
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        std::cout << "\t" << bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
    }
    std::cout << "\n";
}

void processVariant(bcf1_t *rec, bcf_hdr_t *hdr, int n_samples, int *gt_arr, float *ds_arr, 
                   bool has_gt, bool has_ds, bool scaleDosage, double scalingFactor,
                   const std::string &mode) {
    int aa_count = 0, Aa_count = 0, AA_count = 0;
    
    // First pass: count genotypes
    for (int i = 0; i < n_samples; ++i) {
        if (has_gt && !bcf_gt_is_missing(gt_arr[i*2]) && !bcf_gt_is_missing(gt_arr[i*2+1])) {
            int a1 = bcf_gt_allele(gt_arr[i*2]);
            int a2 = bcf_gt_allele(gt_arr[i*2+1]);
            if (a1 == 0 && a2 == 0) aa_count++;
            else if (a1 != a2) Aa_count++;
            else if (a1 > 0 && a2 > 0) AA_count++;
        } else if (has_ds && !std::isnan(ds_arr[i])) {
            int geno = std::round(ds_arr[i]);
            if (geno >= 0 && geno <= 2) {
                if (geno == 0) aa_count++;
                else if (geno == 1) Aa_count++;
                else if (geno == 2) AA_count++;
            }
        }
    }

    if (AA_count == 0) return;

    // Calculate frequencies and dominance values if in dominance mode
    double r = 0, h = 0, a = 0;
    double dom_aa = 0, dom_Aa = 0, dom_AA = 0;
    double minDom = 0, maxDom = 0;

    if (mode == "dominance") {
        r = static_cast<double>(aa_count) / n_samples;
        h = static_cast<double>(Aa_count) / n_samples;
        a = static_cast<double>(AA_count) / n_samples;

        dom_aa = -h * a;
        dom_Aa = 2 * a * r;
        dom_AA = -h * r;

        if (scaleDosage) {
            minDom = std::min(std::min(dom_aa, dom_Aa), dom_AA);
            maxDom = std::max(std::max(dom_aa, dom_Aa), dom_AA);
        }
    }

    // Output variant info
    std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" 
              << (rec->pos + 1) << "\t"
              << ".\t"
              << rec->d.allele[0] << "\t"
              << (rec->n_allele > 1 ? rec->d.allele[1] : ".") << "\t.\t.\t"
              << "AC=" << (Aa_count + 2*AA_count) 
              << ";AN=" << (2*n_samples);

    if (mode == "dominance") {
        std::cout << ";r=" << r 
                  << ";h=" << h
                  << ";a=" << a;

        if (scaleDosage) {
            std::cout << ";localMinDomDosage=" << minDom
                     << ";localMaxDomDosage=" << maxDom;
        }
    }
    
    std::cout << "\tDS";

    // Process each sample
    for (int i = 0; i < n_samples; ++i) {
        double dosage = 0.0;
        
        if (mode == "dominance") {
            if (has_gt && !bcf_gt_is_missing(gt_arr[i*2]) && !bcf_gt_is_missing(gt_arr[i*2+1])) {
                int a1 = bcf_gt_allele(gt_arr[i*2]);
                int a2 = bcf_gt_allele(gt_arr[i*2+1]);
                if (a1 == 0 && a2 == 0) dosage = dom_aa;
                else if (a1 != a2) dosage = dom_Aa;
                else if (a1 > 0 && a2 > 0) dosage = dom_AA;
            } else if (has_ds && !std::isnan(ds_arr[i])) {
                int geno = std::round(ds_arr[i]);
                if (geno == 0) dosage = dom_aa;
                else if (geno == 1) dosage = dom_Aa;
                else if (geno == 2) dosage = dom_AA;
            }

            if (scaleDosage && (maxDom - minDom) > 0) {
                dosage = 2 * ((dosage - minDom) / (maxDom - minDom));
                if (dosage > 2) dosage = 2;
                if (dosage < 0) dosage = 0;
            }
            
            dosage *= scalingFactor;
        } else { // recessive mode
            if (has_gt && !bcf_gt_is_missing(gt_arr[i*2]) && !bcf_gt_is_missing(gt_arr[i*2+1])) {
                int a1 = bcf_gt_allele(gt_arr[i*2]);
                int a2 = bcf_gt_allele(gt_arr[i*2+1]);
                dosage = (a1 > 0 && a2 > 0) ? 1.0 : 0.0;
            } else if (has_ds && !std::isnan(ds_arr[i])) {
                int geno = std::round(ds_arr[i]);
                dosage = (geno == 2) ? 1.0 : 0.0;
            }
        }

        std::cout << "\t" << dosage;
    }
    std::cout << "\n";
}

int main(int argc, char *argv[]) {
    std::string pathInput;
    std::string mode = "dominance";
    double scalingFactor = 1.0;
    bool scaleDosage = false;

    if (!parseArguments(argc, argv, pathInput, mode, scalingFactor, scaleDosage)) return 1;

    htsFile *fp = bcf_open(pathInput.c_str(), "r");
    if (!fp) {
        std::cerr << "Error opening " << pathInput << "\n";
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error reading header\n";
        bcf_close(fp);
        return 1;
    }

    printHeader(hdr, scaleDosage, mode);

    bcf1_t *rec = bcf_init();
    int n_samples = bcf_hdr_nsamples(hdr);
    int *gt_arr = NULL, ngt_arr = 0;
    float *ds_arr = NULL;
    int nds_arr = 0;
    bool has_gt = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT") >= 0;
    bool has_ds = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS") >= 0;

    int *nds_ptr = &nds_arr;
    
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        if (has_gt) bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (has_ds) bcf_get_format_float(hdr, rec, "DS", &ds_arr, nds_ptr);
        processVariant(rec, hdr, n_samples, gt_arr, ds_arr, has_gt, has_ds, 
                      scaleDosage, scalingFactor, mode);
    }

    free(gt_arr);
    free(ds_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return 0;
}
