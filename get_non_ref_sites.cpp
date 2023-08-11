#include <stdio.h>
#include <htslib/vcf.h>

// Structure to store indices of heterozygous genotypes
typedef struct {
    int sample;
    int index;
} heterozygous_t;

int main(int argc, char *argv[]) {
    htsFile *fp;
    bcf_hdr_t *header;
    bcf1_t *rec = bcf_init();
    int variant_counter = 0;
    float maf_threshold = 0.01; // Default MAF threshold

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <in.vcf.gz> [MAF threshold]\n", argv[0]);
        return 1;
    }

    fp = bcf_open(argv[1], "r");
    if (argc > 2)
        maf_threshold = atof(argv[2]); 
    header = bcf_hdr_read(fp);

    printf("Sample\tVariantIndex\tVariant\tGenotype\n");

    // Iterate over all samples and variants
    while (bcf_read(fp, header, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        // We're only interested in bi-allelic sites
        if (rec->n_allele != 2)
            continue;

        float *maf = NULL;
        int maf_len = 0;
        if (bcf_get_info_float(header, rec, "MAF", &maf, &maf_len) <= 0) { 
            // If MAF is not defined, print an error and stop the program
            fprintf(stderr, "Error: MAF not defined for variant at %s:%d\n", bcf_hdr_id2name(header, rec->rid), rec->pos+1);
            return 1;
        }

        // If MAF is higher than the threshold, skip this variant
        if (maf[0] > maf_threshold) {
            free(maf);
            continue;
        }

        int32_t *gt_arr = NULL, ngt_arr = 0, ngt = 0;
        ngt = bcf_get_genotypes(header, rec, &gt_arr, &ngt_arr);

        // Iterate over the samples
        for (int i = 0; i < bcf_hdr_nsamples(header); i++) {
            int32_t *ptr = gt_arr + i*2; // for each sample, we have two allele calls

            if (ptr[0]==bcf_gt_missing || ptr[1]==bcf_gt_missing)
                continue;

            int a1 = bcf_gt_allele(ptr[0]);
            int a2 = bcf_gt_allele(ptr[1]);

            // Homozygous alternative
            if (a1 == 1 && a2 == 1) {
                printf("%s\t%d\t%s:%d:%s:%s\t1|1\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
            }
            // Heterozygous
            else if (a1 != a2) {
                printf("%s\t%d\t%s:%d:%s:%s\t%s\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1], a1 < a2 ? "1|0" : "0|1");
            }
        }
        if (gt_arr) free(gt_arr);
        if (maf) free(maf);
        variant_counter++;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(header);
    bcf_close(fp);

    return 0;
}

