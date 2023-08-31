#include <stdio.h>
#include <zlib.h>
#include <time.h>
#include <string.h>
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
    time_t start_time, current_time;
    time(&start_time);

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <in.vcf/in.vcf.gz/in.bcf> <out.txt.gz> \n", argv[0]);
        fprintf(stderr, "\nNote: Ensure that <in.vcf> has already been filered by MAF/AF.\n");
        return 1;
    }

    const char *suffix = ".gz";
    if (strlen(argv[2]) < strlen(suffix) || strcmp(argv[2] + strlen(argv[2]) - strlen(suffix), suffix) != 0) {
        fprintf(stderr, "Error: The output file must have a '.gz' extension.\n");
        return 1;
    }

    fp = bcf_open(argv[1], "r");
    header = bcf_hdr_read(fp);

    gzFile out = gzopen(argv[2], "wb");

    gzprintf(out, "Sample\tVariantIndex\tVariant\tGenotype\n");

    // Iterate over all samples and variants
    while (bcf_read(fp, header, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        // We're only interested in bi-allelic sites
        if (rec->n_allele != 2)
            continue;

        int32_t *gt_arr = NULL, ngt_arr = 0, ngt = 0;
        ngt = bcf_get_genotypes(header, rec, &gt_arr, &ngt_arr);

        // Iterate over the samples
        for (int i = 0; i < bcf_hdr_nsamples(header); i++) {
            int32_t *ptr = gt_arr + i*2; // for each sample, we have two allele calls

            if (ptr[0]==bcf_gt_missing || ptr[1]==bcf_gt_missing)
                continue;

            int a1 = bcf_gt_allele(ptr[0]);
            int a2 = bcf_gt_allele(ptr[1]);
            if (a1 == 1 && a2 == 1) {
                gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t1|1\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
            }
            else if (a1 != a2) {
                gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t%s\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1], a1 < a2 ? "1|0" : "0|1");
            }
        }
        if (gt_arr) free(gt_arr);
        variant_counter++;
        if (variant_counter % 1000 == 0) {
             time(&current_time);
             double elapsed_time = difftime(current_time, start_time); 
             printf("Processed %d variants in %.2f seconds.\n", variant_counter, elapsed_time);
        }
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(header);
    bcf_close(fp);
    gzclose(out);

    return 0;
}

