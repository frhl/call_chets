#include <iostream>
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
        std::cerr << "\nProgram: chet tools v0.0.1 (get_non_ref_sites)\n" << std::endl;
        std::cerr << "\nUsage: " << argv[0] << " <input vcf/bcf> <output file>" << std::endl;
        std::cerr << "\nDescription:" << std::endl;
        std::cerr << "\tTakes a VCF/BCF with phased genotypes and identifies samples" << std::endl;
        std::cerr << "\twith at least one alternate allele. The result is a file in" << std::endl;
        std::cerr << "\tthe format of sample, variant, genotype. Note: this may" << std::endl;
        std::cerr << "\ttake a very long time if common variants are not excluded " << std::endl;
        std::cerr << "\tusing 'bcftools view --max-af 0.01 in.vcf -Oz -o out.vcf.gz'" << std::endl;
        std::cerr << "\nOptions:" << std::endl;
        std::cerr << " <input vcf> \tInput VCF/BCF with phased genotypes." << std::endl;       
        std::cerr << " <output file> \tOutput file (with .gz extension\n" << std::endl;       
	std::cerr << "Example:" << std::endl;
	std::cerr << " ./get_non_ref_sites.o test/trio.vcf trio.sites.txt.gz\n" << std::endl;
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

    float *pp_arr = NULL;
    int npp_arr = 0;
    int pp_id = bcf_hdr_id2int(header, BCF_DT_ID, "PP");
    bool has_pp = bcf_hdr_idinfo_exists(header, BCF_HL_FMT, pp_id);

    if (has_pp) {
        gzprintf(out, "Sample\tVariantIndex\tVariant\tGenotype\tPP\n");
    } else {
        gzprintf(out, "Sample\tVariantIndex\tVariant\tGenotype\n");
    }

    // Iterate over all samples and variants
    while (bcf_read(fp, header, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        // We're only interested in bi-allelic sites
        if (rec->n_allele != 2)
            continue;

        int32_t *gt_arr = NULL, ngt_arr = 0, ngt = 0;
        ngt = bcf_get_genotypes(header, rec, &gt_arr, &ngt_arr);

	// extract PP
        if (has_pp) {
             bcf_get_format_float(header, rec, "PP", &pp_arr, &npp_arr);
         }

        // Iterate over the samples
        for (int i = 0; i < bcf_hdr_nsamples(header); i++) {
            int32_t *ptr = gt_arr + i*2; // for each sample, we have two allele calls

            if (ptr[0]==bcf_gt_missing || ptr[1]==bcf_gt_missing)
                continue;

            int a1 = bcf_gt_allele(ptr[0]);
            int a2 = bcf_gt_allele(ptr[1]);
            if (a1 == 1 && a2 == 1) {
		if (has_pp){
		    gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t1|1\t\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
		} else {
		    gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t1|1\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
        	}    
	    }
            else if (a1 != a2) {
                if (has_pp) {
			gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t%s\t%.2f\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1], a1 < a2 ? "1|0" : "0|1", pp_arr[i]);
		} else {
			gzprintf(out, "%s\t%d\t%s:%d:%s:%s\t%s\n", header->samples[i], variant_counter+1, bcf_hdr_id2name(header, rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1], a1 < a2 ? "1|0" : "0|1");
		}   
	}
        }
        if (gt_arr) free(gt_arr);
        variant_counter++;
        if (variant_counter % 10 == 0) {
             time(&current_time);
             double elapsed_time = difftime(current_time, start_time); 
             printf("\rWorking.. (Processed %d variants in %.2f seconds).", variant_counter, elapsed_time);
        }
    }
    printf(" Done!\n");

    bcf_destroy(rec);
    bcf_hdr_destroy(header);
    bcf_close(fp);
    gzclose(out);

    return 0;
}

