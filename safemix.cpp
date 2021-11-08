#include "portable_rand.h"

#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "zlib.h"

#include <stdio.h>
#include <math.h>
#include <unistd.h>

const char *GIT_DIFF_FULL =
#include "gitdiff.txt"
;

#define MIN(a, b) (((a)<(b)) ? (a) : (b))
#define MAX(a, b) (((a)>(b)) ? (a) : (b))

bool ispowerof2(int n) {
    return 0 == (n & (n-1));
}

double umistr2prob(const char *str, uint32_t pos, uint32_t isize, uint32_t randseed) {
    const char *umistr1 = strchr(str, '#');
    const char *umistr = (NULL == umistr1 ? str : umistr1);
    uint32_t k1 = __ac_Wang_hash(__ac_X31_hash_string(umistr) ^ randseed);
    uint32_t k2 = __ac_Wang_hash(pos);
    uint32_t k3 = __ac_Wang_hash(isize);
    return (double)((k1+k2+k3)&0xffffff) / 0x1000000;
}

double qname2prob(const char *str, uint32_t randseed) {
    uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(str) ^ randseed);
    return (double)(k&0xffffff) / 0x1000000;
}

double capped(double val, double maxval = 1.0) {
    return ((val > maxval) ? maxval : val);
}


typedef struct {
    const double d = 1.5;
    const double e = 1.5;
    const double f = 0.05;
    const double i = 10.0;
    const double j = 10.0;
    const uint32_t r = 13;
    const uint32_t s = 23;
} arg_default_vals_t;

const arg_default_vals_t arg_default_vals;

void help(int argc, char **argv) {
    fprintf(stderr, "Program %s version %s (%s)\n", argv[0], COMMIT_VERSION, COMMIT_DIFF_SH);
    fprintf(stderr, "  This program mixes two bam files and is aware of the molecular-barcodes (also known as unique molecular identifiers (UMIs))\n");
    
    fprintf(stderr, "Usage: %s -o <OUTPUT-PREFIX> -a <tumor-INPUT-BAM> -b <normal-INPUT-BAM>\n", argv[0]);
    fprintf(stderr, "Optional parameters:\n");
    fprintf(stderr, "  -d <tumor-umi-size> average number of reads in a UMI family in the <TUMOR-INPUT-BAM> file [default to %f]\n", arg_default_vals.d); 
    fprintf(stderr, "  -e <normal-umi-size> average number of reads in a UMI family in the <NORMAL-INPUT-BAM> file [default to %f]\n", arg_default_vals.e);
    fprintf(stderr, "  -f <tumor-fraction> the fraction of DNA that comes from tumor [default to %f]\n", arg_default_vals.f); 
    fprintf(stderr, "  -i <tumor-initial-quantity> initial quantity of DNA in ng sequenced in the <TUMOR-INPUT-BAM> file [default to %f]\n", arg_default_vals.i);
    fprintf(stderr, "  -j <normal-initial-quantity> initial quantity of DNA in ng sequenced in the <NORMAL-INPUT-BAM> file [default to %f]\n",  arg_default_vals.j);
    fprintf(stderr, "  -r <random-seed-for-initial-quantity> random seed used to select the UMI from the initial quantity of DNA [default to %d]\n", arg_default_vals.r);
    fprintf(stderr, "  -s <random-seed-for-umi-size>\n random seed used to select the reads in each UMI [default to %d]\n",  arg_default_vals.s);
    fprintf(stderr, "  -U <use-only-umi>\n set the program to use only UMIs for identifying read families (discard read start and end positions) [default to unset]\n");
    
    fprintf(stderr, "Note:\n");
    fprintf(stderr, "<tumor-INPUT-BAM> and <normal-INPUT-BAM> both have to be sorted and indexed.\n");
    fprintf(stderr, "To detect UMI, this prgram first checks for the MI tag in each alignment record in <tumor-INPUT-BAM> and <normal-INPUT-BAM>. If the MI tag is absent, then the program checks for the string after the number-hash-pound sign (#) in the read name (QNAME).\n");
    fprintf(stderr, "<OUTPUT-PREFIX> is appended by the string literals \".tumor.bam\" and \".normal.bam\" (without the double quotes) to generate the tumor and normal simulated BAM filenames, respectively. The tumor and normal bam files have to be merged (e.g., by samtools merge) to simulate the sequenced sample. \n");
    exit(-1);
}

typedef struct {
    char *filename;
    double allele_frac;
    double init_qty_frac;
    double read_fam_frac;
} subsample_info_t;

int 
main(int argc, char **argv) {
    int flags, opt, option_index;
    char *tbam = NULL;
    char *nbam = NULL;
    char *outpref = NULL;
    double tosd = arg_default_vals.d; // tumor-over-sequencing-depth
    double nosd = arg_default_vals.e; // normal-over-sequencing-depth
    double tiq = arg_default_vals.i; // tumor-initial-quantity of DNA
    double niq = arg_default_vals.j; // normal-initial-quantity of DNA
    double defallelefrac = arg_default_vals.f;
    uint32_t randseed1 = arg_default_vals.r;
    uint32_t randseed2 = arg_default_vals.s;
    int use_only_umi = 0;
    
    while ((opt = getopt(argc, argv, "a:b:d:e:f:i:j:o:r:s:U")) != -1) {
        switch (opt) {
            case 'a': tbam = optarg; break;
            case 'b': nbam = optarg; break;
            case 'd': tosd = atof(optarg); break;
            case 'e': nosd = atof(optarg); break;
            case 'f': defallelefrac = atof(optarg); break;
            case 'i': tiq = atof(optarg); break;
            case 'j': niq = atof(optarg); break;
            case 'o': outpref = optarg; break;
            case 'r': randseed1 = atoi(optarg); break;
            case 's': randseed2 = atoi(optarg); break;
            case 'U': use_only_umi = 1; break;
            default: help(argc, argv);
        }
    }
    if (NULL == tbam || NULL == nbam || NULL == outpref) {
        help(argc, argv);
    }
    if (0 != randseed1) {
        portable_srand(randseed1);
        randseed1 = (uint32_t)portable_rand();
    }
    if (0 != randseed2) {
        portable_srand(randseed2);
        randseed2 = (uint32_t)portable_rand();
    }
    
    double tiqfrac = niq / tiq;
    double tosdfrac = nosd / tosd;
    
    fprintf(stderr, "%s\n=== version ===\n%s\n%s\n%s\n", argv[0], COMMIT_VERSION, COMMIT_DIFF_SH, GIT_DIFF_FULL);
    
    subsample_info_t t_subsample_info = { tbam, defallelefrac,           tiqfrac,     tosdfrac };
    subsample_info_t n_subsample_info = { nbam, 1.0 - defallelefrac, 1.0/tiqfrac, 1.0/tosdfrac };
    
    const double t_umi_draw_prob = capped(t_subsample_info.allele_frac) * capped(t_subsample_info.init_qty_frac);
    const double n_umi_draw_prob = capped(n_subsample_info.allele_frac) * capped(n_subsample_info.init_qty_frac);
    const double umi_draw_prob_mult = 1.0 / MIN(1.0, MAX(t_umi_draw_prob, n_umi_draw_prob));
    
    char *outbam = (char*)malloc(strlen(outpref) + 13);
    for (int i = 0; i < 2; i++) {
        strcpy(outbam, outpref);
        strcat(outbam, ((0 == i) ? ".tumor.bam" : ".normal.bam"));
        samFile *outbam_fp = sam_open(outbam, "wb");
        
        subsample_info_t subsample_info = ((0 == i) ? t_subsample_info : n_subsample_info);
        samFile *bam_fp = sam_open(subsample_info.filename, "r");
        sam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);
        bam1_t *bam_rec = bam_init1();
        
        int write_ret = sam_hdr_write(outbam_fp, bam_hdr);
        if (write_ret < 0) {
            fprintf(stderr, "Failed to write the SAM header to the file %s\n", outbam);
            abort();
        }
        
        const double umi_draw_prob = umi_draw_prob_mult * ((0 == i) ? t_umi_draw_prob : n_umi_draw_prob);
        const double read_draw_given_umi_prob = capped(subsample_info.read_fam_frac);
        
        while (sam_read1(bam_fp, bam_hdr, bam_rec) >= 0) {
            if (0 != (bam_rec->core.flag & 0x900)) { continue; }
            const auto *bam_aux_data = bam_aux_get(bam_rec, "MI"); // this tag is reserved (https://samtools.github.io/hts-specs/SAMtags.pdf)
            const char *umistr = ((bam_aux_data != NULL) ? bam_aux2Z(bam_aux_data) : bam_get_qname(bam_rec));
            const auto abegin = (use_only_umi ? (0) : MIN(bam_rec->core.pos, bam_rec->core.mpos));
            const auto aisize = (use_only_umi ? (0) : abs(bam_rec->core.isize));
            const double prob1 = umistr2prob(umistr, abegin, aisize, randseed1);
            const double prob2 = qname2prob(bam_get_qname(bam_rec), randseed2);
            if (prob1 < umi_draw_prob && prob2 < read_draw_given_umi_prob) {
                int write_ret = sam_write1(outbam_fp, bam_hdr, bam_rec);
                if (write_ret < 0) {
                    fprintf(stderr, "Failed to write the read %s at tid %d pos %d to the file %s\n", bam_get_qname(bam_rec), bam_rec->core.tid, bam_rec->core.pos, outbam);
                    abort();
                }
            }
        }
        bam_destroy1(bam_rec);
        bam_hdr_destroy(bam_hdr);
        sam_close(bam_fp);
        sam_close(outbam_fp);
    }
    free(outbam);
}

