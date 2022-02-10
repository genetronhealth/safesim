#include "portable_rand.h"
#include "version.h"

#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "zlib.h"

#include <deque>
#include <string>
#include <vector>

#include <stdio.h>
#include <math.h>
#include <unistd.h>

/**
>>> This is the pseudocode for the simulator
Input: coordinate-sorted bam and coordinate-sorted vcf
Output: R1 and R2 fastq files with (almost surely) different read orders.
    Please note that fastq-tools can sort these two fastq files by read name.
for each read in the sorted bam
    while the current variant is before the begin position of this read: iterate to the next variant
    calculate the hash value of the UMI
    for each variant between the begin and end positions of the read:
        if the hash value divided by the max hash value is less than the given variant allele fraction:
            modify this bam read with this variant
    flush out this bam read
**/

const char *GIT_DIFF_FULL =
#include "gitdiff.txt"
;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

const double DEFAULT_ALLELE_FRAC = 0.1;
const int DEFAULT_SNV_BQ_PHRED = -1;
const int DEFAULT_INS_BQ_PHRED = 30;
const double DEFAULT_POWER_LAW_EXPONENT = -1.0; // 3.0;
const double DEFAULT_LOGNORMAL_DISP = 15.0;

const uint32_t DEFAULT_RANDSEED = 13;
const uint32_t DEFAULT_NITERS1 = 50;
const uint32_t DEFAULT_NITERS2 = 500;

const char *ACGT = "ACGT";
const char *TAG_FA = "FA";

bool ispowerof2(int n) {
    return 0 == (n & (n-1));
}

double phred2prob(int8_t phred) {
    return pow(10.0, -(phred / 10.0));
}
int prob2phred(double prob) {
    return -(int)round(10.0 / log(10.0) * log(prob));
}

uint32_t hashes2hash(std::vector<uint32_t>  hashes) {
    uint32_t ret = 0;
    for (uint32_t hash : hashes) {
        ret = __ac_Wang_hash(hash ^ ret);
    }
    return ret;
}

double umistr2prob(uint32_t &umihash, uint32_t randseed, uint32_t begpos, uint32_t endpos, const char *str) {
    
    const char *umistr1 = strchr(str, '#');
    const char *umistr = (NULL == umistr1 ? str : (umistr1 + 1));
    
    std::vector<uint32_t>  hashes;
    hashes.reserve(4 + 1);
    hashes.push_back(randseed);
    hashes.push_back(begpos);
    hashes.push_back(endpos);
    size_t umi_strlen = strlen(umistr);
    
    if ((umi_strlen % 2 == 1) && (umistr[(umi_strlen - 1) / 2] == '+') && umi_strlen <= 16 * 2 - 3) {
        char alpha[16] = {0};
        char beta[16] = {0};
        strncpy(alpha, umistr, (umi_strlen - 1) / 2);
        strncpy(beta, umistr + (umi_strlen - 1) / 2 + 1, (umi_strlen - 1) / 2);
        const char *abmin = ((strcmp(alpha, beta) <= 0) ? alpha : beta);
        const char *abmax = ((strcmp(alpha, beta) >= 0) ? alpha : beta);
        hashes.push_back(__ac_X31_hash_string(abmin));
        hashes.push_back(__ac_X31_hash_string(abmax));
        //fprintf(stderr, "The read %s has one duplex UMI with length %d for UMI %s\n", umistr, umi_strlen, umistr);
    } else {
        hashes.push_back(__ac_X31_hash_string(umistr));
        //fprintf(stderr, "The read %s has  NO duplex UMI with length %d for UMI %s\n", umistr, umi_strlen, umistr);
    }
    uint32_t k = hashes2hash(hashes);
    umihash = k;
    return (double)(k&0xffffff) / 0x1000000;
}

double qnameqpos2prob(uint32_t &hash, uint32_t randseed, const char *qname, int qpos) {
    std::vector<uint32_t>  hashes;
    hashes.reserve(3);
    hashes.push_back(randseed);
    hashes.push_back(__ac_X31_hash_string(qname));
    hashes.push_back(qpos);
    uint32_t k = hashes2hash(hashes);
    hash = k;
    return (double)(k&0xffffff) / 0x1000000;
}

double 
allelefrac_powlaw_transform(
        double allelefrac,
        uint32_t tid,
        uint32_t rpos,
        uint32_t samplehash1,
        uint32_t samplehash2,
        double exponent) {
    std::vector<uint32_t> ks1;
    ks1.push_back(samplehash1);
    ks1.push_back(tid);
    ks1.push_back(rpos);
    uint32_t k1 = hashes2hash(ks1);
    std::vector<uint32_t> ks2;
    ks2.push_back(samplehash2);
    ks2.push_back(tid);
    ks2.push_back(rpos);
    uint32_t k2 = hashes2hash(ks2);
    double altfrac =         allelefrac * pow((double)(k1 & 0xffffff) / (double)0x1000000, 1.0 / exponent);
    double reffrac = (1.0 - allelefrac) * pow((double)(k2 & 0xffffff) / (double)0x1000000, 1.0 / exponent);
    double odds_ratio = ((altfrac + 1e-9) / (reffrac + 1e-9));
    return odds_ratio / (1.0 + odds_ratio);
}

int sgn(double x) {
    if (0 == x) { return 0; }
    if (x > 0) { return  1; }
    if (x < 0) { return -1; }
}

double 
allelefrac_lognormal_transform(
        double allelefrac,
        uint32_t tid,
        uint32_t rpos,
        uint32_t samplehash1,
        uint32_t samplehash2,
        double lnsigma) {
    std::vector<uint32_t> ks1;
    ks1.push_back(samplehash1);
    ks1.push_back(tid);
    ks1.push_back(rpos);
    uint32_t k1 = hashes2hash(ks1);
    std::vector<uint32_t> ks2;
    ks2.push_back(samplehash2);
    ks2.push_back(tid);
    ks2.push_back(rpos);
    uint32_t k2 = hashes2hash(ks2);
    // if norm-z-score(rv) is log(2), then its lognormal rv doubles
    //   at rv=log(2), std-norm       density-val is exp(-1/2 * ((log(2) - 0) / 1)**2)
    //   at rv=log(2), norm(0, stdev) density-val is exp(-1/2 * ((log(2) / stdev - 0) / 1)**2), which is frac
    //   frac = exp(-1/2 * ((log(2)- 0) / stdev)**2)
    //   stdev = sqrt(log(frac) / (-1/2)) / log(2)
    
    // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    const double two_pi = 3.14159265358979323846 * 2.0;
    double mu = 0;
    double u1 = (double)(k1 & 0xffffff) / (double)0x1000000 + (0.5 / (double)(0x1000000));
    double u2 = (double)(k2 & 0xffffff) / (double)0x1000000 + (0.5 / (double)(0x1000000));
    auto mag = lnsigma * sqrt(-2.0 * log(u1));
    auto z0  = mag * cos(two_pi * u2) + mu;
    auto z1  = mag * sin(two_pi * u2) + mu;
    double odds_ratio = (allelefrac + 1e-9) / (1.0 - allelefrac + 1e-9) * exp(z0); // ((altfrac + 1e-9) / (reffrac + 1e-9));
    return odds_ratio / (1.0 + odds_ratio);
}

struct _RevComplement {
    char data[256];
    _RevComplement() {
        for (int i = 0; i < 256; i++) {
            data[i] = (char)i;
        }
        data['A'] = 'T';
        data['T'] = 'A';
        data['C'] = 'G';
        data['G'] = 'C';
        data['a'] = 't';
        data['t'] = 'a';
        data['c'] = 'g';
        data['g'] = 'c';
    }
};

const _RevComplement THE_REV_COMPLEMENT;

void complement(std::string & astring) {
    for (int i = 0; i < astring.size(); i++) {
        auto a = astring[i];
        auto b = THE_REV_COMPLEMENT.data[a];
        astring[i] = b;
    }
}

void reverse(std::string & astring) {
    size_t sz = astring.size(); 
    for (int i = 0; i < sz / 2; i++) {
        auto tmp = astring[i];
        astring[i] = astring[sz-1-i];
        astring[sz-1-i] = tmp;
    }
}

int is_var1_before_var2(int tid1, int pos1, int tid2, int pos2) {
    return (tid1 < tid2) || (tid1 == tid2 && pos1 < pos2);
}

int bamrec_write_fastq_raw(const bam1_t *aln, gzFile &outfile) {
    std::string seq;
    std::string qual;
    
    gzprintf(outfile, "@%s\n", bam_get_qname(aln));
    
    for (int i = 0; i < aln->core.l_qseq; i++) {
        char c = seq_nt16_str[bam_seqi(bam_get_seq(aln), i)];
        seq.push_back(c);
        assert(c > ' ' || !fprintf(stderr, "The read with qname %s is invalid with qlen = %d!\n", bam_get_qname(aln), aln->core.l_qseq));
    }
    if ((aln->core.flag & 0x10)) { reverse(seq); complement(seq); }
    gzprintf(outfile, "%s\n+\n", seq.c_str());
    for (int i = 0; i < aln->core.l_qseq; i++) {
        char c = bam_get_qual(aln)[i];
        qual.push_back(c+33);
    }
    if ((aln->core.flag & 0x10)) { reverse(qual); }
    gzprintf(outfile, "%s\n", qual.c_str());
}


int bamrec_write_fastq(const bam1_t *aln, std::string &seq, std::string &qual, gzFile &outfile) {
    
    assert(qual.size() == seq.size());
    for (int i = 0; i < qual.size(); i++) {
        qual[i] = (char)(qual[i] + 33);
    }
    gzprintf(outfile, "@%s\n", bam_get_qname(aln));
    if ((aln->core.flag & 0x10)) { reverse(seq); complement(seq); }
    for (int i = 0; i < seq.size(); i++) {
        assert(seq[i] > ' ' || !fprintf(stderr, "The read with qname %s is invalid!\n", bam_get_qname(aln)));
    }
    gzprintf(outfile, "%s\n+\n", seq.c_str());
    if ((aln->core.flag & 0x10)) { reverse(qual); }
    gzprintf(outfile, "%s\n", qual.c_str());
}

void help(int argc, char **argv) {
    fprintf(stderr, "Program %s version %s (%s)\n", argv[0], FULL_VERSION, COMMIT_DIFF_SH);
    fprintf(stderr, "  This is a NGS variant simulator that is aware of the molecular-barcodes (also known as unique molecular identifiers (UMIs))\n");
    
    fprintf(stderr, "Usage: %s -b <INPUT-BAM> -v <INPUT-VCF> -1 <OUTPUT-R1-FASTQ> -2 <OUTPUT-R2-FASTQ.gz> -0 <OUTPUT-UNPAIRED-FASTQ.GZ>\n", argv[0]);
    fprintf(stderr, "Optional parameters:\n");
    fprintf(stderr, " -f Fraction of variant allele (FA) to simulate. "
            "This value is overriden by the INFO/FA tag (specified by the -F command-line parameter) in the INPUT-VCF. "
            "Please note that INFO/FA must be defined the header of INPUT-VCF to be effective. "
            "Otherwise, the value defined by -f is used in the simulation [default to %f].\n", DEFAULT_ALLELE_FRAC);
    fprintf(stderr, " -p The power-law exponent simulating the over-dispersion of allele fractions in NGS [default to %f] (https://doi.org/10.1093/bib/bbab458). Negative value means that no over-dispersion is simulated. \n", DEFAULT_POWER_LAW_EXPONENT);
    fprintf(stderr, " -q the log-normal over-dispersion parameter in Phred scale [default to %f] (https://doi.org/10.1093/bib/bbab458). Negative value means that no over-dispersion is simulated. \n", DEFAULT_LOGNORMAL_DISP);
    fprintf(stderr, " -s The random seed used to simulate allele fractions from read names labeled with UMIs [default to %u].\n", DEFAULT_RANDSEED);

    
    fprintf(stderr, " -x Phred-scale sequencing error rates of simulated SNV variants "
            "where -2 means zero error and -1 means using sequencer BQ [default to %d].\n", DEFAULT_SNV_BQ_PHRED);
    
    fprintf(stderr, " -i The base quality of the inserted bases in the simulated insertion variants. "
            "[default to %d].\n", DEFAULT_INS_BQ_PHRED);
    fprintf(stderr, " -A The number of reads used to generate the randomness for simulating the nominator of the allele fraction used with the -p cmd-line param [default to %u].\n", DEFAULT_NITERS1);
    fprintf(stderr, " -B The number of reads used to generate the randomness for simulating the denominator of the allele fraction used with the -p cmd-line param [default to %u].\n", DEFAULT_NITERS2);
    fprintf(stderr, " -C The random seed used to simulate basecalling error [default to %u].\n", DEFAULT_RANDSEED);
    fprintf(stderr, " -F allele fraction TAG in the VCF file. " "[default to FA].\n");
    fprintf(stderr, " -S sample name used for the -F command-line parameter. "
                    "The special values NULL pointer, empty-string, and INFO mean using the INFO column instead of the FORMAT column." "[default to NULL pointer].\n");
    
    fprintf(stderr, "Note:\n");
    fprintf(stderr, "Reads in <OUTPUT-R1-FASTQ> and <OUTPUT-R2-FASTQ> are not in the same order, so these output FASTQ files have to be sorted using a tool such as fastq-sort before being aligned again, as most aligners such as BWA and Bowtie2 require reads in the R1 and R2 files to be in the same order (This is VERY IMPORTANT!).\n");
    fprintf(stderr, "<INPUT-BAM> and <INPUT-VCF> both have to be sorted and indexed.\n");
    fprintf(stderr, "To detect UMI, this prgram first checks for the MI tag in each alignment record in <INPUT-BAM>. If the MI tag is absent, then the program checks for the string after the number-hash-pound sign (#) in the read name (QNAME).\n");
    fprintf(stderr, "Each variant record in the INPUT-VCF needs to have only one variant, it cannot be multiallelic.\n");
    fprintf(stderr, "Currently, the simulation of insertion/deletion variants causes longer/shorter-than-expected lengths of read template sequences due to preservation of alignment start and end positions on the reference genome.\n");
    fprintf(stderr, "The symbol '#' denotes the start of UMI sequence so that any string before the '#' symbol is discarded.\n");
    fprintf(stderr, "The symbol '+' in a UMI sequence means that the UMI is a duplex, so the substrings before/after the '+' symbol are respectively the alpha/beta tags.\n");
    fprintf(stderr, "The BAM tag MI has special meaning as mentioned in the BAM file format specification. "
           "Therefore, for each BAM record, this program first searches for the MI tag. If the MI tag is not found, then this program uses the read name QNAME as the string containing UMI sequence\n");

    exit(-1);
}

int 
main(int argc, char **argv) {
    
    int flags, opt;
    char *inbam = NULL;
    char *invcf = NULL;
    char *r0outfq = NULL;
    char *r1outfq = NULL;
    char *r2outfq = NULL;
    double defallelefrac = DEFAULT_ALLELE_FRAC;
    int snv_bq_phred = DEFAULT_SNV_BQ_PHRED;
    int ins_bq_phred = DEFAULT_INS_BQ_PHRED;
    uint32_t randseed = DEFAULT_RANDSEED;
    uint32_t rand_niters1 = DEFAULT_NITERS1;
    uint32_t rand_niters2 = DEFAULT_NITERS2;
    uint32_t randseed_basecall = DEFAULT_RANDSEED;
    const char *tagFA = TAG_FA;
    const char *tagsample = NULL;
    bool is_always_log = false;
    double powerlaw_exponent = DEFAULT_POWER_LAW_EXPONENT;
    double lognormal_disp = DEFAULT_LOGNORMAL_DISP;
    while ((opt = getopt(argc, argv, "b:v:1:2:0:f:i:p:q:s:x:A:B:F:S:L")) != -1) {
        switch (opt) {
            case 'b': inbam = optarg; break;
            case 's': randseed = atoi(optarg); break;
            case 'A': rand_niters1 = atoi(optarg); break;
            case 'B': rand_niters2 = atoi(optarg); break;
            case 'C': randseed_basecall = atoi(optarg); break;
            case 'v': invcf = optarg; break;
            case '0': r0outfq = optarg; break;
            case '1': r1outfq = optarg; break;
            case '2': r2outfq = optarg; break;
            case 'f': defallelefrac = atof(optarg); break;
            case 'i': ins_bq_phred = atof(optarg); break;
            case 'p': powerlaw_exponent = atof(optarg); break;
            case 'q': lognormal_disp = atof(optarg); break;
            case 'x': snv_bq_phred = atof(optarg); break;
            case 'F': tagFA = optarg; break;
            case 'L': is_always_log = true; break;
            case 'S': tagsample = optarg; break;
            default: help(argc, argv);
        }
    }
    if (NULL == inbam || NULL == invcf) {
        fprintf(stderr, "The input BAM and VCF filenames have to be specified on the command line\n");
        help(argc, argv);
    }
    if ((NULL == r0outfq) && (NULL == r1outfq) && (NULL == r2outfq)) {
        fprintf(stderr, "At least one output FASTQ file has to be specified on the command line\n");
        help(argc, argv);
    }
    double lnfrac = pow(10.0, -lognormal_disp / 10.0);
    double lnsigma = log(2.0) / sqrt(log(lnfrac) / (-1.0/2.0));
    lnsigma = lnsigma / sqrt(2.0); // transform obs-to-obs var to obs-to-exp var
    
    if (0 != randseed) {
        randseed  = portable_int2randint(randseed , 1);
    }
    if (0 != randseed_basecall) {
        randseed_basecall = portable_int2randint(randseed_basecall, 4);
    }
    const bool is_FA_from_INFO = ((tagsample == NULL) || (0 == strlen(tagsample)) || !strcmp("INFO", tagsample));
    fprintf(stderr, "%s\n=== version ===\n%s\n%s\n%s\n", argv[0], FULL_VERSION, COMMIT_DIFF_SH, GIT_DIFF_FULL);
    fprintf(stderr, "lnsigma = %f\n", lnsigma);
    
    int64_t num_kept_reads = 0;
    int64_t num_kept_snv = 0;
    int64_t num_kept_mnv = 0;
    int64_t num_kept_ins = 0;
    int64_t num_kept_del = 0;
    int64_t num_skip_reads = 0;
    int64_t num_skip_cmatches = 0;
    
    htsFile *vcf_fp = vcf_open(invcf, "r");
    bcf_hdr_t *vcf_hdr = bcf_hdr_read(vcf_fp);
    int tag_sample_idx = bcf_hdr_nsamples(vcf_hdr) - 1; 
    for (int sidx = 0; sidx < bcf_hdr_nsamples(vcf_hdr); sidx++) {
        if ((tagsample != NULL) && !strcmp(tagsample, vcf_hdr->samples[sidx])) {
            tag_sample_idx = sidx; 
        }
    }
    bcf1_t *vcf_rec = bcf_init();
    vcf_rec->rid = -1;
    vcf_rec->pos = 0;
    samFile *bam_fp = sam_open(inbam, "r");
    sam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);
    bam1_t *bam_rec1 = bam_init1();
    
    gzFile r0file = gzopen(r0outfq, "wb1");
    gzFile r1file = gzopen(r1outfq, "wb1");
    gzFile r2file = gzopen(r2outfq, "wb1");
    
    std::deque<bcf1_t*> vcf_list;
    
    float *bcffloats = NULL; // [256] = {0}; 
    int vcf_read_ret = 0;
    
    std::string newseq, newqual;
    
    samFile *bam_fp2 = sam_open(inbam, "r");
    sam_hdr_t *bam_hdr2 = sam_hdr_read(bam_fp2);
    uint64_t read_cnt = 0;
    uint32_t samplehash1 = 0;
    uint32_t samplehash2 = 0;
    while (sam_read1(bam_fp2, bam_hdr, bam_rec1) >= 0 && read_cnt < MAX(rand_niters1, rand_niters2)) {
        if (read_cnt < rand_niters1) { 
            samplehash1 += __ac_X31_hash_string(bam_get_qname(bam_rec1));
        }
        if (read_cnt < rand_niters2) { 
            samplehash2 += __ac_X31_hash_string(bam_get_qname(bam_rec1));
        }
        read_cnt += 1;
    }
    fprintf(stderr, "The value samplehash1 and samplehash2 are %u and %u, read_cnt = %u\n", samplehash1, samplehash2, read_cnt);
    bam_hdr_destroy(bam_hdr2);
    sam_close(bam_fp2);
    
    while (sam_read1(bam_fp, bam_hdr, bam_rec1) >= 0) {
        const bam1_t *bam_rec = bam_rec1;
        if ((0 != (bam_rec->core.flag & 0x900)) && (0 == (bam_rec->core.flag & 0x4))) { continue; }
        auto &outfile = ((bam_rec->core.flag & 0x40) ? r1file : ((bam_rec->core.flag & 0x80) ? r2file : r0file));
        const char *outfname = ((bam_rec->core.flag & 0x40) ? r1outfq : ((bam_rec->core.flag & 0x80) ? r2outfq : r0outfq));
        if ((0 != (bam_rec->core.flag & 0x4))) {
            int write_ret = bamrec_write_fastq_raw(bam_rec1, outfile);
            continue;
        }
        while (1) {
            if (vcf_read_ret != -1) {
                fprintf(stderr, "The variant at tid %d pos %d is before the read at tid %d pos %d, readname = %s\n", 
                    vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_rec->core.pos, bam_get_qname(bam_rec));
                vcf_read_ret = vcf_read(vcf_fp, vcf_hdr, vcf_rec); // skip this variant
                fprintf(stderr, "The new prep variant is at tid %d pos %d\n", 
                    vcf_rec->rid, vcf_rec->pos);
            }
            if (vcf_read_ret < 0) { break; }
            if (!is_var1_before_var2(vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_rec->core.pos)) {
                vcf_list.push_back(bcf_dup(vcf_rec));  
                break;
            }
        }
        while (is_var1_before_var2(vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_endpos(bam_rec))) {
            if (vcf_read_ret != -1) {
                fprintf(stderr, "The variant at tid %d pos %d is before the read at tid %d endpos %d, readname = %s\n", 
                    vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_endpos(bam_rec), bam_get_qname(bam_rec));
                vcf_read_ret = vcf_read(vcf_fp, vcf_hdr, vcf_rec); // get this variant
                fprintf(stderr, "The new pushed variant is at tid %d pos %d\n", 
                    vcf_rec->rid, vcf_rec->pos);
            }
            if (vcf_read_ret < 0) { break; }
            vcf_list.push_back(bcf_dup(vcf_rec));
        }
        while (vcf_list.size() > 0 && is_var1_before_var2(vcf_list.front()->rid, vcf_list.front()->pos, bam_rec->core.tid, bam_rec->core.pos)) {
            fprintf(stderr, "The variant at tid %d pos %d is destroyed\n", 
                    vcf_list.front()->rid, vcf_list.front()->pos);
            bcf_destroy(vcf_list.front());
            vcf_list.pop_front();
        }
        
        const auto *seq = bam_get_seq(bam_rec);
        const auto *qual = bam_get_qual(bam_rec); 
        newseq.clear();
        newqual.clear();
        
        uint32_t umihash = 0;
        const auto *bam_aux_data = bam_aux_get(bam_rec, "MI"); // this tag is reserved (https://samtools.github.io/hts-specs/SAMtags.pdf)
        const char *umistr = ((bam_aux_data != NULL) ? bam_aux2Z(bam_aux_data) : bam_get_qname(bam_rec));
        
        // TODO: this assumes 100% duplex forming efficiency, which is not what happens in practice. 
        // TODO: introduce another parameter to simulate the efficiency of duplex formation?
        uint32_t begpos = MIN(bam_rec->core.pos, bam_rec->core.mpos);
        const double mutprob = umistr2prob(umihash, randseed, begpos, begpos + abs(bam_rec->core.isize), umistr);
        if (vcf_list.size() > 0) {
            int qpos = 0;
            int rpos = bam_rec->core.pos;
            auto vcf_rec_it = vcf_list.begin();
            int vcf_list_idx = 0;
            for (int i = 0; i < bam_rec->core.n_cigar; i++) {
                const auto c = bam_get_cigar(bam_rec)[i];
                const auto cigar_op = bam_cigar_op(c);
                const auto cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    auto cigar_oplen1 = cigar_oplen;
                    for (int j = 0; j < cigar_oplen1; j++) {
                        while (vcf_rec_it != vcf_list.end() && is_var1_before_var2((*vcf_rec_it)->rid, (*vcf_rec_it)->pos, bam_rec->core.tid, rpos)) {
                            vcf_rec_it++;
                            vcf_list_idx++;
                        }
                        if (vcf_rec_it != vcf_list.end() && rpos == (*vcf_rec_it)->pos) {
                            auto vcf_rec_it_end = vcf_rec_it;
                            auto vcf_list_idx_end = vcf_list_idx;
                            while (vcf_rec_it_end != vcf_list.end() && ((*vcf_rec_it_end)->rid == (*vcf_rec_it)->rid && (*vcf_rec_it_end)->pos == (*vcf_rec_it)->pos)) {
                                vcf_rec_it_end++;
                                vcf_list_idx_end++;
                            }
                            double allelefrac = (double)0;
                            bool is_mutated = false;
for (auto vcf_rec_it2 = vcf_rec_it; vcf_rec_it2 != vcf_rec_it_end; vcf_rec_it2++) {
                            const auto & vcf_rec = *vcf_rec_it2;
                            bcf_unpack(vcf_rec, BCF_UN_ALL);
                            int ndst_val = 0;
                            int valsize = 0;
                            if (is_FA_from_INFO) {
                                valsize = bcf_get_info_float(vcf_hdr, vcf_rec, tagFA, &bcffloats, &ndst_val);
                                allelefrac += (valsize > 0 ? bcffloats[valsize - 1] : defallelefrac);
                            } else {
                                valsize = bcf_get_format_float(vcf_hdr, vcf_rec, tagFA, &bcffloats, &ndst_val);
                                allelefrac += ((valsize > 0 && valsize == bcf_hdr_nsamples(vcf_hdr)) ? bcffloats[tag_sample_idx] : defallelefrac);
                            }
                            
                            double allelefrac2 = allelefrac;
                            if (powerlaw_exponent > 0) {
                                allelefrac2 = allelefrac_powlaw_transform(
                                    allelefrac,
                                    (uint32_t)(vcf_rec->rid),
                                    (uint32_t)(rpos),
                                    (uint32_t)samplehash1,
                                    (uint32_t)samplehash2,
                                    powerlaw_exponent);
                            }
                            double allelefrac3 = allelefrac2;
                            if (lognormal_disp > 0) {
                                allelefrac3 = allelefrac_lognormal_transform(
                                    allelefrac2,
                                    (uint32_t)(vcf_rec->rid),
                                    (uint32_t)(rpos),
                                    (uint32_t)samplehash1,
                                    (uint32_t)samplehash2,
                                    lnsigma);
                            }
                            if (mutprob <= allelefrac3) {
                                const char *newref = vcf_rec->d.allele[0];
                                const char *newalt = vcf_rec->d.allele[1];
                                if (1 == strlen(newref) && 1 == strlen(newalt)) {
                                    uint32_t hash = 0;
                                    double randprob = qnameqpos2prob(hash, randseed_basecall, bam_get_qname(bam_rec), qpos);
                                    const char base = ((-2 == snv_bq_phred)
                                            || (prob2phred(randprob) < (-1 == snv_bq_phred ? qual[qpos] : snv_bq_phred)) 
                                            ? (newalt[0]) : (ACGT[hash % 4]));
                                    newseq.push_back(base);
                                    newqual.push_back(qual[qpos]);
                                    num_kept_snv++;
                                    if (ispowerof2(num_kept_snv) || is_always_log) { 
                                        fprintf(stderr, "The read with name %s is spiked with the snv-variant at tid %d pos %d, FAs = %f,%f,%f\n", 
                                                bam_get_qname(bam_rec), vcf_rec->rid, vcf_rec->pos, allelefrac, allelefrac2, allelefrac3);
                                    }
                                } else if (strlen(newref) == strlen(newalt)) {
                                    fprintf(stderr, "Warning: the MNV at tid %d pos %d is decomposed into SNV and only the first SNV is simulated\n", 
                                            bam_rec->core.tid, bam_rec->core.pos);
                                    uint32_t hash = 0;
                                    double randprob = qnameqpos2prob(hash, randseed_basecall, bam_get_qname(bam_rec), qpos);
                                    const char base = ((-2 == snv_bq_phred)
                                            || (prob2phred(randprob) < (-1 == snv_bq_phred ? qual[qpos] : snv_bq_phred)) 
                                            ? (newalt[0]) : (ACGT[hash % 4]));
                                    newseq.push_back(base);
                                    newqual.push_back(qual[qpos]);
                                    num_kept_mnv++;
                                } else if (strlen(newref) == 1 && strlen(newalt) >  1) {
                                    for (int k = 0; k < strlen(newalt); k++) { 
                                        newseq.push_back(newalt[k]);
                                    }
                                    newqual.push_back(qual[qpos]);
                                    for (int k = 1; k < strlen(newalt); k++) { 
                                        newqual.push_back((char)(ins_bq_phred)); 
                                    }
                                    num_kept_ins++;
                                    if (ispowerof2(num_kept_ins) || is_always_log) { fprintf(stderr, "The read with name %s is spiked with the ins-variant at tid %d pos %d\n", 
                                                bam_get_qname(bam_rec), vcf_rec->rid, vcf_rec->pos); }
                                } else if (strlen(newref) >  1 && strlen(newalt) == 1) {
                                    if (strlen(newref) + j < cigar_oplen1) {
                                        const char nuc = seq_nt16_str[bam_seqi(seq, qpos)];
                                        newseq.push_back(nuc);
                                        newqual.push_back(qual[qpos]);
                                        j += strlen(newref) - 1;
                                        qpos += strlen(newref) - 1;
                                        rpos += strlen(newref) - 1;
                                        num_kept_del++;
                                        if (ispowerof2(num_kept_del) || is_always_log) { fprintf(stderr, "The read with name %s is spiked with the del-variant at tid %d pos %d\n", 
                                                bam_get_qname(bam_rec), vcf_rec->rid, vcf_rec->pos); }
                                    } else {
                                        newseq.push_back(seq_nt16_str[bam_seqi(seq, qpos)]);
                                        newqual.push_back(qual[qpos]);
                                    }
                                } else {
                                    fprintf(stderr, "The variant at tid %d pos %d failed to be processed!\n", bam_rec->core.tid, bam_rec->core.pos);
                                }
                                num_kept_reads++;
                                is_mutated = true;
                                break;
                            } else {
                                if (ispowerof2(num_skip_reads) || is_always_log) {
                                    fprintf(stderr, "The read with name %s is not affected by the variant at tid %d pos %d\n", 
                                    bam_get_qname(bam_rec), vcf_rec->rid, vcf_rec->pos); 
                                }
                            }
}
                            if (!is_mutated) {
                                const char nuc = seq_nt16_str[bam_seqi(seq, qpos)];
                                newseq.push_back(nuc);
                                newqual.push_back(qual[qpos]);
                                num_skip_reads++;
                            }
                        } else {
                            const char nuc = seq_nt16_str[bam_seqi(seq, qpos)];
                            newseq.push_back(nuc);
                            newqual.push_back(qual[qpos]);
                            num_skip_cmatches++;
                        }
                        qpos++;
                        rpos++;
                    }
                } else if ((cigar_op == BAM_CINS) || (cigar_op == BAM_CSOFT_CLIP)) {
                    for (int j = 0; j < cigar_oplen; j++) {
                        const char nuc = seq_nt16_str[bam_seqi(seq, qpos)];
                        newseq.push_back(seq_nt16_str[bam_seqi(seq, qpos)]);
                        newqual.push_back(qual[qpos]);
                        qpos++;
                    }
                } else if (cigar_op == BAM_CDEL) {
                    rpos += cigar_oplen;
                } else if (cigar_op == BAM_CHARD_CLIP) {
                    // fall through
                } else {
                    fprintf(stderr, "The cigar code %d is invalid at tid %d pos %d for read %s!\n", 
                            cigar_op, bam_rec->core.tid, bam_rec->core.pos, bam_get_qname(bam_rec));
                    abort();
                }
            }
            if (outfname != NULL) { bamrec_write_fastq(bam_rec, newseq, newqual, outfile); }
        } else {
            if (outfname != NULL) { bamrec_write_fastq_raw(bam_rec, outfile); } 
        }
    }
    if (bcffloats != NULL) {
        free(bcffloats);
    }
    
    bam_destroy1(bam_rec1);
    bam_hdr_destroy(bam_hdr);
    sam_close(bam_fp);
    
    bcf_destroy1(vcf_rec);
    bcf_hdr_destroy(vcf_hdr);
    vcf_close(vcf_fp);
    
    gzclose(r0file);
    gzclose(r1file);
    gzclose(r2file);
    
    fprintf(stderr, "In total: kept %ld read support, skipped %ld read support"
            ", and skipped %ld no-variant CMATCH cigars.\n", num_kept_reads, num_skip_reads, num_skip_cmatches);
    fprintf(stderr, "Kept %d snv read support\n", num_kept_snv);
    fprintf(stderr, "Kept %d mnv read support\n", num_kept_mnv);
    fprintf(stderr, "Kept %d insertion read support\n", num_kept_ins);
    fprintf(stderr, "Kept %d deletion read support\n", num_kept_del);
}

