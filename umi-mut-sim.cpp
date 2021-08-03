#include <deque>
#include <string>
#include <vector>

#include <stdio.h>

#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

double umistr2prob(const char *str) {
    const char *umistr1 = strchr(str, '#');
    const char *umistr = (NULL == umistr1 ? str : umistr1);
    uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(umistr) ^ 0);
    return (double)(k&0xffffff) / 0x1000000;
}

struct _RevComplement {
    char data[128];
    _RevComplement() {
        for (int i = 0; i < 128; i++) {
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
    for (int i = 0; i < astring.size() / 2; i++) {
        auto tmp = astring[i];
        astring[i] = astring[astring.size()-1-i];
        astring[astring.size()-1-i] = tmp;
    }
}

int is_var1_before_var2(int tid1, int pos1, int tid2, int pos2) {
    return (tid1 < tid2) || (tid1 == tid2 && pos1 < pos2);
}

int bamrec_write_fastq(bam1_t *aln, FILE *r1file, FILE *r2file) {
    FILE *outfile = ((aln->core.flag & 0x40) ? r1file : ((aln->core.flag & 0x80) ? r2file : r1file));
    std::string seq;
    std::string qual;
    
    fprintf(outfile, "@%s\n", bam_get_qname(aln));
    for (int i = 0; i < aln->core.l_qseq; i++) {
        char c = seq_nt16_str[bam_seqi(bam_get_seq(aln), i)];
        seq.push_back(c);
        //sprintf(seq_str, "%c", c);
    }
    if ((aln->core.flag & 0x10)) { reverse(seq); complement(seq); }
    fprintf(outfile, "%s\n+\n", seq.c_str());
    for (int i = 0; i < aln->core.l_qseq; i++) {
        char c = bam_get_qual(aln)[i];
        qual.push_back(c+33);
        // sprintf(qual_str, "%c", (char)(32+(int)c));
    }
    if ((aln->core.flag & 0x10)) { reverse(qual); }
    fprintf(outfile, "%s\n", qual.c_str());
}

#if 0 // dedup
read vcf
for each read in the sorted bam
    while current variant is before the begin position of this read: iterate to the next variant
    if the vcf variant is before the end position of this read:
        calculate the hash value of the UMI
        if the hash value divided by the max hash value is less than the variant allele fraction:
            modify this bam read with this variant
            // modify this bam read with the all the overlapping variants that are inclusively after this variant.
    push this bam read to the pile of bam reads
    for each read2 in the pile of bam reads:
        if read2 ends before the pos of this read:
            flush out read2
    // flush out this bam read
for each read2 in the pile of bam reads: flush out read2
#endif

int 
main(int argc, char **argv) {

    int64_t num_kept_reads = 0;
    int64_t num_skip_reads = 0;
    int64_t num_skip_cmatches = 0;

    htsFile *vcf_fp = vcf_open(argv[2], "r");
    bcf_hdr_t *vcf_hdr = bcf_hdr_read(vcf_fp);
    bcf1_t *vcf_rec = bcf_init();
    vcf_rec->rid = -1;
    vcf_rec->pos = 0;
    samFile *bam_fp = sam_open(argv[1], "r");
    sam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);
    bam1_t *bam_rec = bam_init1();
    
    FILE *r1file = fopen(argv[3], "w");
    FILE *r2file = fopen(argv[4], "w");
    
    std::deque<bcf1_t*> vcf_list;
    std::vector<bam1_t*> bam_list;
    
    float bcffloats[256] = {0};
    
    int vcf_read_ret = 0;
    while (sam_read1(bam_fp, bam_hdr, bam_rec) >= 0) {
        if (0 != (bam_rec->core.flag & 0x900)) { continue; }
        while (is_var1_before_var2(vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_rec->core.pos)) {
            if (vcf_read_ret != -1) {
                fprintf(stderr, "The variant at tid %d pos %d is before the read at tid %d pos %d, readname = %s\n", 
                    vcf_rec->rid, vcf_rec->pos, bam_rec->core.tid, bam_rec->core.pos, bam_get_qname(bam_rec));
                vcf_read_ret = vcf_read(vcf_fp, vcf_hdr, vcf_rec); // skip this variant
                fprintf(stderr, "The new prep variant is at tid %d pos %d\n", 
                    vcf_rec->rid, vcf_rec->pos);
            }
            if (vcf_read_ret < 0) { break; }
            
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
            bcf_unpack(vcf_rec, BCF_UN_ALL);
            vcf_list.push_back(bcf_dup(vcf_rec));
        }
        while (vcf_list.size() > 0 && is_var1_before_var2(vcf_list.front()->rid, vcf_list.front()->pos, bam_rec->core.tid, bam_rec->core.pos)) {
            fprintf(stderr, "The variant at tid %d pos %d is destroyed\n", 
                    vcf_list.front()->rid, vcf_list.front()->pos);
            bcf_destroy(vcf_list.front());
            vcf_list.pop_front();
        }
        
        const double mutprob = umistr2prob(bam_get_qname(bam_rec));
        if (vcf_list.size() > 0) {
            int qpos = 0;
            int rpos = bam_rec->core.pos;
            auto vcf_rec_it = vcf_list.begin();
            for (int i = 0; i < bam_rec->core.n_cigar; i++) {
                const auto c = bam_get_cigar(bam_rec)[i];
                const auto cigar_op = bam_cigar_op(c);
                const auto cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    while (vcf_rec_it != vcf_list.end() && is_var1_before_var2((*vcf_rec_it)->rid, (*vcf_rec_it)->pos, bam_rec->core.tid, rpos)) {
                        vcf_rec_it++;
                    }
                    if (vcf_rec_it != vcf_list.end() && rpos <= (*vcf_rec_it)->pos && (*vcf_rec_it)->pos < rpos + cigar_oplen) {
                        const auto & vcf_rec = *vcf_rec_it;
                        int ndst_val = 0;
                        int valsize = bcf_get_format_int32(vcf_hdr, vcf_rec, "FA", &bcffloats, &ndst_val);
                        const double allelefrac = (ndst_val > 0 ? bcffloats[ndst_val - 1] : 0.1);
                        if (mutprob <= allelefrac) {
                            // SNV only for now
                            auto *seq = bam_get_seq(bam_rec);
                            bcf_unpack(vcf_rec, BCF_UN_ALL);
                            const char newalt = vcf_rec->d.allele[1][0];
                            //fprintf(stderr, "Will pick the read %s at tid %d pos %d for mutation at tid %d pos %d to the base %c\n", 
                            //        bam_get_qname(bam_rec), bam_rec->core.tid, bam_rec->core.pos, vcf_rec->rid, vcf_rec->pos, newalt);
                            num_kept_reads++;
                            bam_set_seqi(seq, qpos + vcf_rec->pos - rpos, seq_nt16_table[newalt]);
                        } else {
                            //fprintf(stderr, "Will skip the read %s at tid %d pos %d for mutation at tid %d pos %d\n", 
                            //        bam_get_qname(bam_rec), bam_rec->core.tid, bam_rec->core.pos, vcf_rec->rid, vcf_rec->pos);
                            num_skip_reads++;
                        }
                    } else {
                        num_skip_cmatches++;
                        //fprintf(stderr, "Will no go to the read %s at tid %d pos %d for mutation at tid %d pos %d\n", 
                        //            bam_get_qname(bam_rec), bam_rec->core.tid, bam_rec->core.pos, vcf_rec->rid, vcf_rec->pos);
                    }
                    qpos += cigar_oplen;
                    rpos += cigar_oplen;
                } else if (cigar_op == BAM_CINS) {
                    qpos += cigar_oplen;
                } else if (cigar_op == BAM_CDEL) {
                    rpos += cigar_oplen;
                } else if (cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
                    qpos += cigar_oplen;
                } else {
                    fprintf(stderr, "The cigar code %d is invalid at tid %d pos %d for read %s!\n", 
                            cigar_op, bam_rec->core.tid, bam_rec->core.pos, bam_get_qname(bam_rec));
                    abort();
                }
            }
        }
        bamrec_write_fastq(bam_rec, r1file, r2file);
    }
    
    bam_destroy1(bam_rec);
    bam_hdr_destroy(bam_hdr);
    sam_close(bam_fp);
    
    bcf_destroy1(vcf_rec);
    bcf_hdr_destroy(vcf_hdr);
    vcf_close(vcf_fp);

    fprintf(stderr, "In total: kept %d reads and skipped %d reads, and skipped %d no-variant CMATCH cigars.\n", num_kept_reads, num_skip_reads, num_skip_cmatches);
}

