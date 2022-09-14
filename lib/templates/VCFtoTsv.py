# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
from dataclasses import field
import os
from stat import FILE_ATTRIBUTE_HIDDEN
import sys

tsv = sys.argv[1]

# open the vcf
with open (tsv, "r") as fh:
    print (f'Gene\tVariant type\tChrom\tPos\tRSID\tREF\tALT\tFilter\tDP\tVD\tAF\tMMQ\tGERMQ\tMBQ\tStrand bias\tF1R2\tF2R1\tROQ\tMFRL\tMPOS\tRPA\tRU\tSTRQ\tECNT\tTLOD\tCONTQ\tNALOD\tNCOUNT\tNLOD\tOCM\tSEQQ\tSTR\tSTRANDQ\tClass\tRefSeq\tHGVS_c\tHGVS_p\tANN\tGT')
    i = 1
    for line in fh:
        if not line.startswith("ANN") and not line.startswith("."):
            fields = line.strip().split('\t')
            gene = fields[0]
            variant_type =fields[1]
            chrom = fields[2]
            pos = fields[3]
            rsid = fields[4]
            ref = fields[5]
            alt = fields[6]
            filter = fields[7]
            depth = fields[8]
            var_depth = fields[9]
            af = float(fields[10])
            mmq = fields[11]
            germq = fields[12]
            mbq = fields[13]
            sb_table = fields[14]
            f1r2 = fields[15]
            f2r1 = fields[16]
            roq = fields[17]
            mfrl = fields[18]
            mpos = fields[19]
            rpa = fields[20]
            ru = fields[21]
            strq = fields[22]
            ecnt = fields[23]
            tlod = fields[24]
            contq = fields[25]
            nalod = fields[26]
            ncount = fields[27]
            nlod = fields[28]
            ocm = fields[29]
            seqq = fields[30]
            str = fields[31]
            strandq = fields[32]
            var_class = fields [33]
            refseq_id = fields[34]
            hgvs_c = fields[35]
            hgvs_p = fields[36]
            ann = fields[37]
            gt = fields[38]

            print (f'{gene}\t{variant_type}\t{chrom}\t {pos}\t{rsid}\t{ref}\t{alt}\t{filter}\t{depth}\t{var_depth}\t{af}\t{mmq}\t{germq}\t{mbq}\t{sb_table}\t{f1r2}\t{f2r1}\t{roq}\t{mfrl}\t{mpos}\t{rpa}\t{ru}\t{strq}\t{ecnt}\t{tlod}\t{contq}\t{nalod}\t{ncount}\t{nlod}\t{ocm}\t{seqq}\t{str}\t{strandq}\t{var_class}\t{refseq_id}\t{hgvs_c}\t{hgvs_p}\t{ann}\t{gt}')

