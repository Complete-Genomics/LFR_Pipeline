"""
1/ get pos of tag (6bp intergenic BC seq, adapter)
2/ for softclip, filter reads with partial map, 
match on left/right of softclip,
output seq with >N seq mapped with softclip
usage:
python BCgDNA_filterCIGAR.py --bam_file Align/noall3.sort.bam
python ${src} --bam_r2_type ${bam_r2_type} > Align/${bam_r2_type}.sort.r2.bam.pos.txt
"""

import argparse
import pysam
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("--bam_r2_type", type=str, help='nobc1')
parser.add_argument("--bam_file", type=str)
parser.add_argument("--map_cutoff", type=int)
parser.add_argument("--softclip_cutoff", type=int)
parser.add_argument("--match_pos_softclip", type=str, help='left of right')

args = parser.parse_args() 
bam_file = args.bam_file
bam_r2_type = args.bam_r2_type
map_cutoff = args.map_cutoff
softclip_cutoff = args.softclip_cutoff
match_pos_softclip = args.match_pos_softclip

def filter_cigar(bam_file):
    outfile_name = bam_file.split(".bam")[0]+"_filtered.bam"
    outfile=pysam.AlignmentFile(outfile_name, "wb", template=samfile)
    filtered_reads = []
    bf = pysam.AlignmentFile(bam_file, 'rb')
    for read in bf:
        cigar = read.cigartuples
        m_flag, s_flag, pos_tag = False, False, False


        ## match on left/right of softclip
        cigar_pos =[c[0] for c in cigar]
        m_pos = cigar_pos.index(0)
        s_pos = cigar_pos.index(4)
        if match_pos_softclip=='left':
            if m_pos<s_pos:
                pos_tag = True
        elif match_pos_softclip=='right':
            if m_pos>s_pos:
                pos_tag = True
        else:
            print('match_pos_softclip error')

        for c in cigar:
            ## match
            if c[0]==0 and c[1]>=map_cutoff:
                m_flag=True
            ## softclip
            if c[0]==4 and c[1]>=soft_cutoff:
                s_flag=True
        if pos_tag and m_flag and s_flag:
            outfile.write(read)
    outfile.close()


def filter_pos(infile):
    readid, cigar, seq = [],[], []
    with open(infile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            readid.append(info[0])
            cigar.append(info[5])
            seq.append(info[9])

    ## anchor 6bp flanking region
    flank = ['TCTGCG', 'CCTTCC', 'TAGCATTTGTCTTCCTAAGA']
    
    for i in range(len(seq)):
        idx1, idx2, idx3=-2,-2,-2
        s = seq[i]

        bc3 = readid[i][-11:-1] 
        try:
            idx1= s.index(flank[0])
        except: pass
        try:
            idx2= s.index(flank[1])
        except: pass
        try:
            idx3= s.index(flank[2])
        except: pass
        print('on {}, pos of TCTGCG={}, CCTTCC={}, adapter_head20bp={}'.format(s, idx1+1, idx2+1, idx3+1))

if __name__ == "__main__":

    bam_r2_txt='Align/{}.sort.r2.bam.txt'.format(bam_r2_type)
    filter_pos(bam_r2_txt)

