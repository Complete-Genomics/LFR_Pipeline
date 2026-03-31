#coding=utf-8
from __future__ import print_function 
import sys
import argparse
example_text = '''example:
    python this.py -f xx.fa -m max_fa
'''
parser = argparse.ArgumentParser(description="To get the longest cdna sequence.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--fafile','-f',type=str,help="xx.fa file",required= True,metavar='')
parser.add_argument('--max_fa', '-m',type= str,help="output file",required= True,metavar='')
args = parser.parse_args()

fafile=args.fafile
max_fa=args.max_fa

def get_all_fa(fa_file):
    seq = {}
    name1=''
    for line in fa_file:
        if line.startswith('>'):
            name1 = line.strip()
            seq[name1] = ''
        else:
            seq[name1] += line
    return seq

def get_max_fa(seq_dict):
    maxseq = {}
    d=0
    for k1,v1 in seq_dict.items(): # k1 = name ; v1 =sequence
        d +=1
        if d == 1:
            maxseq[k1] = v1
            maxlen = len(v1)
            k_last = k1
        elif len(v1) > maxlen:
            maxseq[k1] = v1
            del maxseq[k_last]
            k_last = k1
    return maxseq

if fafile:    
    with open(fafile,'r') as fa_file:
        seq_dict = get_all_fa(fa_file)
        max_seq = get_max_fa(seq_dict)
        with open(max_fa,'w') as output:
            for name, seq in max_seq.items():
                print(name, seq, sep='\n',file = output)
else:
    print(f'denovo failed on {fafile}')
