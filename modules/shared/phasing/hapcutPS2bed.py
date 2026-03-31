"""
convert hapcut2 phased block to .bed file
usage:
hapblock_dir=Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/
python $src --input_dir ${hapblock_dir}
"""
import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", type=str, required=True)
# parser.add_argument("--output_dir", type=str, required=False)
args = parser.parse_args()
_dir = args.input_dir
# output_dir = args.output_dir

def hapcutPS2bed(hapblock_file, out_bed):
    f=open(out_bed, "w")
    blocklist = []
    with open(hapblock_file, 'r') as hbf:
        ps = []
        for line in hbf:
            if len(line) < 3: # empty line
                continue
            if 'BLOCK' in line:
                blocklist.append([])
                continue

            el = line.strip().split('\t')
            if len(el) < 3: # not enough elements to have a haplotype
                continue

            if len(el) >= 5:
                chrom = el[3]
                pos = int(el[4])-1
            blocklist[-1].append((chrom, pos))

    for blk in blocklist:
        first_pos  = -1
        last_pos   = -1

        for chrom, pos in blk:
            if first_pos == -1:
                first_pos = pos
            last_pos = pos
        f.write(f'{chrom}\t{first_pos}\t{last_pos}\n')

    f.close()
    return 

if __name__ == "__main__":
    # _dir='Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/'
    # _dir='./V350158928_L01_L02_V350159512_L03/Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/'
    # name='V350158928_L01_L02_V350159512_L03'
    cmd=f'mkdir -p {_dir}/bed'
    subprocess.call(cmd, shell=True)

    for i in range(1,23):
        hapblock_file=f'{_dir}/data_hapblock_chr{i}'
        out_bed=f'{_dir}/bed/data_hapblock_chr{i}.bed'
        hapcutPS2bed(hapblock_file, out_bed)
    try:
        i='X'
        hapblock_file=f'{_dir}/data_hapblock_chr{i}'
        out_bed=f'{_dir}/bed/data_hapblock_chr{i}.bed'
        hapcutPS2bed(hapblock_file, out_bed)
    except:
        print('XY sample, no chrX phased') 

    cmd = f'cat {_dir}/bed/*.bed > {_dir}/data_hapblock.bed'
    subprocess.call(cmd, shell=True)
