import sys
import os
import matplotlib
matplotlib.use('Agg')

prog_version = '1.0'
prog_date = '2022-08'
usage = '''

     Version %s by Christian Villarosa  %s
     Usage: %s read_file1 [read_file2] [-r Reference.fa] [-o output/dir/path]

''' % (prog_version, prog_date, os.path.basename(sys.argv[0]))

def float_check(float_string):
    try:
        float(float_string)
        return True
    except:
        return False

class GC_Bias(object):
    only_mapped = True

    def __init__(self, data_files, reference, outdir, r1_len, r2_len):
        self.data_files = data_files
        self.reference = reference
        self.outdir = outdir
        print(self.data_files, self.reference, self.outdir)
        if r1_len==0:
            self.read_lengths = [r2_len]
        else:
            self.read_lengths = [r1_len, r2_len]
        return

    def get_fastq_read(self):
        self.rf.readline()
        read = self.rf.readline().strip('\r\n').upper()
        self.rf.readline()
        self.rf.readline()
        return read

    def get_sam_read(self, line=''):
        if not line:
            line = self.rf.readline().strip('\r\n')
            if not line:
                return ''
        info = line.split()
        unmapped = int(info[1]) & 4 and self.only_mapped
        while unmapped:
            line = self.rf.readline().strip('\r\n')
            if not line:
                return ''
            info = line.split()
            unmapped = int(info[1]) & 4 and self.only_mapped
        read = info[9].upper()
        # if info[0][-2] == '/':
        #     self.pe_tag = info[0][-1]
        if int(info[1]) & 64:
            self.pe_tag = '1'
        elif int(info[1]) & 128:
            self.pe_tag = '2'
        return read

    def load_data_file(self, data_fp):
        if data_fp.endswith('.gz'):
            import gzip
            self.file_type = 'fastq'
            self.rf = gzip.open(data_fp, 'rb')
            self.queued_read = self.get_fastq_read()
        else:
            self.file_type = 'sam'
            self.pe_tag = False
            self.rf = open(data_fp, 'r')
            line = self.rf.readline().strip('\r\n')
            while line.startswith('@'):
                line = self.rf.readline().strip('\r\n')
            self.queued_read = self.get_sam_read(line)
        # self.read_lengths.append(len(self.queued_read))
        print("read_lengths: {}".format(self.read_lengths))
        return

    def get_next_read(self):
        read = self.queued_read
        if self.file_type == 'fastq':
            self.queued_read = self.get_fastq_read()
        else:
            if self.pe_tag:
                self.rix = 0 if self.pe_tag == '1' else 1
            self.queued_read = self.get_sam_read()
            if self.pe_tag == '2' and len(self.read_lengths) < 2:
                self.read_lengths.append(len(self.queued_read))
        return read

    def count_GC(self):
        gc_counts = {0:{}, 1:{}}
        for self.rix, data_file in enumerate(self.data_files):
            self.load_data_file(data_file)
            while True:
                read = self.get_next_read()
                if not read: break
                gc_count = read.count('G') + read.count('C')
                if gc_count not in gc_counts[self.rix]:
                    gc_counts[self.rix][gc_count] = 0
                gc_counts[self.rix][gc_count] += 1
        return gc_counts

    def count2bins_(self, gc_counts_, read_length):
        # bin = 1 GC count = counts index
        counts = []
        for rl in range(read_length + 1):
            count = gc_counts_.pop(rl, 0)
            counts.append(count)
        return counts

    def count2bins(self, gc_counts):
        return [self.count2bins_(gc_counts[k], self.read_lengths[k]) for k in sorted(gc_counts.keys()) if gc_counts[k]]

    def get_ref_percs(self, read_length):
        # reference GC distribution for read lengths 20 through 500 have been precalculated
        ref_bin_fn = self.reference
        if os.path.exists(ref_bin_fn):
            ref_bin_fp = ref_bin_fn
        else:
            ref_bin_fp = '..\%s' % ref_bin_fn
        # first line in table corresponds to read_length of rl_start (currently 20)
        rl_start = 20
        with open(ref_bin_fp, 'r') as ref_bin_f:
            while True:
                line = ref_bin_f.readline()
                if not line: return []
                if read_length == rl_start: return [float(rp) for rp in line.split(',') if float_check(rp)]
                rl_start += 1
        return []

    def normalize_coverage(self, read_perc, ref_perc):
        if read_perc == 0.0 and ref_perc == 0.0:
            return 1
        elif ref_perc == 0.0:
            return 3
        elif read_perc == 'NA' or ref_perc == 'NA':
            return 0.0
        else:
            norm_cov = round(read_perc / ref_perc, 2)
            return min(3, norm_cov)

    def calculate_bias(self):
        precision = 3
        self.gc_table = [['GC Bin', 'GC Percent (%)', 'Read #', 'Count', 'Percent of Total (%)', 'Ref. Percent (%)', 'Normalized Coverage']]

        self.gc_percs_list = []
        self.read_percs_list = []
        self.ref_percs_list = []
        self.norm_covs_list = []
        for rix, gc_bins in enumerate(self.gc_bins_list):
            read_length = self.read_lengths[rix]
            ref_percs = self.get_ref_percs(read_length)
            self.ref_percs_list.append(ref_percs)

            read_count = sum(gc_bins)

            gc_percs = []
            read_percs = []
            norm_covs = []
            for gc_value in range(read_length + 1):
                bin_count = gc_bins[gc_value]
                gc_perc = round(float(gc_value) / float(read_length) * 100.0, precision)
                gc_percs.append(gc_perc)
                
                read_perc = round(float(bin_count) / float(read_count) * 100.0, precision) if read_count else 0.00
                read_percs.append(read_perc)
                
                ref_perc = round(ref_percs[gc_value], precision)
                normalize_gc_coveraged = self.normalize_coverage(read_perc, ref_perc)
                norm_covs.append(normalize_gc_coveraged)
                self.gc_table.append([gc_value, gc_perc, rix + 1, bin_count, read_perc, ref_perc, normalize_gc_coveraged])
            self.gc_percs_list.append(gc_percs)
            self.read_percs_list.append(read_percs)
            self.norm_covs_list.append(norm_covs)
        if len(self.read_lengths) ==2 and self.read_lengths[0] == self.read_lengths[1]:
            read_length = self.read_lengths[0]
            self.read_lengths.append(read_length)
            
            ref_percs = self.get_ref_percs(read_length)
            self.ref_percs_list.append(ref_percs)

            gc_bins = [r1 + r2 for r1, r2 in zip(self.gc_bins_list[0], self.gc_bins_list[1])]
            read_count = sum(gc_bins)

            gc_percs = []
            read_percs = []
            norm_covs = []
            for gc_value in range(read_length + 1):
                bin_count = gc_bins[gc_value]
                gc_perc = round(float(gc_value) / float(read_length) * 100.0, precision)
                gc_percs.append(gc_perc)
                
                read_perc = round(float(bin_count) / float(read_count) * 100.0, precision) if read_count else 0.00
                read_percs.append(read_perc)
                
                ref_perc = round(ref_percs[gc_value], precision)
                normalize_gc_coveraged = self.normalize_coverage(read_perc, ref_perc)
                norm_covs.append(normalize_gc_coveraged)
                self.gc_table.append([gc_value, gc_perc, 'Both', bin_count, read_perc, ref_perc, normalize_gc_coveraged])
            self.gc_percs_list.append(gc_percs)
            self.read_percs_list.append(read_percs)
            self.norm_covs_list.append(norm_covs)
            
        outfile = self.outdir + 'GC_Table.csv'
        with open(outfile, 'w') as of:
            for row in self.gc_table:
                of.write(','.join(map(str, row)) + '\n')
        return

    def plot(self):
        import matplotlib.pyplot as plt

        temp_tag = 'Read%(read_num)s_' if len(self.gc_bins_list) == 2 else ''

        rix = 1
        for read_percs, ref_percs, gc_percs, norm_covs, read_length in zip(self.read_percs_list, self.ref_percs_list, self.gc_percs_list,
                                                                           self.norm_covs_list, self.read_lengths):
            read_tag = temp_tag % {'read_num': rix} if rix < 3 else 'Total_'
            
            fig, ax = plt.subplots()
            plt.plot(list(range(read_length + 1)), read_percs, label='Read', color='r')
            plt.plot(list(range(read_length + 1)), ref_percs, label='Reference', color='b')
            ax.grid(linestyle=':')
            plt.xlabel('GC Count (bp)')
            plt.ylabel('Frequency (%)')
            plt.legend()
            png_path = self.outdir + '%sGC_Distribution.png' % read_tag
            fig.savefig(png_path)
            
            fig, ax = plt.subplots()
            plt.plot(gc_percs, norm_covs, label='Read', color='r')
            ax.grid(linestyle=':')
            plt.xlabel('GC Percent (%)')
            plt.ylabel('Normalized Coverage (x)')
            plt.legend()
            png_path = self.outdir + '%sNormalized_GC_Coverage.png' % read_tag
            fig.savefig(png_path)

            if rix == 3:
                read_tag = 'PE_'
                fig, ax = plt.subplots()
                plt.plot(list(range(read_length + 1)), self.read_percs_list[0], label='Read1', color='r')
                plt.plot(list(range(read_length + 1)), ref_percs, label='Reference', color='b')
                plt.plot(list(range(read_length + 1)), self.read_percs_list[1], label='Read2', color='y')
                ax.grid(linestyle=':')
                plt.xlabel('GC Count (bp)')
                plt.ylabel('Frequency (%)')
                plt.legend()
                png_path =self.outdir + '%sGC_Distribution.png' % read_tag
                fig.savefig(png_path)
                
                fig, ax = plt.subplots()
                plt.plot(gc_percs, self.norm_covs_list[0], label='Read1', color='r')
                plt.plot(gc_percs, self.norm_covs_list[1], label='Read2', color='y')
                ax.grid(linestyle=':')
                plt.xlabel('GC Percent (%)')
                plt.ylabel('Normalized Coverage (x)')
                plt.legend()
                png_path = self.outdir + '%sNormalized_GC_Coverage.png' % read_tag
                fig.savefig(png_path)

            rix += 1
    
        return

    def run(self):
        gc_counts = self.count_GC()
        print(gc_counts)
        self.gc_bins_list = self.count2bins(gc_counts)
        self.calculate_bias()
        self.plot()
        return
            

def main():
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage)
    ArgParser.add_argument("-r", "--ref", action="store", dest="reference", default="Ecoli", help="Use specific reference. [%(default)s]")
    ArgParser.add_argument("-o", "--out", action="store", dest="out",
                           default="", 
                           help="desired directory for output files. [cwd]")
    ArgParser.add_argument("-r1", "--r1_len", action="store", dest="r1_len", default=40, help="read1 length", type=int)
    ArgParser.add_argument("-r2", "--r2_len", action="store", dest="r2_len", default=100, help="read2 length", type=int)
    (para, args) = ArgParser.parse_known_args()


    if not args:
        ArgParser.print_help()
        print("\nERROR: Read file(s) must be specified!", file=sys.stderr)
        sys.exit(1)
    else:
        data_files = args
        reference = para.reference
        outdir = para.out
        r1_len = para.r1_len
        r2_len = para.r2_len

    gcb = GC_Bias(data_files, reference, outdir, r1_len, r2_len)
    gcb.run()


if __name__ == '__main__':
    main()
