#!/usr/bin/env perl
## Merged stLFR barcode splitter
## Replaces:
##   split_barcode_stLFR_42_reads_noBC.pl    (--swap none  --output_mode separate)
##   split_barcode_stLFR_42_reads_ac_swap.pl (--swap ac    --output_mode single)
##   split_barcode_stLFR_42_reads_gc_swap.pl (--swap gc    --output_mode single)
##   split_barcode_PEXXX_42_reads_BCgDNA_fullreads.pl (--swap none --output_mode stratified)
##
## Usage:
##   perl split_barcode_stLFR.pl \
##     --barcode barcode.list \
##     --r1 read_1.fq.gz --r2 read_2.fq.gz \
##     --read_len 100 --output data/split_read \
##     --bc_start 1 --gdna_start 43 \
##     --additional_bc_start 0 --additional_bc_len 0 \
##     --gdna_start_r1 1 --read_len_r1 150 \
##     --adapter_len 0 \
##     --swap none \
##     --output_mode single

use strict;
use Getopt::Long;

my ($barcode_file, $r1_file, $r2_file);
my ($read_len, $output_prefix);
my ($bc_start, $gdna_start, $additional_bc_start, $additional_bc_len);
my ($gdna_start_r1, $read_len_r1);
my $adapter_len  = 0;
my $swap         = 'none';   # none | ac | gc
my $output_mode  = 'single'; # single | separate | stratified

GetOptions(
    'barcode=s'             => \$barcode_file,
    'r1=s'                  => \$r1_file,
    'r2=s'                  => \$r2_file,
    'read_len=i'            => \$read_len,
    'output=s'              => \$output_prefix,
    'bc_start=i'            => \$bc_start,
    'gdna_start=i'          => \$gdna_start,
    'additional_bc_start=i' => \$additional_bc_start,
    'additional_bc_len=i'   => \$additional_bc_len,
    'gdna_start_r1=i'       => \$gdna_start_r1,
    'read_len_r1=i'         => \$read_len_r1,
    'adapter_len=i'         => \$adapter_len,
    'swap=s'                => \$swap,
    'output_mode=s'         => \$output_mode,
) or die "Error parsing arguments\n";

die "Missing --barcode\n"  unless $barcode_file;
die "Missing --r1\n"       unless $r1_file;
die "Missing --r2\n"       unless $r2_file;
die "Invalid --swap: $swap (use none|ac|gc)\n"
    unless $swap =~ /^(none|ac|gc)$/;
die "Invalid --output_mode: $output_mode (use single|separate|stratified)\n"
    unless $output_mode =~ /^(single|separate|stratified)$/;

# Convert 1-based positions to 0-based indices
my $bc_start_idx            = $bc_start - 1;
my $gdna_start_idx          = $gdna_start - 1;
my $additional_bc_start_idx = $additional_bc_start - 1;
my $gdna_start_idx_r1       = $gdna_start_r1 - 1;

my ($n1, $n2, $n3, $n4, $n5) = (10, 6, 10, 6, 10);
my %barcode_hash;
my %barcode_hash_string;
my ($bc1_cnt, $bc2_cnt, $bc3_cnt) = (0, 0, 0);
my $hash_string;

# Load string barcode hash (for BC:Z: tag) — only needed for single/separate modes
if ($output_mode ne 'stratified') {
    open IN, $barcode_file or die "can't open $barcode_file";
    while (<IN>) {
        my @line    = split;
        my @barcode = split(//, $line[0]);
        my $id2     = $line[0];
        for my $num (0..9) {
            for my $base ('A','G','C','T') {
                my @mis = @barcode;
                $mis[$num] = $base;
                $barcode_hash_string{join('', @mis)} = $id2;
            }
        }
    }
    close IN;
}

# Load integer barcode hash (for BX:Z: tag)
open IN, $barcode_file or die "can't open $barcode_file";
my $n = 0;
while (<IN>) {
    chomp;
    $n++;
    my @line    = split;
    my @barcode = split(//, $line[0]);
    my $id      = $line[1];
    for my $num (0..9) {
        for my $base ('A','G','C','T') {
            my @mis = @barcode;
            $mis[$num] = $base;
            $barcode_hash{join('', @mis)} = $id;
        }
    }
}
close IN;

my $barcode_types = $n * $n * $n;
my $barcode_each  = $n;

# Open input
open IN1, "gzip -dc $r1_file |" or die "cannot open read1";
open IN2, "gzip -dc $r2_file |" or die "cannot open read2";
open OUT1, "| gzip > $output_prefix.1.fq.gz" or die "Can't write OUT1";
open OUT2, "| gzip > $output_prefix.2.fq.gz" or die "Can't write OUT2";

if ($output_mode eq 'separate') {
    open NOBC_OUT1, "| gzip > data/noBC.1.fq.gz" or die "Can't write NOBC_OUT1";
    open NOBC_OUT2, "| gzip > data/noBC.2.fq.gz" or die "Can't write NOBC_OUT2";
} elsif ($output_mode eq 'stratified') {
    open OUT1_nobc2bc3, "| gzip > ${output_prefix}_nobc2bc3.1.fq.gz" or die "Can't write";
    open OUT2_nobc2bc3, "| gzip > ${output_prefix}_nobc2bc3.2.fq.gz" or die "Can't write";
    open OUT1_noall3,   "| gzip > ${output_prefix}_noall3.1.fq.gz"   or die "Can't write";
    open OUT2_noall3,   "| gzip > ${output_prefix}_noall3.2.fq.gz"   or die "Can't write";
    open OUT1_nobc3,    "| gzip > ${output_prefix}_nobc3.1.fq.gz"    or die "Can't write";
    open OUT2_nobc3,    "| gzip > ${output_prefix}_nobc3.2.fq.gz"    or die "Can't write";
    open OUT1_nobc1,    "| gzip > ${output_prefix}_nobc1.1.fq.gz"    or die "Can't write";
    open OUT2_nobc1,    "| gzip > ${output_prefix}_nobc1.2.fq.gz"    or die "Can't write";
}

$n = 0;
my $reads_num = 0;
my $progress  = 0;
my %index_hash;
my %index_hash_reverse;
my $split_barcode_num = 0;
my $T;
my $id_str;
my @line;
my @Read_num;
$Read_num[0] = 0;
my $split_reads_num = 0;
my ($additional_bc, $seq_r1, $qual_r1);

while (<IN2>) {
    chomp;
    $_ =~ tr/ACac/CAca/ if ($swap eq 'ac' || $swap eq 'gc');
    @line = split;
    $n++;

    if ($n % 4 == 1) {
        $reads_num++;
        my @A = split(/\//, $line[0]);
        $id_str = $A[0];
        if ($reads_num % 1000000 == 1) {
            print "reads processed $progress (M) reads ...\n";
            $progress++;
        }
    }

    if ($n % 4 == 2) {
        my $read   = substr($line[0], $gdna_start_idx, $read_len);
        my $b1     = substr($line[0], $bc_start_idx, $n1);
        my $b2     = substr($line[0], $bc_start_idx+$n1+$n2, $n3);
        my $b3     = substr($line[0], $bc_start_idx+$n1+$n2+$n3+$n4, $n5);
        my $bc_str = substr($line[0], $bc_start_idx, $n1+$n2+$n3+$n4+$n5);

        $additional_bc = ($additional_bc_len > 0)
            ? substr($line[0], $additional_bc_start_idx, $additional_bc_len) : '';

        my $bc1_ok = exists $barcode_hash{$b1};
        my $bc2_ok = exists $barcode_hash{$b2};
        my $bc3_ok = exists $barcode_hash{$b3};
        $bc1_cnt++ if $bc1_ok;
        $bc2_cnt++ if $bc2_ok;
        $bc3_cnt++ if $bc3_ok;

        if ($bc1_ok && $bc2_ok && $bc3_ok) {
            # --- All 3 barcodes valid ---
            my $hash = $barcode_hash{$b1}."_".$barcode_hash{$b2}."_".$barcode_hash{$b3};
            unless (exists $index_hash{$hash}) {
                $split_barcode_num++;
                $index_hash{$hash} = $split_barcode_num;
                $index_hash_reverse{$split_barcode_num} = $hash;
                $Read_num[$split_barcode_num] = 0;
            }
            $split_reads_num++;
            $Read_num[$index_hash{$hash}]++;

            $hash_string = $barcode_hash_string{$b1}.$barcode_hash_string{$b2}.$barcode_hash_string{$b3}
                if $output_mode ne 'stratified';

            if ($read_len_r1 > 0) {
                $T = <IN1>; chomp($T);
                # R1 header
                if ($output_mode eq 'stratified') {
                    print OUT1 "$id_str\#$hash\#\/1\tBX:Z:$hash\n";
                } elsif ($additional_bc_len > 0) {
                    print OUT1 "$id_str\#$hash\/1\tBX:Z:$hash\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
                } else {
                    print OUT1 "$id_str\#$hash\/1\tBX:Z:$hash\tBC:Z:$hash_string\n";
                }
                # R1 seq with optional swap
                $T = <IN1>; chomp($T);
                $T =~ tr/ACac/CAca/ if $swap eq 'ac';
                $T =~ tr/GCgc/CGcg/ if $swap eq 'gc';
                $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                print OUT1 "$seq_r1\n";
                $T = <IN1>; chomp($T); print OUT1 "$T\n";
                $T = <IN1>; chomp($T);
                $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                print OUT1 "$qual_r1\n";
            }
            # R2 header
            if ($output_mode eq 'stratified') {
                if ($additional_bc_len > 0) {
                    print OUT2 "$id_str\#${hash}_${additional_bc}\#\/2\tBX:Z:$hash\n";
                } else {
                    print OUT2 "$id_str\#$hash\#\/2\tBX:Z:$hash\n";
                }
            } elsif ($additional_bc_len > 0) {
                print OUT2 "$id_str\#$hash\/2\tBX:Z:$hash\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
            } else {
                print OUT2 "$id_str\#$hash\/2\tBX:Z:$hash\tBC:Z:$hash_string\n";
            }
            print OUT2 "$read\n";
            $T = <IN2>; $n++; chomp($T); print OUT2 "$T\n";
            $T = <IN2>; $n++; chomp($T);
            my $qual = substr($T, $gdna_start_idx, $read_len);
            print OUT2 "$qual\n";

        } else {
            # --- One or more barcodes invalid ---
            $Read_num[0]++;

            if ($output_mode eq 'single') {
                $hash_string = $b1.$b2.$b3;
                if ($read_len_r1 > 0) {
                    $T = <IN1>; chomp($T);
                    if ($additional_bc_len > 0) {
                        print OUT1 "$id_str\#0_0_0\/1\tBX:Z:0_0_0\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
                    } else {
                        print OUT1 "$id_str\#0_0_0\/1\tBX:Z:0_0_0\tBC:Z:$hash_string\n";
                    }
                    $T = <IN1>; chomp($T);
                    $T =~ tr/ACac/CAca/ if $swap eq 'ac';
                    $T =~ tr/GCgc/CGcg/ if $swap eq 'gc';
                    $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print OUT1 "$seq_r1\n";
                    $T = <IN1>; chomp($T); print OUT1 "$T\n";
                    $T = <IN1>; chomp($T);
                    $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print OUT1 "$qual_r1\n";
                }
                if ($additional_bc_len > 0) {
                    print OUT2 "$id_str\#0_0_0\/2\tBX:Z:0_0_0\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
                } else {
                    print OUT2 "$id_str\#0_0_0\/2\tBX:Z:0_0_0\tBC:Z:$hash_string\n";
                }
                print OUT2 "$read\n";
                $T = <IN2>; $n++; chomp($T); print OUT2 "$T\n";
                $T = <IN2>; $n++; chomp($T);
                my $qual = substr($T, $gdna_start_idx, $read_len);
                print OUT2 "$qual\n";

            } elsif ($output_mode eq 'separate') {
                $hash_string = $b1.$b2.$b3;
                if ($read_len_r1 > 0) {
                    $T = <IN1>; chomp($T);
                    if ($additional_bc_len > 0) {
                        print NOBC_OUT1 "$id_str\#0_0_0\/1\tBX:Z:0_0_0\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
                    } else {
                        print NOBC_OUT1 "$id_str\#0_0_0\/1\tBX:Z:0_0_0\tBC:Z:$hash_string\n";
                    }
                    $T = <IN1>; chomp($T);
                    $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print NOBC_OUT1 "$seq_r1\n";
                    $T = <IN1>; chomp($T); print NOBC_OUT1 "$T\n";
                    $T = <IN1>; chomp($T);
                    $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print NOBC_OUT1 "$qual_r1\n";
                }
                if ($additional_bc_len > 0) {
                    print NOBC_OUT2 "$id_str\#0_0_0\/2\tBX:Z:0_0_0\tBC:Z:$hash_string\tMI:Z:$additional_bc\n";
                } else {
                    print NOBC_OUT2 "$id_str\#0_0_0\/2\tBX:Z:0_0_0\tBC:Z:$hash_string\n";
                }
                print NOBC_OUT2 "$read\n";
                $T = <IN2>; $n++; chomp($T); print NOBC_OUT2 "$T\n";
                $T = <IN2>; $n++; chomp($T);
                my $qual = substr($T, $gdna_start_idx, $read_len);
                print NOBC_OUT2 "$qual\n";

            } elsif ($output_mode eq 'stratified') {
                my $hash     = $barcode_hash{$b1}."_".$barcode_hash{$b2}."_".$barcode_hash{$b3};
                my $full_len = $n1+$n2+$n3+$n4+$n5+$adapter_len+$read_len;
                my $full_read = substr($line[0], $bc_start_idx, $full_len);

                my ($fh1, $fh2);
                if     ($bc1_ok && !$bc2_ok && !$bc3_ok) { ($fh1,$fh2) = (\*OUT1_nobc2bc3, \*OUT2_nobc2bc3); }
                elsif  (!$bc1_ok && !$bc2_ok && !$bc3_ok){ ($fh1,$fh2) = (\*OUT1_noall3,   \*OUT2_noall3);   }
                elsif  (!$bc1_ok && $bc2_ok && $bc3_ok)  { ($fh1,$fh2) = (\*OUT1_nobc1,    \*OUT2_nobc1);    }
                elsif  ($bc1_ok && $bc2_ok && !$bc3_ok)  { ($fh1,$fh2) = (\*OUT1_nobc3,    \*OUT2_nobc3);    }
                else {
                    # other partial combinations: consume and discard
                    if ($read_len_r1 > 0) { for (1..4) { $T = <IN1>; } }
                    $T = <IN2>; $n++; $T = <IN2>; $n++;
                    next;
                }

                if ($read_len_r1 > 0) {
                    $T = <IN1>; chomp($T);
                    print $fh1 "$id_str\#$bc_str\#\/1\tBX:Z:$hash\n";
                    $T = <IN1>; chomp($T);
                    $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print $fh1 "$seq_r1\n";
                    $T = <IN1>; chomp($T); print $fh1 "$T\n";
                    $T = <IN1>; chomp($T);
                    $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
                    print $fh1 "$qual_r1\n";
                }
                if ($additional_bc_len > 0) {
                    print $fh2 "$id_str\#${bc_str}_${additional_bc}\#\/2\tBX:Z:$hash\n";
                } else {
                    print $fh2 "$id_str\#$bc_str\#\/2\tBX:Z:$hash\n";
                }
                print $fh2 "$full_read\n";
                $T = <IN2>; $n++; chomp($T); print $fh2 "$T\n";
                $T = <IN2>; $n++; chomp($T);
                my $qual = substr($T, $bc_start_idx, $full_len);
                print $fh2 "$qual\n";
            }
        }
    }
}

close IN1; close IN2;
close OUT1; close OUT2;
if ($output_mode eq 'separate') {
    close NOBC_OUT1; close NOBC_OUT2;
} elsif ($output_mode eq 'stratified') {
    close OUT1_nobc2bc3; close OUT2_nobc2bc3;
    close OUT1_noall3;   close OUT2_noall3;
    close OUT1_nobc3;    close OUT2_nobc3;
    close OUT1_nobc1;    close OUT2_nobc1;
}

open OUT3, ">split_stat_read1.log" or die "Can't write file";
print OUT3 "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
my $r;
$r = 100 * $split_barcode_num / $barcode_types;
print OUT3 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 * $split_reads_num / $reads_num;
print OUT3 "Reads_pair_num  = $reads_num \n";
print OUT3 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for (my $i = 1; $i <= $split_barcode_num; $i++) {
    print OUT3 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}
close OUT3;

open OUT4, ">bc.individual.count.log" or die "Can't write file";
if ($read_len_r1 == 0) {
    print OUT4 "$reads_num\n";
} else {
    print OUT4 $reads_num * 2, "\n";
}
print OUT4 "bc1_cnt=$bc1_cnt \n";
print OUT4 "bc2_cnt=$bc2_cnt \n";
print OUT4 "bc3_cnt=$bc3_cnt \n";
print OUT4 "$bc1_cnt\t$bc2_cnt\t$bc3_cnt\n";
print OUT4 $bc1_cnt/$reads_num, "\t", $bc2_cnt/$reads_num, "\t", $bc3_cnt/$reads_num, "\n";
close OUT4;

print "all done!\n";
