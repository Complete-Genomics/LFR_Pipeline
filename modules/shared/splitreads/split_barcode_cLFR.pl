#!/usr/bin/env perl
use strict;
use Getopt::Long;
# allow pBC in r1
# Merged from split_barcode_cLFR_randomBC_tag_groupBC.pl and ...rc.pl
# Use --reverse_complement flag to apply RC transform on extracted barcodes

my ($r1_file, $r2_file, $output_prefix);
my $read_len      = 0;
my $cbc_len       = 0;
my $bc_len_redundant = 0;
my $bc_start      = 0;
my $gdna_start    = 0;
my $additional_bc_start  = 0;
my $additional_bc_len    = 0;
my $gdna_start_r1 = 0;
my $read_len_r1   = 0;
my $additional_bc_len_r1 = 0;
my $reverse_complement   = 0;

GetOptions(
    'r1=s'                   => \$r1_file,
    'r2=s'                   => \$r2_file,
    'read_len=i'             => \$read_len,
    'output=s'               => \$output_prefix,
    'cbc_len=i'              => \$cbc_len,
    'bc_len_redundant=i'     => \$bc_len_redundant,
    'bc_start=i'             => \$bc_start,
    'gdna_start=i'           => \$gdna_start,
    'additional_bc_start=i'  => \$additional_bc_start,
    'additional_bc_len=i'    => \$additional_bc_len,
    'gdna_start_r1=i'        => \$gdna_start_r1,
    'read_len_r1=i'          => \$read_len_r1,
    'additional_bc_len_r1=i' => \$additional_bc_len_r1,
    'reverse_complement'     => \$reverse_complement,
) or die "Error in command line arguments\n";

die "Usage: perl split_barcode_cLFR.pl --r1 <r1.fq.gz> --r2 <r2.fq.gz> --read_len <N> --output <prefix> "
  . "--cbc_len <N> --bc_len_redundant <N> --bc_start <N> --gdna_start <N> "
  . "--additional_bc_start <N> --additional_bc_len <N> --gdna_start_r1 <N> "
  . "--read_len_r1 <N> --additional_bc_len_r1 <N> [--reverse_complement]\n"
  unless ($r1_file && $r2_file && $output_prefix);

# Convert from 1-based to 0-based indices
my $bc_start_idx             = $bc_start - 1;
my $gdna_start_idx           = $gdna_start - 1;
my $additional_bc_start_idx  = $additional_bc_start - 1;
my $gdna_start_idx_r1        = $gdna_start_r1 - 1;

my %bc_dict;

open IN1, "gzip -dc $r1_file |" or die "cannot open read1";
open IN2, "gzip -dc $r2_file |" or die "cannot open read2";
open OUT1, "| gzip > $output_prefix.1.fq.gz" or die "Can't write file";
open OUT2, "| gzip > $output_prefix.2.fq.gz" or die "Can't write file";

my $n = 0;
my $reads_num;
my $progress;
my %index_hash;
my %index_hash_reverse;
my $split_barcode_num;
my $T;
my $id;
my @line;
my @Read_num;
$Read_num[0] = 0;
my $split_reads_num;
my $bc;
my $bx;
my $seq;
my $additional_bc;
my $seq_r1;
my $qual_r1;

while(<IN2>){
  chomp;
  @line = split;
  $n++;
  if($n % 4 == 1){
    $reads_num++;
    my @A = split(/\//, $line[0]);
    $id = $A[0];
    if($reads_num % 1000000 == 1){
      print "reads processed $progress (M) reads ...\n";
      $progress++;
    }
  }

  if($n % 4 == 2){
    if($reverse_complement){
      my $bc_tmp = substr($line[0], $bc_start_idx, $bc_len_redundant);
      $bc = reverse $bc_tmp;
      $bc =~ tr/ATGCatgc/TACGtacg/;
      my $bx_tmp = substr($line[0], $bc_start_idx, $cbc_len);
      $bx = reverse $bx_tmp;
      $bx =~ tr/ATGCatgc/TACGtacg/;
    } else {
      $bc = substr($line[0], $bc_start_idx, $bc_len_redundant);
      $bx = substr($line[0], $bc_start_idx, $cbc_len);
    }
    $seq = substr($line[0], $gdna_start_idx, $read_len);

    push @{$bc_dict{$bx}}, $n;

    if(!(exists $index_hash{$bx})){
      $split_barcode_num++;
      $index_hash{$bx} = $split_barcode_num;
      $index_hash_reverse{$split_barcode_num} = $bx;
      $Read_num[$index_hash{$bx}] = 0;
    }
    $split_reads_num++;
    $Read_num[$index_hash{$bx}]++;

    ## for PE
    if($read_len_r1 > 0){
      $T = <IN1>; chomp($T);

      # pBC on r2
      if($additional_bc_len > 0){
        $additional_bc = substr($line[0], $additional_bc_start_idx, $additional_bc_len);
        print OUT1 "$id\#$bc\#$additional_bc\/1\tBX:Z:$bx\n";
        print OUT2 "$id\#$bc\#$additional_bc\/2\tBX:Z:$bx\n";
      }

      $T = <IN1>; chomp($T);
      $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);

      # pBC on r1
      if($additional_bc_len_r1 > 0){
        $additional_bc = substr($T, $additional_bc_start_idx, $additional_bc_len_r1);
        print OUT1 "$id\#$bc\#$additional_bc\/1\tBX:Z:$bx\n";
        print OUT2 "$id\#$bc\#$additional_bc\/2\tBX:Z:$bx\n";
      }
      # no pBC
      if($additional_bc_len == 0 && $additional_bc_len_r1 == 0){
        print OUT1 "$id\#$bc\/1\tBX:Z:$bx\n";
        print OUT2 "$id\#$bc\/2\tBX:Z:$bx\n";
      }

      print OUT1 "$seq_r1\n";
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n";
      $T = <IN1>; chomp($T);
      $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
      print OUT1 "$qual_r1\n";

      print OUT2 "$seq\n";
      $T = <IN2>; $n++; chomp($T);
      print OUT2 "$T\n";
      $T = <IN2>; $n++; chomp($T);
      my $qual = substr($T, $gdna_start_idx, $read_len);
      print OUT2 "$qual\n";

    } else {
      # SE, pBC on r2
      if($additional_bc_len > 0){
        $additional_bc = substr($line[0], $additional_bc_start_idx, $additional_bc_len);
        print OUT2 "$id\#$bc\#$additional_bc\/2\tBX:Z:$bx\n";
      } else {
        print OUT2 "$id\#$bc\/2\tBX:Z:$bx\n";
      }
      print OUT2 "$seq\n";
      $T = <IN2>; $n++; chomp($T);
      print OUT2 "$T\n";
      $T = <IN2>; $n++; chomp($T);
      my $qual = substr($T, $gdna_start_idx, $read_len);
      print OUT2 "$qual\n";
    }
  }
}

close IN1;
close IN2;
close OUT1;
close OUT2;

my $datestring = localtime();
print "writing to split.fq.gz done at $datestring\n";

my $barcode_types = 4**$cbc_len;

open OUT3, ">split_stat_read1.log" or die "Can't write file";
print OUT3 "Barcode_types = 4** $cbc_len = $barcode_types\n";
my $r;
$r = 100 * $split_barcode_num / $barcode_types;
print OUT3 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 * $split_barcode_num / $reads_num;
print OUT3 "Reads_pair_num  = $reads_num \n";
print OUT3 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for(my $i = 1; $i <= $split_barcode_num; $i++){
  print OUT3 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}
close OUT3;
$datestring = localtime();
print "writing to split_stat_read1.log done at $datestring\n";

### check diversity
open OUT4, ">bc.diversity.count.log" or die "Can't write OUT4";
my $bc_diversity_split = $split_barcode_num / $reads_num;
if($read_len_r1 == 0){
  print OUT4 "$reads_num\n";
} else {
  my $reads_num_pe = $reads_num * 2;
  print OUT4 "$reads_num_pe\n";
}
print OUT4 "$bc_diversity_split\n";
close OUT4;

print "all done!\n";
