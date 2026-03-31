#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);

die "Usage:perl $0 <in.picard.bam> </home/tangjb/BIN/samtools-0.1.18/samtools>" unless @ARGV==2;
my $bam = shift;
my $Samtools = shift;
my $clean_reads=0;
my $clean_bases=0;
# my $mapped_reads_with2nd=0;
my $mapped_reads=0;
my $mapped_bases=0;
my $mapping_rate=0;
my $dup_reads=0;
my $dup_rate=0;
my $mismatch_bases=0;
my $mismatch_rate=0;
my $uniq_reads=0;
my $uniq_bases=0;
my $uniq_rate=0;
# my $unmapped_reads=0;
# my $unmapped_rate=0;
# open BAM,"$Samtools view -f 0x1 $bam | " or die;

open BAM,"$Samtools view -X $bam | " or die;

while(<BAM>)
{
		chomp;
		my $info=$_;
        if($info=~/^@/){next;}
		my @info=split /\t/,$info;
		$clean_reads++;
		$clean_bases+=length($info[9]);
		unless($info[1]=~/u/) # unmapped
		{
			if($info[4] >= 1) # mapq 0
			{
				# uniq mapped reads
				if($info=~/AS:i:(\d+)/){$uniq_reads++;$uniq_bases+=$1;}
				else{$uniq_reads++;$uniq_bases+=length($info[9]);}
			} # mapped_reads include mapq 0, double count
		  	if($info=~/AS:i:(\d+)/){$mapped_reads++;$mapped_bases+=$1;}
			else{$mapped_reads++;$mapped_bases+=length($info[9]);}
			if($info[1]=~/d/){$dup_reads++;}
			if($info=~/NM:i:(\d+)/){$mismatch_bases+=$1;}
		}
		# $unmapped_reads++;
}
# $mapped_reads=$clean_reads-$unmapped_reads;
$mismatch_rate=$mismatch_bases/$mapped_bases;
$mapping_rate=$mapped_reads/$clean_reads;
## dup rate against clean_reads
$dup_rate=$dup_reads/$mapped_reads;
$uniq_rate=$uniq_reads/$mapped_reads;
# $unmapped_rate=$unmapped_reads/$clean_reads;

my $name=basename($bam);
my $sample=(split /\./,$name)[0];
print "Sample\t$sample\n";
print "Clean reads\t$clean_reads\n";
print "Clean bases(bp)\t$clean_bases\n";
# print "Mapped reads (count multi-alingment)\t$mapped_reads_with2nd\n";
print "Mapped reads\t$mapped_reads\n";
print "Mapped bases(bp)\t$mapped_bases\n";
print "Mapping rate\t",sprintf("%.2f%%",100*$mapping_rate),"\n";
print "Uniq reads\t$uniq_reads\n";
print "Uniq bases(bp)\t$uniq_bases\n";
print "Unique rate\t",sprintf("%.2f%%",100*$uniq_rate),"\n";
print "Duplicate reads\t$dup_reads\n";
print "Duplicate rate\t",sprintf("%.2f%%",100*$dup_rate),"\n";
print "Mismatch bases(bp)\t$mismatch_bases\n";
print "Mismatch rate\t",sprintf("%.2f%%",100*$mismatch_rate),"\n";
# print "unmapped reads\t$unmapped_reads\n";
# print "unmapped rate\t",sprintf("%.2f%%",100*$unmapped_rate),"\n";
close $bam;
