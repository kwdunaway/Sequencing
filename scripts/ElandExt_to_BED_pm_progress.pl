#!/usr/bin/perl -w
use strict; use warnings;

# ElandExt_to_BEDv2.pl          6/6/12 version
#
# This program takes sequences in Eland Extended format (also known as s_#_export.txt files) and 
# produces multiple files BED which can be used analyze the data.  The bed files go out to the 
# 6th column (strand notation) but have null values for columns 4 (name) and 5 (score).  This
# program will also give statistics and separate out nonuniquely mapped reads.
#
# This is a little more flexible than other ElandExt_to_BED versions because all you need for 
# chromosomal notation is chr##  (ex: chr1 or chr13) somewhere in the chrom notation.  Also,
# you can input read sequence length (which is better than the program figuring it out for speed
# purposes.  Finally, you must input which array field the chromosome, position, and strand are
# in just in case they are different from export file to export file (10, 11, 12 for Berkeley files).
#
# Arguments:
#   1) Input file name
#   2) Output file name
#   3) Read length
#   4) Chromosome array number
#   5) Position array number
#   6) Strand array number
#   7) First Char



##################################################
# Command Line Error Checking and I/O Initiation #
##################################################

die "useage: ElandExt_to_BEDv2.pl 
<input file name> 
<output filename>
<read length>
<Chromosome array number>
<Position array number>
<Strand array number>
<First Char>" unless @ARGV == 7;
my $infile = shift(@ARGV);
my $outfile = shift(@ARGV); 
my $readlength = shift(@ARGV); 
my $chr = shift(@ARGV); 
my $pos = shift(@ARGV); 
my $strand = shift(@ARGV); 
my $firstchar = shift(@ARGV); 
if (! -d $outfile)
{ 
     `mkdir $outfile`; #creates dir if one doesn't exist
     if (! -d $outfile) { die "directory ($outfile) does not exist"} #dir still does not exist
}



##################################################
# Sets up global variables, infile, and outfiles #
##################################################

open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read

my @array;
my $QCcount = 0;
my $NMcount = 0;
my $NonUniqcount = 0;
my $Unknowncount = 0;
my @ChrMapcount;
for (my $n = 0; $n < 27; $n++)
{
    $ChrMapcount[$n] = 0;
}
my $totalcount = 0;
my $chrnum = 1;
my $filename; # Outfile for each chromosome in bed
my @chr_out;

#Create Bed File for each Chromosome
for(my $n = 0; $n < 30; $n++)
{
	if($n > 23)
	{
		if($n == 24)
		{
			$filename = $outfile . "/" . $outfile . "_chrX" . ".bed";
			open($chr_out[$n], ">$filename") or die "cannot open $filename outfile";
		}
	}
	else
	{
		$filename = $outfile . "/" . $outfile . "_chr" . ($n + 1) . ".bed";
		open($chr_out[$n], ">$filename") or die "cannot open $filename outfile";
	}
}

for(my $n = 0; $n < 30; $n++)
{
	close($chr_out[$n]);
}
close(IN);


