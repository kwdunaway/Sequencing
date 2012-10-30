#!/usr/bin/perl -w
use strict; use warnings;

# ElandExt_to_BEDv2.pl          6/6/12 version
#
# This program takes sequences in Eland Extended format (also known as s_#_export.txt files) and 
# produces multiple files BED which can be used analyze the data.  The bed files go out to the 
# 6th column (strand notation) but have null values for columns 4 (name) and 5 (score).  This
# program will also give statistics and separate out nonuniquely mapped reads.
#
# This is a little more flexable than other ElandExt_to_BED versions because all you need for 
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
     if (! -d $outfile) { die "directory ($outfile) does not exist"} #creates dir if one doesn't exist
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


my $filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT1, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT2, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT3, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT4, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT5, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT6, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT7, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT8, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT9, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT10, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT11, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT12, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT13, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT14, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT15, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT16, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT17, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT18, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT19, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT20, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT21, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT22, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
$chrnum=$chrnum+1;

$filename = $outfile . "/" . $outfile . "_chr" . $chrnum . ".bed";
open(OUT23, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
 
$filename = $outfile . "/" . $outfile . "_chrX" . ".bed";
open(OUTX, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
  
$filename = $outfile . "/" . $outfile . "_chrY" . ".bed";
open(OUTY, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

$filename = $outfile . "/" . $outfile . "_chrM" . ".bed";
open(OUTM, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
  
$filename = $outfile . "/" . $outfile . "_NM.fq";
open(OUTNM, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

$filename = $outfile . "/" . $outfile . "_NonUnique.fq";
open(OUTNONUNIQ, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

$filename = $outfile . "/" . "Stats_" . $outfile . ".txt";
open(OUTSTAT, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

$filename = $outfile . "/" . "Unknown_" . $outfile . ".txt";
open(UNKNOWN, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written


 

####################################################################################
# Goes through each line and if it is the chromosome you want, prints info to file #
####################################################################################

while (<IN>)
{
    chomp;
    @array = split ("\t", $_);
#    @assignment = split ("", $array[$chr]);
    $totalcount = $totalcount +1;
    if ($array[$chr] eq "QC")
    { 
	print OUTNM "$_" , "\n"; 
	$QCcount = $QCcount + 1;
    }
    elsif ($array[$chr] eq "NM")
    { 	
	print OUTNM "$_" , "\n"; 
	$NMcount = $NMcount + 1;
    }
#    elsif ($assignment[0] ne "$firstchar")
    elsif ($array[$chr] =~ m/:/)
    { 	
	print OUTNONUNIQ "$_" , "\n"; 
	$NonUniqcount = $NonUniqcount +1;
    }
    elsif ($array[$chr] =~ m/chr10/)
    { 	
	print OUT10 "chr10 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[10] = $ChrMapcount[10] + 1;
    }
    elsif ($array[$chr] =~ m/chr11/)
    { 	
	print OUT11 "chr11 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[11] = $ChrMapcount[11] + 1;
    }
    elsif ($array[$chr] =~ m/chr12/)
    { 	
	print OUT12 "chr12 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[12] = $ChrMapcount[12] + 1;
    }
    elsif ($array[$chr] =~ m/chr13/)
    { 	
	print OUT13 "chr13 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[13] = $ChrMapcount[13] + 1;
    }
    elsif ($array[$chr] =~ m/chr14/)
    { 	
	print OUT14 "chr14 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[14] = $ChrMapcount[14] + 1;
    }
    elsif ($array[$chr] =~ m/chr15/)
    { 	
	print OUT15 "chr15 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[15] = $ChrMapcount[15] + 1;
    }
    elsif ($array[$chr] =~ m/chr16/)
    { 	
	print OUT16 "chr16 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[16] = $ChrMapcount[16] + 1;
    }
    elsif ($array[$chr] =~ m/chr17/)
    { 	
	print OUT17 "chr17 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[17] = $ChrMapcount[17] + 1;
    }
    elsif ($array[$chr] =~ m/chr18/)
    { 	
	print OUT18 "chr18 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[18] = $ChrMapcount[18] + 1;
    }
    elsif ($array[$chr] =~ m/chr19/)
    { 	
	print OUT19 "chr19 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[19] = $ChrMapcount[19] + 1;
    }
    elsif ($array[$chr] =~ m/chr20/)
    { 	
	print OUT20 "chr20 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[20] = $ChrMapcount[20] + 1;
    }
    elsif ($array[$chr] =~ m/chr21/)
    { 	
	print OUT21 "chr21 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[21] = $ChrMapcount[21] + 1;
    }
    elsif ($array[$chr] =~ m/chr22/)
    { 	
	print OUT22 "chr22 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[22] = $ChrMapcount[22] + 1;
    }
    elsif ($array[$chr] =~ m/chr23/)
    { 	
	print OUT23 "chr23 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[23] = $ChrMapcount[23] + 1;
    }
    elsif ($array[$chr] =~ m/chr1/)
    { 	
	print OUT1 "chr1 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[1] = $ChrMapcount[1] + 1;
    }
    elsif ($array[$chr] =~ m/chr2/)
    { 	
	print OUT2 "chr2 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[2] = $ChrMapcount[2] + 1;
    }
    elsif ($array[$chr] =~ m/chr3/)
    { 	
	print OUT3 "chr3 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[3] = $ChrMapcount[3] + 1;
    }
    elsif ($array[$chr] =~ m/chr4/)
    { 	
	print OUT4 "chr4 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[4] = $ChrMapcount[4] + 1;
    }
    elsif ($array[$chr] =~ m/chr5/)
    { 	
	print OUT5 "chr5 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[5] = $ChrMapcount[5] + 1;
    }
    elsif ($array[$chr] =~ m/chr6/)
    { 	
	print OUT6 "chr6 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[6] = $ChrMapcount[6] + 1;
    }
    elsif ($array[$chr] =~ m/chr7/)
    { 	
	print OUT7 "chr7 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[7] = $ChrMapcount[7] + 1;
    }
    elsif ($array[$chr] =~ m/chr8/)
    { 	
	print OUT8 "chr8 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[8] = $ChrMapcount[8] + 1;
    }
    elsif ($array[$chr] =~ m/chr9/)
    { 	
	print OUT9 "chr9 \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[9] = $ChrMapcount[9] + 1;
    }
    elsif ($array[$chr] =~ m/chrX/)
    { 	
	print OUTX "chrX \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[24] = $ChrMapcount[24] + 1;
    }
    elsif ($array[$chr] =~ m/chrY/)
    { 	
	print OUTY "chrY \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[25] = $ChrMapcount[25] + 1;
    }
    elsif ($array[$chr] =~ m/chrM/)
    { 	
	print OUTM "chrM \t" , $array[$pos] , "\t" , $array[$pos]+$readlength, "\t", $outfile , "\t", "0", "\t" , $array[$strand] , "\n";
	$ChrMapcount[0] = $ChrMapcount[0] + 1;
	$ChrMapcount[26] = $ChrMapcount[26] + 1;
    }
    else
    {
	print UNKNOWN "$_" , "\n";
	$Unknowncount = $Unknowncount + 1;
    }
}
print "Finished separating data by chromosome \n";



####################################
# Prints Stats to the Stat outfile #
####################################

my $weirdcount = $totalcount - $QCcount - $NMcount - $NonUniqcount - $ChrMapcount[0];
for (my $n = 1; $n < 24; $n++)
{
    if ($ChrMapcount[$n] > 0)
    {
	print OUTSTAT "Number of reads mapped to Chromosome $n is:\t" , $ChrMapcount[$n], "\n";
    }
}
print OUTSTAT "Number of reads mapped to Chromosome X is:\t", $ChrMapcount[24], "\n";
print OUTSTAT "Number of reads mapped to Chromosome Y is:\t", $ChrMapcount[25], "\n";
print OUTSTAT "Number of reads mapped to Mitochondria is:\t", $ChrMapcount[26], "\n \n";

print OUTSTAT "Number of Uniquely mapped reads is:\t", $ChrMapcount[0], "\n";
print OUTSTAT "Number of NonUniquely mapped reads is:\t", $NonUniqcount, "\n";
print OUTSTAT "Number of QC (Low Quality reads) is:\t", $QCcount, "\n";
print OUTSTAT "Number of NM (nonmappable reads) is:\t", $NMcount, "\n";
print OUTSTAT "Number of unknown results (printed in separate file) is:\t", $Unknowncount, "\n";
print OUTSTAT "Number of weird reads (should be same as previous line) is:\t", $weirdcount, "\n";
print OUTSTAT "Number of Total mapped reads is:\t", $totalcount, "\n";



#########################################################
# Sends progress message and closes necessary I/O files #
#########################################################

close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;
close OUT9;
close OUT10;
close OUT11;
close OUT12;
close OUT13;
close OUT14;
close OUT15;
close OUT16;
close OUT17;
close OUT18;
close OUT19;
close OUT20;
close OUT21;
close OUT22;
close OUT23;
close OUTX;
close OUTY;
close OUTM;
close OUTNM;
close OUTNONUNIQ;
close OUTSTAT;



################################
# Sorts all reads within files #
################################

my $filehead = $outfile . "/" . $outfile;

for (my $num = 1; $num < 24; $num++)
{
     if ($ChrMapcount[$num] > 0)
     {
          my $commandinput = "sort -n +1 -2 " . $filehead . "_chr" . $num . ".bed > " . $filehead . "_Chr" . $num . ".bed";
          `$commandinput`;
     }
     my $commandinput = "rm $filehead" . "_chr" . $num . ".bed"; 
     `$commandinput`;
}
my $num = "X";
my $commandinput = "sort -n +1 -2 " . $filehead . "_chr" . $num . ".bed > " . $filehead . "_Chr" . $num . ".bed";
`$commandinput`;
$commandinput = "rm $filehead" . "_chr" . $num . ".bed"; 
`$commandinput`;

$num = "Y";
$commandinput = "sort -n +1 -2 " . $filehead . "_chr" . $num . ".bed > " . $filehead . "_Chr" . $num . ".bed";
`$commandinput`;
$commandinput = "rm $filehead" . "_chr" . $num . ".bed"; 
`$commandinput`;

$num = "M";
$commandinput = "sort -n +1 -2 " . $filehead . "_chr" . $num . ".bed > " . $filehead . "_Chr" . $num . ".bed";
`$commandinput`;
$commandinput = "rm $filehead" . "_chr" . $num . ".bed"; 
`$commandinput`;

$commandinput = "gzip " . $outfile . "/" . $outfile . "_NM.fq";
`$commandinput`;
$commandinput = "gzip " . $outfile . "/" . $outfile . "_NonUnique.fq";
`$commandinput`;
