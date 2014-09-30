#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Last Update Date: 7-8-2014
# Version: 2.0
#
# Takes SAM output from BS_Seeker2 and creates percentage methylation BED files that
# can be uploaded to the UCSC genome browser or further analyzed through StochHMM.
#
# PCR duplicate filter: This script takes the longest read that matches an strand and
# position of the same chromosome. If more than one read are the longest, it only takes
# whichever read came first in the SAM file.
#
# The positions in the resulting percent methylation (permeth) BED files are what you
# would get if you go to the following website. For example, if you go here: 
#    http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chrY:59032572,59032573
# if would return CG. However, when you look at the position on the genome browswer, the
# color will only cover 1 base (the 2nd one).
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following parameters:
    1) Input sorted SAM file
    2) Output file
" unless @ARGV == 2;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output = shift(@ARGV);
open(OUT, ">$output") or die "cannot open $output output file";

my $bscg = 0;
my $bsat = 0;
my $rawcg = 0;
my $rawat = 0;
my $total = 0;
my $milestone = 10000000;

# Check for @ symbols, and then process first line
my $header = <IN>;
while($header =~ /^@/g) {$header = <IN>};
my @line = split("\t", $header);
my @num = ($line[9] =~ /[CG]/g);
my $count = @num;
$bscg += $count;
$total += $count;
@num = ($line[9] =~ /[AT]/g);
$count = @num;
$bsat += $count;
$total += $count;
my @dash = split("_", $line[15]);
@num = ($dash[1] =~ /[CG]/g);
$count = @num;
$rawcg += $count;
@num = ($dash[1] =~ /[AT]/g);
$count = @num;
$rawat += $count;
if($total > $milestone)
{
	print "Total bases has crossed $milestone\n";
	$milestone = $milestone + 10000000;
}

while(<IN>){
	my @line = split("\t", $_);
	my @num = ($line[9] =~ /[CG]/g);
	my $count = @num;
	$bscg += $count;
	$total += $count;
	@num = ($line[9] =~ /[AT]/g);
	$count = @num;
	$bsat += $count;
	$total += $count;
	my @dash = split("_", $line[15]);
	@num = ($dash[1] =~ /[CG]/g);
	$count = @num;
	$rawcg += $count;
	@num = ($dash[1] =~ /[AT]/g);
	$count = @num;
	$rawat += $count;
	if($total > $milestone)
	{
		print "Total bases has crossed $milestone\n";
		$milestone = $milestone + 10000000;
	}
}

print OUT "BS Read Totals:\n\tCG:", $bscg/$total, "\tAT:", $bsat/$total, "\n\n";
print OUT "Raw Read Totals:\n\tCG:", $rawcg/$total, "\tAT:", $rawat/$total;

close IN;
close OUT;


__END__
	my @line = split("\t", $_);
	my @num = ($line[9] =~ /C/g);
	my $count = @num;
	$bscg += $count;
	$total += $count;
	@num = ($line[9] =~ /G/g);
	$count = @num;
	$bscg += $count;
	$total += $count;
	@num = ($line[9] =~ /A/g);
	$count = @num;
	$bsat += $count;
	$total += $count;
	@num = ($line[9] =~ /T/g);
	$count = @num;
	$bsat += $count;
	$total += $count;
	my @dash = split("_", $line[15]);
	@num = ($dash[1] =~ /C/g);
	$count = @num;
	$rawcg += $count;
	@num = ($dash[1] =~ /G/g);
	$count = @num;
	$rawcg += $count;
	@num = ($dash[1] =~ /A/g);
	$count = @num;
	$rawat += $count;
	@num = ($dash[1] =~ /T/g);
	$count = @num;
	$rawat += $count;

