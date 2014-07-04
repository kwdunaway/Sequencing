#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 6-3-2014
# Version: 1.0
#
# Takes multiple bed prefixes and combines 
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "BSSeeker2_permeth.pl needs the following parameters:
    1) Output prefix
    2) Genome (hg19, mm9)
    3) Track name prefix
    4-?) Percent Methylation track prefixes (with chr but not #)
" unless @ARGV > 4;

my $outprefix = shift(@ARGV);
my $genome = shift(@ARGV);
my $bedprefix = shift(@ARGV);
my @infiles = @ARGV;

my $autosome_chromnumber = 0;
if($genome eq "hg19"){$autosome_chromnumber = 22;}
elsif($genome eq "mm9"){$autosome_chromnumber = 19;}
else {die "$genome needs to be hg19 or mm9";}


# Builds the chromosome array for the genome you are looking at
my @chromosomes;
for (my $c = 1; $c <= $autosome_chromnumber; $c++) {push(@chromosomes,$c);}
push(@chromosomes , "X");
push(@chromosomes , "Y");
push(@chromosomes , "M");



#############
# Main Loop #
#############

while(@chromosomes){
	my $chrom = shift(@chromosomes);
	my %Methylation;
	for(my $n = 0; $n < @infiles; $n++){
		my $infile = $infiles[$n] . $chrom . ".bed";
		open(IN, "<$infile") or next;
		print "Adding $infile to hash\n";
		my $headerline = <IN>;
		while(<IN>){
			my @line = split("\t", $_);
			my $start = $line[1];
			if (exists $Methylation{$start}) {
				my $methcell = $line[3];
				my @methcur = split("-", $methcell);
				my @methprev = split("-", $Methylation{$start});
				my $methcount = $methcur[1] + $methprev[1];
				my $methper = (($methcur[0] * $methcur[1]) + ($methprev[0] * $methprev[1])) / $methcount;
				$methper = sprintf("%.2f", $methper);
				$Methylation{$start} = $methper . "-" . $methcount;
			}
			else {$Methylation{$start}=$line[3];}
		}
		close IN;
	}
	my $outfile = $outprefix . "_chr" . $chrom . ".bed";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	print "Writing to $outfile \n";
	print OUT "track name=" , $bedprefix , $chrom, " description=" , $bedprefix, $chrom, " useScore=0 itemRgb=On db=" , $genome , "\n";
	foreach my $start (sort { $a <=> $b }  keys %Methylation) {
		my $end = $start + 1;
		my @methinfo = split("-", $Methylation{$start});
		my $methperc = $methinfo[0];
		my $color = "0,0,0"; #black
		if ($methperc > 0 && $methperc <= .6) {$color = "27,74,210";} #blue
		elsif ($methperc > .6 && $methperc <= .8) {$color = "27,210,57";} #green
		elsif ($methperc > .8) {$color = "210,27,27";} #red
		print OUT "chr" , $chrom , "\t" , $start , "\t" , $end , "\t" , $Methylation{$start}, "\t" , "0\t+\t0\t0\t" , $color , "\n";
	}
	close OUT;
}
