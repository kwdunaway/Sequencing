#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 10/31/2014
#
# This script takes a fasta read file and makes a Bisulfite converted
# fastq file with mutations
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) input fasta file
    2) output fastq file
    3) readlength (ex: 100)
    4) mutation rate (ex: .01)
    5) CG methylation rate (ex: .75)
    6) CH methylation rate (ex: .02)
" unless @ARGV == 6;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
my $readlength = shift(@ARGV);
my $mutrate = shift(@ARGV);
my $CGmethrate = shift(@ARGV);
my $CHmethrate = shift(@ARGV);

my $strandrate = .5;      #rate of finding reverse strand
my $qualline = "";
for(my $i = 0; $i < $readlength; $i++){
		$qualline = $qualline . "I";
}
$readlength = $readlength - 1;
my $ycoord = 1;

#############
# Main Loop #
#############

while (<IN>)
{
	chomp;
	my $name = substr($_ , 1);
	my $seq = <IN>;
	chop($seq);
	$seq = uc$seq;
	my $strand = "forward";
	if (rand() < $strandrate){
		$seq = revcomp($seq);
		$strand = "reverse";
	}
	
	my $seqline = "";
	my $mutations = "";
	my $CGmethsites = "";
	my @sequence = split ("", $seq);
	for(my $n = 0; $n < $readlength; $n++){
		my $base = $sequence[$n];
		if(rand() < $mutrate){
			#Mutate base
			$mutations = $mutations . ($n +1) . $base . ">";
			$base = mutate($base);
			$mutations = $mutations . $base . ",";
		}
		if($base eq "C"){
			if($sequence[$n+1] eq "G"){
				$CGmethsites = $CGmethsites . ($n +1) . ">";
				if(rand() < $CGmethrate) {$seqline = $seqline . "C"; $CGmethsites = $CGmethsites . "C" . ",";}
				else {$seqline = $seqline . "T"; $CGmethsites = $CGmethsites . "T" . ",";}
			}
			else{
				if(rand() < $CHmethrate) {$seqline = $seqline . "C";}
				else {$seqline = $seqline . "T";}
			}
		}
		else{$seqline = $seqline . $base;}
	}
	$seqline = $seqline . $sequence[$readlength];
	print OUT "\@" , $name , "_" ,  $strand , "_" , $mutations , " HS2:306:C1A3MACXX:6:1101:1:$ycoord 1:N:0:\n" , 
		$seqline, "\n" , 
		"+" , $CGmethsites , 
		"\n", $qualline , "\n";
	$ycoord++;
}
close IN;
close OUT;



###############
# Subroutines #
###############

sub revcomp 
{
	my ($seq) = @_;
	$seq = uc($seq); # convert to uppercase
	
	# convert to matching base pair
	$seq =~ tr/(G)(C)(A)(T)(W)(S)(M)(K)(R)(Y)(B)(V)(D)(H)/(C)(G)(T)(A)(S)(W)(K)(M)(Y)(R)(V)(B)(H)(D)/;
	my $rseq = scalar reverse($seq);
	

	return $rseq;
}

sub mutate 
{
	my ($base) = @_;
	my $mutbase = $base;
	my @chars = ("A","T","C","G");
	while($mutbase eq $base){
		$mutbase = $chars[rand @chars];
	}
	return($mutbase)
}