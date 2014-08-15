#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway and Roy Chu
# Email: kwdunaway@ucdavis.edu AND rgchu@ucdavis.edu
# Date: 8-12-2014
#
# This script deals with certain genome browser errors, allowing replacement of 
# the header and gets rid of positions past the reference chromosomes.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: $0
    1) Input File
    2) Output File
    3) Length of reads (window size)
    4) Step size (sliding size)
" unless @ARGV == 4;

my $input = shift(@ARGV);	
open(IN, "<$input") or die "Error: Missing Chromosome, cannot open $input input file";
my $outputfile = shift(@ARGV);
open(OUT, ">$outputfile") or die "Error: cannot open $outputfile output file";
my $windowsize = shift(@ARGV);
my $stepsize = shift(@ARGV);

# Global variables
my $sequence = "";
my $ycoord = 0;
my $xcoord = 0;
my $qualline = "";
for(my $i = 0; $i < $windowsize; $i++){
		$qualline = $qualline . "I";
}

while(<IN>)
{
	chomp;
	if($_ =~ /^>/){
		if(length($sequence) >= $windowsize){
			my $seqline = substr($sequence, 0, $windowsize);
			if($seqline !~ /N/){
				$ycoord++;
				print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $seqline, "\n" , "+\n", $qualline , "\n";			
				my $rcseqline = revcomp($seqline);
				$ycoord++;
				print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $rcseqline, "\n" , "+\n", $qualline , "\n";
			}
			my $lastwindow = $windowsize * -1;
			$seqline = substr($sequence, 0, $lastwindow);
			if($seqline !~ /N/){
				$ycoord++;
				print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $seqline, "\n" , "+\n", $qualline , "\n";
				my $rcseqline = revcomp($seqline);
				$ycoord++;
				print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $rcseqline, "\n" , "+\n", $qualline , "\n";
			}
		}
		$sequence = "";
		$xcoord++;
		$ycoord = 0;
		next;
	}
	$sequence = $sequence . uc($_);
	while(length($sequence) > $windowsize + $stepsize){
		my $seqline = substr($sequence, 0, $windowsize);
		if($seqline !~ /N/){
			$ycoord++;
			print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $seqline, "\n" , "+\n", $qualline , "\n";
			
			my $rcseqline = revcomp($seqline);
			$ycoord++;
			print OUT "\@HS2:306:C1A3MACXX:6:1101:$xcoord:$ycoord 1:N:0:\n" , $rcseqline, "\n" , "+\n", $qualline , "\n";
		}
		$sequence = substr($sequence, $stepsize);
	}
}
close IN;
close OUT;

################################################################################################
#                 Subroutines     
################################################################################################

################################################################################################
#    Complement of DNA Sequence 
#     Input: DNA Sequence      
#    Output: Complementary DNA Sequence      
#   					
#    Note: translates ambiguous codons into complementary ambiguous codons      
################################################################################################

sub revcomp 
{
	my ($seq) = @_;
	$seq = uc($seq); # convert to uppercase
	
	# convert to matching base pair
	$seq =~ tr/(G)(C)(A)(T)(W)(S)(M)(K)(R)(Y)(B)(V)(D)(H)/(C)(G)(T)(A)(S)(W)(K)(M)(Y)(R)(V)(B)(H)(D)/;
	my $rseq = scalar reverse($seq);
	

	return $rseq;
}
