#!/usr/bin/perl 
use strict; use warnings;
#use Scalar::Util;
# qw(looks_like_number);


###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 1-4-2013
# Script Name: SPKMWIG_visualize.pl
#
# Takes a SPKM file and creates a WIG file.
#
# Arguments:
#    (see below)
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: SPKMWIG_visualize.pl 
    1) Input prefix FPKM WIG files [leave out chrom # and .wig]
    2) Output prefix FPKM wig files (they will be zipped and uploadable to genome browser)
    3) Window Size
    4) Step Size
    5) Color (RRR,GGG,BBB format)
    6) Track name (do NOT include window or step size, they will be added to description)
" unless @ARGV == 6;
my $inputprefix = shift(@ARGV);
my $outprefix = shift(@ARGV);
my $windowsize = shift(@ARGV);
my $stepsize = shift(@ARGV);
my $color = shift(@ARGV);
my $tracknameprefix = shift(@ARGV);

my @Chr;             # array that contains all the the names of the mouse chromosomes
for (my $n = 1; $n< 20; $n++)
{
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");



#############################################
# Visualize and print OUT FPKM from file(s) #
#############################################

while(@Chr)
{
	my $chr = shift(@Chr);
	print "Now Visualizing: Chr$chr\n";
	my $inputWIG = $inputprefix . $chr . ".wig";
	open(IN, "<$inputWIG") or die "cannot open $inputWIG IN infile";
	my $outfile = $outprefix . "_Chr" . $chr . ".wig";	
	open(OUT, ">$outfile") or die "cannot open $outfile OUT outfile";
	my $trackname = $tracknameprefix . "_Chr" . $chr;
	my $description = $tracknameprefix . "_Window:" . $windowsize . "_Step:" . $stepsize . "_Chr" . $chr;
    print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $trackname, "\" description=\"", $description, "\" color=", $color, "\n";
	print OUT "variableStep chrom=chr", $chr," span=", $stepsize ,"\n";

	my %HeightHash;
	my $currentpos = 0;
	my $startpos = 0;
	while(<IN>){
		chomp;
		my @line = split("\t",$_);

		# If a header line, go to next line
		if (! /^[0-9]/){ next;}
#		if(! looks_like_number($line[0])) 
		
		# Change current position
		$currentpos = $line[0];
		
		# If a certain distance through, print to outfile and flush hash
		while($startpos + $windowsize < $currentpos){
			if(keys(%HeightHash) == 0){
				$startpos = $currentpos;
				last;
			}
			my $height=0;	
			for(my $n = $startpos; $n < ($startpos + $windowsize); $n++){
				if(exists $HeightHash{$n}){	
					$height = $height + $HeightHash{$n};
					if($n < $startpos + $stepsize){
						delete $HeightHash{$n};
					}
				}
			}
			$height = $height / $windowsize;
			$height = sprintf("%.2f", $height);
			my $printpos = int($startpos + ($windowsize/2));
			if ($height != 0){
				print OUT $printpos , "\t" , $height , "\n";
			}
			$startpos = $startpos + $stepsize;
		}
		
		# Add position and height to Hash
		$HeightHash{$line[0]}=$line[1];
		if($currentpos % 10000 == 0) {print "currentpos = $currentpos \n";}
	}

	while(%HeightHash){
		my $height=0;	
		for(my $n = $startpos; $n < ($startpos + $windowsize); $n++){
			if(exists $HeightHash{$n}){	
				$height = $height + $HeightHash{$n};
				if($n < $startpos + $stepsize){
					delete $HeightHash{$n};
				}
			}
		}
		$height = $height / $windowsize;
		$height = sprintf("%.2f", $height);
		my $printpos = int($startpos + ($windowsize/2));
		if ($height != 0){
			print OUT $printpos , "\t" , $height , "\n";
		}
		$startpos = $startpos + $stepsize;
	}

	close IN;
	close OUT;	
	print "Finished making $outfile, now zipping it.\n";
	my $commandline = "gzip $outfile";
	`$commandline`;
}


