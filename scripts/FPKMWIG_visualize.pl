#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-17-2012
# Script Name: FPKMWIG_visualize.pl
#
# Combines two FPKM files given multipliers for them
#
# Arguments:
#    1) Input prefix for first FPKM WIG files [leave out chrom # and .wig]
#    6) Color (RRR,GGG,BBB format)
#    7) Track name
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: FPKMWIG_visualize.pl 
    1) Input prefix FPKM WIG files [leave out chrom # and .wig]
    2) Output prefix FPKM wig files (they will be zipped and uploadable to genome browser)
    3) Step Size
    4) Color (RRR,GGG,BBB format)
    5) Track name
" unless @ARGV == 5;
my $inputprefix = shift(@ARGV);
my $outprefix = shift(@ARGV);
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
	# Shave off first two header lines
	my $lin = <IN>;
	$lin = <IN>;
	my $outfile = $outprefix . "_Chr" . $chr . ".wig";	
	open(OUT, ">$outfile") or die "cannot open $outfile OUT outfile";
	my $trackname = $tracknameprefix . "_Chr" . $chr;
	
    print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $trackname, "\" description=\"", $trackname, "\" color=", $color, "\n";
	print OUT "variableStep chrom=chr", $chr," span=", $stepsize ,"\n";

	my %HeightHash;
	$lin = <IN>;
	my @firstline = split("\t",substr($lin, 0, -1));
	my $currentpos = $firstline[0];
	my $startpos = $firstline[0];
	$HeightHash{$firstline[0]}=$firstline[1];
	while(<IN>){
		chomp;
		my @line = split("\t",$_);
		$HeightHash{$line[0]}=$line[1];
		$currentpos = $line[0];
		if ($startpos + $stepsize <= $currentpos){
			my $height=0;
			for(my $n = $startpos; $n < ($startpos + $stepsize); $n++){
				if(exists $HeightHash{$n}){	
					$height = $height + $HeightHash{$n};
					delete $HeightHash{$n};
				}
			}
			$height = $height / $stepsize;
			$height = sprintf("%.2f", $height);
			print OUT $startpos , "\t" , $height , "\n";
			$startpos = $currentpos;
		}
	}
	my $height=0;
	for(my $n = $startpos; $n < ($startpos + $stepsize); $n++){
		if(exists $HeightHash{$n}){	
			$height = $height + $HeightHash{$n};
		}
	}
	$height = $height / $stepsize;
	$height = sprintf("%.5f", $height);
	print OUT $startpos , "\t" , $height , "\n";
	close IN;
	close OUT;	
	print "Finished making $outfile, now zipping it.\n";
	my $commandline = "gzip $outfile";
	`$commandline`;
}

