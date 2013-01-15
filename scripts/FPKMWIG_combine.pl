#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-17-2012
# Script Name: FPKMWIG_combine.pl
#
# Combines two FPKM files given multipliers for them
#
# Arguments:
#    1) Input prefix for first FPKM WIG files [leave out chrom # and .wig]
#    2) Multiplier for first Input file
#    3) Input prefix for second FPKM WIG files [leave out chrom # and .wig]
#    4) Multiplier for second Input file
#    5) Output prefix FPKM combined wig file
#    6) Color (RRR,GGG,BBB format)
#    7) Track name
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: FPKMWIG_combine.pl 
    1) Input prefix for first FPKM WIG files [leave out chrom # and .wig]
    2) Multiplier for first Input file
    3) Input prefix for second FPKM WIG files [leave out chrom # and .wig]
    4) Multiplier for second Input file
    5) Output prefix FPKM combined wig file
    6) Color (RRR,GGG,BBB format)
    7) Track name
" unless @ARGV == 7;
my $inputfirstprefix = shift(@ARGV);
my $inputfirstmultiplier = shift(@ARGV);
my $inputsecondprefix = shift(@ARGV);
my $inputsecondmultiplier = shift(@ARGV);
my $outprefix = shift(@ARGV);
my $color = shift(@ARGV);
my $tracknameprefix = shift(@ARGV);

my @Chr;             # array that contains all the the names of the mouse chromosomes
for (my $n = 1; $n< 20; $n++){
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");

#push(@Chr, "7");

#############################################
# Calculate and print OUT FPKM from file(s) #
#############################################

while(@Chr)
{
	my $chr = shift(@Chr);
	print "Now Combining: Chr$chr\n";
	my $inputfirstWIG = $inputfirstprefix . $chr . ".wig";
	open(INFIRST, "<$inputfirstWIG") or die "cannot open $inputfirstWIG INFIRST infile";
	# Shave off first two header lines
	my $lin = <INFIRST>;	$lin = <INFIRST>;
	my $inputsecondWIG = $inputsecondprefix . $chr . ".wig";
	open(INSECOND, "<$inputsecondWIG") or die "cannot open $inputsecondWIG INSECOND infile";
	# Shave off first two header lines
	$lin = <INSECOND>;	$lin = <INSECOND>; 
	my $outfile = $outprefix . "_Chr" . $chr . ".wig";	
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	my $trackname = $tracknameprefix . "_Chr" . $chr;
	
    print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $trackname, "\" description=\"", $trackname, "\" color=", $color, "\n";
	print OUT "variableStep chrom=chr", $chr," span=1\n";

	my $inloop = 1;
	my $firstpos = 0;
	my $firstheight = 0;
	my $secondpos = 0;
	my $secondheight = 0;
	my $printedlinecount = 0;
	
	while($inloop)	{
		my $line = <INFIRST>;
		my @linefirst = split("\t",substr($line, 0, -1));
#		print "infirst line:\t" , $linefirst[0] , "\t" , $linefirst[1] , "\n";
		if(exists $linefirst[1]){
			$firstpos = $linefirst[0];
			$firstheight = $linefirst[1];
		}
		else{
			print "End of infirst file. Last position was $firstpos \n";
			$firstpos = 999999999999;
			$firstheight = 999999999999;
		}
		while($secondpos < $firstpos){
			$line = <INSECOND>;
			my @linesecond = split("\t",substr($line, 0, -1));
#			print "insecond line:\t" , $linesecond[0] , "\t" , $linesecond[1] , "\n";
			if(exists $linesecond[1]){
				$secondpos = $linesecond[0];
				$secondheight = $linesecond[1];
				if($secondpos < $firstpos){
					my $position = $secondpos;
					my $height = $secondheight * $inputsecondmultiplier;
#						if ($height > 0) {print $height , "\n";} 
					if($height != 0){
#							if ($height > 0) {print $height , "\t printed", "\n";} 
						$printedlinecount++;
						if ($printedlinecount % 2500000 == 0){print "Printed ", $printedlinecount , " lines\n"; }
						print OUT $position , "\t" , $height , "\n";
					}
				}
			}
			else{
				print "End of insecond file. Last position was $secondpos \n";
				$secondpos = 999999999999;
				$secondheight = 999999999999;
			}
		}
		if ($firstpos == 999999999999 && $secondpos == 999999999999){ 
			$inloop = 0; 
			last;
		}
		else{
			if($firstpos == $secondpos){
				my $position = $firstpos;
				my $height = ($firstheight * $inputfirstmultiplier) + ($secondheight * $inputsecondmultiplier);
#					if ($height > 0) {print $height , "\n";} 
				if($height != 0){
#						if ($height > 0) {print $height , "\t printed", "\n";} 
					$printedlinecount++;
					if ($printedlinecount % 2500000 == 0){print "Printed ", $printedlinecount , " lines\n"; }
					print OUT $position , "\t" , $height , "\n";
				}
			}
			else{
				my $position = $firstpos;
				my $height = $firstheight * $inputfirstmultiplier;
#					if ($height > 0) {print $height , "\n";} 
				if($height != 0){
#						if ($height > 0) {print $height , "\t printed", "\n";} 
					$printedlinecount++;
					if ($printedlinecount % 2500000 == 0){print "Printed ", $printedlinecount , " lines\n"; }
					print OUT $position , "\t" , $height , "\n";
				}
			}
		}
	}	
	close INFIRST;
	close INSECOND;
	close OUT;
}
