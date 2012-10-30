#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-17-2012
# Script Name: VarStepWIG_to_FPKMWIG.pl
#
# Calculates FPKM from raw WIG files using read length and read count (in millions)
#
# Arguments:
#    1) Input VarStepWIG file (format: position [tab] height)
#    2) Output FPKM wig file 
#    3) Read Length
#    4) Read Count (in millions)
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: RPKM.pl 
    1) Input VarStepWIG file prefix (without #.wig) (format: position [tab] height)
    2) Output FPKM wig file 
    3) Read Length
    4) Read Count (in millions)
" unless @ARGV == 4;
my $inputWIGprefix = shift(@ARGV);
my $outprefix = shift(@ARGV);
my $readlength = shift(@ARGV);
my $readcount = shift(@ARGV);

my @Chr;             # array that contains all the the names of the mouse chromosomes
for (my $n = 1; $n< 20; $n++)
{
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");



#############################################
# Calculate and print OUT FPKM from file(s) #
#############################################

while(@Chr)
{
	my $chr = shift(@Chr);
	print "Now Converting: Chr$chr\n";
	my $inputWIG = $inputWIGprefix . $chr . ".wig";
	my $outfile = $outprefix . "_Chr" . $chr . ".wig";
	open(INWIG, "<$inputWIG") or die "cannot open $inputWIG INWIG infile";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	my $lin = <INWIG>;
	print OUT $lin;
	$lin = <INWIG>;
	print OUT $lin;
	while(<INWIG>){
		chomp;
	    my @line = split ("\t", $_);
	    if(exists $line[1]){
	    	my $FPKM = ($line[1] * 1000) / ($readlength * $readcount);
		$FPKM = sprintf("%.4f",$FPKM);
		print OUT $line[0], "\t" , $FPKM,"\n";
		}
	}
	close OUT;
	close INWIG;
}
