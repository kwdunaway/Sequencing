#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 6-6-2012
# Script Name: extend_read_length.pl
#
# Takes a folder of bed files and creates a new one with read length extended directionally
# based on arguments input on command line.
#
# Arguments:
#    1) Input Bed File prefix (ex: DY_Chr)
#    2) Output Bed prefix without Chr (ex: NewDY will make NewDY_Chr*.bed files) 
#    3) Read Length extension (ex: 97)
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: extend_read_length.pl 
    1) Input Bed Dir
    2) Output Bed Dir 
    3) Read Length extension (ex: 97)
" unless @ARGV == 3;

my $InDir = shift(@ARGV);
my $OutDir = shift(@ARGV);
my $readlengthextension = shift(@ARGV);

my @Bedfiles;

my $dir = $InDir;
opendir(DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {

    # We only want files
    next unless (-f "$dir/$file");

    # Use a regular expression to find files ending in .txt
    next unless ($file =~ m/\.bed$/);

#    print "$file\n";
    push (@Bedfiles,$file);
}
closedir(DIR);



###################
# Main Input loop #
###################

while(@Bedfiles){
	my $infilename = $InDir . "/" . $Bedfiles[0];
	my $outfilename = $OutDir . "/" . $Bedfiles[0];
	open(IN, "<$infilename") or die "cannot open $infilename IN infile";
	open(OUT, ">$outfilename") or die "cannot open $outfilename OUT outfile";

#	my $chr = shift(@Chr);
#	my $inputfile = $inputbedprefix . $chr . ".bed";
#	open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";
#	my $outfile = $outputbedprefix . "_Chr" . $chr . ".bed";
#	open(OUT, ">$outfile") or die "cannot open $outfile outfile";

#	print "Processing $inputfile \n";
	
	print "Starting: " , $Bedfiles[0] , "\n";
	while(<IN>){
		chomp;
    	my @line = split ("\t", $_);
		my $start = $line[1];
		my $end = $line[2];
		if($line[5] eq "+"){
			my $end = $line[2]+$readlengthextension;
			print OUT $line[0],"\t",$line[1],"\t",$end,"\t",$line[3],"\t",$line[4],"\t",$line[5],"\n";
		}
		elsif($line[5] eq "-"){
			my $start = $line[1]-$readlengthextension;
			print OUT $line[0],"\t",$start,"\t",$line[2],"\t",$line[3],"\t",$line[4],"\t",$line[5],"\n";
		}
		else {die "$line[5] does not equal + or - \n";}
	}
	shift(@Bedfiles);
}

