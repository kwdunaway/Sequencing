#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: Line1_analysis.pl
# Version: 1.0
# Last Updated: March 3, 2014
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the input pattern.  Then, it puts those reads in a separate file.  Also, the
# script quantifies methylation of these sequences across three potential 
# methylation sites.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
	1) Sequence (with ambiguous symbols: 
    2) Stats outfile
    3) Reads that match suffix (adds to the end of every input file for output)
    4-?) Input fastq file
" unless @ARGV > 3;

my $in_sequence = shift(@ARGV);
my $stats_filename = shift(@ARGV);
open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";
my $out_suffix = shift(@ARGV);
my @Infiles = @ARGV;

W	A	T
S	C	G
M	A	C		
K	G	T
R	A	G	
Y	C	T
##############################################
# Convert in_sequence to the 4 possibilities #
##############################################

# Forward sequence BS converted on same (forward) strand
my $seq_fs = ;

# Forward sequence BS converted on opposite (reverse) strand
my $seq_fo = ;

# Reverse sequence BS converted on same (reverse) strand
my $seq_rs = ;

# Reverse sequence BS converted on opposite (forward) strand
my $seq_ro = ;


#################
# In Files Loop #
#################

while(@Infiles){
	my $infile = shift(@Infiles);
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "\nStarting file: $infile\nNumber or Reads matched so far:\t";
	my $output_filename = $infile . $out_suffix;
	open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


###############################
# Initialization of Variables #
###############################

	my $meth1 = 0;
	my $meth2 = 0;
	my $meth3 = 0;
	my $meth4 = 0;
	my $unmeth1 = 0;
	my $unmeth2 = 0;
	my $unmeth3 = 0;
	my $unmeth4 = 0;
	my $SNPG = 0;
	my $SNPT = 0;
	my $LINE1 = 0;

#############################
# Main Search Function Loop #
#############################

	while (<IN>)
	{
	    my $ID = $_;
   		my $seq = <IN>;
	    my $third = <IN>;
	    my $quality = <IN>;
	    if($seq =~ /(TT([CT])GTGGTG([CT])GT([CT])GTTTTTTAA([GT])T([CT])GGTT)/){
			if($2 eq "C") {$meth1++;} elsif($2 eq "T") {$unmeth1++;} else{die "Site1: $2";}
			if($3 eq "C") {$meth2++;} elsif($3 eq "T") {$unmeth2++;} else{die "Site2: $3";}
			if($4 eq "C") {$meth3++;} elsif($4 eq "T") {$unmeth3++;} else{die "Site3: $4";}
			if($6 eq "C") {$meth4++;} elsif($6 eq "T") {$unmeth4++;} else{die "Site4: $6";}
			if($5 eq "G") {$SNPG++;} elsif($5 eq "T") {$SNPT++;} else{die "A/T: $5";}
	    	$LINE1++;
	    	if($LINE1 % 100 == 0) {print $LINE1 , "\n";}
			print OUT $ID , $seq , $third , $quality;
    	}
	}

#######################
# Print Stats to File #
#######################

	my $perc1 = ($meth1 / $LINE1) * 100; 
	my $perc2 = ($meth2 / $LINE1) * 100; 
	my $perc3 = ($meth3 / $LINE1) * 100; 
	my $perc4 = ($meth4 / $LINE1) * 100; 
	print STATS "$infile\n";
	print STATS "Total lines: $LINE1\n";
	print STATS "Sites\tMethylated\tUnmethylated\tPercentage\n";
	print STATS "Site1\t$meth1\t$unmeth1\t$perc1\n";
	print STATS "Site2\t$meth2\t$unmeth2\t$perc2\n";
	print STATS "Site3\t$meth3\t$unmeth3\t$perc3\n";
	print STATS "Site4\t$meth4\t$unmeth4\t$perc4\n";
	print STATS "G or T\t$SNPG\t$SNPT\n\n";
	
	close IN;
	close OUT;
}    

close STATS;

__END__
Search string:

TTYGT
GGTGY
GTYGT
TTTTT
A[A(t?)]GTY
GGTTT

