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
	1) Infile
#    2) Matched outfile
" unless @ARGV == 1;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
#my $output_filename = shift(@ARGV);
#open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


while (<IN>)
{
	chomp;
    my $seq = $_;
    $seq =~ tr/a-z/A-Z/;
#    print $seq;
#    if($seq =~ /[CT][CT]CGTGGTGCG[CT]CG[CT][CT][CT][CT][CT][CT]AA[GT][CT]CGG[CT][CT]/){
    if($seq =~ m/CGTGGTGCG/){
			print $seq , "\n";
			my @tmp = split("CGTGGTGCG", $seq);
			print $tmp[0] , "         " , $tmp[1] , "\n";
 	}
}
close IN;

__END__
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

