#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: Line1_analysis.pl
# Version: 1.2
# Last Updated: March 4, 2014
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the Line1 pattern.  Then, it puts those reads in a separate file.  Also, the
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
    1) Stats outfile (also makes a read file with counts named seqcount_statsoutfile)
    2) Outsuffix (to output reads that have the desired sequence)
    3-?) Input fastq file
" unless @ARGV > 2;

my $stats_filename = shift(@ARGV);
open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";
my $seqcount = "seqcount_" . $stats_filename;
open(SEQOUT, ">$seqcount") or die "cannot open $seqcount outfile";
my $out_suffix = shift(@ARGV);
my @Infiles = @ARGV;

#my $output_filename = shift(@ARGV);
#open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


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
	
	my %Reads;

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
			
			if (exists $Reads{$seq}){
				$Reads{$seq}++;
			}
			else{
				if($2 eq "C") {$meth1++;} elsif($2 eq "T") {$unmeth1++;} else{die "Site1: $2";}
				if($3 eq "C") {$meth2++;} elsif($3 eq "T") {$unmeth2++;} else{die "Site2: $3";}
				if($4 eq "C") {$meth3++;} elsif($4 eq "T") {$unmeth3++;} else{die "Site3: $4";}
				if($6 eq "C") {$meth4++;} elsif($6 eq "T") {$unmeth4++;} else{die "Site4: $6";}
				if($5 eq "G") {$SNPG++;} elsif($5 eq "T") {$SNPT++;} else{die "A/T: $5";}
				$Reads{$seq} = 1;
				$LINE1++;
	    		if($LINE1 % 20 == 0) {print $LINE1 , "\n";}
			}
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
	my $perc5 = ($SNPG / $LINE1) * 100; 
	print STATS "$infile\n";
	print STATS "Total lines: $LINE1\n";
	print STATS "Sites\tMethylated\tUnmethylated\tPercentage\n";
	print STATS "Site1\t$meth1\t$unmeth1\t$perc1\n";
	print STATS "Site2\t$meth2\t$unmeth2\t$perc2\n";
	print STATS "Site3\t$meth3\t$unmeth3\t$perc3\n";
	print STATS "Site4\t$meth4\t$unmeth4\t$perc4\n";
	print STATS "G or T\t$SNPG\t$SNPT\t$perc5\n\n";

	print SEQOUT "$infile\n";
	foreach my $key (sort { $Reads{$b} <=> $Reads{$a} } keys %Reads) {
        print SEQOUT $Reads{$key}, "\t" , $key , ;
    }
    print SEQOUT "\n\n";

}    


__END__
Search string:

TTYGT
GGTGY
GTYGT
TTTTT
A[A(t?)]GTY
GGTTT

