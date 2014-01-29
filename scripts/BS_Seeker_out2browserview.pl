#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: BS_Seeker_out2browserview.pl
# Version: 1.0
# Last Updated: 6/14/2013
#
# This is a wrapper that puts has minimal input and automates output. It is NOT
#   memory efficient.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

# Global paths
my $BS_sort_path = "/home/kwdunaway/perl_script/BSSeeker_sort.pl";
my $permeth_path = "/home/kwdunaway/perl_script/BSSeeker_percentmeth_BED_oneside_v3.pl";

# I/O check
die "This script needs the following arguments:
    1) Aligned BSOUT file
    2) Sorted output folder
    3) Permeth output folder
    4) genome (hg18 or mm9)
    5) zipoutput? (Y or N)
" unless @ARGV == 5;

my $alignedfile = shift(@ARGV);
my $sortedout = shift(@ARGV);
my $permethout = shift(@ARGV);
my $genome = shift(@ARGV);
my $zipout = shift(@ARGV);

die "Paramater 5 (zipout) needs to be Y or N, not $zipout" unless ($zipout eq "Y" | $zipout eq "y" | $zipout eq "N" | $zipout eq "n");

my %Chroms;
if($genome eq "hg18"){
	%Chroms = ('0001' => "chr1",
	           '0002' => "chr2",
	           '0003' => "chr3",
	           '0004' => "chr4",
	           '0005' => "chr5",
	           '0006' => "chr6",
	           '0007' => "chr7",
	           '0008' => "chr8",
	           '0009' => "chr9",
	           '0010' => "chr10",
	           '0011' => "chr11",
	           '0012' => "chr12",
	           '0013' => "chr13",
	           '0014' => "chr14",
	           '0015' => "chr15",
	           '0016' => "chr16",
	           '0017' => "chr17",
	           '0018' => "chr18",
	           '0019' => "chr19",
	           '0020' => "chr20",
	           '0021' => "chr21",
	           '0022' => "chr22",
	           '0023' => "chrX",
	           '0024' => "chrY",
	           '0025' => "chrM",);
}
elsif($genome eq "mm9"){
	%Chroms = ('0011' => "chr1",
	           '0012' => "chr2",
	           '0013' => "chr3",
	           '0014' => "chr4",
	           '0015' => "chr5",
	           '0016' => "chr6",
	           '0017' => "chr7",
	           '0018' => "chr8",
	           '0019' => "chr9",
	           '0001' => "chr10",
	           '0002' => "chr11",
	           '0003' => "chr12",
	           '0004' => "chr13",
	           '0005' => "chr14",
	           '0006' => "chr15",
	           '0007' => "chr16",
	           '0008' => "chr17",
	           '0009' => "chr18",
	           '0010' => "chr19",
	           '0021' => "chrX",
	           '0022' => "chrY",
	           '0020' => "chrM",);
}
else{die "$genome is not hg18 or mm9";}

################
# Main Section #
################
my $commandline = "mkdir " . $sortedout;
print $commandline , "\n";
`$commandline`;
$commandline = "mkdir " . $permethout;
print $commandline , "\n";
`$commandline`;

foreach my $key (sort keys %Chroms) {
	my $chrom = $Chroms{$key};
	
	my $soout = $sortedout . "/" . $sortedout . "_" . $chrom . ".txt";
	$commandline = "perl " . $BS_sort_path . " " . $alignedfile . " " . $soout . " " . $key;
	print $commandline , "\n";
	`$commandline`;
	
	my $perout = $permethout . "/" . $permethout . "_" . $chrom . ".bed";
	$commandline = "perl " . $permeth_path . " " . $soout . " " . $perout . " " . $chrom;
	print $commandline , "\n";
	`$commandline`;
	
	if($zipout eq "Y" | $zipout eq "y") {
		$commandline = "gzip " . $perout;
		print $commandline , "\n";
		`$commandline`;
	}
}

