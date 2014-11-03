#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: BS_Seeker_out2browserview_v6_BS2.pl
# Version: 6.0
# Last Updated: 6/17/2013
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
my $BS_sort_path = "/data/scratch/programs/perl_script/BSSeeker2_sort.pl";
my $permeth_path = "/data/scratch/programs/perl_script/BSSeeker_percentmeth_BED_v4.pl";

# I/O check
die "This script needs the following arguments:
    1) Aligned BSOUT file
    2) Sorted output folder
    3) Permeth output folder
    4) genome (hg18, hg19, rn4, or mm9)
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
elsif($genome eq "hg19"){
	%Chroms = ('chr1' => "chr1",
	           'chr2' => "chr2",
	           'chr3' => "chr3",
	           'chr4' => "chr4",
	           'chr5' => "chr5",
	           'chr6' => "chr6",
	           'chr7' => "chr7",
	           'chr8' => "chr8",
	           'chr9' => "chr9",
	           'chr10' => "chr10",
	           'chr11' => "chr11",
	           'chr12' => "chr12",
	           'chr13' => "chr13",
	           'chr14' => "chr14",
	           'chr15' => "chr15",
	           'chr16' => "chr16",
	           'chr17' => "chr17",
	           'chr18' => "chr18",
	           'chr19' => "chr19",
	           'chr20' => "chr20",
	           'chr21' => "chr21",
	           'chr22' => "chr22",
	           'chrX' => "chrX",
	           'chrY' => "chrY",
	           'chrM' => "chrM",);
}
elsif($genome eq "rn4"){
	%Chroms = ('chr1' => "chr1",
	           'chr2' => "chr2",
	           'chr3' => "chr3",
	           'chr4' => "chr4",
	           'chr5' => "chr5",
	           'chr6' => "chr6",
	           'chr7' => "chr7",
	           'chr8' => "chr8",
	           'chr9' => "chr9",
	           'chr10' => "chr10",
	           'chr11' => "chr11",
	           'chr12' => "chr12",
	           'chr13' => "chr13",
	           'chr14' => "chr14",
	           'chr15' => "chr15",
	           'chr16' => "chr16",
	           'chr17' => "chr17",
	           'chr18' => "chr18",
	           'chr19' => "chr19",
	           'chr20' => "chr20",
	           'chrX' => "chrX",
	           'chrY' => "chrY",
	           'chrM' => "chrM",);
}
else{die "$genome is not hg18, hg19, rn4, or mm9";}

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
	$commandline = "perl " . $permeth_path . " " . $soout . " " . $perout . " " . $chrom . " " . $genome . " CG combined";
	print $commandline , "\n";
	`$commandline`;
	
	if($zipout eq "Y" | $zipout eq "y") {
		$commandline = "gzip " . $perout;
		print $commandline , "\n";
		`$commandline`;
	}
}

