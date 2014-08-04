#!/usr/bin/perl
use strict;
use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 7-25-2013
# Script Name: Remove_dup_readsv2.pl
#
# This script removes any duplicate reads from bed files in a bed dir.
#
# Arguments: See Below
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0.pl needs the following parameters:
    1) Input bed dir
    2) Output bed dir
    3) Max reads per duplication
" unless @ARGV == 3;

my $InDir = shift(@ARGV);
my $OutDir = shift(@ARGV);
my $duplimit = shift(@ARGV);
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

while(@Bedfiles){
	my $infilename = $InDir . "/" . $Bedfiles[0];
	my $outfilename = $OutDir . "/" . $Bedfiles[0];
	open(IN, "<$infilename") or die "cannot open $infilename IN infile";
	open(OUT, ">$outfilename") or die "cannot open $outfilename OUT outfile";
	my $lastline = <IN>;
	chop($lastline);
	my $dupnumber = 0;
	while(<IN>){
		chomp;
		if($lastline eq $_){
			$dupnumber++;
			if ($dupnumber < $duplimit){
				print OUT $lastline , "\n";
			}			
		}
		else{
			print OUT $lastline , "\n";
			$dupnumber = 0;
		}
		$lastline = $_;
	}
	print OUT $lastline;
	close IN;
	close OUT;
	shift(@Bedfiles);
}

