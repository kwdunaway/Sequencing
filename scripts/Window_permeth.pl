#!/usr/bin/perl 
use strict; use warnings; 
use POSIX 'floor';

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu and rgchu@ucdavis.edu
# Date: 7-21-2014 (happy birthday)
#
# This script takes windows (user defined size)
#
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 usage, needs the following parameters: 
    1) Output table file
    2) NONE or CpG island GTF (or bed) file to mask. If no masking, put NONE
    3) Window size
    4) Min # of CpGs per window (otherwise prints NA)
    5) Min # of reads per CpG counted
    6) Min number of files have info
    7,9+) Permeth prefix (leave off chr#.bed)
    8,10+) Name of experiments in output file
" unless @ARGV > 7;

my $outputname = shift(@ARGV);
open(OUT, ">$outputname") or die "Error: cannot open $outputname OUT outfile";
my $CPGinputname = shift(@ARGV);
open(CPG, "<$CPGinputname") or die "Error: cannot open $CPGinputname GTF infile";

#my $windowsize = 10000;
#my $mincpg = 20;
#my $mincoverage = 3;
#my $minfiles = 10;
my $windowsize = shift(@ARGV);
my $mincpg = shift(@ARGV);
my $mincoverage = shift(@ARGV);
my $minfiles = shift(@ARGV);


my @Permethfiles; 
my @Permethnames; 

while(@ARGV){
	push(@Permethfiles, shift(@ARGV));
	push(@Permethnames, shift(@ARGV));
}
my $commandline = "";




###############
# CPG Islands #
###############

my %CPG;
while(<CPG>){
	my @line = split("\t",$_);
	if ($line[0] =~ /chr/){} else {next;}
	if ($line[0] =~ /_/){next;}	
	$CPG{$line[0]}{$line[1]} = $line[2];
}
close CPG;
print "CPGs loaded, analyzing chromosomes:\n";
foreach my $key (sort keys %CPG) {
     print "$key" , "\n";
}
print "\n";

#############
# Main Loop #
#############

my @samplenames;
print OUT "chr" , "\t" , "start" , "\t" , "end";
for (my $n = 0; $n < @Permethnames; $n++){
	print OUT "\t" , $Permethnames[$n];
}
print OUT "\n";

foreach my $key (sort keys %CPG) {
#my $key = "chr15";
	my %MethCpG;
	my @Permeth = @Permethfiles;
	my @Names = @Permethnames;
	while (@Permeth){
		my $inprefix = shift(@Permeth);
		my $sampname = shift(@Names);
		my $infile = $inprefix . $key . ".bed";
		open(IN, "<$infile") or die "Error: cannot open $infile infile";
		print "Analyzing $sampname $key \n";
		while(<IN>){
			my @line = split("\t",$_);
			my @methinfo = split("-",$line[3]);
			if ($methinfo[1] < $mincoverage) {next;}

			my $start = $line[1];
			my $newkey = floor($start/$windowsize);
			if(defined $MethCpG{$sampname}{$newkey}{"sum"}) {
				$MethCpG{$sampname}{$newkey}{"sum"} = $MethCpG{$sampname}{$newkey}{"sum"} + $methinfo[0];
				$MethCpG{$sampname}{$newkey}{"CPGcount"}++;
				if ($MethCpG{$sampname}{$newkey}{"CPGcount"} == $mincpg){
					if(defined $MethCpG{"count"}{$newkey}) {$MethCpG{"count"}{$newkey}++;}
					else {$MethCpG{"count"}{$newkey} = 1;}
				}
			}
			else{
				$MethCpG{$sampname}{$newkey}{"sum"} = $methinfo[0];
				$MethCpG{$sampname}{$newkey}{"CPGcount"} = 1;
			}
		}
		close IN;		
	}
	
	@Permeth = @Permethfiles;
	@Names = @Permethnames;
	#print to outfile
	print "Printing $key \n\n";
	foreach my $countkey (sort { $a <=> $b } keys( %{$MethCpG{"count"}} ) ){
		if ($MethCpG{"count"}{$countkey} < $minfiles) {next;}
		my $start = $countkey * $windowsize;
		my $end = $start + $windowsize;
		print OUT $key , "\t" , $start , "\t" , $end;
		for (my $n = 0; $n < @Names; $n++){
			if(defined $MethCpG{$Names[$n]}{$countkey}){
				if ( $MethCpG{$Names[$n]}{$countkey}{"CPGcount"} >= $mincpg) {
					my $avgmeth = sprintf ("%.3f", ($MethCpG{$Names[$n]}{$countkey}{"sum"} / $MethCpG{$Names[$n]}{$countkey}{"CPGcount"}));
					print OUT "\t" , $avgmeth;
				}
				else{
					print OUT "\t" , "NA";
				}
			}
			else{
				print OUT "\t" , "NA";
			}
		}
		print OUT "\n";
	}	
}

close OUT;

__END__

#		$commandline = "sed 1d $infile > tmppermeth.bed";
#		`$commandline`;
#		print $commandline , "\n";
#		$commandline = "/bin/bedtools subtract -a tmppermeth.bed -b $CPGinputname > masked.bed";
#		print $commandline , "\n";

