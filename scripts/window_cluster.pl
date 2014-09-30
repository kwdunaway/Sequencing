#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 
#
# Sliding window that returns results
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) In table
    2) Out clustered bed file
    3) Out hyper bed file
    4) Out hypo bed file
    5) trackname
    6) column of p values (ex: 15)
    7) column of hyper/hypo calls (ex:16)
    8) p-value cutoff (ex: .05)
    9) windowsize (ex: 17 or 12)
    10) minimum number of significant calls in window (ex: 5 or 7)
" unless @ARGV == 10;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";
my $output_hyper = shift(@ARGV);
open(HYPER, ">$output_hyper") or die "cannot open $output_hyper outfile";
my $output_hypo = shift(@ARGV);
open(HYPO, ">$output_hypo") or die "cannot open $output_hypo outfile";
my $trackname = shift(@ARGV);
my $pvaluecol = shift(@ARGV);
my $dircol = shift(@ARGV);
my $pval = shift(@ARGV);
my $windowsize = shift(@ARGV);
my $minsig = shift(@ARGV);

###################
# Another section #
###################

my @window;
my @windowstarts;
my @windowends;
my $chr = "";
my $runningtot = 0;
my %PosHash;

#remove header
<IN>;
print OUT "track name=" , $trackname , " description=" , $trackname , " itemRgb=\"On\"\n";
print HYPER "track name=" , $trackname , "_hyper description=" , $trackname , " color=255,0,0\n";
print HYPO "track name=" , $trackname , "_hypo description=" , $trackname , " color=0,0,255\n";

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    if($line[$pvaluecol]=~ /D/){next;}
    
    if($chr ne $line[0]){
    	@window = ();
    	@windowstarts = ();
    	@windowends = ();
    	$chr = $line[0];
    	%PosHash = ();
    	$runningtot = 0;
    	print "\n\n" , $chr , "\n";
    }
    
    if($line[$pvaluecol] > $pval){push(@window,0);}
    elsif($line[$dircol] eq "hyper"){
    	push(@window,1); 
    	$runningtot++;
    	print HYPER $line[0] , "\t" , $line[1] , "\t" , $line[2] , "\t" , $line[$dircol] , "\n";
    }
    elsif($line[$dircol] eq "hypo"){
    	push(@window,-1); 
    	$runningtot--;
    	print HYPO $line[0] , "\t" , $line[1] , "\t" , $line[2] , "\t" , $line[$dircol] , "\n";
    }
    else{ die "$line[$dircol] not hyper or hypo";}
    push(@windowstarts,$line[1]);
    push(@windowends,$line[2]);
    if(int(@window) < $windowsize) {next;}
    
    if($runningtot >= $minsig){
		print $runningtot , "\t";
		for(my $n=0;$n<$windowsize;$n++){
			if($window[$n] > 0) {$PosHash{$windowstarts[$n]} = $windowends[$n]; print "\t" , $windowstarts[$n];}
		}
		print "\n";
    }
    elsif($runningtot <= (-1*$minsig)){
		print $runningtot , "\t";
		for(my $n=0;$n<$windowsize;$n++){
			if($window[$n] < 0) {$PosHash{$windowstarts[$n]} = $windowends[$n]; print "\t" , $windowstarts[$n];}
		}
		print "\n";
    }
    elsif (keys %PosHash){
    	print $runningtot , "\n";
       	my $lowest_key = 0;
    	my $highest_value = 0;
    	my $numkeys = 0;
    	foreach my $start (sort { $PosHash{$a} <=> $PosHash{$b} } keys %PosHash) {
    		if($lowest_key == 0){ $lowest_key = $start;}
    		$highest_value = $PosHash{$start};
    		$numkeys++;
    	}
    	my $direction = "hyper";
    	if ($runningtot < 0) { $direction = "hypo";}

		print OUT $chr , "\t" , $lowest_key , "\t" , $highest_value , "\t" , $direction , "-" , $numkeys, "\t0\t+\t" , $lowest_key , "\t" , $highest_value , "\t";
		if($direction eq "hyper") {print OUT "255,0,0\n";}
		else{print OUT "0,0,255\n";}
		%PosHash = ();
    }
    
    $runningtot = $runningtot - shift(@window);
    shift(@windowstarts);
    shift(@windowends);
}
