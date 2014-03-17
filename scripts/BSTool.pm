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
    2) Min cutoff of reads for primer sets
    2-?) Input fastq file
" unless @ARGV > 2;

my $stats_filename = shift(@ARGV);	
my $cutoff = shift(@ARGV);
my @Infiles = @ARGV;

open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";
my $seqcount = "seqcount_" . $stats_filename;
open(SEQOUT, ">$seqcount") or die "cannot open $seqcount outfile";
#my $out_suffix = shift(@ARGV);

#my $output_filename = shift(@ARGV);
#open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";

#################
# In Files Loop #
#################

while(@Infiles){
	my $infile = shift(@Infiles);
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "\nStarting file: $infile\n";
#	print "Number or Reads matched so far:\t";

#	my $output_filename = $infile . $out_suffix;
#	open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";

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
	my %FlankBefore;
#	my %FlankAfter;
#	my %FlankBoth;

#############################
# Main Search Function Loop #
#############################

	while (<IN>)
	{
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
	    if($seq =~ /(TT([CT])GTGGTG([CT])GT([CT])GTTTTTTAA([GT])T([CT])GGTT)/){
			if (exists $Reads{$seq}){
				$Reads{$seq}{"copy_number"}++;
			}
			else{
				# Initialize hash variables
				$Reads{$seq}{"copy_number"} = 1;
				$Reads{$seq}{"methylation"} = 0;
				$Reads{$seq}{"site1"}=0;
				$Reads{$seq}{"site2"}=0;
				$Reads{$seq}{"site3"}=0;
				$Reads{$seq}{"site4"}=0;
				
				# Get site methylation information data
				#Variable that counts and bins read methylation, 
				my $read_meth = 0;
				if($2 eq "C") {$meth1++; $read_meth++; $Reads{$seq}{"site1"}=1;} elsif($2 eq "T") {$unmeth1++;} else{die "Site1: $2";}
				if($3 eq "C") {$meth2++; $read_meth++; $Reads{$seq}{"site2"}=1;} elsif($3 eq "T") {$unmeth2++;} else{die "Site2: $3";}
				if($4 eq "C") {$meth3++; $read_meth++; $Reads{$seq}{"site3"}=1;} elsif($4 eq "T") {$unmeth3++;} else{die "Site3: $4";}
				if($6 eq "C") {$meth4++; $read_meth++; $Reads{$seq}{"site4"}=1;} elsif($6 eq "T") {$unmeth4++;} else{die "Site4: $6";}
				if($5 eq "G") {$SNPG++;} elsif($5 eq "T") {$SNPT++;} else{die "A/T: $5";}
				$Reads{$seq}{"methylation"} = $read_meth;

				# Get flanking region information
				my @Flank = split(/TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA[GT]T[CT]GGTT/, $seq, 2);
#				print $seq ,"\n"; 				
#				print $Flank[0], "                             " , $Flank[1], "\n";
				my $FlankUp = substr($Flank[0], -30);
#				my $FlankDown = substr($Flank[1], 0, 20);
#				print "Up:\t$FlankUp\nDown:\t$FlankDown\n\n";
				if(length($FlankUp) == 30){
					$FlankUp = substr($FlankUp, 0, 20);
					if (exists $FlankBefore{$read_meth}{$FlankUp}){$FlankBefore{$read_meth}{$FlankUp}++;} else {$FlankBefore{$read_meth}{$FlankUp}=1;}
				}
#				if(length($FlankDown) == 20){
#					if (exists $FlankAfter{$read_meth}{$FlankDown}){$FlankAfter{$read_meth}{$FlankDown}++;} else {$FlankAfter{$read_meth}{$FlankDown}=1;}
#				}
#				if (length($FlankUp) + length($FlankDown) == 40){
#					my $flankb = $FlankUp . "  " . $FlankDown;
#					if (exists $FlankBoth{$read_meth}{$flankb}) {$FlankBoth{$read_meth}{$flankb}++;} else {$FlankBoth{$read_meth}{$flankb}=1;}				
#				}
				
				$LINE1++;
#	    		if($LINE1 % 20 == 0) {print $LINE1 , "\n";}
			}
#			print OUT $ID , $seq , "\n", $third , $quality;
    	}
		else { die "$seq does not match up";}
	}
	
	
	for(my $read_meth = 0; $read_meth < 5; $read_meth++){
#		last;
#		next;
#		delete $FlankBefore{$read_meth}{"GGATTTTTTGAGTTAGGTGT"};

		# Upstream flanking region
		my $first = 0;	
		foreach my $key (sort { $FlankBefore{$read_meth}{$b} <=> $FlankBefore{$read_meth}{$a} } keys %{ $FlankBefore{$read_meth} }) {
			if($first eq "0") { 
				$first = $key; 
				print "Methylation:\t" , $read_meth, "\n";
				print $key , "\t" , $FlankBefore{$read_meth}{$key}, "\n";
			}
			elsif($FlankBefore{$read_meth}{$key} >= $cutoff){
				my $s1 = $first;
				my $s2 = $key;	
				my @s1 = split //,$s1;
				my @s2 = split //,$s2;
				my $i = 0;
				foreach  (@s1) {
				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
				    else {print " ";}
			    	$i++;
				}
	       		print "\t" , $FlankBefore{$read_meth}{$key}, "\n";
		    }
	    }
    
    	# Downstream flanking region
#		$first = 0;
#		foreach my $key (sort { $FlankAfter{$read_meth}{$b} <=> $FlankAfter{$read_meth}{$a} } keys $FlankAfter{$read_meth}) {
#			if($first eq "0") { 
#				$first = $key; 
#				print $key , "\t" , $FlankAfter{$read_meth}{$key}, "\n";
#			}
#			elsif($FlankAfter{$read_meth}{$key} >= $cutoff){
#				my $s1 = $first;
#				my $s2 = $key;
#				my @s1 = split //,$s1;
#				my @s2 = split //,$s2;
#				my $i = 0;
#				foreach  (@s1) {
#				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
#				    else {print " ";}
#			    	$i++;
#				}
#	       		print "\t" , $FlankAfter{$read_meth}{$key}, "\n";
#		    }
#	    }

    	# Both Upstream and Downstream flanking region
#		$first = 0;
#		foreach my $key (sort { $FlankBoth{$read_meth}{$b} <=> $FlankBoth{$read_meth}{$a} } keys $FlankBoth{$read_meth}) {
#			if($first eq "0") { 
#				$first = $key; 
#				print $key , "\t" , $FlankBoth{$read_meth}{$key}, "\n";
#			}
#			elsif($FlankBoth{$read_meth}{$key} >= $cutoff){
#				my $s1 = $first;
#				my $s2 = $key;
#				my @s1 = split //,$s1;
#				my @s2 = split //,$s2;
#				my $i = 0;
#				foreach  (@s1) {
#				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
#				    else {print " ";}
#			    	$i++;
#				}
#	       		print "\t" , $FlankBoth{$read_meth}{$key}, "\n";
#		    }
#	    }

#		print "\n\n";
	}    
	
	my $printseq = "GCGTATTACGAGATTATATT";
#	$printseq = "GGATTTTTTGAGTTAGGTGT";
	print "Stats for: $printseq \n";
	
	my $runningmeth = 0;
	my $runningcount = 0;
	for(my $read_meth = 0; $read_meth < 5; $read_meth++){
		my $methperc = $read_meth / 4;
		print "methylation:\t" , $methperc , "\t" , $FlankBefore{$read_meth}{$printseq}, "\n";
		$runningmeth = $runningmeth + $methperc * $FlankBefore{$read_meth}{$printseq};
		$runningcount = $runningcount + $FlankBefore{$read_meth}{$printseq};
	}
	print "Total:\t", $runningmeth / $runningcount , "\t" , $runningcount , "\n";


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

	print SEQOUT "$infile\n", "Sequence\tCopy Number\tMethylation\tSite 1\tSite 2\tSite 3\tSite 4\n";
	foreach my $key (sort { $Reads{$b}{"methylation"} <=> $Reads{$a}{"methylation"} } keys %Reads) {
        print SEQOUT $key , "\t" , $Reads{$key}{"copy_number"}, "\t" , $Reads{$key}{"methylation"}, "\t", $Reads{$key}{"site1"}, "\t", $Reads{$key}{"site2"}, "\t", $Reads{$key}{"site3"}, "\t", $Reads{$key}{"site4"}, "\n";
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

