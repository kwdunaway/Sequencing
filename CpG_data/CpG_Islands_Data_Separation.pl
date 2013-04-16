#!/usr/bin/perl
BEGIN {push @INC, "/home/kwdunaway/perl_script";}
use strict; use warnings;

###############################################################################################
# Author: Roy Chu
# Email: rgchu@ucdavis.edu
# Date: 3/5/2013
# Script Name: CpG_Islands_Data_Separation.pl
#
# This script separates CpG Island data into seven files:
#	A) Data with Sole CpG Islands (Not part of archipelagos)(Shores are 2000 bp)
#		1) CpG Islands only (no shores)
#		2) CpG Islands with shores 
#		3) CpG Shores of Islands only
#	B) Data with CpG Archipelagos (CpG islands close enough to have overlapping shores)
#		4) CpG Islands of Archipelagos only
#		5) CpG Archipelagos with Outer Shores (Shores bordering start and end of archipelago)
#		6) CpG Archipelagos without Outer Shores (Start and end at first and last CpG island)
#		7) CpG Shores of Archipelagos only
#
# Arguments: See Below
#
################################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "CpG_Islands_Data_Separation.pl needs the following parameters:
    1) Input CpG GTF file
    2) Folder path to contain data (e.g. /home/user/Folder/)(No spaces in path)(Will make folder if nonexistent)
    3) File prefix for created CpG data files (name for new files)
" unless @ARGV == 3;

my $InputGTFfile = shift(@ARGV);
my $ExperimentTopDir = shift(@ARGV);
my $FilePrefix = shift(@ARGV);

# Adjusting Folder Names
$ExperimentTopDir = $ExperimentTopDir . "\/" unless $ExperimentTopDir =~ m/\/$/; # Ensures ending "/"

# Making folder to contain output
print "CpG Islands Data Separation: Checking/making $ExperimentTopDir directory in current directory\n";

if (! -d $ExperimentTopDir) 
{ 
	`mkdir $ExperimentTopDir`;
}

# File Name Creation
my $islandsfilename = $ExperimentTopDir . $FilePrefix . "_Sole_CpG_Islands.gtf";
my $islandsandshoresfilename = $ExperimentTopDir . $FilePrefix . "_Sole_CpG_Islands_And_Shores.gtf";
my $shoresfilename = $ExperimentTopDir . $FilePrefix . "_Sole_CpG_Shores.gtf";
my $archipelagoislandsfilename = $ExperimentTopDir . $FilePrefix . "_CpG_Archipelagos_Islands.gtf";
my $archipelagowithshoresfilename = $ExperimentTopDir . $FilePrefix . "_CpG_Archipelagos_With_Outer_Shores.gtf";
my $archipelagowithoutshoresfilename = $ExperimentTopDir . $FilePrefix . "_CpG_Archipelagos_Without_Outer_Shores.gtf";
my $archipelagoshoresfilename = $ExperimentTopDir . $FilePrefix . "_CpG_Archipelagos_Shores.gtf";

####################
# Global Variables #
####################

my %GTFHash; # Hash to store all the variables in every line

# cpgcount - 1 is total reads in specified file
#		0) CpG Islands only (no shores)
#		1) CpG Islands with shores 
#		2) CpG Shores of Islands only
#		3) CpG Islands of Archipelagos only
#		4) CpG Archipelagos with Outer Shores (Shores bordering start and end of archipelago)
#		5) CpG Archipelagos without Outer Shores (Start and end at first and last CpG island)
#		6) CpG Shores of Archipelagos only
my @cpgcount = (1, 1, 1, 1, 1, 1, 1);

my $archstart = 0;
my $archend = 0;
my $archlinestart = 0;
my $archlineend = -1;
my $linenum = 1;

# Opening Files        
open(IN, "<$InputGTFfile") or die "Error: CpG Islands Data Separation: cannot open $InputGTFfile GTF infile";
open(ISLANDS, ">$islandsfilename") or die "Error: CpG Islands Data Separation: cannot open $islandsfilename outfile";
open(ISLANDSANDSHORES, ">$islandsandshoresfilename") or die "Error: CpG Islands Data Separation: cannot open $islandsandshoresfilename outfile";
open(SHORES, ">$shoresfilename") or die "Error: CpG Islands Data Separation: cannot open $shoresfilename outfile";
open(ARCHIPELAGOISLANDS, ">$archipelagoislandsfilename") or die "Error: CpG Islands Data Separation: cannot open $archipelagoislandsfilename outfile";
open(ARCHIPELAGOWITHSHORES, ">$archipelagowithshoresfilename") or die "Error: CpG Islands Data Separation: cannot open $archipelagowithshoresfilename outfile";
open(ARCHIPELAGOWITHOUTSHORES, ">$archipelagowithoutshoresfilename") or die "Error: CpG Islands Data Separation: cannot open $archipelagowithoutshoresfilename outfile";
open(ARCHIPELAGOSHORES, ">$archipelagoshoresfilename") or die "Error: CpG Islands Data Separation: cannot open $archipelagoshoresfilename outfile";

####################################################################
#          Processing GTF file into multiple output files          #
####################################################################

print "CpG Islands Data Separation: Processing $InputGTFfile and splitting to output files\n";
# Process data for first line
$_ = <IN>;
chomp($_);
my @line = split("\t",$_);
for(my $n = 0; $n < 7; $n++)
{
	push(@{$GTFHash{$linenum}},$line[$n]);
}
$linenum++;

while(<IN>){
	chomp;
	my @line = split("\t",$_);
	for(my $n = 0; $n < 7; $n++)
	{
		push(@{$GTFHash{$linenum}},$line[$n]);
		#print($GTFHash{$linenum}[$n], "\t");
	}

	if($GTFHash{$linenum}[3] < $GTFHash{$linenum-1}[4] + 4000) # If previous island within 4000 bp, then part of archipelago
	{
		if($GTFHash{$linenum-1}[0] eq $GTFHash{$linenum}[0]) # Make sure same chromosome
		{
			if($archstart == 0) # Checks for start of new archipelago or addition to existing
			{
				$archstart = $GTFHash{$linenum-1}[3]; # Only initializes start if new
				$archlinestart = $linenum-1;
			}
			$archend = $GTFHash{$linenum}[4]; # End of archipelago is current end
			$archlineend = $linenum;
		}
		else # Not in archipelago or end of current archipelago
		{
			if($archlinestart <= $archlineend) # If end of an archipelago, print to archipelago
			{
				# Entire Archipelago - starting and ending at the first and last island
				if ($cpgcount[5] > 1) {print ARCHIPELAGOWITHOUTSHORES "\n";}
				++$cpgcount[5];
				print ARCHIPELAGOWITHOUTSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Without_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[5]-1 , "\t" , $archstart , "\t" , $archend , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

				# Entire Archipelago - starting and ending at the outer shores
				if ($cpgcount[4] > 1) {print ARCHIPELAGOWITHSHORES "\n";}
				++$cpgcount[4];
				print ARCHIPELAGOWITHSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_With_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1 , "\t" , $archstart - 2000 , "\t" , $archend + 2000, "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

				# Archipelago Shores Only - CpG number will be that of the start of the archipelago's
				# Outer Shore at the start
				if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
				++$cpgcount[6];
				print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_1" , "\t" , $archstart-2000 , "\t" , $archstart , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

				# Inner Shores
				for(my $k = 0; $k < $archlineend - $archlinestart; $k++)
				{
					if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
					++$cpgcount[6];
					print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $k+2 , "\t" , $GTFHash{$archlinestart+$k}[4] , "\t" , $GTFHash{$archlinestart+$k+1}[3] , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];
				}

				# Outer Shore at the end
				if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
				++$cpgcount[6];
				print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $archlineend-$archlinestart+2 , "\t" , $archend , "\t" , $archend+2000 , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

				# Archipelago Islands Only
				for(my $j = 0; $j <= $archlineend - $archlinestart; $j++)
				{
					if ($cpgcount[3] > 1) {print ARCHIPELAGOISLANDS "\n";}
					++$cpgcount[3];
					print ARCHIPELAGOISLANDS $GTFHash{$archlinestart+$j}[0] , "\t" , "hg18_CpG_Archipelagos_Islands" , "\t" , $GTFHash{$archlinestart+$j}[2] , "\t" , $GTFHash{$archlinestart+$j}[3] , "\t" , $GTFHash{$archlinestart+$j}[4] , "\t" , $GTFHash{$archlinestart+$j}[5] , "\t" , $GTFHash{$archlinestart+$j}[6];
				}
	
				# Since the archipelago data has been printed, the variables can be reset
				$archlinestart = 0; # Reset start and end lines for archipelago
				$archlineend = -1;
				$archstart = 0; # Reset start coordinates of archipelago
			}
			else # Read is an island
			{
				if ($cpgcount[0] > 1) {print ISLANDS "\n";}
				++$cpgcount[0];
				print ISLANDS $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Islands" , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

				if ($cpgcount[1] > 1) {print ISLANDSANDSHORES "\n";}
				++$cpgcount[1];
				print ISLANDSANDSHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Islands_And_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

				# Print Left Shore(1) and Right Shore(2)
				if ($cpgcount[2] > 1) {print SHORES "\n";}
				++$cpgcount[2];
				print SHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "_1" , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

				if ($cpgcount[2] > 1) {print SHORES "\n";}
				++$cpgcount[2];
				print SHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "_2" , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];
			}
		}
	}
	else # Not in archipelago or end of current archipelago
	{
		if($archlinestart <= $archlineend) # If end of an archipelago, print to archipelago
		{
			# Entire Archipelago - starting and ending at the first and last island
			if ($cpgcount[5] > 1) {print ARCHIPELAGOWITHOUTSHORES "\n";}
			++$cpgcount[5];
			print ARCHIPELAGOWITHOUTSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Without_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[5]-1 , "\t" , $archstart , "\t" , $archend , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

			# Entire Archipelago - starting and ending at the outer shores
			if ($cpgcount[4] > 1) {print ARCHIPELAGOWITHSHORES "\n";}
			++$cpgcount[4];
			print ARCHIPELAGOWITHSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_With_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1 , "\t" , $archstart - 2000 , "\t" , $archend + 2000, "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

			# Archipelago Shores Only - CpG number will be that of the start of the archipelago's
			# Outer Shore at the start
			if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
			++$cpgcount[6];
			print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_1" , "\t" , $archstart-2000 , "\t" , $archstart , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

			# Inner Shores
			for(my $k = 0; $k < $archlineend - $archlinestart; $k++)
			{
				if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
				++$cpgcount[6];
				print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $k+2 , "\t" , $GTFHash{$archlinestart+$k}[4] , "\t" , $GTFHash{$archlinestart+$k+1}[3] , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];
			}

			# Outer Shore at the end
			if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
			++$cpgcount[6];
			print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $archlineend-$archlinestart+2 , "\t" , $archend , "\t" , $archend+2000 , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

			# Archipelago Islands Only
			for(my $j = 0; $j <= $archlineend - $archlinestart; $j++)
			{
				if ($cpgcount[3] > 1) {print ARCHIPELAGOISLANDS "\n";}
				++$cpgcount[3];
				print ARCHIPELAGOISLANDS $GTFHash{$archlinestart+$j}[0] , "\t" , "hg18_CpG_Archipelagos_Islands" , "\t" , $GTFHash{$archlinestart+$j}[2] , "\t" , $GTFHash{$archlinestart+$j}[3] , "\t" , $GTFHash{$archlinestart+$j}[4] , "\t" , $GTFHash{$archlinestart+$j}[5] , "\t" , $GTFHash{$archlinestart+$j}[6];
			}
	
			# Since the archipelago data has been printed, the variables can be reset
			$archlinestart = 0; # Reset start and end lines for archipelago
			$archlineend = -1;
			$archstart = 0; # Reset start coordinates of archipelago
		}
		else # Read is an island
		{
			if ($cpgcount[0] > 1) {print ISLANDS "\n";}
			++$cpgcount[0];
			print ISLANDS $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Islands" , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

			if ($cpgcount[1] > 1) {print ISLANDSANDSHORES "\n";}
			++$cpgcount[1];
			print ISLANDSANDSHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Islands_And_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

			# Print Left Shore(1) and Right Shore(2)
			if ($cpgcount[2] > 1) {print SHORES "\n";}
			++$cpgcount[2];
			print SHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "_1" , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

			if ($cpgcount[2] > 1) {print SHORES "\n";}
			++$cpgcount[2];
			print SHORES $GTFHash{$linenum-1}[0] , "\t" , "hg18_Sole_CpG_Shores" , "\t" , $GTFHash{$linenum-1}[2] , "_2" , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];
		}
	}
	$linenum++;
}

# Print information for last line of input GTF file
if($archlinestart <= $archlineend) # If end of an archipelago, print to archipelago
{
	# Entire Archipelago - starting and ending at the first and last island
	if ($cpgcount[5] > 1) {print ARCHIPELAGOWITHOUTSHORES "\n";}
	++$cpgcount[5];
	print ARCHIPELAGOWITHOUTSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Without_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[5]-1 , "\t" , $archstart , "\t" , $archend , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

	# Entire Archipelago - starting and ending at the outer shores
	if ($cpgcount[4] > 1) {print ARCHIPELAGOWITHSHORES "\n";}
	++$cpgcount[4];
	print ARCHIPELAGOWITHSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_With_Outer_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1 , "\t" , $archstart - 2000 , "\t" , $archend + 2000, "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

	# Archipelago Shores Only - CpG number will be that of the start of the archipelago's
	# Outer Shore at the start
	if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
	++$cpgcount[6];
	print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_1" , "\t" , $archstart-2000 , "\t" , $archstart , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

	# Inner Shores
	for(my $k = 0; $k < $archlineend - $archlinestart; $k++)
	{
		if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
		++$cpgcount[6];
		print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $k+2 , "\t" , $GTFHash{$archlinestart+$k}[4] , "\t" , $GTFHash{$archlinestart+$k+1}[3] , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];
	}

	# Outer Shore at the end
	if ($cpgcount[6] > 1) {print ARCHIPELAGOSHORES "\n";}
	++$cpgcount[6];
	print ARCHIPELAGOSHORES $GTFHash{$archlinestart}[0] , "\t" , "hg18_CpG_Archipelagos_Shores" , "\t" , "CpG_Archipelago_" , $cpgcount[4]-1, "_" , $archlineend-$archlinestart+2 , "\t" , $archend , "\t" , $archend+2000 , "\t" , $GTFHash{$archlinestart}[5] , "\t" , $GTFHash{$archlinestart}[6];

	# Archipelago Islands Only
	for(my $j = 0; $j <= $archlineend - $archlinestart; $j++)
	{
		if ($cpgcount[3] > 1) {print ARCHIPELAGOISLANDS "\n";}
		++$cpgcount[3];
		print ARCHIPELAGOISLANDS $GTFHash{$archlinestart+$j}[0] , "\t" , "hg18_CpG_Archipelagos_Islands" , "\t" , $GTFHash{$archlinestart+$j}[2] , "\t" , $GTFHash{$archlinestart+$j}[3] , "\t" , $GTFHash{$archlinestart+$j}[4] , "\t" , $GTFHash{$archlinestart+$j}[5] , "\t" , $GTFHash{$archlinestart+$j}[6];
	}
	
	# Since the archipelago data has been printed, the variables can be reset
	$archlinestart = 0; # Reset start and end lines for archipelago
	$archlineend = -1;
	$archstart = 0; # Reset start coordinates of archipelago
}
else # Read is an island
{
	if ($cpgcount[0] > 1) {print ISLANDS "\n";}
	++$cpgcount[0];
	print ISLANDS $GTFHash{$linenum-1}[0] , "\t" , $GTFHash{$linenum-1}[1] , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

	if ($cpgcount[1] > 1) {print ISLANDSANDSHORES "\n";}
	++$cpgcount[1];
	print ISLANDSANDSHORES $GTFHash{$linenum-1}[0] , "\t" , $GTFHash{$linenum-1}[1] , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

	# Print Left Shore(1) and Right Shore(2)
	if ($cpgcount[2] > 1) {print SHORES "\n";}
	++$cpgcount[2];
	print SHORES $GTFHash{$linenum-1}[0] , "\t" , $GTFHash{$linenum-1}[1] , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[3]-2000 , "\t" , $GTFHash{$linenum-1}[3] , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];

	if ($cpgcount[2] > 1) {print SHORES "\n";}
	++$cpgcount[2];
	print SHORES $GTFHash{$linenum-1}[0] , "\t" , $GTFHash{$linenum-1}[1] , "\t" , $GTFHash{$linenum-1}[2] , "\t" , $GTFHash{$linenum-1}[4] , "\t" , $GTFHash{$linenum-1}[4]+2000 , "\t" , $GTFHash{$linenum-1}[5] , "\t" , $GTFHash{$linenum-1}[6];
}

# Closing Files
close IN;
close ISLANDS;
close ISLANDSANDSHORES;
close SHORES;
close ARCHIPELAGOISLANDS;
close ARCHIPELAGOWITHSHORES;
close ARCHIPELAGOWITHOUTSHORES;
close ARCHIPELAGOSHORES;

print "CpG Islands Data Separation: Finished with output\n";
print "\nStatistics: Total reads in each file\n
  0) CpG Islands only (no shores) - ", $cpgcount[0]-1, "\n",
"  1) CpG Islands with shores - ", $cpgcount[1]-1, "\n",
"  2) CpG Shores of Islands only - ", $cpgcount[2]-1, "\n",
"  3) CpG Islands of Archipelagos only - ", $cpgcount[3]-1, "\n",
"  4) CpG Archipelagos with Outer Shores (Shores bordering start and end of archipelago) - ", $cpgcount[4]-1, "\n",
"  5) CpG Archipelagos without Outer Shores (Start and end at first and last CpG island) - ", $cpgcount[5]-1, "\n",
"  6) CpG Shores of Archipelagos only - ", $cpgcount[6]-1, "\n";
