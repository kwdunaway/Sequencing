#!/usr/bin/perl
use warnings; use strict;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Script Name: Combine_HMM_counts.pl
# Version: 1.1
# Last Updated: 5/5/2014
#
# This script takes multiple models and combines their counts into one model. 
# This is useful for training a model on multiple samples.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Output Model File Name
    2) First Model File (all information except counts are copied from this file)
    3-?) Other Model Files
" unless @ARGV > 2;

my $outmodel_file = shift(@ARGV);
open my $outmodel, ">$outmodel_file" or die "Can't write to $outmodel_file\n";
my @model_infiles = @ARGV;


##############################################
# Read First Model for format and basic info #
##############################################

my @entries;		# Holds all numbers, $entries[row][col]
while (@model_infiles){
	my $model_infile = shift(@model_infiles);
	open(IN, "<$model_infile") or die "cannot open $model_infile infile";
	my $row = 0;

	if (!@model_infiles) { # If this is the last file, use it as template
		while(<IN>){
			chomp;
			if(substr($_, 0, 1) =~ /^\d/){ #Number line
				my $i = 0; # Column
				my @line = split ("\t", $_);
				foreach my $entry (@line) {
					$entries[$row][$i] += $entry;
					print $outmodel "$entries[$row][$i]\t";	
					$i++;
				}
				print $outmodel "\n";
				$row++;
			}
			else{ #Text line
				print $outmodel "$_\n";
			}
		}
	}
	else { # This is not the last file, only do addition
		while(<IN>){
			chomp;
			if(substr($_, 0, 1) =~ /^\d/){ #Number line
				my $i = 0;
				my @line = split ("\t", $_);
				foreach my $entry (@line) {
					$entries[$row][$i] += $entry - 1;	
					$i++;
				}
				$row++;
			}
			else{ #Text line
			}
		}
	}
	close IN;
}

close $outmodel;
