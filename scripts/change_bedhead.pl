#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Script Name: change_bedhead.pl
# Version: 0.1
# Last Updated: 5-15-2013
#
# This script reads every bed file in the input directory and changes "PercMethylation" 
# in the track name and "PercentMethylation" in the description to a new name.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input folder
    2) NewID
" unless @ARGV == 2;

my $inputfolder = shift(@ARGV);
my $newid = shift(@ARGV);
# Read all .bed files
# Open file
# Rename 1st row
#
# Was:
# track name=PercMethylationchr1 description=PercentMethylationchr1 useScore=0 itemRgb=On db=hg18
#
# Now:
# track name=$NewID_chr1 description=$NewID_chr1 useScore=0 itemRgb=On db=hg18



###################
# Another section #
###################

change_bedhead($inputfolder,$newid);

sub change_bedhead
{
	# Input
	my ($inputfolder, $newid) = @_;

	my $filedir = $inputfolder;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /\.bed/, @files; # Takes only the bed files
	foreach my $infile (@files) {

		# File I/O
		open(IN, "<$infile") or die "Error: change_bedhead.pl: cannot open $inputfolder infile";
		my $outfile = "temp_chr.bed";
		open(OUT, ">$outfile") or die "Error: change_bedhead.pl: cannot open $outfile outfile";
		
		# Extracting Track Name and Description
		my $firstline = <IN>;
		$firstline =~ /name=PercMethylation(\w+)[ ]/;
		my $trackname = $1;
		$firstline =~ /description=PercentMethylation(\w+)[ ]/;
		my $description = $1;

		# The following method of printing is to ensure that there is no space at the end of the line
		my @line = split(" ",$firstline); # Split line by spaces
		print OUT $line[0]; # Print first field
		foreach (@line) { # For every field
        		if($_ eq $line[0]){ # Do not print first field again
			}
			elsif($_ =~ /name=/) { # If field is part of track name, replace
				print OUT " name=", $newid, $trackname;
			}
			elsif($_ =~ /description=/) { # If field is part of description, replace
				print OUT " description=", $newid, $description;
			}
			else { # Otherwise, print original field
				print OUT " ", $_;
			}
		}
		print OUT "\n";
		while (<IN>) { print OUT $_; }
		close IN;
		close OUT;
		`rm $infile`; # Remove original file
		`mv $outfile $infile`; # NewID file becomes the original file

	}
}

__END__
my $fline = chop($firstline);
my @line = split ("PercMethylation", $fline);
my @line2 = split ("PercentMethylation", $line[1]);

