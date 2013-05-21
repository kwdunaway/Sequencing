#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 7-31-2012
# Script Name: bisulfite_stats.pl
#
#
#
# Arguments:
#    See Below
#
################################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input file
" unless @ARGV == 1;

my $inputfile = shift(@ARGV);

####################
#    Processing    #
####################

bisulfite_stats($inputfile);

sub bisulfite_stats
{
	# Input
	my ($inputfile) = @_;

	# File I/O
	open(IN, "<$inputfile") or die "Error: Bisulfite Stats: cannot open $inputfile infile";
	(my $filename = $inputfile) =~ s{\.[^.]+$}{}; # Find file name
	my ($extension) = $inputfile =~ /(\.[^.]+)$/; # Find file extension
	my $outfile = $filename . "_stats" . $extension; # Output file is Filename_stats.ext
	open(OUT, ">$outfile") or die "Error: Bisulfite Stats: cannot open $outfile outfile";

	# Hash{Number of Sites}->[Number of Reads, Number of Sites]

	my %CGhash; 
	my %CHGhash;
	my %CHHhash;

	while(<IN>){
		chomp;
		my @line = split("\t",$_);
		my $sites_line = $line[6];
		my $CGsites = 0;
		my $num_methylated = 1;

		$CGsites++ while ($sites_line =~ m/[xX]/g);
		$num_methylated++ while ($sites_line =~ m/[X]/g);
		$CGhash{$CGsites}[0]++;
		$CGhash{$CGsites}[$num_methylated]++;
	}

	foreach my $sites (keys %CGhash) {
		print $sites, " ", $CGhash{$sites}[0], " ";
		for (my $i =0; $i < $sites; $i++) {
			print $CGhash{$sites}[$i], ",";
		}
		print "\n";
	}

	close IN;
	close OUT;
}

__END__
		
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
