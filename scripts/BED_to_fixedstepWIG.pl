#!/usr/bin/perl 
#use strict; use warnings;
#BED_to_fixedstepWIG.pl

##################################################
# Command Line Error Checking and I/O Initiation #
##################################################
die "useage: $0 
  1) In BED file name 
  2) Out WIG file name
  3) Name of WIG Track 
  4) WIG Track color (format: RRR,GGG,BBB)" unless @ARGV == 4;
my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read (must be .bed)
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile"; #opens output file to write to (.wig)
my $name = shift(@ARGV);
my $color = shift(@ARGV);

# Prints the head of the track (necessary for genome browser to read file properly) (customizable through terminal)
print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $name, "\" description=\"", $name, "\" color=", $color, "\n";

###################################################
# Global variables needed and their explainations #
###################################################

my %PosVal;              # hash that contains the peak heights, kept small
my $position = 0;        # position of window
my @line;                # temp array used to retreive the information of each line
my $chrom = "notyet";               # chromosome of current line
my $startread = 0;       # start of read at current line
my $endread = 0;         # end of read at current line (does not include this position)



while (<IN>)
{
####################################################################################
# Grabs the information from each line and assigns it to the appropriate variables #
####################################################################################

    chomp;
    @line = split ("\t", $_);
	# Gets rid of any header lines (lines without at least 3 cells in array)
    if (exists $line[2]){}
    else {next;}
	# gets rid of all lines before 0 point
    if ($startread < 0)  {next;}


    my $chr = $line[0];            # makes $chrom have the chromosome information
    
    if ($chr ne $chrom){
		if ($chrom eq "notyet") { 
		    $position = $startread;	
	    	$chrom = $chr;
		}
		else {
			while ($position < $endread) {
			    if (exists $PosVal{$position}) {
			        print OUT $PosVal{$position}, "\n";
					delete $PosVal{$position};
			        $position = $position + 1;
    			}
			    else  {      # skips gaps by jumping to the next position and starting a new fixed step line
			        $position = $startread;
			        print OUT "fixedStep chrom=", $chrom," start=",$startread," step=1\n";
			    }
			}		
		}
		print OUT $chr , "\t" , $chrom , "\n";
    }
    
    $chrom = $chr;
    $startread = $line[1];        # sets the value for the start of the read to $startread
    $endread = $line[2];          # sets the value for the end of the read to $endread
    
	# Makes the program run faster by skipping the first gap
	
	# adds height to PosVal for the sequence found
	my $addcount = $startread;
	while ($addcount < $endread) {
	    if (exists $PosVal{$addcount}) {$PosVal{$addcount} = $PosVal{$addcount} + 1;}
	    else {$PosVal{$addcount} = 1;}
		$addcount = $addcount +1;
	}


##################################################
# Prints data to a fixed step wiggle file (.wig) #
##################################################

	# Since you should NEVER have a read's start before the current read's start, 
	#  this will print all positions until that point
	while ($position < $startread) {
	    if (exists $PosVal{$position}) {
			print OUT $PosVal{$position}, "\n";
			delete $PosVal{$position};
			$position = $position + 1;
	    }
	    else        # skips gaps by jumping to the next position and starting a new fixed step line
	    {
			$position = $startread;
			print OUT "fixedStep chrom=", $chrom," start=",$startread," step=1\n";
	    }
	}
#    }
}


##############################################################
# Prints data past the window to $outfile (ONE LAST TIME!!!) #
##############################################################

#  This will print the rest of the reads
while ($position < $endread) 
{
    if (exists $PosVal{$position})
    {
        print OUT $PosVal{$position}, "\n";
	delete $PosVal{$position};
        $position = $position + 1;
    }
    else        # skips gaps by jumping to the next position and starting a new fixed step line
    {
        $position = $startread;
        print OUT "fixedStep chrom=", $chrom," start=",$startread," step=1\n";
    }
}

close IN;
close OUT;
