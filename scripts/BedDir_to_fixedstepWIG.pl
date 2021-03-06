#!/usr/bin/perl 
#use strict; use warnings;
#BedDir_to_fixedstepWIG.pl



##################################################
# Command Line Error Checking and I/O Initiation #
##################################################
die "useage: BED_to_fixedstepWIG.pl 
<Infile name> 
<WIG Track color (format: RRR,GGG,BBB)> 
<Track Name>
[Note: the output filename will be the same as input except .wig]
****** BE IN SAME DIR AS INFILE TO GET CORRECT WIG FILE NAMES ******" unless @ARGV == 3;
my $infileroot = shift(@ARGV);
my $color = shift(@ARGV);
my $trackname = shift(@ARGV);


###################################################
# Global variables needed and their explainations #
###################################################

my @Chr;                                               # array that contains all the the names of the chromosomes
push(@Chr, "1");
for (my $n = 10; $n< 20; $n++)
{
    push(@Chr, $n);
}
for (my $n = 2; $n< 10; $n++)
{
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");



####################################################################################
# Grabs the information from each line and assigns it to the appropriate variables #
####################################################################################

while(@Chr)
{
    my $infile = $infileroot . $Chr[0] . ".bed";
    open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read (must be .bed)
    my @infileroot = split(".bed",$infile);
    my $outfile = $infileroot . $Chr[0] . ".wig";
    open(OUT, ">$outfile") or die "cannot open $outfile outfile"; #opens output file to write to (.wig)
    # Prints the head of the track (necessary for genome browser to read file properly) (customizable through terminal)
    my $name = $trackname . "_chr_" . $Chr[0];
    print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $name, "\" description=\"", $name, "\" color=", $color, "\n";

    my %PosVal;              # hash that contains the peak heights, kept small
    my $position = 0;        # position of window
    my $chromcheck;          # variable that checks each chrom to make sure they are the same chrom
    my @line;                # temp array used to retreive the information of each line
    my $chrom;               # chromosome of current line
    my $startread = 0;       # start of read at current line
    my $endread = 0;         # end of read at current line (does not include this position)

    while (<IN>)
    {
	chomp;
	@line = split ("\t", $_);
	$chrom = $line[0];            # makes $chrom have the chromosome information
	$startread = $line[1];        # sets the value for the start of the read to $startread
	$endread = $line[2];          # sets the value for the end of the read to $endread

	if ($startread > -1)  # gets rid of all lines before 0 point
	{
#    print "$chrom \t $startread \t $endread \n";
	    # Makes the program run faster by skipping the first gap
	    if ($position == 0)
	    { 
		$position = $startread;	
		$chromcheck = $chrom;
		print OUT "fixedStep chrom=", $chrom," start=",$startread," step=1\n";
	    }
	    
	    # ends program if you have different chromosomes in your data
	    die "You have different chromosomes in this bed file (Program ended before completion)/n" unless $chrom == $chromcheck;
	    
	    # adds height to PosVal for the sequence found
	    my $addcount = $startread;
	    while ($addcount < $endread)
	    {
		if (exists $PosVal{$addcount})
		{
		    $PosVal{$addcount} = $PosVal{$addcount} + 1;
#	        print "$PosVal{$addcount} \n";
		}
		else
		{
		    $PosVal{$addcount} = 1;
		}
#           print "$addcount \t",$PosVal{$addcount}, "\n";
		$addcount = $addcount +1;
	    }
	    
	    
	    
##################################################
# Prints data to a fixed step wiggle file (.wig) #
##################################################

	# Since you should NEVER have a read's start before the current read's start, 
	#  this will print all positions until that point
	    while ($position < $startread) 
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
	}
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
    print "Finished with Chromosome $Chr[0] \n";
    shift(@Chr);
    my $commandline = "gzip " . $outfile;
    `$commandline`;
}
