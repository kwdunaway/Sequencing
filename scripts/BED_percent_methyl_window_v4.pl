#!/usr/bin/perl -w

# Diane Schroeder, 8/22/13 version

# This version will allow multiple BED files (for the same sample).  It will
# also be multi-chromosome friendly (up to 39 chromosomes).  
# Importantly, it will allow adding additional
# data to the output file later.  To do that, windows with not enough data
# will be given "NA" instead of being omitted.

#####################################################################

# This program takes a UCSC Genome Browser BED file with percent methylation
# for CpG sites and runs non-overlapping windows across it, collecting the
# average percent methylation. 

# Output file will have: 
#  1) the start coordinate of each window 
#  2) the average percent methylation over the window.

# Arguments:
#  1) BED input files, using # as a wildcard for chromosome number
#  2) output file
#  3) window size (will be non-overlapping, so make it small enough)
#  4) minimum number of covered CpG sites in the window 
#  5) use X and Y chromosomes?  (Y|N)  (default is N)
#  6) tissue name (.methyl will be added to end, for column info)

use POSIX;

unless ($ARGV[5]) {

    die "Arguments:\n\tDNA methylation BED input files, # as chromosome wildcard\n\toutput file\n\twindow size\n\tminimun number of (covered) CpG sites in window\n\tuse X and Y chromosomes?  (Y|N)\n\ttissue name\n\n";
}

$singlefile = 0;

if ($ARGV[0] !~ m/#/) {  #no wildcard, just use file

    $ARGV[0] =~ m/chr(\w{1,2})/;

    push @chromlist, $1;

    $singlefile = 1;
}

elsif ($ARGV[4] eq "Y") { #wildcard, with X and Y

    @chromlist = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "X", "Y");
}
else {  #wildcard, no X or Y

  @chromlist = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39");
}  

open (OUTPUT, ">$ARGV[1]") or die "Couldn't create $ARGV[1]\n";

unless ($window_size = $ARGV[2]) {
    $window_size = 20000;
}

unless ($min_CpGs = $ARGV[3]) {
    $min_CpGs = 500;
}

$tissue_name = $ARGV[5] . ".methyl";

print OUTPUT "start\t$tissue_name\n";

###########################################################################

foreach $chromnum (@chromlist) {

    $filename = $ARGV[0];

    unless ($singlefile) {

	$filename =~ s/#/$chromnum/;
    }

    open (INPUT, "<$filename") or next;

    $line = <INPUT>;  #get rid of header line
    
    $line = <INPUT>;  #first line of data.  Start windowing from here
    
    @fields = split /\t/, $line;
    
    $window_start = $fields[1];   #position in chromosome
    
    $window_end = $window_start + $window_size - 1;
    
    $num_CpGs = 1;
    
	@tmethfields = split /-/, $fields[3];

    $methyl_sum = $tmethfields[0];  #methylation at this site
    
    $finished_file = 0;
    
    while (!$finished_file) { #go through file
	
	while () { #go through window
	    
	    if ($line = <INPUT>) { #get next line
		
		@fields = split /\t/, $line;
		@methfields = split /-/, $fields[3];

		$this_methyl = $methfields[0];
		
		$this_position = $fields[1];
		
		if ($this_position <= $window_end) {
		    
		    $num_CpGs++;
		    
		    $methyl_sum += $this_methyl;
		}
		else {  #end of window
		    last;
		}
	    }
	    
	    else { #no more lines in file
		
		$finished_file = 1;
		
		#we'll just leave this partial window out of the output
		
		last;
	    }
	    
	} #end while (done going through window)
	
	unless($finished_file) {  #print to file and start new window

	    $output_start = "chr$chromnum:$window_start";
	    
	    if ($num_CpGs >= $min_CpGs) {
		
		$ave_methyl = $methyl_sum / $num_CpGs;
		
		print OUTPUT "$output_start\t$ave_methyl\n";
	    }
	    else {

		print OUTPUT "$output_start\tNA\n";
	    }
	    
	    $num_CpGs = 1;
	    
	    $methyl_sum = $this_methyl;
	    
	    $window_start = $window_end + 1;
	    
	    $window_end = $window_start + $window_size - 1;
	}
	
    } #end while (done going through file)

} #end foreach (going through all files)

