#!/usr/bin/perl -w

# Diane Schroeder, 9/10/12 version

# This program takes the output from BED_percent_methyl_window_v3.pl or this
# program and adds an extra column for a new sample.  

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
#  5) tissue name (.methyl will be added to end, for column info)
#  6) input file name

use POSIX;

unless ($ARGV[5]) {

    die "Arguments:\n\tDNA methylation BED input files, # as chromosome wildcard\n\toutput file\n\twindow size\n\tminimun number of (covered) CpG sites in window\n\ttissue name\n\tinput file name\n\n";
}

$file_base = $ARGV[0];

open (OUTPUT, ">$ARGV[1]") or die "Couldn't create $ARGV[1]\n";

unless ($window_size = $ARGV[2]) {
    $window_size = 20000;
}

unless ($min_CpGs = $ARGV[3]) {
    $min_CpGs = 500;
}

open (INPUT, "<$ARGV[5]") or die "Couldn't open $ARGV[5]\n";

$header = <INPUT>;
chomp $header;

$tissue_name = $ARGV[4] . ".methyl";

print OUTPUT "$header\t$tissue_name\n";


###########################################################################

$this_chrom = "NA";

while ($input_line = <INPUT>) {

    chomp $input_line;

    @fields = split /\t/, $input_line;

    $coord = shift @fields; 

    ($next_chrom, $window_start) = split /:/, $coord;

    $window_end = $window_start + $window_size - 1;

    $methyls = join "\t", @fields; 

    if ($this_chrom ne $next_chrom) {  #get new BED file

	$this_chrom = $next_chrom;
	$this_chrom =~ m/chr(\w+)/;
	$chrom_num = $1;

	$filename = $file_base;
	$filename =~ s/#/$chrom_num/;

	open (BED, "<$filename") or die "Couldn't open $filename\n";

	$line = <BED>;  #get rid of header line

	$line = <BED>;  #get first data line
	@fields = split /\t/, $line;
	$this_position = $fields[1];

	while ($this_position < $window_start) { #get to first window

	    $line = <BED>;  #get first data line
	    @fields = split /\t/, $line;
	    $this_position = $fields[1];
	}
    }

    #now get methylation

    $num_CpGs = 0;
    $methyl_sum = 0;
    
    $finished_file = 0;
    
    while (!$finished_file) { #go through window
		    
	if ($line = <BED>) { #get next line
	    
	    @fields = split /\t/, $line;
	    
	    $this_methyl = $fields[3];
	    
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
	
	if ($num_CpGs >= $min_CpGs) {
	    
	    $ave_methyl = $methyl_sum / $num_CpGs;
	    
	    print OUTPUT "$coord\t$methyls\t$ave_methyl\n";
	}
	else {
	    
	    print OUTPUT "$coord\t$methyls\tNA\n";
	}
	
	$num_CpGs = 1;
	
	$methyl_sum = $this_methyl;
    }
    
} #end while (going through input file)

