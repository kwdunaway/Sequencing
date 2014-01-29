#!/usr/bin/perl -w

# Diane Schroeder, 4/25/11 version

# This version removes the requirement for the read length to be 
# put in at the command line and allows for the -a arugment in
# BSSeeker, which crops off the end of the
# read if it contains adapter sequence.  This version will adjust
# read lengths for each read individually.

#----------------------------------------------------------------

# As of 3/29/11 version, errors causing a small number of CpG sites to 
# have >100% methylation has been corrected.  This occurred in special
# cases where reads on the reverse strand ended in a methylated C
# (causing methylation to be incremented but not coverage beings it
# didn't span the adjacent G).  However, this program still won't behave
# properly if, for example, a sole read is in the above situation.  In
# that case, the program will pretend the methylated C at the tail of
# the reverse orientation read didn't happen and will be considered
# "not covered".  So you basically loose a small amount of your data at
# the very end of a handful of reads.

#------------------------------------------------------------------

#  This program generates a BED file for the UCSC Genome Browser,
#  showing percent methylation.  For now, this program will only run properly
#  on a file with data from a single chromosome!

#  For now, this program will also only run properly if only one sequence
#  in the file has a given start site (which you probably don't want anyway
#  because that means you might have clones).  Use BSSeeker_sort.pl to do 
#  this.

#  For now, this program will only generate data for mCG methylation, not
#  mCHG or mCHH.  Those will be ignored for now.

#  Note that for BS_Seeker output:
#     x = unmeth CG, X = meth CG
#     y = unmeth CHG, Y = meth CHG
#     z = unmeth CHH, Z = meth CHH

#  For colors, 
#    No methylation = black (0,0,0)
#    0 > % >= 60     blue 
#    60 > % >= 80    green
#    80 > % >= 100   red

#  For now, the BED file will convert all reads to the + orientation since right now
#  we don't have good sequence coverage to get good stats for both strands.  For example,
#  in a CG pair, the C position will be methylated on the + strand and the G position
#  will be methylated on the opposite strand (relative to the + strand).  This will work
#  fine as long as we can assume CpG methylation is always symmetric.  Which it mostly is.

# Arguments:
#   1) name of sorted BS Seeker output file
#   2) UCSC Genome Browser BED output file for percent methylation
#   3) chromsome name for BED file

my $debug = 0;  #set to 1 if want to debug a position, 0 otherwise
my $test_min_pos = 1572560;  #set to desired chromosome interval
my $test_max_pos = 1572660;

my $test_region = 0;  #don't set.  Used to find regions for detailed output

unless ($ARGV[2]) {

    print "Arguments:\n\tSorted BS Seeker output file\n\tBED percent methylation output file\n\tchromosome name (ex. chr7)\n\n";
    die;
}

open (INPUT, "<$ARGV[0]")
    or die "Couldn't open $ARGV[0]\n";

open (OUTPUT, ">$ARGV[1]") or die "Couldn't open $ARGV[1]\n";

$chrom = $ARGV[2];

print OUTPUT "track name=PercMethylation$chrom description=PercentMethylation$chrom useScore=0 itemRgb=On db=hg18\n";


# We'll implement a queue (push to add elements to end, shift to take
#  elements off the front of the array) of the lines in the input file,
#  taking only as many as we need.  Within the queue loop, a while
#  statement will add lines to the queue until the start of sequence
#  read is after the end of the sequence shifted off.  The sequences
#  in the queue will then be the window for creating a coverage estimate.

# For this program, we'll actually have three queues to maintain.  The second
#  queue will have the list of methylated cytosines cooresponding to the
#  sequences in the first queue.  The third queue will have the list of
#  all CpG's, methylated or not (to find totally unmethylated sites) 
#  corresponding to the sequences in the first queue.  

# Here we go....  First initialize our queue and get chromosome info

$line = <INPUT>;
@fields = split /\s+/, $line;

$orient = substr($fields[2], 0, 1); 

$meth_line = $fields[6];

$meth_line =~ s/\s//;

$seq_length = length ($meth_line);

$fields[3] =~ m/(\d+)(\+|-)0*(\d+)/;

$start_pos = $3;
$end_pos = $start_pos + $seq_length - 1;

push @queue, $start_pos;
push @end_queue, $end_pos;

#### Get CpGs from first line

if ($meth_line =~ m/x/i) { #there are CpGs

    undef @CG_temp_array;
    undef @meth_temp_array;

    while ($meth_line =~ m/\G[-yz]*(x)/gi) { #getting coordinates of CpGs

	$x_matched = $1;

	if ($orient eq "+") {

	    $temp = $start_pos -1 + pos $meth_line; 
	}

	else {  #if on - strand, need to reverse positions

	    $temp = pos $meth_line;
	    $temp = ($seq_length - $temp) -1 + $start_pos; 
	}

	if ($x_matched eq "x") {

	    push @CG_temp_array, $temp;
	}

	elsif ($x_matched eq "X") {

	    push @CG_temp_array, $temp;
	    push @meth_temp_array, $temp;
	}
    }

    $pos = join ':', @CG_temp_array;

    push @CG_queue, $pos;


    if (@meth_temp_array) { #were methylated CpGs

	$meth_pos = join ':', @meth_temp_array;
    
	push @meth_queue, $meth_pos;
    }

    else {
	
	push @meth_queue, "0";
    }

} #end if (were any CpGs in line)

else {

    push @meth_queue, "0";   #this will be false, no methylation

    push @CG_queue, "0";   #no CpGs either
}

if ($debug) {
    
    if ( ($start_pos <= $test_max_pos) && (($start_pos + 76) >= $test_min_pos)) {

	$test_region = 1;

	print "Added read start position $start_pos to queue ----------------\n";

	for ($i = 0; $i <= $#meth_queue; $i++) {

	    print "   queue : $queue[$i]    meth_queue : $meth_queue[$i]\n";
	}
    }
}

$last_assessed_pos = 1;  #keep track of the positions we've written to file

$last_start_pos_added = 0;  #for keeping track if new sequences need to be added to queue


### Now do the rest of the lines in the file


while ($first_start = shift @queue) {

    $first_end = shift @end_queue;

    if ($debug && $test_region) {

	print "first_start = $first_start\n";
	print "first_end = $first_end\n";
    }

    $last_possible_start = $first_end;

    # Now populate the queue with new sequences in the file that
    #   are within the window (start site less than or equal to
    #   $last_possible_start)

    unless ($last_start_pos_added > $last_possible_start) { #need to add more sequences
	                                                    #  to queue

	while ($line = <INPUT>) {  #will quit when run out of reads
	    
	    @fields = split /\s+/, $line;
	    $orient = substr($fields[2], 0, 1); 
	    $meth_line = $fields[6];
	    $fields[3] =~ m/(\d+)(\+|-)0*(\d+)/;
	    $start_pos = $3;

	    $meth_line =~ s/\s//;
	    $seq_length = length($meth_line);

	    $end_pos = $start_pos + $seq_length - 1;

	    push @queue, $start_pos;
	    push @end_queue, $end_pos;
	    
	    if ($meth_line =~ m/x/i) { #there are CpGs
		
		undef @CG_temp_array;
		undef @meth_temp_array;
		
		while ($meth_line =~ m/\G[-yz]*(x)/gi) { #getting coordinates of CpGs
		    
		    $x_matched = $1;
		    
		    if ($orient eq "+") {
			
			$temp = $start_pos -1 + pos $meth_line; 
		    }
		    
		    else {  #if on - strand, need to reverse positions
			
			$temp = pos $meth_line;
			$temp = ($seq_length - $temp) -1 + $start_pos; 
		    }
		    
		    if ($x_matched eq "x") {
			
			push @CG_temp_array, $temp;
		    }
		    
		    elsif ($x_matched eq "X") {
			
			push @CG_temp_array, $temp;
			push @meth_temp_array, $temp;
		    }
		}
		
		$pos = join ':', @CG_temp_array;
		
		push @CG_queue, $pos;
		
		
		if (@meth_temp_array) { #were methylated CpGs
		    
		    $meth_pos = join ':', @meth_temp_array;
		    
		    push @meth_queue, $meth_pos;
		}
		
		else {
		    
		    push @meth_queue, "0";
		}
		
	    } #end if (were any CpGs in line)
	    
	    else {
		
		push @meth_queue, "0";   #this will be false, no methylation
		
		push @CG_queue, "0";   #no CpGs either
	    }
	    

	    if ($debug) {
		
		if ( ($start_pos <= $test_max_pos) && (($start_pos + 76) >= $test_min_pos)) {

		    $test_region = 1;
		    
		    print "\nAdded read start position $start_pos to queue, last assessed position = $last_assessed_pos, last possible start for queue = $last_possible_start\n";
		    
		    print "   queue : $first_start    meth_queue : $meth_queue[0]\n";
		    
		    for ($i = 1; $i <= $#meth_queue; $i++) {
			
		    print "   queue : $queue[$i-1]    meth_queue : $meth_queue[$i]\n";
		    }
		}
		else {

		    $test_region = 0;
		}
	    } #end if (debug)
	    
	    
	    if ($start_pos > $last_possible_start) {  #last sequence in queue
		#will be outside of range
		last;  #exit while loop
	    }
	}

	$last_start_pos_added = $start_pos;

    } #end unless (don't need to add more sequences to queue)

    if ($debug && $test_region) {
	print "\nFinished filling queue\n";
    }


    #################################################################################
    
    # Now compute coverage starting at $last_assessed_pos until
    #   the end of the sequence ($last_possible_start)

    if ($last_assessed_pos < $first_start) {  #in case no coverage before this

	$last_assessed_pos = $first_start;

    }

    while ($last_assessed_pos <= $last_possible_start) { #go until end of seq

	if ($debug && $test_region) {

	    print "\nlast_assessed_pos = $last_assessed_pos\n";
	    print "\nfirst_start = $first_start\n";
	}

	$next_start = 0;

	$coverage_count = 0;  #number of sequences in this interval

	#First count shifted sequence

	$coverage_count++;

	if ($queue[0]) {  #There are sequences in queue, meaning we haven't run out in file

	    for ($i = 0; $i <= $#queue; $i++) {  #go through sequences in queue
		
		$this_start = $queue[$i];
		
		$next_start = $this_start;  #use this later to find next interval
		
		if ($this_start <= $last_assessed_pos) { #within interval
		    
		    $coverage_count++;
		}
		else {
		    
		    last;  #left the interval
		}
	    }
	    
	} #end if (not at end of queue and input file)

	else {  #near the end

	    $next_start = $last_possible_start + 1;
	}

	if ($debug && $test_region) {

	    print "last_possible_start = $last_possible_start\n";
	    print "next_start = $next_start\n";
	}

	$old_assessed_pos = $last_assessed_pos;

	if ($next_start <= $last_possible_start) { #more coverage intervals
	    #  for this shifted seq
	    
	    if (($queue[$#queue] == $next_start) &&
		($old_assessed_pos == $next_start))   { #ran out of data in input file
		
		$last_assessed_pos = $last_possible_start + 1;
		
	    }
	    
	    else {  #more coverage intervals for this shifted sequence
		
		$last_assessed_pos = $next_start;
		
	    }
	}
	
	else { #next_start might be far away

	    $last_assessed_pos = $last_possible_start + 1;
	}
	

	# Now we will go through the methylated residues within this interval
	#   First, create a hash of all methylated positions in these
	#   sequences.  The value will be the counts

	undef %meth_hash;  #get rid of old values
	
	for ($j = 0; $j <= $#meth_queue; $j++) {  #going through sequences

	    $meth_string = $meth_queue[$j];

	    if ($meth_string) {  #this would be 0 (false) if no methylation

		@meths = split /:/, $meth_string;     # /
		
		foreach $position (@meths) {  #going through positions in sequence
		    
		    if (exists $meth_hash{$position}) {
			
			$meth_hash{$position}++;
		    }
		    
		    else {
			
			$meth_hash{$position} = 1;
		    }
		    
		} #end foreach (going through methyl positions in sequence)
	    } #end if (the are methyl positions in sequence)
        } #end for (going through sequences in queue)
	

	# We will also keep a hash of all CpG positions using positions from @CG_queue 

	undef %CG_hash;

	$special_methyl_case = 0;  #tells us if last CpG needs more coverage due to special case

	for ($k = 0; $k <= $#CG_queue; $k++) {  #going through sequences

	    $CG_string = $CG_queue[$k];

	    if ($CG_string) {  #this would be 0 (false) if no CpGs

		@CGs = split /:/, $CG_string;     # /
		    
		foreach $position (@CGs) {  #going through positions in sequence
		    
		    unless (exists $CG_hash{$position}) {
			
			$CG_hash{$position} = 1;  #useless value
		    }
		    
		}

		# Now for dealing with the special case when a read is in the reverse
		#  strand and has a methylated C at its last position.  This will be
		#  converted to a methylated C on the forward strand, which unfortunately
		#  is BEFORE the read's start site.  The code snippet below will allow
		#  it to be counted normally, at least in the case where it would normally
		#  cause an inflated % methylation because it didn't contribute to coverage
		#  but did to methylation

		if ($k > 0) {

		    $this_read_start = $queue[$k-1];  #one sequence already shifted off
                                                      #  queue
		    
		    if (($this_read_start == $next_start) &&  #read starts the next interval
			($CGs[$#CGs] == $this_read_start - 1)) {  #methylation "before" read 
			                                       # start (because in opposite
			                                       # orientation)

			$special_methyl_case = 1;

			$special_methyl_pos = $CGs[$#CGs];  #position needing higher coverage
			
		    }
		}
	    }
        } #end for (going through CG queue)

        
        # Now take the hash of CpG positions, find which are in
        #  the interval, and write percent methylation to file
	
        foreach $position (sort {$a <=> $b} keys %CG_hash) {
	    
            if (($position >= $old_assessed_pos) && ($position < $last_assessed_pos )) {  #in our interval

		$percent = 0;

		if (($special_methyl_case) && ($special_methyl_pos == $position)) {

		    $coverage_count++;
		}

		
		if (exists($meth_hash{$position})) {  #position was methylated

		    $percent = $meth_hash{$position} / $coverage_count;
		    
		}

		if ($percent > 1.00) {  #we have problems....

		    print "WARNING: chromosome = $chrom, position = $position, methylation = $meth_hash{$position}, coverage = $coverage_count, percent = $percent\n\n";
		}

		$percent = sprintf("%.2f",$percent);
		$color = 0;
		
		
		if ($debug) {
		    
		    if ( ($position <= $test_max_pos) && ($position >= $test_min_pos)) {
			print "\nposition = $position\ncoverage count = $coverage_count\npercent = $percent\nnext start = $next_start\nlast assessed position = $last_assessed_pos\nlast possible start = $last_possible_start\nold assessed position = $old_assessed_pos\n";
		    }
		}
		
		
		if ($percent == 0) {
		    $color = "0,0,0";
		}
		
		#if ($percent > 0 && $percent <= .6) {
		elsif ($percent > 0 && $percent <= .6) {
		    $color = "27,74,210";   #blue
		    
		}
		elsif ($percent > .6 && $percent <= .8) {
		    $color = "27,210,57";   #green
		    
		}
		elsif ($percent > .8) {
		    $color = "210,27,27";   #red
		    
		}
		
		$temp = $position - 1;  #note that we need to subtract positions by one

		if ($debug) {  #we'll add extra information to the percent data

		    $percent = $percent . "-" . $coverage_count;

		}
		
		if ($color) { #in case you don't want to output 0% methylation
		    
		    print OUTPUT "$chrom\t$temp\t$position\t$percent\t0\t+\t0\t0\t$color\n";
		}
            }
	    
        } #end foreach (going through methylated positions)
	
    }  #end while (going through shifted sequence, going through intervals)
	
    shift @meth_queue;
    
    shift @CG_queue;
    
} #end while (going through queue)
    


