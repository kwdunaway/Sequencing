#!/usr/bin/perl 
use strict; use warnings;
use Math::BigFloat;
use Math::BigInt;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 9-17-14
#
# Gets the probability of finding 
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) p value (ex: .05)
    2) window size (ex: 10)
    3) number of times sampled (ex: 87107)
" unless @ARGV == 3;

my $p = shift(@ARGV);
my $q = 1-$p;
my $window = shift(@ARGV);
my $sampnum = shift(@ARGV);


###################
# Another section #
###################

my $result = $q**$window;
my $total = $result;
print "$result\t$total\n";
$result = $window * ($q**($window - 1)) * $p;
$total += $result;
print "$result\t$total\n";
$result = ($window*($window-1)/2) * ($q**($window - 2)) * $p**2;
$total += $result;
print "$result\t$total\n";
print $sampnum*(1-$total) , "\n";

my $x = 10;
my $y = 2;
my $tmp = $x->bnok($y);
print $tmp , "\n";

__END__
for($n = 0; $n<$window; $n+){
	

}