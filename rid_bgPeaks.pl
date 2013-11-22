#!/usr/bin/perl
use strict;
use warnings;

die "
enter a bg peaks file followed by background 30mers file
" unless @ARGV;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

open (CHR, $file1);
my @peaks = <CHR>;
close CHR;
my ($pc, $d, $d1, $ps, $pe) = split ("\t", $peaks[0]);
$ps -= 200;
$pe += 200;

open (IN, $file2);

LINE: while (<IN>){
    my ($chr, $start, $stop) = split;

    if ($chr eq $pc and $start >= $ps && $start <= $pe){
	next;
    }
    if ($chr eq $pc &&  $start > $pe){
	shift @peaks;
	unless ($peaks[0]){
	    last LINE;
	}
	($pc, $d, $d1, $ps, $pe) = split ("\t", $peaks[0]);
	$ps -= 200;
	$pe += 200;
	if ($chr eq $pc and $start >= $ps && $start <= $pe){
	    next;
	} else {
	    print "$chr\t$start\t$stop\n";
	}
    } else {
	print "$chr\t$start\t$stop\n";
	
    }
}
	
close IN;
