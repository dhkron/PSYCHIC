#!/usr/bin/perl

use strict;

MAIN : {

    my $prev_state = 0;
    my $prev_start;
    my $prev_end;
    while (my $line = <>) {
	chomp $line;
	my ($chr, $start, $end, $state) = split(/\t/,$line);
	if ($prev_state == 0) {
	    unless ($state == 3) {
		next;
	    }
	}
	if ($state != $prev_state) {
	    if (($state == 3) && ($prev_state == 0)) {
		print $chr . "\t" . $start . "\t";
	    }
	    if (($state == 3) && ($prev_state == 1)) {
		print $prev_end . "\n" . $chr . "\t" . $start . "\t";
	    }

	}
	$prev_state = $state;
	$prev_start = $start;
	$prev_end = $end;
    }
    if ($prev_state == 1) {
	print $prev_end . "\n";
    }
}
