#!/usr/bin/perl

use strict;

MAIN : {

    my ($genome_size_file, $chr) = @ARGV;
    if ((not defined $genome_size_file) ||
	(not defined $chr)) {
	die ("Usage: ./hmm-state_caller.pl <*fai file> <chr>\n");
    }

	my $genome_size;
	open(FILE,$genome_size_file);
	while (my $line = <FILE>) {
		chomp $line;
		my ($chr, $size) = split(/\t/,$line);
		$genome_size->{$chr} = $size;
	}
	close(FILE);

    my $prev_state;
    my $prev_start;
    my $prev_end;
    while (my $line = <STDIN>) {
	chomp $line;
	my ($chr, $start, $end, $state) = split(/\s/,$line);
	if ($state != $prev_state) {
	    if ($prev_state == 2) {
		print $chr . "\t" . $start . "\t";
	    }
	    if ($state == 2) {
		print $prev_end . "\t" . $prev_state . "\n";
	    }
	    if (($state != 2) && ($prev_state != 2)) {
		print $prev_end . "\t" . $prev_state . "\n" . $chr . "\t" . $start . "\t";
	    }
	}
	$prev_state = $state;
	$prev_start = $start;
	$prev_end = $end;
    }

        if ($prev_state != 2)
	{print $genome_size->{$chr} . "\t" . $prev_state . "\n";}
    
}
