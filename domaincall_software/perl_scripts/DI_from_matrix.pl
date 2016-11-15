#!/usr/bin/perl
use strict;

MAIN : {

    my ($matrix, $bin_size, $window_size, $genome_size_file) = @ARGV;
    if ((not defined $matrix) ||
	(not defined $bin_size) ||
	(not defined $window_size)) {
	die ("Usage: ./DI_from_matrix.pl <matrix> <bin size> <window size> <genome size file>\n");
     }

    my $genome_size;
    open(FILE,$genome_size_file);
    while (my $line = <FILE>) {
	chomp $line;
	my ($chr, $size) = split(/\t/,$line);
	$genome_size->{$chr} = $size;
    }
    close(FILE);

    my $bins = $window_size/$bin_size;

    open(FILE,$matrix);
    my @array = <FILE>;
    close(FILE);

    my $array_size = scalar(@array);
    
    my $starter;
    
    for (my $i = 0; $i < $array_size; $i++) {
	my $line = $array[$i];
	chomp $line;
	my ($chr, $start, $end, @row) = split(/\t/,$line);
	my $tally = 0;
	foreach my $value (@row) {
	    $tally += $value;
	}
	if ($tally > 0) {
	    $starter = $i;
	    last;
	}
    }



    for (my $i = $starter; $i < $array_size; $i++) {
	my $line = $array[$i];
        chomp $line;
        my ($chr, $start, $end, @row) = split(/\t/,$line);
	my $A = 0;
	my $B = 0;
	for (my $z = $i - $bins; $z < $i; $z++) {
	    unless (($z < $starter) || ($z >= $array_size)) {
		$A += $row[$z];
	    }
	}
	for (my $z = $i + 1; $z <= $i + $bins; $z++) {
	    unless (($z < $starter) || ($z >= $array_size)) {
                $B += $row[$z];
            }
	}
	my $E = ($A + $B)/2;
	my $DI;
	if (($E == 0) || ($A == $B)) {
	    $DI = 0;
	} else {
	    $DI = (($B - $A)/abs($B - $A))*((($A - $E)**2)/$E + (($B - $E)**2)/$E);
	}
	if ($end < $genome_size->{$chr}) {
	    my $tmp = $chr; $tmp =~ s/chr//g;
	    print $tmp . "\t" . $start . "\t" . $end . "\t" . $DI . "\n";
	} else {
	    my $tmp = $chr; $tmp =~ s/chr//g;
	    print $tmp . "\t" . $start . "\t" . $genome_size->{$chr} . "\t" . $DI . "\n";
	}
    }

}
