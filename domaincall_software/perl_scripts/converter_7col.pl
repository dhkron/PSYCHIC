#!/usr/bin/perl
use strict;

MAIN : {

    my $spot = 6;
    my $lines;
    my $sum;
    my @array;
    my $i = 0;    

    foreach my $line (<>) {
	chomp $line;
        $array[$i] = $line;
	my @row = split(/\t/,$line);
	if ($row[$spot - 1] == 1) {
	    $lines->{1}++;
	    $sum->{1} += $row[3];
	}
	if ($row[$spot - 1] == 2) {
	    $lines->{2}++;
	    $sum->{2} += $row[3];
	}
	if ($row[$spot - 1] == 3) {
	    $lines->{3}++;
	    $sum->{3} += $row[3];
	}

        $i++;
    }

    my $mean;
    my $ref;
    foreach my $key (sort {$a <=> $b} keys %$lines) {
	$mean->{$key} = $sum->{$key}/$lines->{$key};
	$ref->{$mean->{$key}} = $key;
    }

    my @order;
    foreach my $key (sort {$a <=> $b} keys %$ref) {
	push(@order,$key);
    }

    my @sorted = sort {$a <=> $b} @order;

    my $change;
    my $state_id;
    for (my $i = 1; $i <= 3; $i++) {
	$change->{$ref->{$sorted[$i - 1]}} = $i;
	$state_id->{$i} = $ref->{$sorted[$i - 1]};
    }

    foreach my $line (@array) {
	chomp $line;
	my @row = split(/\t/,$line);
	my $digit = int($row[$spot - 1]);
	print(join("\t",
		   "chr" . $row[0],
		   $row[1],
		   $row[2],
		   $change->{$digit},
		   $row[5 + $state_id->{"1"}],
		   $row[5 + $state_id->{"2"}],
		   $row[5 + $state_id->{"3"}]) . "\n");
    }

}
