#!/usr/bin/perl

use strict;

MAIN : {

    my ($file, $boundry) = @ARGV;
    if ((not defined $file) ||
	(not defined $boundry)) {
	die ("Usage: ./file_ends_cleaner.pl <hmm output file> <hmm Input DI file>\n");
    }

    open(FILE,$file);
    my @array = <FILE>;
    close(FILE);
    
    open(FILE,$boundry);
    my @sizer = <FILE>;
    close(FILE);
    
    my $size = scalar(@sizer);
    
    for (my $i = 2; $i < $size + 2; $i++) {
	my $other = $sizer[$i - 2];
	chomp $other;
	my @o_row = split(/\t/,$other);
	my $line = $array[$i];
	$line =~ s/^\s+//g;
	$line =~ s/\s+/\t/g;
	chomp $line;
	my @row = split(/\t/,$line);
	print(join("\t",@o_row,@row) . "\n");
    }
	
}
