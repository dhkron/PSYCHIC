#!/usr/bin/perl
use strict;

# combines the files
if (@ARGV != 4)
{ print "Please provide filename count_thr prob_thr binsize\n";exit;}
my $infile  = $ARGV[0]; open(IN,"<$infile"); open (OUT,">tmp");

my $flag = 0; my $prev_state; my $start; my $end; my $prob; my $chr; my $prob_str;
while (my $line = <IN>)
{ chomp $line;
  my @liner = split('\s',$line);

  # initial setup
  if ($flag == 0)
  { $chr = $liner[0]; $prev_state = $liner[3]; $start = $liner[1]; $end = $liner[2]; $prob = $liner[3+$prev_state]; $prob_str = $liner[3+$prev_state].','; $flag = 1; next;}

  # no transition
  if ($liner[3] == $prev_state)
  { $chr = $liner[0]; $prev_state = $liner[3]; $end = $liner[2]; $prob *= $liner[3+$prev_state]; $prob_str .= $liner[3+$prev_state].',';}
  else
  { my $dan = ($end-$start+1)/$ARGV[3]; my @danner = split(',',$prob_str); my $new = median_fn(@danner);
    print OUT "$chr\t$start\t$end\t$prev_state\t$dan\t$prob\t$new\n";
    $chr = $liner[0]; $prev_state = $liner[3]; $start = $liner[1]; $end = $liner[2]; $prob = $liner[3+$prev_state]; $prob_str = $liner[3+$prev_state].',';
  }
}
my $dan = ($end-$start+1)/$ARGV[3]; my @danner = split(',',$prob_str); my $new = median_fn(@danner);
print OUT "$chr\t$start\t$end\t$prev_state\t$dan\t$prob\t$new\n";
close OUT; close IN;

sub median_fn
{ my @array = @_; my $count = scalar @array;
  @array = sort { $a <=> $b } @array; 
  if ($count % 2) { 
      return $array[int($count/2)]; 
  } else { 
      return ($array[$count/2] + $array[$count/2 - 1]) / 2;}
}

open (IN,"<tmp");
my @file = <IN>; close IN; chomp @file; system("rm tmp");

for (my $i = 0; $i < scalar @file; $i++)
{ chomp $file[$i];
  my @liner = split('\t',$file[$i]);
  $liner[4] = int($liner[4]); # to balance out the number of bins

  # first line -- if first three bins are different than the next bins, we change these to next bins, since hmm starts shaky   
  if ($i == 0)
  {  if ($liner[4] <= 3) { my $next = $file[$i+1]; chomp $next; my @ne = split('\t',$next); print "$liner[0]\t$liner[1]\t$liner[2]\t$ne[3]\n";}
     else                { print "$liner[0]\t$liner[1]\t$liner[2]\t$liner[3]\n";}
  } # last line -- similar to first line
  elsif ($i == ((scalar @file) - 1))
  { if ($liner[4] <= 3) { my $pre = $file[$i-1]; chomp $pre; my @pr = split('\t',$pre); print "$liner[0]\t$liner[1]\t$liner[2]\t$pr[3]\n";}
    else                { print "$liner[0]\t$liner[1]\t$liner[2]\t$liner[3]\n";}
  }
  else
  {   # bin count > thr_bin_count
      if ($liner[4] >= $ARGV[1]) {print "$liner[0]\t$liner[1]\t$liner[2]\t$liner[3]\n";}
      else 
      { # even if bin count is small, median prob are good
        if ($liner[6] >= $ARGV[2]) {print "$liner[0]\t$liner[1]\t$liner[2]\t$liner[3]\n";}
        # small bin and median are bad, so look at previous and next state --> pr==ne? pr; 2.
        else 
	{ my $previous = $file[$i-1]; chomp $previous;
          my $next     = $file[$i+1]; chomp $next;
          
          my @pr = split('\t',$previous); my @ne = split('\t',$next);
          if ($pr[3] == $ne[3]) {print "$liner[0]\t$liner[1]\t$liner[2]\t$pr[3]\n";}
          else {print "$liner[0]\t$liner[1]\t$liner[2]\t2\n";}       
         }
      }      
  }
}
