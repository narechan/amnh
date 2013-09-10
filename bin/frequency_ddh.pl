#!/usr/bin/perl

# Purpose: Accept a single column data file and generate sequence stats

# Options:
#          -i input file
#          -m is the bin maximum
#          -s is the bin step size

# use the following modules 
use Getopt::Std;

# get inputs                                                                    
my %opts = ();
getopts ('i:m:s:', \%opts);
my $infile = $opts{'i'};
my $max   = $opts{'m'};
my $step   = $opts{'s'};

# create the bin datastuc
my $bins = {};
my $counter = 0;
until ($counter > $max){
    $bins->{$counter}->{'YES'} = 0;
    $bins->{$counter}->{'NO'} = 0;
    $bins->{$counter}->{'NA'} = 0;
    $counter += $step;
}

open (FILE, "$infile");
while (my $line = <FILE>){
    chomp $line;
    my ($ddh, $cov) = split (/\t/, $line);
    
    foreach my $bin (sort {$a <=> $b} keys %$bins){
	if (($cov <= $bin) and ($cov >= $bin - $step)){
	    if ($ddh eq 'YES'){
		$bins->{$bin}->{'YES'}++;
	    }
	    if ($ddh eq 'NO'){
                $bins->{$bin}->{'NO'}++;
            }
	    if ($ddh eq 'N/A'){
                $bins->{$bin}->{'NA'}++;
            }

	    last;
	}
	else {
	    next;
	}
    }
}

foreach my $bin (sort {$a <=> $b} keys %$bins){
    print "$bin\t$bins->{$bin}->{'YES'}\t$bins->{$bin}->{'NO'}\t$bins->{$bin}->{'NA'}\n";
}

