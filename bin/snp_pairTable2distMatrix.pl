#!/usr/bin/perl

use Getopt::Std;
use List::Util qw(shuffle);

my %opts = ();
getopts ('i:d:r:', \%opts);
my $infile    = $opts{'i'};
my $datafile = $opts{'d'};
my $randoms  = $opts{'r'};

# parse strain file with st/cc info and the correct order
my @order;
my $sts = {};
open (D, "$datafile");
while (my $line = <D>){
    chomp $line;
    my ($strain, $crap, $st, $cc) = split (/\t/, $line);
    push (@order, $strain);
    push (@{$sts->{$st}}, $strain);
}
close (O);

# parse the snp table                                                                                      
my $matrix = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my ($ref, $qry, $snps) = split (/\t/, $line);
    $matrix->{$ref}->{$qry} = $snps;
}
close (F);


# subsample taxa from sts
my $samps = {};
if ($randoms){
    foreach my $st (keys %$sts){
	my @stsample = (shuffle(@{$sts->{$st}}))[0..$randoms-1];
	foreach my $sts (@stsample){
	    $samps->{$sts} = 1;
	}
    }
}
else {
    foreach my $sts (@order){
	$samps->{$sts} = 1;
    }
}

my @header;
#push (@header, "\t");
foreach my $r (@order){
    (push (@header, $r)) if (exists ($samps->{$r}));
}
my $header = join "\t", @header;
print "$header\n";

foreach my $r (@order){
    next unless (exists ($samps->{$r}));
    print "$r\t";
    my @string;
    foreach my $q (@order){
	next unless (exists ($samps->{$q}));
	push (@string, $matrix->{$r}->{$q});
    }
    my $string = join "\t", @string;
    print "$string\n";
}
