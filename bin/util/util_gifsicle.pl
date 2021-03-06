#!/usr/bin/perl

use Getopt::Std;
use List::MoreUtils qw(part);

my %opts = ();
getopts ('d:o:', \%opts);
my $dir  = $opts{'d'};
my $outdir = $opts{'o'};

opendir (D, "$dir");
my @resdirs = sort(readdir (D));
shift @resdirs;
shift @resdirs;
closedir (D);

my @sorted = sort { $a <=> $b } @resdirs;

my $i = 0;
my @part = part { int( $i++ / 1000 ) } @sorted;

my $counter = 0;
my @counter;
foreach my $part (@part){
    $counter++;
    my @string;
    foreach my $p (@$part){
	my $string = $dir . "/" . $p;
	push (@string, $string);
    }
    my $frames = join " ", @string;
    `gifsicle -d 10 -D bg $frames > $outdir/$counter.final.gif`;
    push (@counter, "$outdir/$counter.final.gif");
}

my $finalframes = join " ", @counter;
`gifsicle -d 10 -D bg $finalframes > $outdir/final.gif`;
