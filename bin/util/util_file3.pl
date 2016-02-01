#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:m:o:l:s:', \%opts);
my $mlstfile  = $opts{'m'};
my $lookupfile = $opts{'l'};
my $stfile = $opts{'s'};
my $indir   = $opts{'i'};
my $outdir = $opts{'o'};

my $mlst = {};
open (M, "$mlstfile");
while (my $line = <M>){
    chomp $line;
    my @line = split (/\,/, $line);
    my $acc = shift @line;
    my $st  = pop @line;
    $mlst->{$st}->{$acc} = 1;
}
close (M);
    
my $lookup = {};
open (L, "$lookupfile");
while (my $line = <L>){
    chomp $line;
    my ($strain, $acc) = split (/\t/, $line);
    $lookup->{$acc} = $strain;
}
close (L);

open (F, "$stfile");
while (my $line = <F>){
    chomp $line;
    `mkdir -p $outdir/$line`;

    foreach my $acc (sort keys %{$mlst->{$line}}){
	`cp $indir/$lookup->{$acc}/tblastn.nuc.seq $outdir/$line/$lookup->{$acc}.seq`;
    }
}
close (F);

opendir (D, "$outdir");
my @ccs = sort (readdir (D));
shift @ccs;
shift @ccs;

foreach my $cc (@ccs){
    `cat $outdir/$cc/* > $outdir/$cc.fa`;
    `rm -rf $outdir/$cc`;
}

