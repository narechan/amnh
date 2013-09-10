#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:a:b:c:', \%opts);
my $indir    = $opts{'i'};
my $afile = $opts{'a'};
my $bfile = $opts{'b'};
my $cfile = $opts{'c'};

opendir (D, "$indir");
my @exptfiles = sort (readdir(D));
shift @exptfiles;
shift @exptfiles;
closedir (E);

my $readprc = {};
foreach my $file (@exptfiles){
    open (F, "$indir/$file");
    while (my $line = <F>){
	chomp $line;
	my @line = split (/\t/, $line);
	shift @line;
	my $newline = join "\t", @line;
	$readprc->{$file} = $newline;
    }
    close (F);
}

my $genomeprc = {};
open (A, "$afile");
while (my $line = <A>){
    chomp $line;
    my ($comp, $rest) = split (/\t/, $line, 2);
    $genomeprc->{$comp} = $rest;
}
close (A);

my $readdist = {};
my $genomedist = {};
open (B, "$bfile");
while (my $line = <B>){
    chomp $line;
    my ($comp, $g, $r) = split (/\t/, $line);
    $genomedist->{$comp} = $g;
    $readdist->{$comp} = $r;
}
close (B);

open (C, "$cfile");
while (my $line = <C>){
    chomp $line;
    my ($comp, $ddh) = split (/\t/, $line);
    print "$comp\t$ddh\t$genomeprc->{$comp}\t$genomedist->{$comp}\t$readprc->{$comp}\t$readdist->{$comp}\n";
}
close (C);
