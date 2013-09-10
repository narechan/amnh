#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('a:b:n:', \%opts);
my $afile  = $opts{'a'};
my $bfile  = $opts{'b'};
my $node   = $opts{'n'};

my $a = {};
open (A, "$afile");
while (my $line = <A>){
    chomp $line;
    my @a = split (/\t/, $line);
    $a->{$a[0]}->{$a[2]} = $a[10];
}
close (A);

my $b = {};
open (B, "$bfile");
while (my $line = <B>){
    chomp $line;
    my @b = split (/\t/, $line);
    $b->{$b[0]} = $b[1];
}
close (B);

foreach my $gene (keys %$b){
    my $pbs = $a->{$gene}->{$node};
    
    if ($a->{$gene}){
	print "$gene\t$node\t$pbs\t$b->{$gene}\n";
    }
    else {
	print "$gene\t$node\tNOINF\t$b->{$gene}\n";
    }
}
