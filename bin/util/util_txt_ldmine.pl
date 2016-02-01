#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:i:x:y:', \%opts);
my $ldfile  = $opts{'l'};
my $ldfile2 = $opts{'y'};
my $infile  = $opts{'i'};
my $ildfile = $opts{'x'};

my $ldstrings = {};
my $lds = {};
open (L, "$ldfile");
while (my $line = <L>){
    chomp $line;
    my ($comp, $ld) = split (/\t/, $line);
    $ldstrings->{$comp} = $ld;
    my ($one, $two) = split (/v/, $comp);
    $lds->{$one}->{$two} = $ld;
    $lds->{$two}->{$one} = $ld;
}
close (L);

my $ldstrings2 = {};
open (Y, "$ldfile2");
while (my $line = <Y>){
    chomp $line;
    my ($comp, $ld) = split (/\t/, $line);
    $ldstrings2->{$comp} = $ld;
}
close (Y);

#=head
open (X, "$ildfile");
while (my $line = <X>){
    my ($comp, $pval, $ild) = split (/\t/, $line);
    if (($pval <= 0.05) and (($ldstrings->{$comp} == 0) or ($ldstrings2->{$comp} == 0))){
	print "$comp\t$pval\t$ldstrings->{$comp}\t$ldstrings2->{$comp}\n";
    }
}
close (X);
#=cut

=head
my $qs = {};
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    $qs->{$line} = 1;
}
close (I);

foreach my $id1 (sort keys %$qs){
    foreach my $id2 (sort keys %$qs){
	next if ($id2 eq $id1);
	print "$id1\t$id2\t$lds->{$id1}->{$id2}\n";
    }
    delete $qs->{$id1};
}
=cut
