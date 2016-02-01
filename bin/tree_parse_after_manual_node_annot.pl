#!/usr/bin/perl -w

# -r is the raw nodes file from tree_parse.pl
# -a is the annotated nodes file from tree_parse.pl (eliminate all nodes you're not interested in)
#    everything else will become the core.

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::TreeIO;

my ($rawfile, $annotfile);
GetOptions(
	   'r|rawfile=s'      => \$rawfile,
	   'a|annotfile=s'    => \$annotfile,
	   );

my $taxa = {};
open (R, "$rawfile");
while (my $line = <R>){
    chomp $line;
    my ($index, $genes, $support, $count) = split (/\t/, $line);
    my @genes = split (/\,/, $genes);
    foreach my $gene (@genes){
	$taxa->{$gene} = 1;
    }
}
close (R);

my $counter = 0;
open (A, "$annotfile");
while (my $line = <A>){
    chomp $line;
    $counter++;
    my ($index, $genes, $support, $count) = split (/\t/, $line);
    my @genes = split (/\,/, $genes);
    foreach my $gene (@genes){
        delete $taxa->{$gene};
    }
    print "$counter\t$genes\t$support\t$count\n";
}
close (A);

$counter++;
my @join;
my $corecounter = 0;
foreach my $tax (keys %$taxa){
    $corecounter++;
    push (@join, $tax);
}
my $join = join ",", @join;
#print "$counter\t$join\tNA\t$corecounter\n";
print "$counter\t$join\t100\t$corecounter\n";
