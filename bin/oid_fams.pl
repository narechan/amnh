#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('o:l:', \%opts);
my $oid  = $opts{'o'};
my $list = $opts{'l'};

my @species;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    push (@species, $line);
}
close (L);

opendir (E, "$oid/data");
my @families = readdir (E);
shift @families;
shift @families;
closedir (E);

my %data;
my $count;
foreach my $family (@families){
    next if ($family eq "singlets");
    my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$oid/data/$family/FAMILY");

    my $counter = 0;
    while (my $sequence_obj = $seqin -> next_seq()){
	$counter++;
	my $id       = $sequence_obj -> display_id();
	my ($species, $acc) = split (/\#/, $id);
	$data{$family}{$species}++;
    }
    $count{$family} = $counter;

}

foreach my $family (sort keys %data){
    print "$family\t";
#    foreach my $species (keys %{$data{$family}}){
    foreach my $species (@species){
	my $spcnt = $data{$family}{$species};
	my $sppct = sprintf ("%.2f", ($spcnt / $count{$family}));
	print "$sppct\t";
    }
    print "\n";
}
