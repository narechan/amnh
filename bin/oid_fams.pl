#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('o:', \%opts);
my $oid  = $opts{'o'};

opendir (E, "$oid");
my @families = readdir (E);
shift @families;
shift @families;
closedir (E);

my $partdist = {};
foreach my $family (@families){
    my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$oid/$family/FAMILY");

    my $sp = {};
    my $counter = 0;
    while (my $sequence_obj = $seqin -> next_seq()){
	$counter++;
	my $i       = $sequence_obj -> display_id();
	my ($id, $acc) = split (/\#/, $i);
	    
	if (($id eq "PE001") or ($id eq "PE002") or ($id eq "PE003") or ($id eq "PE004") or ($id eq "PE005")){
	    $partdist->{$family}->{'low'}++;
	}
	elsif (($id eq "PE006") or ($id eq "PE007") or ($id eq "PE008") or ($id eq "PE009") or ($id eq "PE010")){
	    $partdist->{$family}->{'high'}++;
	}
	else {
	    $partdist->{$family}->{'other'}++;
	}
    }

}
foreach my $gene (sort keys %$partdist){                          
    my @cats;                                                      
    foreach my $cat ("low", "high", "other"){          
        if ($partdist->{$gene}->{$cat}){                          
            push (@cats, $partdist->{$gene}->{$cat});             
        }                                                        
        else {                                                     
            push (@cats, 0);                                      
        }                                                         
    }                                                              
    my $cats = join "\t", @cats;                                   
    print "$gene\t$cats\n";      
}
