#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:i:c:', \%opts);
my $list  = $opts{'l'}; #tuberculist
my $in    = $opts{'i'}; #snpeff txt output file (can be concatenate of many)
my $cutoff = $opts{'c'};

my $listdata = {};
open (F, "$list");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    $listdata->{$data[0]} = $data[3];
}
close (F);

# search for dup ref coords and determine whether it's an overlapping gene (keep)
# or a het where both alleles are diff from the ref (discard)

my $pos = {};
open (I, "$in");
while (my $line = <I>){
    chomp $line;
    my @data = split (/\t/, $line);
    if ($data[4] eq "SNP"){
	$pos->{$data[1]}->{$data[3]} = 1;
    }
}
close (I);

my $skippos = {};
foreach my $p (sort keys %$pos){
    my $counter = 0;
    foreach my $a (sort keys %{$pos->{$p}}){
	$counter++;
    }
    if ($counter > 1){
	$skippos->{$p} = 1;
    }
}

# both cannot be simeltaneously true for first position
open (I, "$in");
while (my $line = <I>){
    chomp $line;
    my @data = split (/\t/, $line);
    unless  ($data[4] eq "SNP"){ #snps only
	next;
    }
    unless ($data[5] eq "Hom"){ #homozygous snps only (getting rid of hets where one allele is the ref)
	next;
    }
    unless ($data[6] >= $cutoff){ #above quality cutoff
	next;
    }
    if (exists ($skippos->{$data[1]})){ #getting rid if hets different from the ref
	next;
    }
    else {
	print "$in\t";
	print "$line\t";
	if (exists ($listdata->{$data[10]})){
	    print "$listdata->{$data[10]}\n";
	}
	else {
	    print "NO FUNC CAT\n";
	}
    }
}
close (I);
