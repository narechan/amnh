#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('d:', \%opts);
my $dir  = $opts{'d'};

opendir (D, "$dir");
my @resdirs = sort(readdir (D));
shift @resdirs;
shift @resdirs;
closedir (D);

my $data = {};
foreach my $resdir (@resdirs){

    # store the expt's data
    open (F, "$dir/$resdir/summary");
    while (my $line = <F>){
        chomp $line;
	
	my @array = split (/\t/, $line);
	my $cov = $array[2];
	my $refcov = $array[7];
	$data->{$cov}->{$resdir} = $refcov;
    }
    close (S);
}


# print out the sorted data
my %orgs;
foreach my $cov (sort {$a <=> $b} keys (%$data)){
    print "$cov\t";
    foreach my $org (sort keys %{$data->{$cov}}){
	$orgs{$org} = 1;
	print "$data->{$cov}->{$org}\t";
    }
    print "\n";
}
my $orgs = join "\t", sort keys %orgs;
print "\t$orgs\n";
