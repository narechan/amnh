#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $opticsfinalfile, $indir); 
# where indir is the directory that contains the fastas generated from randTrees_seq-gen.pl
# and opticsfinal is the parsed optics classification file from the final frame
GetOptions(
    'h|help'          => \$help,
    'o|opticsfinal=s' => \$opticsfinalfile,
    'i|indir=s'       => \$indir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $opticsfinalfile){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

####MAIN####                                                                         

# parse the optics file
my $opticsclass = {};
open (O, "$opticsfinalfile");
while (my $line = <O>){
    chomp $line;
    my ($cluster, $gene) = split (/\t/, $line);
    $opticsclass->{$gene} = $cluster;
}
close (O);

# parse the true classes coded in the fasta names
# and calculate rand index relative to the optics
# flocking classifications

# parse the flocking info
opendir (I, "$indir");
my $files = {};
my @files = sort (readdir (I));
shift @files;
shift @files;
closedir (I);

foreach my $file (@files){
    $files->{$file} = 1;
}

# calculate the rand index
my $tp = 0;
my $tn = 0;
my $fp = 0;
my $fn = 0;
foreach my $f1 (keys %$files){
    foreach my $f2 (keys %$files){
	next if ($f1 eq $f2);

	# get data for the true clusters
#	my ($cluster1, $gene1) = split (/\_/, $f1);
#	my ($cluster2, $gene2) = split (/\_/, $f2);
	my ($acc1, $cluster1, $rate1) = split (/\_/, $f1);
	my ($acc2, $cluster2, $rate2) = split (/\_/, $f2); 
	
	# get data for the calculated clusters
	my $bcluster1 = $opticsclass->{$f1};
	my $bcluster2 = $opticsclass->{$f2};

	# contingencies
	if (($cluster1 == $cluster2) and ($bcluster1 == $bcluster2)){
	    $tp++;
	}
	elsif (($cluster1 != $cluster2) and ($bcluster1 != $bcluster2)){
	    $tn++;
	}
	elsif (($cluster1 != $cluster2) and ($bcluster1 == $bcluster2)){
	    $fp++;
	}
	elsif (($cluster1 == $cluster2) and ($bcluster1 != $bcluster2)){
            $fn++;
        }
	else{
	    print STDERR "Unknown contingency\n";
	    die;
	}
	print STDERR "$f1\t$f2\t$tp\t$tn\t$fp\t$fn\n";
    }
    delete $files->{$f1};
}

my $randindex = ($tp + $tn) / ($tp + $tn + $fp + $fn);
my $jaccard   = $tp / ($tp + $fp + $fn);

my $tpr;
unless (($tp + $fn) == 0){
   $tpr = $tp / ($tp + $fn);
}
else{
    $tpr = "NA";
}

my $tnr;
unless (($tn + $fp) == 0){
    $tnr = $tn / ($tn + $fp);
}
else{
    $tnr = "NA";
}

my $fpr;
unless (($fp + $tn) == 0){
    $fpr = $fp / ($fp + $tn);
}
else{
    $fpr = "NA";
}

my $fnr;
unless (($fn + $tp) == 0){
    $fnr = $fn / ($fn + $tp);
}
else {
    $fnr = "NA";
}

#print "$tpr\t$tnr\t$fpr\t$fnr\t$jaccard\t$randindex\n";
print "$tpr\t$tnr\t$jaccard\t$randindex\n";
