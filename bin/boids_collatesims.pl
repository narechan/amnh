#!/usr/bin/perl

# use the following modules 
use Getopt::Std;
use Statistics::Descriptive;

# get inputs                                                                    
my %opts = ();
getopts ('i:r:', \%opts);
my $indir = $opts{'i'}; #contains subdirs corresponding to the number of trees modeled
my $reps  = $opts{'r'}; #the number of cf reps

opendir (I, "$indir");
my @treedirs = sort (readdir (I));
shift @treedirs;
shift @treedirs;
closedir (I);

foreach my $treedir (@treedirs){
    my @radmeans;
    my @jaccards;
    my @tprs;
    my @tnrs;
    my @fprs;
    my @fnrs;
    for (my $i = 1; $i <= $reps; $i++){
	open (F, "$indir/$treedir/$i.optics.radmean");
	while (my $line = <F>){
	    my ($tpr, $tnr, $fpr, $fnr, $j, $randindex) = split (/\t/, $line);
	    push (@tprs, $tpr);
	    push (@tnrs, $tnr);
	    push (@fprs, $fpr);
	    push (@fnrs, $fnr);
	    push (@jaccards, $j);
	    push (@radmeans, $randindex);
	}
    }
    
    my $radmean = ave (\@radmeans);
    my $jmean = ave (\@jaccards);
    my $tp = ave (\@tprs);
    my $tn = ave (\@tnrs);
    my $fp = ave (\@fprs);
    my $fn = ave (\@fnrs);

    print "$treedir\t$radmean\t$jmean\t$tp\t$tn\t$fp\t$fn\n";
}

sub ave {
    my $array = shift;
    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@$array);
    my $mean  = $statobj->mean();
    return ($mean);
}



    
