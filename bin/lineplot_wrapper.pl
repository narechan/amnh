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

my @treedirs = ("t1.tree", "t2.tree", "t3.tree", "t4.tree", "t5.tree", "t6.tree");

foreach my $resdir (@resdirs){
    foreach my $treedir (@treedirs){
	`R --slave --args $dir/$resdir/summaries/$treedir/coords.out < ~apurva/bin/lineplot.R`;
	`mv Rplots.pdf $dir/$resdir/$treedir/summaries/$resdir.$treedir.pdf`;
    }
}
