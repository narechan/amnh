#!/usr/bin/perl

# Purpose: accept a super-directory of sub-matrix summaries and collate stats into table.
#          each directory must contain one subdir for every topology tested, and topologies
#          must be consistent across every sub-matrix directory.
# Options:
#          -i input dir

# use the following modules 
use Getopt::Std;

# get inputs                                                                    
my %opts = ();
getopts ('i:', \%opts);
my $indir = $opts{'i'};

# store the matrix dirs
opendir (DIRS, "$indir");
my @dirs = sort (readdir (DIRS));
shift @dirs;
shift @dirs;
closedir (DIRS);

my $data = {};
foreach my $dir (@dirs){

    # store the topology dirs
    opendir (D, "$indir/$dir");
    my @tdirs = sort (readdir (D));
    shift @tdirs;
    shift @tdirs;
    closedir (D);

    foreach my $tdir (@tdirs){
	open (F, "$indir/$dir/$tdir/integration.out");
	while (my $line = <F>){
	    chomp $line;
	    my ($node, $auc, $aucloess) = split (/\t/, $line);
	    ($data->{$dir}->{$tdir} = $auc) if ($node == 0); #take out the CFI AUCs
	}
	close (F);
    }

}
    
foreach my $m (sort keys %$data){
    my $datastring;
    foreach my $d (sort keys %{$data->{$m}}){
	$datastring .= $data->{$m}->{$d} . "\t";
    }
    chop $datastring;
    print "$m\t$datastring\n";
}
