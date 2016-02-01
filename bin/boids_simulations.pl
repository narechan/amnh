#!/usr/bin/perl -w

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Statistics::Descriptive;

my ($help, $treesamp, $ldconfig, $cfconfig, $cfreps, $procs, $outdir, $fracmissing, $rates);

# hard coded in spots for simulation (especially check the cf iterations!)

GetOptions(
    'h|help'          => \$help,
    't|treesamp=s'    => \$treesamp,
    'l|ldconfig=s'    => \$ldconfig,
    'c|cfconfig=s'    => \$cfconfig,
    'o|outdir=s'      => \$outdir,
    'r|cfreps=s'      => \$cfreps,
    'p|procs=s'       => \$procs,
    'm|missing=s'     => \$fracmissing,
    'x|rates=s'       => \$rates,
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir/cf`;

#####MAIN#####

# run the simulations script
print STDERR "$treesamp sims\n";
#if ($fracmissing){
#    `randTrees_seq-gen.pl -t 10 -x 100 -g 100 -l 100 -y 100 -z $treesamp -m $fracmissing -o $outdir/data &> $outdir/stderr.sim`;
#}
#elsif ($rates){
    `randTrees_seq-gen-partrates.pl -t 10 -x $treesamp -m $fracmissing -g 100 -l 100 -r $rates -o $outdir/data &> $outdir/stderr.sim`;
#}
#else {
#    print "no sim method!\n";
#    die;
#}

# convert to nexus
print STDERR "$treesamp convert\n";
`converter_anytonexus_lite.pl -i $outdir/data/fastas -o $outdir/data`;

# calc the lds (inherently parallel -p hardcoded)
print STDERR "$treesamp lds\n";
#`ld_massivePW.pl -c $ldconfig -p 20 -m $outdir/data/nexus -o $outdir/data_lds &> $outdir/stderr.data_lds`;
`ld_massivePW.pl -c $ldconfig -p $procs -m $outdir/data/nexus -o $outdir/data_lds &> $outdir/stderr.data_lds`;

# parse the lds
print STDERR "$treesamp ldparse\n";
`ld_massivePW_parse.pl -m $outdir/data/nexus -o $outdir/data_lds`;

# MDS with kmeans and optics of the collpased 2D space
`R --slave --args $outdir/data_lds/matrix.header < ~apurva/amnh/bin/mds.R > $outdir/mds.coords`;
`mv Rplots.pdf $outdir/mds.coords.pdf`;
`boids_optics.pl -m 10 -x 0.1 -i $outdir/mds.coords -j ~apurva/packages/elki.jar -o $outdir > $outdir/mds.clusters.optics`;
`boids_kmeans.pl -i $outdir/mds.coords -k $treesamp -j ~apurva/packages/elki.jar -o $outdir > $outdir/mds.clusters.kmeans`;
`boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/mds.clusters.optics 1>$outdir/mds.clusters.optics.radmean 2>$outdir/stderr.mds.clusters.optics.radmean`;
`boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/mds.clusters.kmeans 1>$outdir/mds.clusters.kmeans.radmean 2>$outdir/stderr.mds.clusters.kmeans.radmean`;
#`R --slave --args $outdir/mds.coords $treesamp < ~apurva/amnh/bin/kmeans.R > $outdir/mds.clusters.kmeans`;
#`mv Rplots.pdf $outdir/mds.clusters.kmeans.pdf`;

# kmedoids (PAM) directly on distances
`R --slave --args $outdir/data_lds/matrix.header $treesamp < ~apurva/amnh/bin/kmedoids.R > $outdir/pam.clusters`;
`mv Rplots.pdf $outdir/pam.clusters.pdf`;
`cut -f1 $outdir/pam.clusters > c1`;
`cut -f2 $outdir/pam.clusters > c2`;
`paste c2 c1 > $outdir/pam.clusters`;
`rm c1 c2`;
`boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/pam.clusters 1>$outdir/pam.clusters.radmean 2>$outdir/stderr.pam.clusters.radmean`;


# hierarchical clustering on euclidian transformation of the distances
`R --slave --args $outdir/data_lds/matrix.header $treesamp average < ~apurva/amnh/bin/hclust.R > $outdir/hclust.clusters`;
`mv Rplots.pdf $outdir/hclust.clusters.pdf`;
`cut -f1 $outdir/hclust.clusters > c1`;
`cut -f2 $outdir/hclust.clusters > c2`;
`paste c2 c1 > $outdir/hclust.clusters`;
`rm c1 c2`;
`boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/hclust.clusters 1>$outdir/hclust.clusters.radmean 2>$outdir/stderr.hclust.clusters.radmean`;

### fork the clusterflock processes
my $pm = Parallel::ForkManager->new($procs);
for (my $i = 1; $i <= $cfreps; $i++){
    $pm->start and next;

    print STDERR "$treesamp cflocks $i\n";
    `clusterflock.pl -i $outdir/data/fastas -c $cfconfig -l $outdir/data_lds/table2 -s all -b 1 -o $outdir/cf/$i-cf &> $outdir/cf/stderr.$i-cf`;

#    `clusterflock.pl -i $outdir/data/fastas -c $cfconfig -l $outdir/data_lds/table2 -s all -b 1 -d -x -o $outdir/cf/$i-cf &> $outdir/cf/stderr.$i-cf`;

    print STDERR "$treesamp drawing $i\n";
    `boids_gnuplot.pl -i $outdir/cf/$i-cf/logs/500.log -d 50 -o $outdir/cf/$i-cf/images/500.gif`;

    print STDERR "$treesamp optics $i\n";
    `boids_optics.pl -i $outdir/cf/$i-cf/logs/500.log -m 10 -x 0.1 -j ~apurva/packages/elki.jar -o $outdir/cf/$i-cf/flocks > $outdir/cf/$i-cf/flocks/500.optics`;

    print STDERR "$treesamp kmeans $i\n";
    `boids_kmeans.pl -i $outdir/cf/$i-cf/logs/500.log -k $treesamp -j ~apurva/packages/elki.jar -o $outdir/cf/$i-cf/flocks > $outdir/cf/$i-cf/flocks/500.kmeans`;
    
    print STDERR "$treesamp optics rand means $i\n";
    `boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/cf/$i-cf/flocks/500.optics 1>$outdir/cf/$i.optics.radmean 2>$outdir/cf/stderr.$i-optics-radmean`;
    
    print STDERR "$treesamp kmeans rand means $i\n";
    `boids_clustersimanalysis.pl -i $outdir/data/fastas -o $outdir/cf/$i-cf/flocks/500.kmeans 1>$outdir/cf/$i.kmeans.radmean 2>$outdir/cf/stderr.$i-kmeans-radmean`;

    $pm->finish;
}
$pm->wait_all_children;

# get averages of simulation metrics across all cf reps
my @ktpr;
my @ktnr;
my @kjaccard;
my @krandindex;

my @otpr;
my @otnr;
my @ojaccard;
my @orandindex;

for (my $i = 1; $i <= $cfreps; $i++){
    open (K, "$outdir/cf/$i.kmeans.radmean");
    while (my $line = <K>){
	chomp $line;
	my @line = split (/\t/, $line);
	push (@ktpr, $line[0]);
	push (@ktnr, $line[1]);
	push (@kjaccard, $line[2]);
	push (@krandindex, $line[3]);
    }
    close (K);
}

for (my $i = 1; $i <= $cfreps; $i++){
    open (K, "$outdir/cf/$i.optics.radmean");
    while (my $line = <K>){
        chomp $line;
        my @line = split (/\t/, $line);
        push (@otpr, $line[0]);
        push (@otnr, $line[1]);
        push (@ojaccard, $line[2]);
        push (@orandindex, $line[3]);
    }
    close (K);
}

my $ktpr = avg (\@ktpr);
my $ktnr = avg (\@ktnr);
my $kjaccard = avg (\@kjaccard);
my $krandindex = avg (\@krandindex);
open (KCF, ">$outdir/cf.clusters.kmeans.radmean");
print KCF "$ktpr\t$ktnr\t$kjaccard\t$krandindex\n";
close (KCF);

my $otpr = avg (\@otpr);
my $otnr = avg (\@otnr);
my $ojaccard = avg (\@ojaccard);
my $orandindex = avg (\@orandindex);
open (OCF, ">$outdir/cf.clusters.optics.radmean");
print OCF "$otpr\t$otnr\t$ojaccard\t$orandindex\n";
close (OCF);

`paste $outdir/cf.clusters.kmeans.radmean $outdir/cf.clusters.optics.radmean $outdir/mds.clusters.kmeans.radmean $outdir/mds.clusters.optics.radmean > $outdir/cf.mds.final.stats`;

### subs ###

sub avg {
    my $array = shift;
    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@$array);
    my $mean = $statobj->mean();
    return ($mean);
}
