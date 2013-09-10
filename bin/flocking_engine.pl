#!/usr/bin/perl -w

=head1 NAME

flocking_engine.pl

=head1 SYNOPSIS

  flocking_engine.pl -- 
    Implements Reynolds' flocking algorithm

Options:

 --help        Show brief help and exit
 --fasta       Is your fasta file of sequences (boids)
 --config      Is the configuration for your simulation
 --outdir      Is your output dir
 --procs       Is the number of forks to run

=head1 DESCRIPTION

Given the configuration, simulate flocking behavior and 
derive clusters.

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2012 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SeqIO;

my ($help, $fasta, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# count the number of boids we have
=head
my $boids = {};
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    $boids->{$id} = $seq;
}
#my $boidcount = keys %$boids;
=cut

my $boidcount = 3;

# set the initial random positions
my $boidtracker = {};
for (my $nboid = 0; $nboid < $boidcount; $nboid++){
    $boidtracker->{$nboid}->{'xcoord'} = rand($conf->{'DIMMENSIONS'});
    $boidtracker->{$nboid}->{'ycoord'} = rand($conf->{'DIMMENSIONS'}); 
    $boidtracker->{$nboid}->{'course'} = rand(1);
    $boidtracker->{$nboid}->{'velocity'} = rand(1);
}

# iterate the boids and simulate the motion
for (my $it = 0; $it < $conf->{'ITERATIONS'}; $it++){
    for (my $nboid = 0; $nboid < $boidcount; $nboid++){
	separation ($boidtracker, $nboid, $boidcount, $conf);
	cohesion($Nboid);
	alignment($Nboid);
    }
   
    for($Nboid=0; $Nboid < $Number_of_boids; $Nboid++)
    {
	$locationX[$Nboid] =  $locationX[$Nboid] + cos($course[$Nboid]) * $velocity[$Nboid] * $dTime;
	$locationY[$Nboid] =  $locationY[$Nboid]    + sin($course[$Nboid]) * $velocity[$Nboid] * $dTime;
     

	print $locationX[$Nboid] ;
	print " ===";
	print $locationY[$Nboid] ;
	print " ===";
	print $velocity[$Nboid] ;
	print " ===";
	print $course[$Nboid] ;

	print "\n";

    }

    print "\n******************************************\n"; 
}
=cut




#####SUBS#####


sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
	$line =~s/\s//g;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }

    return (\%config);
    close (F);
}

=head
sub dist{
    my ($b1, $b2) = shift;
    return abs(sqrt(($b1->{'xcoord'} - $b2->{'xcoord'})**2 + ($b1->{'ycoord'} - $b2->{'ycoord'})**2));
}


sub cohesion ($)
{

    my ($bcurr)=shift;
    $X=0;
    $Y=0;
    $size=0;
    for($boid=0; $boid < $Number_of_boids; $boid++)
    {
     if ($bcurr != $boid)
     {
        if (dist($bcurr, $boid) < $max_distance)
        {
            $X=$X+$locationX[$boid];
            $Y=$Y+$locationY[$boid];
            $size=$size+1;

        }
    }
 }






      if ($size != 0 )
      {
          $X=$X/($size);
          $Y=$Y/($size);

          $locationX[$bcurr] = $locationX[$bcurr] + ($X-$locationX[$bcurr])*$alpha;
          $locationY[$bcurr] = $locationY[$bcurr] + ($Y-$locationY[$bcurr])*$alpha;
      }
}


sub separation{
    my $boidtracker = shift;
    my $bcurr       = shift;
    my $boidcount   = shift;
    my $conf        = shift;

    my $X    = 0;
    my $Y    = 0;
    my $influencers = 0;

    for (my $boid = 0; $boid < $boidcount; $boid++){
	if ($boid != $bcurr){
	    if (dist ($boidtracker->{$bcurr}, $boidtracker->{$boid}) < $conf->{'RADIUS'}){
		$X = $X + $boidtracker->{$bcurr}->{'xcoord'} - $boidtracker->{$boid}->{'xcoord'};
		$Y = $Y + $boidtracker->{$bcurr}->{'ycoord'} - $boidtracker->{$boid}->{'ycoord'};
		$influencers++;
	    }
	}
	
    }
    if ($influencers > 0){
	$locationX[$bcurr] = $locationX[$bcurr] + $alpha1*$X/$size;
          $locationY[$bcurr] = $locationY[$bcurr] + $alpha1*$Y/$size;
      }

}

sub alignment ($)
{
    my ($bcurr) = shift;

    $dCourse=0;
    $dvelocity=0;

    $size=0;
    for($boid=0; $boid < $Number_of_boids; $boid++)
{
     if ($bcurr != $boid)
     {
        if (dist($bcurr, $boid) < $max_distance)
        {
            $dCourse  = $dCourse + $course[$boid] - $course[$bcurr];
            $dvelocity  = $dvelocity + $velocity[$boid] - $velocity[$bcurr];
            $size=$size+1;
        }

    }
}



if ($size != 0)
{
    $course[$bcurr]=$course[$bcurr] + $dCourse / $size;
    $velocity[$bcurr]=$velocity[$bcurr]+$dvelocity /  $size;
}

}
=cut
