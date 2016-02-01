#!/usr/bin/perl
=head1 NAME                                                                       
                                                                                  
boids_engine_continuous_forked.pl
                                                                                  
=head1 SYNOPSIS                                                                   
                                                                                  
  boids_engine.pl --                                                           
    Implements Reynolds' flocking algorithm for the evolutionary 
    analysis of gene partitions without cutoff. Vectors are weighted
    by distance.
                                                                                  
Options:                                                                          
                                                                                  
 --help        Show brief help and exit                                           
 --matrix      Is your nexus file of sequences (boids are partitions)
 --config      Is the configuration for your flocking simulation                           
 --outdir      Is your output dir                                                 
 --procs       Is the number of forks to run                                      
 --list        Is a list of precomputed ILD stats
 --highlights  Is a list file of partitions you want to highlight in the images (optional)
 --chars       Is a list of charsets (either this or matrix must be invoked)
 --perception  Is the perceptive ability of each boid (heuristic for number of checks within a grid)
               If you want to do all checks, input 'all'; otherwise input the maximum number
 --binlattice  Is set to 1 if you want to do spatial hashing to save computes or 2 otherwise

=head1 DESCRIPTION                                                                
                                                                                  
Given the configuration, simulate flocking behavior and                           
derive clusters without estimating cluster number or supplying a cutoff 
for flocking.                                                                  

Note: either matrix or chars must be specified

Note: reduce parrallelism of imagemagick if running multiple simualtions on 
a single multicore machine:

export MAGICK_THREAD_LIMIT=1
                                                                                  
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

# TODO: 
# add randomness beyond initial position to avoid determinism and assure through search of the space?
# add back in the ILD calcs in forked fashion
# remove bioperl dependencies and input 1 fasta file per boid
# clean up required configuration
# make gnuplot and optics optional or per every X

use strict;
use warnings;
use POSIX;
use List::MoreUtils qw(part);
use List::Util qw(shuffle);
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SeqIO;
use Ild;

my ($help, $matrix, $config, $outdir, $procs, $ildlist, $highlights, $charfile, $perceive, $binlattice);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrix,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    'l|list=s'	      => \$ildlist,
    'h|highlights=s'  => \$highlights,
    'y|chars=s'       => \$charfile,
    's|perceive=s'    => \$perceive,
    'b|binlattice=s'  => \$binlattice,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrix, $outdir, $config, $procs, $perceive, $binlattice){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/images`;
`mkdir -p $outdir/logs`;
`mkdir -p $outdir/flocks`;
`mkdir -p $outdir/posvel`;
`mkdir -p $outdir/ild/cmds`;
`mkdir -p $outdir/ild/logs`;
`mkdir -p $outdir/ild/nxs`;

####MAIN####

# get the config
my $conf = parse_config ($config);

# parse out the boids
my $charsets;
my $partitions;
my $lengths;
my $ildobj;
if ($matrix){
    $ildobj = Ild->new;
    $ildobj->load_aln ($matrix);
    $ildobj->load_config ($config);
    $charsets = $ildobj->get_charsets;
    ($partitions, $lengths) 
	= $ildobj->store_alignment_seq($matrix, $charsets);
}
elsif ($charfile){
    open (C, "$charfile");
    while (my $line = <C>){
	chomp $line;
	$charsets->{$line} = 1;
    }
}

# assign each partition a boid array slot
my $boidcount = 0;
my $boidmap = {};
my $partmap = {};
foreach my $part (sort keys %$charsets){
    $boidcount++;
    my $boidindex = $boidcount - 1;
    $boidmap->{$boidindex} = $part;
    $partmap->{$part} = $boidindex;
}

# if available store the precomputed ILD list
my $ilds = {};
if ($ildlist){
    open (L, "$ildlist");
    while (my $line = <L>){
	chomp $line;
	my ($pair, $pval, $treelen) = split (/\t/, $line);
	my ($b1, $b2) = split (/v/, $pair);
	my $i1 = $partmap->{$b1};
	my $i2 = $partmap->{$b2};
	$ilds->{$i1}->{$i2} = $pval;
	$ilds->{$i2}->{$i1} = $pval;
    }
    close (L);
}

# set the initital random conditions in the datastruc
# arrayofarrays struc contains (x_location, y_location, x_velocity, y_velocity)
open (LOGZERO, ">$outdir/posvel/0.posvel");
for (my $nboid = 0; $nboid < $boidcount; $nboid++){

    # seed the velocities with potentially negative random numbers
    my $xv;
    my $yv;
    my $xvsign = int(rand(2));
    my $xvmag  = rand($conf->{'INIT_VELOCITY'});
    my $yvsign = int(rand(2));
    my $yvmag  = rand($conf->{'INIT_VELOCITY'});
    if ($xvsign == 0){
	$xv = -1 * $xvmag;
    }
    else {
	$xv = $xvmag;
    }
    if ($yvsign == 0){
        $yv = -1 * $yvmag;
    }
    else {
        $yv = $yvmag;
    }
    my $xp = rand($conf->{'DIMMENSIONS'});
    my $yp = rand($conf->{'DIMMENSIONS'});

    print LOGZERO "$nboid\t$xp\t$yp\t$xv\t$yv\n";
}
close (LOGZERO);

# begin iterations (simulation frames)
my @frames;
open (CLOG, ">$outdir/flocks.log");
my $icounter = 0;
my $icachecounter = 0;
my $globalildlog = {};
my $localtime;
for (my $it = 1; $it <= $conf->{'ITERATIONS'}; $it++){
    $localtime = localtime;
    print STDERR "$it start $localtime\n";

    # read in the boid locations and velocities for this frame
    my @boidList;
    my $reader = $it - 1;
    open (BLIST, "$outdir/posvel/$reader.posvel");
    while (my $line = <BLIST>){
	chomp $line;
	my ($nb, $xp, $yp, $xv, $yv) = split (/\t/, $line);
	$boidList[$nb] = [ $xp,
			   $yp,
			   $xv,
			   $yv
	    ];
    }
    close (BLIST);
	
    # set up the grid for spatial hash lookups
    # take cell size to be the radius size for now
    # and width should be a whole number (Dimmensons divisible by radius)
    my $boidGrid  = {};
    my $boidIndex = {};
    my $width;
    if ($binlattice == 1){
	$width = $conf->{'DIMMENSIONS'} / $conf->{'RADIUS'};
	for (my $b = 0; $b < $boidcount; $b++){
	    
	    my $cell  = (floor ($boidList[$b]->[0] / $conf->{'RADIUS'})) +
		(floor ($boidList[$b]->[1] / $conf->{'RADIUS'}) * $width);
	    
	    push @{$boidGrid->{$cell}}, $b;
	    $boidIndex->{$b} = $cell;
	}
    }
    else {
	for (my $b = 0; $b < $boidcount; $b++){
	    $boidIndex->{$b} = 1;
	}
    }

    # this loop calculates potential nabes for a given boid 
    # using the spatial hash and proceeds to ILDs
    open (IL, ">$outdir/logs/$it.ildaccess");
    open (LOGIT, ">$outdir/logs/$it.log");
    open (LOG, ">$outdir/posvel/$it.posvel");
    my $boidCullTot = 0;
    for (my $counter = 0; $counter < $boidcount; $counter++){

	# find all the cells around the cell of interest
	# cell of interest given by middle value in both loops
	# or discontinue the bin-lattice routine
	my @boidNabe;
	if ($binlattice == 1){
	    my @boidCells;
	    foreach my $i ($boidList[$counter]->[0] - $conf->{'RADIUS'}, $boidList[$counter]->[0], $boidList[$counter]->[0] + $conf->{'RADIUS'}){
		next if (($i < 0) or ($i > $conf->{'DIMMENSIONS'}));
		foreach my $j ($boidList[$counter]->[1] - $conf->{'RADIUS'}, $boidList[$counter]->[1], $boidList[$counter]->[1] + $conf->{'RADIUS'}){
		    next if (($j < 0) or ($j > $conf->{'DIMMENSIONS'}));
		    
		    my $ncell  = (floor ($i / $conf->{'RADIUS'})) +
			(floor ($j / $conf->{'RADIUS'}) * $width);
		    
		    push (@boidCells, $ncell);
		}
	    }
	    
	    # dump all the boids from the primary cells and the nabe cells into an array
	    foreach my $bc (@boidCells){
		if ($boidGrid->{$bc}){
		    push (@boidNabe, @{$boidGrid->{$bc}});
		}
		else{
		    next;
		}
	    }
	}
	else{
	    @boidNabe = keys %$boidIndex;
	}

	# cull the array randomly to account for selective perception / awareness
	my @boidCull;
	if ($perceive eq "all"){
	    @boidCull = @boidNabe;
	}
	else{
	    my $boidNabeCnt = @boidNabe;
	    if ($perceive <= $boidNabeCnt){
		@boidCull = (shuffle(@boidNabe))[0..$perceive-1];
	    }
	    else {
		@boidCull = @boidNabe;
	    }
	}

	my $boidCull = @boidCull;
	$boidCullTot += $boidCull;

	# discount self and do full dist calcs only on culled nabe boids
	my $boidstocheck = {};
	foreach my $counter2 (@boidCull){
	    if($counter2 != $counter){
		if (distance (\@boidList, $counter2, $counter) < $conf->{'RADIUS'}){
#		    print STDERR "$counter\t$counter2\n";
		    # do an ILD test unless ILDs were precalculated
		    # or if this one was already encountered
		    my $party = $boidmap->{$counter} . "v" . $boidmap->{$counter2};
		    my $counterstring = $counter . ":" . $counter2;

		    # global ild logging
		    my @party;
		    push (@party, $boidmap->{$counter}, $boidmap->{$counter2});
		    my @sortedparty = sort @party;
		    my $partystring = join "v", @sortedparty;

		    if (exists ($globalildlog->{$partystring})){
			$globalildlog->{$partystring}++;
		    }
		    else {
			print IL "$partystring\n";
			$globalildlog->{$partystring}++;
		    }

		    # record that the job has been seen but that 
		    # it's data must be retrieved since in nabe
		    $icachecounter++;
		    $boidstocheck->{$counter2} = 1;
		}
		else {
		    next;
		}
	    }
	    else {
		next;
	    }
	}

	# three velocity vectors for the three rules                                                         
        my @v1 = (0.0, 0.0);
        my @v2 = (0.0, 0.0);
        my @v3 = (0.0, 0.0);
	
	# call each rule for boid behavior and pass the ILDs to weight behavior
	@v1 = cohesion   (\@boidList, 
			  $counter, 
			  $boidcount, 
			  $conf, 
			  $boidstocheck,
			  $ilds);
	@v3 = alignment  (\@boidList,
			  $counter,
			  $boidcount,
			  $conf,
			  $boidstocheck,
			  $ilds);
	@v2 = separation (\@boidList, 
			  $counter, 
			  $boidcount, 
			  $conf, 
			  $boidstocheck,
			  $ilds);

        # update boid velocity vector                                                                         
	$boidList[$counter]->[2] += $v1[0] + $v2[0] + $v3[0];
	$boidList[$counter]->[3] += $v1[1] + $v2[1] + $v3[1];
	
	# limit the velocities in case too high                                                            
        if (abs($boidList[$counter]->[2]) > $conf->{'VELOCITY_LIMIT'}){
            $boidList[$counter]->[2] =
                (abs($boidList[$counter]->[2]) / $boidList[$counter]->[2]) * $conf->{'VELOCITY_LIMIT'};
        }
        if (abs($boidList[$counter]->[3]) > $conf->{'VELOCITY_LIMIT'}){
            $boidList[$counter]->[3] =
                (abs($boidList[$counter]->[3]) / $boidList[$counter]->[3]) * $conf->{'VELOCITY_LIMIT'};
        }
	
	# get boid location given the new vectors                           
	my $x = $boidList[$counter]->[0] + $boidList[$counter]->[2];
	my $y = $boidList[$counter]->[1] + $boidList[$counter]->[3];
	
        # steer boids off walls
	if ($x > $conf->{'DIMMENSIONS'} - $conf->{'BOUNDARY'}){
	    $boidList[$counter]->[2] = $conf->{'INIT_VELOCITY'} * -1;
	}
	if ($x < $conf->{'BOUNDARY'}){
	    $boidList[$counter]->[2] = $conf->{'INIT_VELOCITY'};
	}
	if ($y > $conf->{'DIMMENSIONS'} - $conf->{'BOUNDARY'}){
            $boidList[$counter]->[3] = $conf->{'INIT_VELOCITY'} * -1;
        }
        if ($y < $conf->{'BOUNDARY'}){
            $boidList[$counter]->[3] = $conf->{'INIT_VELOCITY'};
        }

	# update the locations in the datastruc
	$boidList[$counter]->[0] = $x;
	$boidList[$counter]->[1] = $y;
	    	
	# logs
	print LOG "$counter\t";
	print LOG "$boidList[$counter]->[0]\t$boidList[$counter]->[1]\t";
	print LOG "$boidList[$counter]->[2]\t$boidList[$counter]->[3]\n";
	print LOGIT "$boidmap->{$counter}\t$boidList[$counter]->[0]\t$boidList[$counter]->[1]\n";
	
    }
    close (IL);
    close (LOGIT);
    close (LOG);
=head    
    # print the frame for this iteration and convert to gif
    print STDERR "$it GNUPLOT\n";
    if ($highlights){
	`boids_gnuplot.pl -i $outdir/logs/$it.log -h $highlights -o $outdir/images/$it.gif -d $conf->{'DIMMENSIONS'}`;
    }
    else {
	`boids_gnuplot.pl -i $outdir/logs/$it.log -o $outdir/images/$it.gif -d $conf->{'DIMMENSIONS'}`;
    }
    push (@frames, "$outdir/images/$it.gif");
    
    # find the flocks using optics
    print STDERR "$it OPTICS\n";
    `boids_optics.pl -i $outdir/logs/$it.log -m $conf->{'MINPTS'} -x $conf->{'XI'} -j $conf->{'JARFILE'} -o $outdir/flocks > $outdir/flocks/$it.optics`;
    
    open (OPTICS, "$outdir/flocks/$it.optics");
    my $flocks = {};
    while (my $line = <OPTICS>){
	chomp $line;
	my ($num, $acc) = split (/\t/, $line);
	$flocks->{$num}++;
    }
    close (OPTICS);
    
    my @flocks = keys %$flocks;
    my $flockcount = @flocks;
    print CLOG "$flockcount\n";
=cut
    my @gilds = keys %$globalildlog;
    my $gilds = @gilds;
    
    $localtime = localtime;
    my $boidCullRatio = $boidCullTot / $boidcount;

    # $gilds is the number of uniq ild lookups for the entire simulation
    # $icachecounter is the total number of ild lookups for the simulation
    # $boidCullTot is the number of pairwise comps that need to be sampled in this iteration post heuristics
    # $boidCullRation is the tot divided by the number of boids
    print STDERR "ILDTEST\t$it\t$gilds\t$icounter\t$icachecounter\t$boidCullTot\t$boidCullRatio\n";
    print STDERR "$it end $localtime\n";
}
close (CLOG);
=head
# generate the final animated gif using gifsicle
my $framecount = @frames;
if ($framecount <= 1000){
    my $frames = join " ", @frames;
    `gifsicle -d 10 -D bg $frames > $outdir/final.gif`;
}
else {
    my $i = 0;
    my @part = part { int( $i++ / 1000 ) } @frames;
    
    my $gifcounter = 0;
    my @gifcounter;
    foreach my $part (@part){
	$gifcounter++;
	my $frames = join " ", @$part;
	`gifsicle -d 10 -D bg $frames > $outdir/$gifcounter.final.gif`;
	push (@gifcounter, "$outdir/$gifcounter.final.gif");
    }
    my $finalframes = join " ", @gifcounter;
    `gifsicle -d 10 -D bg $finalframes > $outdir/final.gif`;
}
=cut
####SUBS####

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

# cohesion: center of mass tendency
sub cohesion{
    my $boidList  = shift;
    my $ctr1      = shift;
    my $boidcount = shift;
    my $config    = shift;
    my $boidstocheck = shift;
    my $ilds      = shift;

    my @tmpV = (0.0, 0.0);
    
    # loop thru boids and sum positions
    my $inradcount = 0;
    my @cohweights;
    foreach my $ctr2 (keys %$boidstocheck){
	push (@cohweights, $ilds->{$ctr1}->{$ctr2});
	$tmpV[0] = $tmpV[0] + $boidList->[$ctr2]->[0];
	$tmpV[1] = $tmpV[1] + $boidList->[$ctr2]->[1];
	$inradcount++;
    }

    # divide sum by in-radius boids only if there 
    # are boids within the radius
    if ($inradcount > 0){
	$tmpV[0] /= $inradcount;
	$tmpV[1] /= $inradcount;
    }
    
    # find the cummulative weight of all these cohesions
    # given the ILDs
    my $counter = 0;
    my $sum = 0;
    foreach my $cw (@cohweights){
	$counter++;
#	$cw = ($cw + 0.01) * 100;
	$sum += $cw;
    }
    
    my $cohesion_factor;
    if ($counter > 0){
	$cohesion_factor = $sum / $counter;
    }
    else {
	$cohesion_factor = 1;
	$tmpV[0] = $boidList->[$ctr1]->[0]; #so the cohesion vec reduces to zero
	$tmpV[1] = $boidList->[$ctr1]->[1]; #when there's nothing in the radius
    }
=head
    # choose steering directions randomly
    my $rand1 = int (rand(2));
    my $rand2 = int (rand(2));
    print STDERR "$rand1\t$rand2\n";

    if ($rand1 == 0){
	$tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / (10 ** $cohesion_factor);                      
    }
    else {
	$tmpV[0] = ($tmpV[0] + $boidList->[$ctr1]->[0]) / (10 ** $cohesion_factor);
    }
    if ($rand2 == 0){
        $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / (10 ** $cohesion_factor);
    }
    else {
        $tmpV[1] = ($tmpV[1] + $boidList->[$ctr1]->[1]) / (10 ** $cohesion_factor);
    }
=cut

    # calculate and return 'center of mass' for boids
    # decreasing cohesion factor increase the cohesion
#    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / $config->{'COHESION_FACTOR'};
#    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / $config->{'COHESION_FACTOR'};
#    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / (10 ** $cohesion_factor);
#    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / (10 ** $cohesion_factor);
#    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / (2 ** $cohesion_factor);
#    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / (2 ** $cohesion_factor);
    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) * ((1 - $cohesion_factor) / $config->{'COHESION_FACTOR'});
    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) * ((1 - $cohesion_factor) / $config->{'COHESION_FACTOR'});
#    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / 10;
#    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / 10;

    # return the updated vector
    return @tmpV;
}

# separation: maintain a minimum small distance between boids
sub separation{
    my $boidList  = shift;
    my $ctr1      = shift;
    my $boidcount = shift;
    my $config    = shift;
    my $boidstocheck = shift;
    my $ilds      = shift;

    my @tmpV = (0.0, 0.0);
    
    # repel only those between radius and separation factor a value proprtional to distance
    foreach my $ctr2 (keys %$boidstocheck){
	if ((distance ($boidList, $ctr2, $ctr1) <= $config->{'RADIUS'}) and (distance ($boidList, $ctr2, $ctr1) > $config->{'SEPARATION_DISTANCE'})){
#	    $tmpV[0] = $tmpV[0] - (($boidList->[$ctr2]->[0] - $boidList->[$ctr1]->[0]) * $ilds->{$ctr1}->{$ctr2});
#            $tmpV[1] = $tmpV[1] - (($boidList->[$ctr2]->[1] - $boidList->[$ctr1]->[1]) * $ilds->{$ctr1}->{$ctr2});
#	    $tmpV[0] = $tmpV[0] - (($boidList->[$ctr2]->[0] - $boidList->[$ctr1]->[0]) * (5 * $ilds->{$ctr1}->{$ctr2}));
#            $tmpV[1] = $tmpV[1] - (($boidList->[$ctr2]->[1] - $boidList->[$ctr1]->[1]) * (5 * $ilds->{$ctr1}->{$ctr2}));
	    $tmpV[0] = $tmpV[0] - (($boidList->[$ctr2]->[0] - $boidList->[$ctr1]->[0]) * ($config->{'REPEL_FACTOR'} * $ilds->{$ctr1}->{$ctr2}));
            $tmpV[1] = $tmpV[1] - (($boidList->[$ctr2]->[1] - $boidList->[$ctr1]->[1]) * ($config->{'REPEL_FACTOR'} * $ilds->{$ctr1}->{$ctr2}));


	}
	elsif (distance ($boidList, $ctr2, $ctr1) <= $config->{'SEPARATION_DISTANCE'}){
	    $tmpV[0] = $tmpV[0] - ($boidList->[$ctr2]->[0] - $boidList->[$ctr1]->[0]);
	    $tmpV[1] = $tmpV[1] - ($boidList->[$ctr2]->[1] - $boidList->[$ctr1]->[1]);
	}
	else{
	    next;
	}
    }
    
    # if boids repel, double the separation distance
#    foreach my $ctr3 (keys %$repelboids){
#        if (distance ($boidList, $ctr3, $ctr1) < $config->{'REPEL_FACTOR'}){
#            $tmpV[0] = $tmpV[0] - ($boidList->[$ctr3]->[0] - $boidList->[$ctr1]->[0]);
#            $tmpV[1] = $tmpV[1] - ($boidList->[$ctr3]->[1] - $boidList->[$ctr1]->[1]);
#        }
#    }
    
    # return the updated vector
    return @tmpV;
}

# alignment: match velocities
sub alignment{
    my $boidList  = shift;
    my $ctr1      = shift;
    my $boidcount = shift;
    my $config    = shift;
    my $boidstocheck = shift;
    my $ilds      = shift;

    my @pvJ = (0.0, 0.0);

    # loop through boids and sum their velocities
    my $inradcount = 0;
    my @aliweights;
    foreach my $ctr2 (keys %$boidstocheck){
	push (@aliweights, $ilds->{$ctr1}->{$ctr2});
	$pvJ[0] = $pvJ[0] + $boidList->[$ctr2]->[2];
	$pvJ[1] = $pvJ[1] + $boidList->[$ctr2]->[3];
	$inradcount++;
    }
    
    # avg boid velocities
    if ($inradcount > 0){
	$pvJ[0] /= $inradcount;
	$pvJ[1] /= $inradcount;
    }

    # ILD based cummulative weight of aln factors
    my $counter = 0;
    my $sum = 0;
    foreach my $al (@aliweights){
        $counter++;
#        $al = ($al + 0.01) * 100;
        $sum += $al;
    }

    my $alignment_factor;
    if ($counter > 0){
        $alignment_factor = $sum / $counter;
    }
    else {
        $alignment_factor = 1;
	$pvJ[0] = $boidList->[$ctr1]->[2]; #so the alignment vec reduces to zero                             
        $pvJ[1] = $boidList->[$ctr1]->[3]; #when there's nothing in the radius 
    }

    # calculate and return 'percieved velocity'
    # decreasing the alignment factor increases the alignment effect
#    my @tmpV = (0.0, 0.0);
#    $tmpV[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) / $config->{'ALIGNMENT_FACTOR'};
#    $tmpV[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) / $config->{'ALIGNMENT_FACTOR'};
#    $pvJ[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) / (10 ** $alignment_factor);
#    $pvJ[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) / (10 ** $alignment_factor);
#    $pvJ[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) / (2 ** $alignment_factor);
#    $pvJ[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) / (2 ** $alignment_factor);
    $pvJ[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) * ((1 - $alignment_factor) / $config->{'ALIGNMENT_FACTOR'});
    $pvJ[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) * ((1 - $alignment_factor) / $config->{'ALIGNMENT_FACTOR'});
#    $pvJ[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) / 10;
#    $pvJ[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) / 10;
    
    # return updated velocity vector
    return @pvJ;
}

# returns distance between 2 points
sub distance{
    my $boidList = shift;
    my $ct1      = shift;
    my $ct2      = shift;

    my @distV = (0.0, 0.0);
    
    $distV[0] = $boidList->[$ct1]->[0] - $boidList->[$ct2]->[0];
    $distV[1] = $boidList->[$ct1]->[1] - $boidList->[$ct2]->[1];

    $distV[0] *= $distV[0];
    $distV[1] *= $distV[1];
    
    return sqrt($distV[0] + $distV[1]);
}

# other undeveloped code:
# wrap the coordinates around the edges                                                                    
#       if ($x > $conf->{'DIMMENSIONS'}){                                                                     
#           $x -= $conf->{'DIMMENSIONS'};                                                                    
#       }                                                                                                     
#       if ($x < 0){                                                                                         
#           $x += $conf->{'DIMMENSIONS'};                                                                    
#       }                                                                                                   
#       if ($y > $conf->{'DIMMENSIONS'}){                                                                 
#           $y -= $conf->{'DIMMENSIONS'};                                                                 
#       }                                                                                                 
#       if ($y < 0){                                                                                        
#           $y += $conf->{'DIMMENSIONS'};                                                                   
#       }                                        

#       $boidList[$counter]->[0] += $conf->{'VELOCITY_LIMIT'} * $boidList[$counter]->[2];                  
#       $boidList[$counter]->[1] += $conf->{'VELOCITY_LIMIT'} * $boidList[$counter]->[3];           
