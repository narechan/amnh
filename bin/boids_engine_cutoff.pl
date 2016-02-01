#!/usr/bin/perl
=head1 NAME                                                                       
                                                                                  
boids_engine.pl
                                                                                  
=head1 SYNOPSIS                                                                   
                                                                                  
  boids_engine.pl --                                                           
    Implements Reynolds' flocking algorithm for the evolutionary 
    analysis of gene partitions.
                                                                                  
Options:                                                                          
                                                                                  
 --help        Show brief help and exit                                           
 --matrix      Is your nexus file of sequences (boids are partitions)
 --config      Is the configuration for your flocking simulation                           
 --outdir      Is your output dir                                                 
 --procs       Is the number of forks to run                                      
 --list        Is a list of precomputed ILD stats (optional)                                                  
 --highlights  Is a list file of partitions you want to highlight in the images (optional)
 --threshold   Is the ILD value threshold beyond which only the 
                  separation parameter is in effect. If threshold has a comma, 
                   then repel anything below or above the threshold values input
 --type        Is the type of ILD value being compared: raw or pval
                  note raw only currently works for precomputed LD raw stats
 --chars       Is a list of charsets (either this or matrix must be invoked)

=head1 DESCRIPTION                                                                
                                                                                  
Given the configuration, simulate flocking behavior and                           
derive clusters.                                                                  

Note: either matrix or chars must be specified

Note: the threshold parameter can be two comma separated values if attracting 
    only in the interval.

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


use strict;
use warnings;
use POSIX;
use List::MoreUtils qw(part);
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SeqIO;
use Ild;

my ($help, $matrix, $config, $outdir, $procs, $ildlist, $threshold, $highlights, $type, $charfile);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrix,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    'l|list=s'	      => \$ildlist,
    't|threshold=s'   => \$threshold,
    'h|highlights=s'  => \$highlights,
    'x|type=s'        => \$type,
    'y|chars=s'       => \$charfile,
	   ) or pod2usage;

pod2usage if $help;

#for my $option ($matrix, $outdir, $config, $procs, $threshold){
#    (warn ("Missing a required option\n") and pod2usage)
#        unless ($option);
#}

`mkdir -p $outdir/images`;
`mkdir -p $outdir/logs`;
`mkdir -p $outdir/flocks`;
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
my @boidList;
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

    $boidList[$nboid] = [ $xp,
			  $yp,
			  $xv,
			  $yv
			  ];
    print STDERR "$nboid\t$xp\t$yp\t$xv\t$yv\n";
}

# begin iterations (simulation frames)
my @frames;
open (LOG, ">$outdir/pos_vel.log");
open (CLOG, ">$outdir/flocks.log");
my $icounter = 0;
my $icachecounter = 0;
my $globalildlog = {};
my $localtime;
for (my $it = 0; $it < $conf->{'ITERATIONS'}; $it++){
    $localtime = localtime;
    print STDERR "$it start $localtime\n";

    # iteration specific vars
    my $attractboids = {};
    my $repelboids   = {};
    my $pjobs = {};
    my $pjobscached = {};

    # set up the grid for spatial hash lookups
    # take cell size to be the radius size for now
    # and width should be a whole number (Dimmensons divisible by radius)
    my $boidGrid  = {};
    my $boidIndex = {};
    my $width = $conf->{'DIMMENSIONS'} / $conf->{'RADIUS'};
    for (my $b = 0; $b < $boidcount; $b++){

	my $cell  = (floor ($boidList[$b]->[0] / $conf->{'RADIUS'})) +
	            (floor ($boidList[$b]->[1] / $conf->{'RADIUS'}) * $width);

	push @{$boidGrid->{$cell}}, $b;
	$boidIndex->{$b} = $cell;
    }
    
    # this loop calculates potential nabes for a given boid 
    # using the spatial hash and proceeds to ILDs
    open (IL, ">$outdir/logs/$it.ildaccess");
    for (my $counter = 0; $counter < $boidcount; $counter++){
	my @boidCells;

	# find all the cells around the cell of interest
	# cell of interest given by middle value in both loops
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
	my @boidNabe;
	foreach my $bc (@boidCells){
	    if ($boidGrid->{$bc}){
		push (@boidNabe, @{$boidGrid->{$bc}});
	    }
	    else{
		next;
	    }
	}
	
	# discount self and do full dist calcs only on nabe boids
	foreach my $counter2 (@boidNabe){
	    if($counter2 != $counter){
		if (distance (\@boidList, $counter2, $counter) < $conf->{'RADIUS'}){
		    
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

		    unless ((exists ($ilds->{$counter}->{$counter2})) or (exists ($ilds->{$counter2}->{$counter}))){
			$icounter++;
			$ildobj->generate_nxs_pairwise ($partitions->{$boidmap->{$counter}}, 
							$partitions->{$boidmap->{$counter2}}, 
							$lengths->{$boidmap->{$counter}}, 
							$lengths->{$boidmap->{$counter2}},
							$boidmap->{$counter},
							$boidmap->{$counter2},
							"$outdir/ild/nxs");
			$ildobj->generate_ild_part ($icounter,
						    $boidmap->{$counter},
						    $boidmap->{$counter2},
						    "$outdir/ild/cmds",
						    "$outdir/ild/logs",
						    "$outdir/ild/nxs/$party.nxs");
			$pjobs->{$icounter} = $counterstring;
			
			# cache the ILDs
			$ilds->{$counter}->{$counter2} = 1;                                               
			$ilds->{$counter2}->{$counter} = 1;
		    }
		    else {

			# record that the job has been seen but that 
			# it's data must be retrieved since in nabe
			$icachecounter++;
			$pjobscached->{$icachecounter} = $counterstring;
		    }
		}
		else {
		    next;
		}
	    }
	    else {
		next;
	    }
	}
    }
    close (IL);
    
    # run all the ILDs for this iteration
    my $pm = Parallel::ForkManager->new($procs);
    foreach my $job (keys %$pjobs){
	$pm->start and next;
	$ildobj->run_paup ("$outdir/ild/cmds/$job.nex");
	$pm->finish;
    }
    $pm->wait_all_children;
    
    # parse the ILDs, update the cache with pvals,
    # and find attractors and repellers
    foreach my $job (keys %$pjobs){
	my ($pval, $slen, $ildchars) = $ildobj->parse_ild ("$outdir/ild/logs/$job.log");
	my $string = $pjobs->{$job};
	my ($c1, $c2) = split (/\:/, $string);
	$ilds->{$c1}->{$c2} = $pval;
	$ilds->{$c2}->{$c1} = $pval;

	if ($type eq "pval"){
	    if ($threshold =~m/\,/){
		my ($a, $b) = split (/\,/, $threshold);
		if (($ilds->{$c1}->{$c2} >= $a) and ($ilds->{$c1}->{$c2} <= $b)){
		    $attractboids->{$c1}->{$c2} = 1;
		}
		else{
		    $repelboids->{$c1}->{$c2} = 1;
		}
	    }
	    else{
		if ($ilds->{$c1}->{$c2} > $threshold){
		    $attractboids->{$c1}->{$c2} = 1;
		}
		else {
		    $repelboids->{$c1}->{$c2} = 1;
		}
	    }
	}
	elsif ($type eq "raw"){
	    if ($threshold =~m/\,/){
		my ($a, $b) = split (/\,/, $threshold);
                if (($ilds->{$c1}->{$c2} >= $a) and ($ilds->{$c1}->{$c2} <= $b)){
                    $attractboids->{$c1}->{$c2} = 1;
                }
                else{
                    $repelboids->{$c1}->{$c2} = 1;
                }
	    }
	    else{
		if ($ilds->{$c1}->{$c2} <= $threshold){ # alteration when not using sig; using raw ilds
		    $attractboids->{$c1}->{$c2} = 1;
		}
		else {
		    $repelboids->{$c1}->{$c2} = 1;
		}
	    }
	}
	else {
	    print STDERR "Unrecognized ILD data\n";
	    die;
	}
    }

    # find attractors and repellers among the cached pairs
    foreach my $cachedjob (keys %$pjobscached){
	my $string = $pjobscached->{$cachedjob};
        my ($c1, $c2) = split (/\:/, $string);
	
	my $tild;
	if ($ilds->{$c1}->{$c2}){
	    $tild = $ilds->{$c1}->{$c2};
	}
	else {
	    $tild = $ilds->{$c2}->{$c1};
	}

	if ($type eq "pval"){
	    if ($threshold =~m/\,/){
                my ($a, $b) = split (/\,/, $threshold);
                if (($tild >= $a) and ($tild <= $b)){
                    $attractboids->{$c1}->{$c2} = 1;
                }
                else{
                    $repelboids->{$c1}->{$c2} = 1;
                }
	    }
	    else{
		if ($tild > $threshold){
		    $attractboids->{$c1}->{$c2} = 1;
		}
		else {
		    $repelboids->{$c1}->{$c2} = 1;
		}
	    }
	}
	elsif ($type eq "raw"){
	    if ($threshold =~m/\,/){
		my ($a, $b) = split (/\,/, $threshold);
		if (($tild >= $a) and ($tild <= $b)){
                    $attractboids->{$c1}->{$c2} = 1;
                }
                else{
                    $repelboids->{$c1}->{$c2} = 1;
                }
	    }
	    else{
		if ($tild <= $threshold){ # alteration when not using sig; using raw ilds 
		    $attractboids->{$c1}->{$c2} = 1;
		}
		else {
		    $repelboids->{$c1}->{$c2} = 1;
		}
	    }
	}
	else {
	    print STDERR "Unrecognized ILD data\n";
	    die;
	}
    }
    
    # run through the boids again and update their velocities
    # note that we do not need to recalculate distances
    open (LOGIT, ">$outdir/logs/$it.log");
    for (my $counter = 0; $counter < $boidcount; $counter++){
	
	# three velocity vectors for the three rules                                                         
        my @v1 = (0.0, 0.0);
        my @v2 = (0.0, 0.0);
        my @v3 = (0.0, 0.0);
	
	# call each rule for boid behavior only if the boid has attract/repulse nabes
	# cohesion and alignment only operate on nonsig ILD comparisons
	# separation acts on all but doubles the separation factor for repelled boids
	if ($attractboids->{$counter}){
	    @v1 = cohesion   (\@boidList, 
			      $counter, 
			      $boidcount, 
			      $conf, 
			      $attractboids->{$counter});
	    @v3 = alignment  (\@boidList,
			      $counter,
			      $boidcount,
			      $conf,
			      $attractboids->{$counter});
	}
	if ($attractboids->{$counter} or $repelboids->{$counter}){
	    @v2 = separation (\@boidList, 
			      $counter, 
			      $boidcount, 
			      $conf, 
			      $attractboids->{$counter}, 
			      $repelboids->{$counter});
	}
	
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
	print LOG "$it\t$counter\t";
	print LOG "$boidList[$counter]->[0]\t$boidList[$counter]->[1]\t";
	print LOG "$boidList[$counter]->[2]\t$boidList[$counter]->[3]\n";
	print LOGIT "$boidmap->{$counter}\t$boidList[$counter]->[0]\t$boidList[$counter]->[1]\n";
	
    }
    close (LOGIT);
    
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

    my @gilds = keys %$globalildlog;
    my $gilds = @gilds;
    
    $localtime = localtime;
    print STDERR "ILDTEST\t$it\t$gilds\t$icounter\t$icachecounter\n";
    print STDERR "$it end $localtime\n";
}
close (LOG);
close (CLOG);

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

    my @tmpV = (0.0, 0.0);
    
    # loop thru boids and sum positions
    my $inradcount = 0;
    foreach my $ctr2 (keys %$boidstocheck){
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
    
    # calculate and return 'center of mass' for boids
    # decreasing cohesion factor increase the cohesion
    $tmpV[0] = ($tmpV[0] - $boidList->[$ctr1]->[0]) / $config->{'COHESION_FACTOR'};
    $tmpV[1] = ($tmpV[1] - $boidList->[$ctr1]->[1]) / $config->{'COHESION_FACTOR'};
    
    # return the updated vector
    return @tmpV;
}

# separation: maintain a minimum small distance between boids
sub separation{
    my $boidList  = shift;
    my $ctr1      = shift;
    my $boidcount = shift;
    my $config    = shift;
    my $attractboids = shift;
    my $repelboids = shift;

    my @tmpV = (0.0, 0.0);
    
    # if boids attract, repel only those within separation factor
    foreach my $ctr2 (keys %$attractboids){
	if (distance ($boidList, $ctr2, $ctr1) < $config->{'SEPARATION_FACTOR'}){
	    $tmpV[0] = $tmpV[0] - ($boidList->[$ctr2]->[0] - $boidList->[$ctr1]->[0]);
	    $tmpV[1] = $tmpV[1] - ($boidList->[$ctr2]->[1] - $boidList->[$ctr1]->[1]);
	}
    }
    
    # if boids repel, double the separation distance
    foreach my $ctr3 (keys %$repelboids){
        if (distance ($boidList, $ctr3, $ctr1) < $config->{'REPEL_FACTOR'}){
            $tmpV[0] = $tmpV[0] - ($boidList->[$ctr3]->[0] - $boidList->[$ctr1]->[0]);
            $tmpV[1] = $tmpV[1] - ($boidList->[$ctr3]->[1] - $boidList->[$ctr1]->[1]);
        }
    }
    
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

    my @pvJ = (0.0, 0.0);

    # loop through boids and sum their velocities
    my $inradcount = 0;
    foreach my $ctr2 (keys %$boidstocheck){
	$pvJ[0] = $pvJ[0] + $boidList->[$ctr2]->[2];
	$pvJ[1] = $pvJ[1] + $boidList->[$ctr2]->[3];
	$inradcount++;
    }
    
    # avg boid velocities
    if ($inradcount > 0){
	$pvJ[0] /= $inradcount;
	$pvJ[1] /= $inradcount;
    }
    
    # calculate and return 'percieved velocity'
    # decreasing the alignment factor increases the alignment effect
    my @tmpV = (0.0, 0.0);
    $tmpV[0] = ($pvJ[0] - $boidList->[$ctr1]->[2]) / $config->{'ALIGNMENT_FACTOR'};
    $tmpV[1] = ($pvJ[1] - $boidList->[$ctr1]->[3]) / $config->{'ALIGNMENT_FACTOR'};
  
    # return updated velocity vector
    return @tmpV;
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
