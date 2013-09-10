#!/usr/bin/perl -w

=head1 NAME

treeSupports_parse.pl

=head1 SYNOPSIS

treeSupports_parse.pl

Options:

--outdir is your project directory
--config is your configfile
--tree is your treefile
--matrix is your alignfile
--window       is your window size if doing                                                             
    sildeRule (optional)                                                                        
--motion       is the number of residues the                                                           
    window moves per calculation (5' -> 3')                                                           
    if doing slideRule (optional)
--start        is the start position of slideRule                                                    
    (optional)                                                                                           
--end          is the end position of slideRule                                                           
    (optional)                             

Both the matrix and tree must be in 
nexus format.

Dependencies:

PAUP must be installed and in $PATH
Requires the bioperl libs. 
Requires TreeSupports.pm

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

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;

my ($help, $treefile, $matrixfile, $configfile, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_config ($configfile);
my $root = $supportobj->get_root;

my $config = $supportobj->get_config;
my %tests;
my @tests = split (/\s/, $config->{'TESTS'});
foreach my $test (@tests){
    $tests{$test} = 1;
}


$supportobj->load_aln    ($matrixfile);
$supportobj->load_tree   ($treefile, $root);

# override partitions in the nxs file if slideRule chosen                                       
if ( ($window) and ($motion) and ($start) and ($end) ){
    $supportobj->generate_partitions ($window, $motion, $start, $end);
}

# get the charsets and lineages                                                                  
my $charset  = $supportobj->get_charsets;
my $lineages = $supportobj->get_lineages;

# parse the ndi data
opendir (NDI, "$outdir/ndi/logs");
my @ndilogs = grep (/^.+\..+$/, readdir(NDI));
closedir (NDI);

my $ndimaster = {};
foreach my $ndilog (@ndilogs){
    my ($node, $direction, $stuff) = split (/\./, $ndilog, 3);
    my $ndidata = $supportobj->parse_ndi ("$outdir/ndi/logs/$ndilog");
    ($ndimaster->{$node}->{$direction} = $ndidata) if (%$ndidata);
}

# parse the bs data
opendir (BS, "$outdir/bs/logs");
my @bslogs = grep (/^.+\..+$/, readdir(BS));
closedir (BS);

my $bsmaster = {};
foreach my $bslog (@bslogs){
    my ($charset, $stuff) = split (/\./, $bslog);
    my $bsdata = $supportobj->parse_bs ("$outdir/bs/logs/$bslog");
    ($bsmaster->{$charset} = $bsdata) if (%$bsdata);
}

# parse the pbs data
opendir (PBS, "$outdir/pbs/logs");
my @pbslogs = grep (/^.+\..+$/, readdir(PBS));
closedir (PBS);

my $pbsmaster = {};
foreach my $pbslog (@pbslogs){
    my ($node, $direction, $stuff) = split (/\./, $pbslog, 3);
    my $pbsdata = $supportobj->parse_pbs ("$outdir/pbs/logs/$pbslog");
    ($pbsmaster->{$node}->{$direction} = $pbsdata) if (%$pbsdata);
}

# tabulate output                                                                                        
open (SUM, ">$outdir/summary");
print SUM "Part\tCoords\tNode\tTopo\t";
#print SUM "BSCh\tBSwT\tBSwoT\t";
print SUM "BSwL\tBSwoL\tBS\t";
print SUM "PBSunL\tPBSwL\tPBSwoL\t";
print SUM "PBS\tPBSL\tPHBS\tPHBSL\t";
if (exists ($tests{'ndi'})){
    print SUM "NDIwLall\tNDIwoLall\tNDIwLrm\tNDIwoLrm\tNDI\n";
}
else {
    print SUM "\n";
}

open (BTREE, ">$outdir/bs.trees");
foreach my $char (sort keys %$charset){
    
    # tablulate all the BS tree topologies used
    my $const = {};
    my @inf;
    open (B, "$outdir/bs/cmds/$char.bs_commands.nex");
    while (my $line = <B>){
	chomp $line;
	
	if ($line =~m/^constraints/i){
	    my ($con, $tax) = split (/\=/, $line);
	    $con =~s/constraints//i;
	    $con =~s/\s//g;
	    $tax =~s/\;//;
	    $tax =~s/[\(\)]//g;

	    my @tax = split (/\,/, $tax);
	    $const->{$con} = [@tax];
	}
	
	if ($line =~m/^taxset/i){
	    my ($scrim, $set) = split (/\=/, $line);
	    $set =~s/\;//;
	    @inf = split (/\s/, $set);
	}
    }
    close (B);
    
    foreach my $c (sort keys %$const){
	open (BTREE, ">$outdir/bs/trees/$char-$c.tre");
	print BTREE "$char-$c=((";
	
	my %cache;
	foreach my $tax (@{$const->{$c}}){
	    $cache{$tax} = 1;
	}
	
	my @in;
	my @out;
	foreach my $inf (@inf){
	    if (exists ($cache{$inf})){
		push (@in, $inf);
	    }
	    else {
		push (@out, $inf);
	    }
	}
	
	my $in  = join (",", @in);
	my $out = join (",", @out);
	print BTREE "$in),$out)\n";
	close (BTREE);
    }

    foreach my $constraint (sort {$a<=>$b} keys %$lineages){

	# get the topology test result
	open (TP, "$outdir/bs/prunes/$char.$constraint.prun");
	my @topo;
	while (my $line = <TP>){
	    chomp $line;
	    my @splits = split (/\t/, $line);
	    push (@topo, $splits[0]);
	}
	
	my $topo = join (";", @topo);
	close (TP);
	
	my $signal = 0;
	my $string;
	my $bs;
	my $pbs;
	my $pbsl;
	my $phbs;
	my $phbsl;
	my $ndi;
        if ($bsmaster->{$char}->{$constraint}){
	    $signal++;
	    $string .= $char . "\t" . $charset->{$char} . "\t" . $constraint . "\t" . $topo . "\t";

	    # bs values
	    $string .= $bsmaster->{$char}->{$constraint}->{'with'}[2] . "\t";
	    $string .= $bsmaster->{$char}->{$constraint}->{'without'}[2] . "\t";
            $bs = $bsmaster->{$char}->{$constraint}->{'without'}[2] -
                $bsmaster->{$char}->{$constraint}->{'with'}[2];
	    $string .= $bs . "\t";

	    # pbs values
	    $string .= $pbsmaster->{'unconstrained'}->{'none'}->{$char} . "\t";
	    $string .= $pbsmaster->{$constraint}->{'with'}->{$char} . "\t";
	    $string .= $pbsmaster->{$constraint}->{'without'}->{$char} . "\t";
            $pbs = $pbsmaster->{$constraint}->{'without'}->{$char} -
                $pbsmaster->{'unconstrained'}->{'none'}->{$char};
	    $pbsl = $pbsmaster->{$constraint}->{'without'}->{$char} -
		$pbsmaster->{$constraint}->{'with'}->{$char};
	    $string .= $pbs . "\t" . $pbsl . "\t";
            $phbs = $pbs - $bs;
	    $phbsl = $pbsl - $bs;
	    $string .= $phbs . "\t" . $phbsl . "\t";
	    
	    # ndi values
	    if (exists ($tests{'ndi'})){
		$string .= $ndimaster->{$constraint}->{'with'}->{'all'} . "\t";
		$string .= $ndimaster->{$constraint}->{'without'}->{'all'} . "\t";
		$string .= $ndimaster->{$constraint}->{'with'}->{$char} . "\t";
		$string .= $ndimaster->{$constraint}->{'without'}->{$char} . "\t";
		$ndi = ($ndimaster->{$constraint}->{'without'}->{'all'} - 
			   $ndimaster->{$constraint}->{'with'}->{'all'}) -
			   ($ndimaster->{$constraint}->{'without'}->{$char} - 
			    $ndimaster->{$constraint}->{'with'}->{$char});
		$string .= $ndi . "\n";
	    }
	    else {
		(chop $string) if ($string);
		$string .= "\n";
	    }
        }
        else{
            if ($pbsmaster->{$constraint}){
		$string .= $char . "\t" . $charset->{$char} . "\t" . $constraint . "\t" . $topo . "\t";
		$string .= "-" . "\t" . "-" . "\t" . "0" . "\t";

		# pbs vals
		$string .= $pbsmaster->{'unconstrained'}->{'none'}->{$char} . "\t";
		$string .= $pbsmaster->{$constraint}->{'with'}->{$char} . "\t";
		$string .= $pbsmaster->{$constraint}->{'without'}->{$char} . "\t";

		$pbs = $pbsmaster->{$constraint}->{'without'}->{$char} -
		    $pbsmaster->{'unconstrained'}->{'none'}->{$char};
		$pbsl = $pbsmaster->{$constraint}->{'without'}->{$char} -
		    $pbsmaster->{$constraint}->{'with'}->{$char};
		$string .= $pbs . "\t" . $pbsl . "\t";
		$phbs = $pbs - 0;
		$phbsl = $pbsl - 0;
		$string .= $phbs . "\t" . $phbsl . "\t";
		$signal++;
	    }
	    
	    # ndi values
	    if (exists ($tests{'ndi'})){
		if ($ndimaster->{$constraint}){
		    $string .= $ndimaster->{$constraint}->{'with'}->{'all'} . "\t";
		    $string .= $ndimaster->{$constraint}->{'without'}->{'all'} . "\t";
		    $string .= $ndimaster->{$constraint}->{'with'}->{$char} . "\t";
		    $string .= $ndimaster->{$constraint}->{'without'}->{$char} . "\t";
		    $ndi = ($ndimaster->{$constraint}->{'without'}->{'all'} -
			       $ndimaster->{$constraint}->{'with'}->{'all'}) -
			       ($ndimaster->{$constraint}->{'without'}->{$char} -
				$ndimaster->{$constraint}->{'with'}->{$char});
		    $string .= $ndi . "\n";
		    $signal++;
		}
	    }
	    else {
		(chop $string) if ($string);
		$string .= "\n";
	    }
	}
	if ($signal == 0){
	    if (exists ($tests{'ndi'})){
		print SUM "$char\t$charset->{$char}\t$constraint\t$topo\t-\t-\t0\t-\t-\t-\t0\t0\t0\t0\t-\t-\t-\t-\t0\n";
	    }
	    else{
		print SUM "$char\t$charset->{$char}\t$constraint\t$topo\t-\t-\t0\t-\t-\t-\t0\t0\t0\t0\n";
	    }
	}
	else {
	    print SUM "$string";
	}
	
    }
}
close (SUM);
close (BTREE);

