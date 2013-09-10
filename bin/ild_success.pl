#!/usr/bin/perl -w

=head1 NAME

ild_success.pl

=head1 SYNOPSIS

ild_success.pl -- This program does susccessive ilds according to the following rules:
    PUT IN RULES

Options:

--procs        the number of procs you want to use
--config       is the config file
--matrix       is the data matrix
--outdir       is your output dir for the run
--infile       is the set of partitions you want to test
--rounds       is the number of rounds
--type         is whether you want to go after "high" or "low"
    use those quoted strings
--signifincance  is the pval threshold for whether a pair is incongruent

the matrix must be in
nexus format.

Requires the bioperl libs.
Requires Ild.pm

=head1 DESCRIPTION


=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Ild;
use Parallel::ForkManager;

my ($help, $configfile, $matrixfile, $outdir, $infile, $procs, $rounds, $type, $significance);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    'i|infile=s'       => \$infile,
    'p|procs=s'        => \$procs,	   
    'r|rounds=s'       => \$rounds,
    't|type=s'         => \$type,
    's|significance=s' => \$significance	   
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $matrixfile, $outdir, $infile, $procs, $rounds, $type, $significance){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# create dir structure for the all agnst all results                                                    
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;

# instantiate the object and load required data
my $ildobj = Ild->new;
$ildobj->load_config ($configfile);
$ildobj->load_aln    ($matrixfile);

# create all pairwise ILDs from the infile (loads whole matrix for now)
my $file = {};
my $counter = 0;
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    $file->{$line} = 1;
}
close (F);

foreach my $line1 (sort keys %$file){
    foreach my $line2 (sort keys %$file){
	next if ($line1 eq $line2);
	
	$counter++;
	$ildobj->generate_ild_part ($counter,
				    $line1,
				    $line2,
				    "$outdir/cmds",
				    "$outdir/logs",
				    $matrixfile);
    }
    delete $file->{$line1};
}

# run all ILDS from the infile
print STDERR "Running all against all\n";
run ("$outdir/cmds");

# parse the initial all against all run
opendir (R, "$outdir/logs");
my @logs = sort (readdir(R));
shift @logs;
shift @logs;
closedir (R);

open (AOUT, ">$outdir/summary.0");
my $pvals = {};
my $ilds  = {};
my $charall = {};
my $charsig = {};
foreach my $log (@logs){
    my ($pval, $slen, $ildchars) = $ildobj->parse_ild ("$outdir/logs/$log");
    print AOUT "$ildchars\t$pval\t$slen\n";
    $ildchars =~s/v/\_/;

    my ($char1, $char2) = split (/\_/, $ildchars);
    $charall->{$char1}++;
    $charall->{$char2}++;
    ($charsig->{$char1}++) if ($pval <= $significance);
    ($charsig->{$char2}++) if ($pval <= $significance);
    
    push @{$pvals->{$pval}}, $ildchars;
    push @{$ilds->{$ildchars}}, $pval, $slen;
}
close (AOUT);

# determine the 5 genes with the most extreme incongruences
my %extremes;
open (EOUT, ">$outdir/avg_incong.0");
foreach my $char (sort keys %$charall){
    my $fracsig;

    if ($charsig->{$char}){
        $fracsig = $charsig->{$char} / $charall->{$char};
    }
    else {
        $fracsig = 0;
    }

    $extremes{$char} = $fracsig;
    
    if ($charsig->{$char}){
        print EOUT "$char\t$charsig->{$char}\t$charall->{$char}\t$fracsig\n";
    }   
    else{
	print EOUT "$char\t0\t$charall->{$char}\t$fracsig\n";
    }

}
close (EOUT);

my @extremeselections;
my $extremeselections = {};
my $extcounter = 0;

if ($type eq "low"){
    foreach my $ext (sort {$extremes{$b}<=>$extremes{$a}} keys %extremes){
	$extcounter++;
	push (@extremeselections, $ext);
	$extremeselections->{$ext} = 1;
	last if ($extcounter == 5);
    }
}
elsif ($type eq "high"){
    foreach my $ext (sort {$extremes{$a}<=>$extremes{$b}} keys %extremes){
        $extcounter++;
	push (@extremeselections, $ext);
	$extremeselections->{$ext} = 1;
        last if ($extcounter == 5);
    }
}
else {
    print STDERR "Unrecognized type\n";
    die;
}
    
# cycle through the pvals for the all against all 
# and select pairs using the extremes 
my @selections;
my %existing;
foreach my $extsel (@extremeselections){
    my $pvalsarray = sorter ($pvals, $type);

    my $selection_counter = 0;
    foreach my $pval (@$pvalsarray){
	my @pairs = @{$pvals->{$pval}};
	my $pairs = @pairs;

	while (1){
	    last if ($pairs == 0);

	    my $rand = int(rand(@pairs));
	    my $pair = $pairs[$rand];
	    my @components = split (/\_/, $pair);
	
	    splice (@pairs, $rand, 1);
	    $pairs--;
	
	    # do the extreme checks
	    # makes sure $extsel is one of the pairs
	    unless (($extsel eq $components[0]) or 
		    ($extsel eq $components[1])){
		next;
	    }
	    
	    # makes sure that the pair does not contain two extremes
	    if (exists ($extremeselections->{$components[0]}) and
		exists ($extremeselections->{$components[1]})){ 
		next;
	    }

	    my $signal = 0;
	    foreach my $component (@components){
		if (exists ($existing{$component})){
		    $signal++;
		}
		else {
		    next;
		}
	    }
	    
	    if ($signal == 0){
		push (@selections, $pair);
		$selection_counter++;
		foreach my $comp (@components){
		    $existing{$comp} = 1;
		}
		last;
	    }
	    else {
		next;
	    }
	}
    
	last if ($selection_counter == 1);
    }
}

# print the selections and their order from the all-all
open (A, ">$outdir/selections.0");
my $s = 0;
foreach my $sel (@selections){
    $s++;
    print A "$s\t$sel\t$ilds->{$sel}[0]\t$ilds->{$sel}[1]\n";
}

# ilds for all selections against all original individual genes
# as in, ilds with more than one charset per partition (loads whole matrix)
for (my $i = 1; $i <= $rounds; $i++){
    `mkdir -p $outdir/cmds.$i`;
    `mkdir -p $outdir/logs.$i`;

    # selections defined by all-all for first round;
    # defined by concat criterion for next rounds.
    my $counter = 0;
    foreach my $line1 (@selections){
	my @line1 = split (/_/, $line1);
    
	open (F, "$infile");
	while (my $line2 = <F>){
	    chomp $line2;
	    
	    # check to see if there is charset overlap
	    my $signal = 0;
	    foreach my $ele1 (@line1){
		($signal++) if ($ele1 eq $line2);
	    }
	    next if ($signal > 0);
	    
	    $counter++;
	    $ildobj->generate_ild_part ($counter,
					$line1,
					$line2,
					"$outdir/cmds.$i",
					"$outdir/logs.$i",
					$matrixfile);
	}
    }

    # run all ILDS from the next iteration                                                       
    print STDERR "Running round $i\n";
    run ("$outdir/cmds.$i");

    # parse the next iteration's runs                                                                 
    opendir (R, "$outdir/logs.$i");
    my @logs = sort (readdir(R));
    shift @logs;
    shift @logs;
    closedir (R);

    open (IOUT, ">$outdir/summary.$i");
    my $pvals = {};
    my $ilds  = {};
    my %pairs;
    my %existing;
    foreach my $log (@logs){
	my ($pval, $slen, $ildchars) = $ildobj->parse_ild ("$outdir/logs.$i/$log");
	my ($pair, $gene) = split (/v/, $ildchars);
	my @p = split (/\_/, $pair);
	
	foreach my $p (@p){
	    $existing{$p} = 1;
	}
	$pairs{$pair} = 1;

	push @{$pvals->{$pair}->{$pval}}, $gene;
	push @{$ilds->{$pair}->{$gene}}, $pval, $slen;

	print IOUT "$ildchars\t$pval\t$slen\n";
    }
    close (IOUT);

    # concatenate to existing pairs using the following rules:
    # TK
    @selections = ();
    my @pairs = keys %pairs;
    my $pairs = @pairs;
    until ($pairs == 0){
	my $rand = int(rand(@pairs));
	my $pair = $pairs[$rand];
	
	splice (@pairs, $rand, 1);
	$pairs--;
	
	my $pvalsarray = sorter ($pvals->{$pair}, $type);
	foreach my $pval (@$pvalsarray){
	    my @genes = @{$pvals->{$pair}->{$pval}};
	    my $signal = 0;
	    my $genes = @genes;
	    
	    until ($signal == 1){
		last if ($genes == 0);
		
		my $randgene = int(rand(@genes));
		my $gene = $genes[$randgene];

		splice (@genes, $randgene, 1);
		$genes--;
		
		if (exists ($existing{$gene})){
		    next;
		}
		else {
		    my $selection = $pair . "_" . $gene;
		    push (@selections, $selection);
		    $existing{$gene} = 1;
		    $signal = 1;
		}
	    }
	    last if ($signal == 1);
	}
    }
    
    # print the selections and their order and associated data
    open (I, ">$outdir/selections.$i");
    my $s = 0;
    foreach my $sel (@selections){
	$s++;
	
	my @sel = split (/\_/, $sel);
	my $g = pop @sel;
	my $p = join "_", @sel;
	
	# determine the fraction of that pair's tests that were significant
	my $sigcounter = 0;
	my $totcounter = 0;
	foreach my $gs (sort keys %{$ilds->{$p}}){
	    $totcounter++;
	    ($sigcounter++) if ($ilds->{$p}->{$gs}[0] <= $significance);
	}

	my $sigfrac = $sigcounter / $totcounter;

	print I "$s\t$sel\t$ilds->{$p}->{$g}[0]\t$ilds->{$p}->{$g}[1]\t$sigfrac\n";
    }
    close (I);
}

###SUBS###

sub sorter {
    my $pvals = shift;
    my $type = shift;
    
    my @pvals;
    if ($type eq "low"){
	@pvals = sort {$a <=> $b} keys %$pvals;
    }
    elsif ($type eq "high"){
	@pvals = sort {$b <=> $a} keys %$pvals;
    }
    else {
	print STDERR "Unknown type\n";
	die;
    }

    return (\@pvals);
}


sub run {    
    my $dir = shift;
    opendir (E, "$dir");
    my @expts = grep (/^.+\..+$/, readdir(E));
    closedir (E);
    
    my $pm = Parallel::ForkManager->new($procs);
    foreach my $expt (@expts){
	$pm->start and next;
	warn "ILD for $expt\n";
	$ildobj->run_paup ("$dir/$expt");
	$pm->finish;
    }
    $pm->wait_all_children;
}
