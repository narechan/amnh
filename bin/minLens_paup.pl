#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;
use Parallel::ForkManager;

my ($help, $matrixfile, $configfile, $outdir, $procs, $tssummary);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'p|procs=s'        => \$procs,
    't|tssum=s'        => \$tssummary,
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $outdir, $procs){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/parts`;
`mkdir -p $outdir/logs`;


#####MAIN#####

# instantiate the object and load required data
my $supportobj = TreeSupports->new;
$supportobj->load_config ($configfile);
my $config = $supportobj->get_config;

$supportobj->load_aln    ($matrixfile);
my $charsets = $supportobj->get_charsets;
=head
my $pm = Parallel::ForkManager->new($procs);
foreach my $charset (keys %$charsets){
    $pm->start and next;
    print STDERR "$charset\n";

    open (P, ">$outdir/parts/$charset.nex");
    print P "#NEXUS\n";
    print P "set increase=auto;\n";
    print P "set WarnRoot=no;\n";
    print P "execute $matrixfile;\n";
    print P "log start file=$outdir/logs/$charset replace=yes;\n";
    print P "exclude all;\n";
    print P "include $charsets->{$charset};\n";
    print P "$config->{'TREECOMMAND'};\n";
    print P "pscores /ci=yes;\n";
    print P "log stop;\n";
    print P "quit /warnTsave=no;\n";
    close (P);

    `paup $outdir/parts/$charset.nex`;

    $pm->finish;
}
$pm->wait_all_children;
=cut
opendir (L, "$outdir/logs");
my @logs = sort (readdir (L));
shift @logs;
shift @logs;
closedir (L);

open (ML, ">$outdir/minLens.out");
my $minlens = {};
foreach my $log (@logs){
    open (LF, "$outdir/logs/$log");
    while (my $line = <LF>){
	chomp $line;
	
	if ($line =~m/min\.\spossible\slengths/){
	    $line =~s/\s//g;
	    my ($stuff, $val) = split (/\=/, $line);
	    $minlens->{$log} = $val;
	    print ML "$log\t$val\n";
	}
	else{
	    next;
	}
    }
    close (LF);
}
close (ML);	    

open (TS, "$tssummary");
while (my $line = <TS>){
    chomp $line;
    next if ($line =~m/Part/);
    my @line = split (/\t/, $line);
    my $pbs = $line[10];
    my $pbsl = $line[11];
    my $gene = $line[0];
    
    my $pbsnorm = $pbs / $minlens->{$gene};
    my $pbslnorm = $pbsl / $minlens->{$gene};

    print "$line\t$pbsnorm\t$pbslnorm\n";
}
close (TS);
