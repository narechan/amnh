#!/usr/bin/perl -w

=head1 NAME

ild_genome.pl

--config       is the config file
--matrix       is the data matrix
--outdir       is your output dir for the run
--procs        is the number of processors you want to use
--outdir       is your output directory

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($help, $configfile, $matrixfile, $outdir, $procs);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    'p|procs=s'        => \$procs,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $matrixfile, $outdir, $procs){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# create dir structure for the results
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;

#####MAIN#####

# load the config
my $config = load_config ($configfile);

# load the charsets and find nchar
my ($charsets, $nchar) = load_aln ($matrixfile);

# build all the commands
foreach my $charset (sort keys %$charsets){
    my ($beg, $fin) = split (/\-/, $charsets->{$charset});
    generate_ild_sr ($beg,
		     $fin,
		     "$outdir/cmds",
		     "$outdir/logs",
		     $nchar,
		     $charset,
		     $matrixfile,
	             $config);
}

# run all the commands in parallel using forks
opendir (D, "$outdir/cmds");
my @cmds = sort (readdir (D));
shift @cmds;
shift @cmds;
closedir (D);

my $pm = Parallel::ForkManager->new($procs);
foreach my $cmd (@cmds){
    $pm->start and next;
    `paup $outdir/cmds/$cmd`;
#    last;
    $pm->finish;
}
$pm->wait_all_children;

# parse all the log files and generate report
opendir (L, "$outdir/logs");
my @logs = sort (readdir (L));
shift @logs;
shift @logs;
closedir (L);

foreach my $log (@logs){
    my $pval;
    my $len;

    open (LOG, "$outdir/logs/$log");
    while (my $line = <LOG>){
	chomp $line;
	if ($line =~m/\s*P\s{1}value\s*\=.*\=\s*(.+)/){
	    $pval = $1;
	}
	if ($line =~m/\s*(\d+)\*\s*\d+/){
	    $len  = $1;
	}
    }
    print "$log\t$len\t$pval\n";
}

########SUBS#########
    

sub generate_ild_sr{
    my $start      = shift;
    my $end        = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $nchar      = shift;
    my $partition  = shift;
    my $alignfile  = shift;
    my $config     = shift;
    
    open (ILD, ">$outdir_cmd/$partition.nex");
    
    # print out the header                                                           
    print ILD "#NEXUS\n";
    print ILD "SET increase = auto;\n";
    print ILD "EXECUTE $alignfile;\n";
    print ILD "LOG START file = $outdir_log/$partition.log replace = yes;\n";
    
    # create the partition set                                                       
    print ILD "BEGIN SETS;\n";
    print ILD "CHARSET C$partition = ";
    print ILD "$start-$end;\n";
    print ILD "END;\n";

    # create the experiment                                                          
    print ILD "BEGIN PAUP;\n";

    my $fivestart  = 1;
    my $fiveend    = $start - 1;

    my $threestart = $end + 1;
    my $threeend   = $nchar;

    print ILD "[!==>ILD: C$partition]\n";
    if ($fiveend == 0){
        print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$threestart-$threeend;\n";
    }
    elsif ($threestart > $nchar){
        print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$fivestart-$fiveend;\n";
    }
    else{
        print ILD "CHARPARTITION CharPart$partition = ";
        print ILD "1:C$partition, 2:$fivestart-$fiveend $threestart-$threeend;\n";
    }

    print ILD "$config->{HOMPART} ";
    print ILD "partition = CharPart$partition";

    if ($config->{SEARCH}){
        print ILD " / $config->{SEARCH};\n";
    }
    else{
        print ILD ";\n";
    }
    
    print ILD "[!==>end_record]\n";
    
    print ILD "END;\n";
    print ILD "LOG STOP;\n";
    print ILD "QUIT /warnTsave=no;\n";    
}

sub load_config{
    my $file = shift;

    # get and store key/val pairs in config                                          
    my $config = {};
    open (F, "$file") or
        die "Unable to open config file\n";
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) =
            ($line =~ m/^\s*(\w+)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    return ($config);
}

sub load_aln{
    my $alignfile = shift;

    open (NEX, "$alignfile");
    my $charset = {};
    my $nchar;
    while (my $line = <NEX>){
        chomp $line;

	# mine for nchar
        ($nchar = $1) if ($line =~m/nchar\s*=\s*(\d+)/i);

	# mine for charsets
        if ($line =~m/^charset/i){
            $line =~s/charset//ig;
            $line =~s/\s+//g;
            $line =~s/\;//g;

            my ($partition, $coords) =
                split (/\=/, $line);
            $charset->{$partition} = $coords;
        }
    }
    close (NEX);

    return ($charset, $nchar);
}
