#!/usr/bin/perl -w

# For any program that uses one file input and one dir output, use this
# script to fork any program across all infiles.

# config would look something like this, where [IN] and [OUT] are specified 
# through the options below

#PROGRAM=~apurva/packages/vector_screening/crossmatch/cross_match.manyreads [INDI\
#R]/[INFILE] ~apurva/packages/vector_screening/UniVec_Core -screen 1>[OUT]/[INFIL\
#E].out 2>[OUT]/[INFILE].stderr

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

use Mummer;

my ($help, $procs, $configfile, $indir, $outdir);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'c|config=s'      => \$configfile,
    'o|outdir=s'      => \$outdir,
    'i|indir=s'       => \$indir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $procs, $outdir, $indir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

# load the config
my $config = load_config ($configfile);

# open all the indir files
opendir (D, "$indir");
my @files = sort (readdir (D));
shift @files;
shift @files;
closedir (D);

# execute the program in the config for every file in the indir
my $mummerpm = Parallel::ForkManager->new($procs);
foreach my $file (sort @files){
    $mummerpm->start and next;
    
    my $runner = $config->{'PROGRAM'};
    $runner =~s/\[INDIR\]/$indir/g;
    $runner =~s/\[INFILE\]/$file/g;
    $runner =~s/\[OUT\]/$outdir/g;
    
    print STDERR "$runner\n";
    
    `$runner`;
    
    $mummerpm->finish;
}
$mummerpm->wait_all_children;

####SUBS####

sub load_config{
    my $file = shift;

    # get and store key/val pairs in config                                                              
    my $config = {};
    open (F, "$file") or
        die "Unable to open config file\n";
    while (my $line = <F>){
	chomp $line;
        my ($key, $value) =
            ($line =~ m/^\s*(.*)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    return ($config);
}
