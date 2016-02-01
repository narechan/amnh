#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Ild;

my ($help, $outdir, $its); #where outdir is the directory from flocking that wraps around all the reps
GetOptions(
    'h|help'          => \$help,
    'o|outdir=s'      => \$outdir,
    'i|iterations=s'  => \$its,
    ) or pod2usage;

pod2usage if $help;

for my $option ($outdir){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

####MAIN####                                                                         

# parse the flocking info
opendir (D, "$outdir");
my @repdirs = sort (readdir (D));
shift @repdirs;
shift @repdirs;
closedir (D);
my $repdirs = @repdirs;

# stats by rep
my $repwt = {}; #index by rep
foreach my $repdir (@repdirs){
    print STDERR "rep $repdir\n";
    
    my $count = 0;
    for (my $i = 0; $i < $its; $i++){
	open (I, "$outdir/$repdir/logs/$i.ildaccess");
	while (my $line = <I>){
	    chomp $line;
	    $count++;
	}
	push (@{$repwt->{$repdir}}, $count);
    }
}

my $totwt = {}; #index by iteration                                                                           
my $seen  = {};
for (my $i = 0; $i < $its; $i++){
    print STDERR "it $i\n";
    
    foreach my $repdir (@repdirs){
        open (I, "$outdir/$repdir/logs/$i.ildaccess");
        while (my $line = <I>){
            chomp $line;
            ($totwt->{$i}->{$line}++) unless (exists ($seen->{$line}));
            $seen->{$line}++;
	}
    }
}


# print out the per rep ild comp matrix
open (PR, ">$outdir/ildburden_perrep.out");
foreach my $rep (sort keys %$repwt){
    my $string = join "\t", @{$repwt->{$rep}};
    print PR "$rep\t$string\n";
}
close (PR);

# print out the total ild comp burden
open (TR, ">$outdir/ildburden_acrossreps.out");
foreach my $it (sort {$a <=> $b} keys %$totwt){
    my @array = keys %{$totwt->{$it}};
    my $array = @array;
    print TR "$it\t$array\n";
}
close (TR);
		       
