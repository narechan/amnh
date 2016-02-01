#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Ild;

my ($help, $matrix, $indir, $outdir, $burnstart, $burnend); #where indir is the directory from flocking
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrix,
    'i|indir=s'       => \$indir,
    'o|outdir=s'      => \$outdir,
    's|startburn=s'   => \$burnstart,
    'e|endburn=s'     => \$burnend,
    ) or pod2usage;

pod2usage if $help;

for my $option ($matrix, $outdir, $indir, $burnstart, $burnend){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

`mkdir -p $outdir`;

####MAIN####                                                                         

# create a symbol library (each letter of the alphabet)
my $library = {};
my $c = -1;
my @letters;
for my $letter ("A".."Z"){
    $c++;
    $library->{$c} = $letter;
    push (@letters, $letter);
}
my $symbolstring = join " ", @letters;

# parse the nexus file and store the data                                            
my $ildobj = Ild->new;
$ildobj->load_aln ($matrix);
my $charsets = $ildobj->get_charsets;
my @charsets = keys %$charsets;
my $charsetcnt = @charsets;

# parse the flocking info
opendir (D, "$indir");
my @repdirs = sort (readdir (D));
shift @repdirs;
shift @repdirs;
closedir (D);
my $repdirs = @repdirs;

my $pa = {};
my $apd = {};
foreach my $repdir (@repdirs){
    opendir (R, "$indir/$repdir");
    my @its = sort (readdir (R));
    shift @its;
    shift @its;
    closedir (R);
    
    for (my $i = $burnstart; $i <= $burnend; $i++){
	print STDERR "$repdir\t$i\n";
	
	# if the optics file exists, parse it,
	# otherwise insert missing data for a replicate
	# that failed at this point
	if (-e "$indir/$repdir/flocks/$i.optics"){
	    open (I, "$indir/$repdir/flocks/$i.optics");
	    while (my $line = <I>){
		chomp $line;
		my ($flock, $gene) = split (/\t/, $line);
		if ($flock > 26){
		    print STDERR "Too many flocks: not enough characters!\n";
		    die;
		}
		else{
		    push (@{$pa->{$gene}}, $library->{$flock});
		    push (@{$apd->{$gene}}, $flock);
		}
	    }
	    close (I);
	}
	else {
	    foreach my $charset (@charsets){
		push (@{$pa->{$charset}}, "?");
		push (@{$apd->{$charset}}, "?");
	    }
	}
    }
}

my $ntax = 0;
foreach my $tax (sort keys %$pa){
    $ntax++;
}

open (O, ">$outdir/nexus");
print O "#NEXUS\n";
print O "BEGIN DATA;\n";
print O "DIMENSIONS NTAX=$ntax NCHAR=$repdirs;\n"; #$charsetcnt
print O "FORMAT SYMBOLS=\"$symbolstring\";\n";
print O "MATRIX\n";
foreach my $tax (sort keys %$pa){
    my @paarray = @{$pa->{$tax}};
    my $pastring = join "", @paarray;
    print O "$tax\t$pastring\n";
}
print O ";\n";
print O "END;\n";
close (O);

open (A, ">$outdir/apd");
print A "ST\t";
my @header;
for (my $i = 1; $i <= $repdirs; $i++){
    my $string = "rep" . $i;
    push (@header, $string);
}
my $header = join "\t", @header;
print A "$header\n";

foreach my $tax (sort keys %$apd){
    my @paarray = @{$apd->{$tax}};
    my $pastring = join "\t", @paarray;
    print A "$tax\t$pastring\n";
}
close (A);


