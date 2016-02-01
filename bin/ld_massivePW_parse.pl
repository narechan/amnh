#!/usr/bin/perl -w

=head1 NAME

ild_massivePW_parse.pl

=head1 SYNOPSIS

ild_massivePW_parse.pl

Options:

--outdir       is your output dir for the run
--matrix       is the matrix used for calculation

the matrix must be in
nexus format.

the program deposits table and matrix into the 
root level of $outdir

Requires the bioperl libs.
Requires Lild.pm

=head1 DESCRIPTION

This program does the parse on ild_massivePW.pl

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

my ($help, $outdir, $matrixfile);
GetOptions(
    'h|help'           => \$help,
    'o|outdir=s'       => \$outdir,
    'm|matrix=s'       => \$matrixfile,
    ) or pod2usage;
pod2usage if $help;

for my $option ($outdir, $matrixfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


#####MAIN#####

# instantiate the object and load required data
my $ildobj = Ild->new;
$ildobj->load_aln ($matrixfile);
my $charsets = $ildobj->get_charsets;

# parse all the sg commands
my $sglens = {};

opendir (L, "$outdir/sg_logs");
my @sglogs = sort (readdir (L));
shift @sglogs;
shift @sglogs;
closedir (L);

foreach my $sglog (@sglogs){
#    print STDERR "$sglog\n";

    my ($lens, $numlens) = $ildobj->parse_treelens ("$outdir/sg_logs/$sglog");

    foreach my $part (keys %$lens){
	$sglens->{$part} = $lens->{$part};
    }
}

# parse all the pw commands
my $pwlens = {};

opendir (X, "$outdir/pw_logs");
my @pwlogs = sort (readdir (X));
shift @pwlogs;
shift @pwlogs;
closedir (X);

foreach my $pwlog (@pwlogs){
#    print STDERR "$pwlog\n";

    my ($lens, $numlens) = $ildobj->parse_treelens ("$outdir/pw_logs/$pwlog");

    foreach my $part (keys %$lens){
	my ($part1, $part2) = split (/v/, $part);
	$pwlens->{$part1}->{$part2} = $lens->{$part};
#        $pwlens->{$part} = $lens->{$part};
    }
}

# get the header
my @header;
foreach my $cs (sort keys %$charsets){
    push (@header, $cs);
}
my $header = join "\t", @header;

# calc the ILD stats and bulid the matrix
open (M, ">$outdir/matrix");
open (H, ">$outdir/matrix.header");
print H "\t$header\n";
foreach my $charset1 (sort keys %$charsets){
    my @row;

    foreach my $charset2 (sort keys %$charsets){
	if (exists $pwlens->{$charset1}->{$charset2}){
	    if ($pwlens->{$charset1}->{$charset2} == 0){
		push (@row, 0);
	    }
	    else{
		my $ild =
		    ($pwlens->{$charset1}->{$charset2} - ($sglens->{$charset1} + $sglens->{$charset2})) / $pwlens->{$charset1}->{$charset2};
		push (@row, $ild);
	    }
	}
	elsif (exists $pwlens->{$charset2}->{$charset1}){
	    if ($pwlens->{$charset2}->{$charset1} == 0){
		push (@row, 0);
            }
	    else{
		my $ild =
		    ($pwlens->{$charset2}->{$charset1} - ($sglens->{$charset1} + $sglens->{$charset2})) / $pwlens->{$charset2}->{$charset1};
		push (@row, $ild);
	    }
	}
	elsif ($charset1 eq $charset2){
	    push (@row, 0);
	}
	else{
	    print STDERR "Problem with $charset1 $charset2\n";
	    die;
	}
    }
    
    my $rowstring = join "\t", @row;
    print M "$charset1\t$rowstring\n";
    print H "$charset1\t$rowstring\n";
}
close (M);
close (H);

# calc the ILD stats
open (T, ">$outdir/table");
open (Z, ">$outdir/table2");
foreach my $pw1 (sort keys %$pwlens){
    foreach my $pw2 (sort (keys %{$pwlens->{$pw1}})){
	my $acc = $pw1 . "v" . $pw2;
	if ($pwlens->{$pw1}->{$pw2} == 0){
	    print T "$pw1\t$pw2\t$pwlens->{$pw1}->{$pw2}\t$sglens->{$pw1}\t$sglens->{$pw2}\t0\n";
	    print Z "$acc\t0\n";
	}
	else {
	    my $ild = 
		($pwlens->{$pw1}->{$pw2} - ($sglens->{$pw1} + $sglens->{$pw2})) / $pwlens->{$pw1}->{$pw2};
	    print T "$pw1\t$pw2\t$pwlens->{$pw1}->{$pw2}\t$sglens->{$pw1}\t$sglens->{$pw2}\t$ild\n";
	    print Z "$acc\t$ild\n";
	}
    }
}
close (T);
close (Z);
