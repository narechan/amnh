#!/usr/bin/perl -w

=head1 NAME

atomizer.pl

=head1 SYNOPSIS

  atomizer.pl

Options:

 --help        Show brief help and exit
 --file        Is your input summary file
 --type        Is the type of summary file input

=head1 DESCRIPTION

Given support values from a slideRule run,
average supports across each residue, and atomize
support calcs to the residue level

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

#####TODO:
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);
use Math::Round;

my ($help, $infile, $type);
GetOptions(
    'h|help'          => \$help,
    'f|infile=s'      => \$infile,
    't|type=s'        => \$type,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile, $type){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# parse the TS summary file and store data
my $atoms = {};
open (S, "$infile");
while (my $line = <S>){
    chomp $line;
    next if ($line =~m/^Part\tCoords/); # all headers should start this way
    
#    my ($part, $coords, $node, $topo, $bsch, $bswt, $bswot, 
#	$bswl, $bswol, $bs, $pbsunl, $pbswl, $pbswol, $pbs, 
#	$pbsl, $phbs, $phbsl) = split (/\t/, $line);

    my @data = split (/\t/, $line);
    if ($type eq "ild"){
	my $node = "NA";
	my ($start, $end) = split (/\-/, $data[1]);
	for (my $pos = $start; $pos <= $end; $pos++){                                   
	    push @{$atoms->{$pos}->{$node}->{'2'}}, $data[2];
	}
    }
    elsif ($type eq "lild"){
	my $node = $data[2];
	my ($start, $end) = split (/\-/, $data[1]);
	for (my $pos = $start; $pos <= $end; $pos++){
            if ($data[6]){
		push @{$atoms->{$pos}->{$node}->{'6'}}, $data[6];
	    }
	    else{
		push @{$atoms->{$pos}->{$node}->{'6'}}, 1; # a very nonsig pval!
	    }
        }
    }
    elsif ($type eq "ts"){
	my $node = $data[2];
	my ($start, $end) = split (/\-/, $data[1]);
	for (my $pos = $start; $pos <= $end; $pos++){
	    push @{$atoms->{$pos}->{$node}->{'9'}}, $data[9];
	    push @{$atoms->{$pos}->{$node}->{'13'}}, $data[13];
	    push @{$atoms->{$pos}->{$node}->{'15'}}, $data[15];
	}
    }
    else {
	print STDERR "Unknown summary file type.\n";
	die;
    }

#    for (my $pos = $start; $pos <= $end; $pos++){
#	push @{$atoms->{$pos}->{$node}->{'bs'}}, $bs;
#	push @{$atoms->{$pos}->{$node}->{'pbs'}}, $pbs;
#	push @{$atoms->{$pos}->{$node}->{'phbs'}}, $phbs;
#    }
}
close (S);

foreach my $pos (sort {$a <=> $b} keys %$atoms){
    foreach my $node (sort {$a <=> $b} keys %{$atoms->{$pos}}){
	my $mastavg;
	my $mastraw;
	foreach my $col (sort {$a <=> $b} keys %{$atoms->{$pos}->{$node}}){
	    my @data = @{$atoms->{$pos}->{$node}->{$col}};
	    $mastavg .= nearest (.01, sum(@data) / @data);
	    $mastraw .= join ",", @data;
	    $mastavg .= "\t";
	    $mastraw .= "\t";
	}
	chop $mastavg;
	chop $mastraw;
	
	print "$pos\t$node\t$mastavg\t$mastraw\n";

	
#	my @bs      = @{$atoms->{$pos}->{$node}->{'bs'}};
#	my $bsavg   = nearest (.01, sum(@bs) / @bs);
#	my @pbs     = @{$atoms->{$pos}->{$node}->{'pbs'}};
#	my $pbsavg  = nearest (.01, sum(@pbs) / @pbs);
#	my @phbs    = @{$atoms->{$pos}->{$node}->{'phbs'}};
#	my $phbsavg = nearest (.01, sum(@phbs) / @phbs);
	
#	my $bs = join ",", @bs;
#	my $pbs = join ",", @pbs;
#	my $phbs = join ",", @phbs;
	
#	print "$pos\t$node\t$bsavg\t$pbsavg\t$phbsavg\t";
#	print "$bs\t$pbs\t$phbs\n";
    }
}

    
