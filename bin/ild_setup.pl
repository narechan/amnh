#!/usr/bin/perl -w

=head1 NAME

ild_setup.pl

=head1 SYNOPSIS

ild_setup.pl

Options:

--config       is the config file
--matrix       is the data matrix
--outdir       is your output dir for the run
--window       is your window size if doing 
                LILD sildeRule (optional)
--motion       is the number of residues the                                                               
                window moves per calculation (5' -> 3') 
                if doing LILD slideRule (optional)
--start        is the start position of slideRule; set to 1 by default (optional)                            
--end          is the end position of slideRule; set to nchar by default    (optional)
--file1        is your first file of partitions (optional)                                            
--file2        is your second set of partitions (optional)  
--multiple     is if your infiles have multiple charsets per ild partition (optional)

the matrix must be in
nexus format.

Requires the bioperl libs.
Requires Lild.pm

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

my ($help, $configfile, $matrixfile, $outdir, $window, $motion, $start, $end, $file1, $file2, $multiple);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    'a|file1=s'        => \$file1,
    'b|file2=s'        => \$file2,
    'z|multiple'     => \$multiple,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $matrixfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


# create dir structure for the results
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;

#####MAIN#####

# instantiate the object and load required data
my $ildobj = Ild->new;
$ildobj->load_config ($configfile);
$ildobj->load_aln    ($matrixfile);

# define start as position 1 and and end as nchar
# unless otherwise specified
my $nchar = $ildobj->get_nchar;
($start = 1) unless ($start);
($end   = $nchar) unless ($end);

# get the charsets depending on expt (slide rule or defined charsets)
my $charsets;
if ( ($window) and ($motion) and ($start) and ($end) ){
    $charsets = $ildobj->generate_partitions ($window, $motion, $start, $end);
}
else {
    $charsets = $ildobj->get_charsets;
}

# generate the paup ild cmds
my $counter = 0;

# for slideRule case
if ( ($window) and ($motion) and ($start) and ($end) ){
    foreach my $charset (sort keys %$charsets){
	my ($beg, $fin) = split (/\-/, $charsets->{$charset});
	$ildobj->generate_ild_sr ($beg,
				  $fin,
				  "$outdir/cmds",
				  "$outdir/logs",
				  $end,
				  $charset,
				  $matrixfile);
    }
}

# for the case where pairwise ilds are defined by infiles (no overlaps)
elsif ( ($file1) and ($file2) and !$multiple){
    `mkdir -p $outdir/nxs`;
    my ($partitions, $lengths) = $ildobj->store_alignment_seq($matrixfile, $charsets);

    open (FA, "$file1");
    while (my $line1 = <FA>){
	chomp $line1;
	
	open (FB, "$file2");
	while (my $line2 = <FB>){
	    chomp $line2;

	    $counter++;
	    my $part1 = $partitions->{$line1};
	    my $part2 = $partitions->{$line2};
	    my $len1  = $lengths->{$line1};
	    my $len2  = $lengths->{$line2};

	    my $party = $line1 . "v" . $line2;
	    $ildobj->generate_nxs_pairwise ($part1, $part2, $len1, $len2, $line1, $line2, "$outdir/nxs");
	    $ildobj->generate_ild_part ($counter,
					$line1,
                                        $line2,
                                        "$outdir/cmds",
                                        "$outdir/logs",
                                        "$outdir/nxs/$party.nxs");
	}
	close (FB);
    }
    close (FA);
}
	    
# for the case where you want to do all pairwise ilds from a single file
elsif ( ($file1) and !($file2) ){
    `mkdir -p $outdir/nxs`;
    my ($partitions, $lengths) = $ildobj->store_alignment_seq($matrixfile, $charsets);
    
    my $file = {};
    open (F, "$file1");
    while (my $line = <F>){
	chomp $line;
	$file->{$line} = 1;
    }
    close (F);
    
    foreach my $line1 (sort keys %$file){
        foreach my $line2 (sort keys %$file){
            next if ($line1 eq $line2);
	    
	    $counter++;
	    my $part1 = $partitions->{$line1};
            my $part2 = $partitions->{$line2};
            my $len1  = $lengths->{$line1};
            my $len2  = $lengths->{$line2};

	    my $party = $line1 . "v" . $line2;
            $ildobj->generate_nxs_pairwise ($part1, $part2, $len1, $len2, $line1, $line2, "$outdir/nxs");
            $ildobj->generate_ild_part ($counter,
					$line1,
                                        $line2,
                                        "$outdir/cmds",
                                        "$outdir/logs",
                                        "$outdir/nxs/$party.nxs");
	}
	delete $file->{$line1};
    }
}

# for the case where your infiles may have more than one charset per partition
# this case does not build isolating nexus files so it may be slow to
# load the matrix for large datasets.
# infiles must have multiple charsets separated by underscores!!!
elsif ( ($file1) and ($file2) and ($multiple) ){
    open (FA, "$file1");
    while (my $line1 = <FA>){
        chomp $line1;
	my @line1 = split (/\_/, $line1);
	
        open (FB, "$file2");
        while (my $line2 = <FB>){
            chomp $line2;
	    my @line2 = split (/\_/, $line2);
	    
	    # check to see if there is charset overlap
	    my $signal = 0;
	    foreach my $ele2 (@line2){
		foreach my $ele1 (@line1){
		    ($signal++) if ($ele1 eq $ele2);
		}
	    }
	    next if ($signal > 0);
	    
	    $counter++;
	    $ildobj->generate_ild_part ($counter,
					$line1,
                                        $line2,
                                        "$outdir/cmds",
                                        "$outdir/logs",
					$matrixfile);
        }
	close (FB);
    }
    close (FA);
}


# for the pairwise all-all case for defined charsets
# STILL NEEDS PAIRWISE MATRIX GENERATION!!! running on whole right now!!
else {
    foreach my $charset1 (sort keys %$charsets){
	foreach my $charset2 (sort keys %$charsets){
	    next if ($charset1 eq $charset2);
	    $counter++;
	    $ildobj->generate_ild_part ($counter,
					$charset1,
					$charset2,
					"$outdir/cmds",
					"$outdir/logs",
					$matrixfile);
	}
	delete $charsets->{$charset1};
#	last;
    }
}

