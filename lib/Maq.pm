package Maq;

=head1 NAME                                                                       
                                                                                  
Maq.pm - contains methods for running maq
                                                                                  
=head1 SYNOPSIS                                                                   
                                                                                  
                                                                                  
                                                                                  
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

#####TODO: 

use strict;
use warnings;

sub new {
    my ($class, %args) = @_;

    # create the object with room for                                             
    # a tree, and aln, and the config                                             
    my $self = bless {
        'configuration' => undef,

    }, $class;

    # if any args, store those                                                    
    foreach my $arg (keys %args){
        $self->{$arg} = $args{$arg};
    }

    return $self;
}

sub load_config{
    my $self = shift;
    my $file = shift;

    # get and store key/val pairs in config                                       
    my $config = {};
    open (F, "$file") or
        die "Unable to open config file\n";
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) =
            ($line =~ m/^\s*(.*)\s*=\s*(.*)$/);
	if ($key eq "COVERAGES"){
	    my @covs = split (/\,/, $value);
	    $config->{$key} = \@covs;
	}
	else{
	    $config->{$key} = $value;
	}
    }
    close (F);

    $self->{'configuration'} = $config;
}

sub maq_simulate{
    my $self   = shift;
    my $reads  = shift;
    my $fasta  = shift;
    my $data   = shift;
    my $out    = shift;
    my $stderr_out = shift;
    
    # maq simulate                                                                     
    my $maqsim = "maq simulate";
    ($maqsim .= " $self->{'configuration'}->{'SIMULATE'}") 
	if ($self->{'configuration'}->{'SIMULATE'});
    $maqsim .= " -N $reads";
    $maqsim .= " $out/$reads-1.fastq";
    $maqsim .= " $out/$reads-2.fastq";
    $maqsim .= " $fasta";
    $maqsim .= " $data";
    `$maqsim &> $stderr_out/$reads.sim.stderr`;

}

sub maq_fastq2bfq{
    my $self = shift;
    my $in   = shift;
    my $out  = shift;
    my $stderr = shift;
    
    my $cmd = "maq fastq2bfq";
    $cmd .= " $in";
    $cmd .= " $out";
    `$cmd &> $stderr`;
}

sub maq_match{
    my $self = shift;
    my $reads   = shift;
    my $refname = shift;
    my $indir  = shift;
    my $outdir = shift;
    my $stderr = shift;

    my $cmd = "maq match";
    ($cmd .= " $self->{'configuration'}->{'MATCH'}")
        if ($self->{'configuration'}->{'MATCH'});
    $cmd .= " $outdir/$reads-v-$refname.match";
    $cmd .= " $indir/$refname.bfa";
    $cmd .= " $indir/$reads-1.bfq";
    $cmd .= " $indir/$reads-2.bfq";
    `$cmd &> $stderr/$reads-v-$refname.match.stderr`;
}

sub maq_mapview{
    my $self = shift;
    my $reads   = shift;
    my $refname = shift;
    my $resdir  = shift;
    my $stderr = shift;
    
    my $cmd = "maq mapview";
    $cmd .= " $resdir/$reads-v-$refname.match";
    `$cmd 1>$resdir/$reads-v-$refname.mapview 2>$stderr/$reads-v-$refname.mapview.stderr`;
}


sub maq_pileup{
    my $self = shift;
    my $reads   = shift;
    my $refname = shift;
    my $indir  = shift;
    my $resdir = shift;
    my $stderr = shift;

    my $cmd = "maq pileup";
    ($cmd .= " $self->{'configuration'}->{'PILEUP'}")
        if ($self->{'configuration'}->{'PILEUP'});
    $cmd .= " $indir/$refname.bfa $resdir/$reads-v-$refname.match";
    `$cmd 1>$resdir/$reads-v-$refname.pileup 2>$stderr/$reads-v-$refname.pileup.stderr`;
}

sub maq_fasta2bfa{
    my $self = shift;
    my $reffile  = shift;
    my $out  = shift;
    my $stderr_out = shift;
    
    my $ref_name;
    if ($reffile =~/\//g){
	$reffile =~m/.*\/(.*)$/;
	$ref_name = $1;
    }
    else {
	$ref_name = $reffile;
    }

    my $cmd = "maq fasta2bfa";
    $cmd .= " $reffile";
    $cmd .= " $out/$ref_name.bfa";
    `$cmd &> $stderr_out/$ref_name.bfa.stderr`;

    return ($ref_name);
}

sub get_coverages{
    my $self = shift;
    return $self->{'configuration'}->{'COVERAGES'};
}

sub get_length{
    my $self = shift;
    return $self->{'configuration'}->{'LENGTH'};
}

sub get_datafile{
    my $self = shift;
    return $self->{'configuration'}->{'DATAFILE'};
}


1;
