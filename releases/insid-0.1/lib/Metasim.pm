package Metasim;

=head1 NAME                                                                       
                                                                                  
Metasim.pm - contains methods for running metasim
                                                                                  
=head1 AUTHOR                                                                     
                                                                                  
Apurva Narechania                                                                 
anarechania *a|t* amnh.org                                                        
                                                                                  
=head1 COPYRIGHT                                                                  
                                                                                  
Copyright (c) 2010 American Museum of Natural History                             
                                                                                  
This library is free software;                                                    
you can redistribute it and/or modify                                             
it under the same terms as Perl itself.                                           
    
=cut                                             
    
# ----------------------------------------------------                            

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

sub generate_reads{
    my $self   = shift;
    my $readnum  = shift;
    my $fasta  = shift;
    my $out    = shift;
    my $stderr_out = shift;
    my $counter = shift;
    
    # name the reads with the counter
    # if counter defined, else with the readnum
    my $reads;
    if ($counter){
	$reads = $counter;
    }
    else {
	$reads = $readnum;
    }
    
    # make tmp outdir                                                             
    `mkdir $out/$reads-tmp`;

    # metasim                                                                     
    my $metasim = "MetaSim cmd";
    ($metasim .= " $self->{'configuration'}->{'METASIM'}") 
	if ($self->{'configuration'}->{'METASIM'});
    $metasim .= " -r $readnum";
    $metasim .= " -d $out/$reads-tmp";
    $metasim .= " $fasta";
    `$metasim &> $stderr_out/stderr.$reads`;
#    print STDERR "$metasim\n";
    `mv $out/$reads-tmp/* $out/$reads.fasta`;
    `rm -rf $out/$reads-tmp/`;

}

sub get_coverages{
    my $self = shift;
    return $self->{'configuration'}->{'COVERAGES'};
}

sub get_metasim_params{
    my $self = shift;
    return $self->{'configuration'}->{'METASIM'};
}



1;
