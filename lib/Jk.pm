package Jk;

=head1 NAME

Jk.pm - contains methods for the JK pipeline

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
#####      Document

use strict;
use warnings;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::SeqIO;

our $VERSION = "0.01";

sub new {
    my ($class, %args) = @_;
    
    # create the object with room for
    # a tree, and aln, and the config
    my $self = bless {
	'configuration' => undef,
	'alignment'     => undef,
	
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
	    ($line =~ m/^\s*(\w+)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    $self->{'configuration'} = $config;
}

sub load_aln{
    my $self      = shift;
    my $alignfile = shift;

    open (NEX, "$alignfile");
    my $charset = {};
    my $nexus   = {};
    while (my $line = <NEX>){
        chomp $line;

        # take only first instances of all of these things                                                 
        # header only                                                                                 
        ($nexus->{'nchar'}    = $1) if (($line =~m/nchar\s*=\s*(\d+)/i) and (!$nexus->{'nchar'}));
        ($nexus->{'ntax'}     = $1) if (($line =~m/ntax\s*=\s*(\d+)/i) and (!$nexus->{'ntax'}));
        ($nexus->{'datatype'} = $1) if (($line =~m/datatype\s*=\s*(\w+)/i) and (!$nexus->{'datatype'}));
        ($nexus->{'missing'}  = $1) if (($line =~m/missing\s*=\s*(.{1})/i) and (!$nexus->{'missing'}));
        ($nexus->{'gap'}      = $1) if (($line =~m/gap\s*=\s*(.{1})/i) and (!$nexus->{'gap'}));

        if ($line =~m/outgroup/i){
            $line =~s/outgroup//ig;
            $line =~s/\s+//g;
            $line =~s/\;//g;

            # any instances of more than one outgroup???? <====FIX                                     
            $nexus->{'outgroup'} = $line;
        }

        if ($line =~m/^charset/i){
            $line =~s/charset//ig;
            $line =~s/\s+//g;
            $line =~s/\;//g;

            my ($partition, $coords) =
                split (/\=/, $line);
            $charset->{$partition} = $coords;
        }
        $nexus->{'charset'} = $charset;
    }
    close (NEX);

    $self->{'alignment'} = $nexus;
}


sub get_pctdel{
    my $self = shift;
    return $self->{'configuration'}->{'PCTDELETE'};
}

sub get_tax{
    my $self = shift;
    return $self->{'alignment'}->{'ntax'};
}

sub run_paup{
    my $self   = shift;
    my $exptfile   = shift;
    
    `paup -n $exptfile`;
    
}

sub generate_jk_part{
    my $self       = shift;
    my $pdel      = shift;
    my $outdir_cmd = shift;
    my $outdir_log = shift;
    my $alignfile  = shift;

    open (JK, ">$outdir_cmd/$pdel.nex");

    # print out the header                                                                           
    print JK "#NEXUS\n";
    print JK "SET increase = auto;\n";
    print JK "EXECUTE $alignfile;\n";
    print JK "LOG START file = $outdir_log/$pdel.log replace = yes;\n";

    # create the experiment              
    print JK "BEGIN PAUP;\n";
    print JK "[!==>JK: $pdel]\n";
    print JK "$self->{'configuration'}->{'JACKKNIFE'} ";
    print JK "pctdelete=$pdel";

    if ($self->{'configuration'}->{'SEARCH'}){
        print JK " / $self->{'configuration'}->{'SEARCH'};\n";
    }
    else{
        print JK ";\n";
    }

    print JK "[!==>end_record]\n";

    print JK "END;\n";
    print JK "LOG STOP;\n";
    print JK "QUIT /warnTsave=no;\n";
}

sub parse_jk{
    my $self    = shift;
    my $file    = shift;

    my @freqs;
    open (JK, "$file");
    while (my $line = <JK>){
        chomp $line;

	if ($line =~m/\*/){
	    next if ($line =~m/P A U P/);
	    my ($taxcode, $freq) = split (/\s+/, $line);
	    push (@freqs, $freq);
	}
    }
    close (JK);
    
    return (\@freqs);
}


1;
