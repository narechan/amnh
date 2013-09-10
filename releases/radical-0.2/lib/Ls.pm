package Ls;

=head1 NAME

Ls.pm - contains methods for the LS pipeline

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2011 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

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
	'tree'          => undef,
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
    my $header;
    my $configuration = {};
    open (C, "$file");
    while (my $line = <C>){
	next if ($line =~/^\n/);
	chomp $line;

	if ($line =~m/^\[(.*)\]/){
	    $header = $1;
	    next;
	}

	$line =~s/\s//g;
	my ($key, $value) = split (/\=/, $line);

	$configuration->{$header}->{$key} = $value;
    }
    close (C);

    $self->{'configuration'} = $configuration;
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

sub load_tree{
    my $self     = shift;
    my $treefile = shift;
    my $root     = shift;

    # generate the tree object                                                                         
    my $in_tree  = Bio::TreeIO->new(-file   => "$treefile",
                                    -format => "nexus");
    my $tree     = $in_tree->next_tree;
    
    # re-root if required                                                                              
    if ($root){
        my @nodes = $tree->find_node(-id => $root);
	my $node = $nodes[0];
        $tree->reroot($node);
    }

    # get all nodes and internal nodes                                                                  
    my @allnodes      = $tree->get_nodes;
    my @internalnodes = grep {!$_->is_Leaf} @allnodes;

    # cycle through internal nodes and get each internal                                                
    # node's immediate descendents. Get all leaves from                                                 
    # each immediate descendent and store.                                                           
    my $lineages = {};
    my $internal_counter = 0;
    foreach my $node (@internalnodes){
	$internal_counter++;
        my @lineages = $node->each_Descendent;

        my $lineage_counter = 0;
        foreach my $lineage (@lineages){
            $lineage_counter++;

            my @children;
            if ($lineage->is_Leaf){
		push (@children, $lineage);
            }
	    else {
                my @descendents = $lineage->get_all_Descendents;
                @children = grep{$_->is_Leaf} @descendents;
            }

            foreach my $leaf (@children){
                my $id = $leaf->id;
                push @{$lineages->{$internal_counter}->{$lineage_counter}}, $id;
            }
        }
    }
    $self->{'tree'} = $lineages;

}

sub run_garli{
    my $self = shift;
    my $config = shift;
    `Garli-1.0 $config`;
}

sub get_charsets{
    my $self = shift;
    return $self->{'alignment'}->{'charset'};
}

sub get_nchar{
    my $self = shift;
    return $self->{'alignment'}->{'nchar'};
}

sub get_config{
    my $self = shift;
    return $self->{'configuration'};
}

sub get_leaves{
    my $self = shift;

    my %nodes;
    my $lineages = $self->{'tree'};
    foreach my $node (keys %$lineages){
        foreach my $lineage (keys %{$lineages->{$node}}){
            push @{$nodes{$node}}, @{ $lineages->{$node}->{$lineage} };
        }
    }

    return (\%nodes);
}

sub get_lineages{
    my $self = shift;
    return $self->{'tree'};
}

1;
