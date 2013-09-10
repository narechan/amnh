package Mummer;

=head1 NAME                                                                       
                                                                                  
Mummer.pm - contains methods for running wga using mummer
                                                                                  
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
            ($line =~ m/^\s*(.*)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    $self->{'configuration'} = $config;
}

sub run_nucmer{
    my $self   = shift;
    my $fasta  = shift;
    my $ref    = shift;
    my $out    = shift;
    my $stderr_out = shift;

    # nucmer                                                                      
    my $nucmer = "nucmer";
    ($nucmer .= " $self->{'configuration'}->{'NUCMER'}") 
	if ($self->{'configuration'}->{'NUCMER'});
    $nucmer .= " -p $out/out";
    $nucmer .= " $ref";
    $nucmer .= " $fasta";
    `$nucmer &> $stderr_out/nucmer.stderr`;

    # show-coords                                                                 
    my $showcoords = "show-coords";
    ($showcoords .= " $self->{'configuration'}->{'SHOW-COORDS'}") 
	if ($self->{'configuration'}->{'SHOW-COORDS'});
    $showcoords .= " $out/out.delta";
    $showcoords .= " 1>$out/out.showcoords";
    `$showcoords 2>$stderr_out/showcoords.stderr`;
}

sub run_mummer{
    my $self   = shift;
    my $fasta  = shift;
    my $ref    = shift;
    my $out    = shift;
    my $stderr_out = shift;

    # mummer                                                            
    my $mummer = "mummer";
    ($mummer .= " $self->{'configuration'}->{'MUMMER'}")
        if ($self->{'configuration'}->{'MUMMER'});
    $mummer .= " $ref";
    $mummer .= " $fasta";
    $mummer .= " 1>$out/out.mummer";
    `$mummer 2>$stderr_out/mummer.stderr`;
}


sub parse_nucmer{
    my $self = shift;
    my $file = shift;

    my @refdata;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my @data = split (/\t/, $line);

        push @refdata,
        {
            'chrom' => $data[11],
            'start' => $data[0],
            'end'   => $data[1],
	};
    }
    close (F);

    return (\@refdata);
}

sub parse_mummer{
    my $self = shift;
    my $file = shift;

    my @refdata;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
	next if ($line =~m/^>/);
        my @data = split (/\s+/, $line);

        push @refdata,
        {
            'chrom' => $data[1],
            'start' => $data[2],
            'end'   => $data[2] + $data[4] - 1,
        };
    }
    close (F);

    return (\@refdata);
}

sub genome_coverage {
    my $self    = shift;
    my $refdata = shift;

    # sort ref data by chrom and by start                                                               
    my @refdata = sort {$a->{chrom} cmp $b->{chrom} ||
                            $a->{start} <=> $b->{start}
                    } @$refdata;
    
    # get the first element                                                                        
    my $first  = shift @refdata;
    my $chrom1 = $first->{chrom};
    my $start1 = $first->{start};
    my $end1   = $first->{end};

    # do joins on overlapping hits to get the                                                          
    # non redundant coverage                                                                      
    my $coverage = {};
    my $counter  = 0;
    foreach my $point (@refdata){

        # get data                                                                               
        my $chrom2 = $point->{chrom};
        my $start2 = $point->{start};
        my $end2   = $point->{end};

        if ($chrom1 eq $chrom2){
            if (($start1 > $end2) or ($end1 < $start2)) {
                $counter++;

		# they don't overlap at all                                                         
		$coverage->{$chrom1}->{$counter} = [$start1, $end1];
		$chrom1 = $chrom2;
		$start1 = $start2;
                $end1   = $end2;
            }
            else {
                my @coor = sort {$a <=> $b} ($start1, $start2, $end1, $end2);

                # get the extremes as the new range                                                   
                $start1 = shift @coor;
                $end1   = pop @coor;
            }
        }
        else {

            # put in the last values                                                                
            $coverage->{$chrom1}->{$counter + 1} = [$start1, $end1];

            # redefine for new chrom                                                             
            $counter = 0;
            $chrom1 = $point->{chrom};
            $start1 = $point->{start};
            $end1   = $point->{end};
        }

    }

    # put in the last values                                                                       
    $coverage->{$chrom1}->{$counter + 1} = [$start1, $end1];
    return ($coverage);

}

1;
