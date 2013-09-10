#!/usr/bin/perl -w

=head1 NAME

rast_submit.pl

=head1 SYNOPSIS

  rast_submit.pl -- 
              

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta directory or fasta file (script can handle both)
 --name        Is org name
 --taxid       Is the taxonomy id
 --domain      Is the domain (archae, bacteria)
 --outdir      Is your output dir
 --username    RAST username
 --password    RAST password

=head1 DESCRIPTION

Submit RAST annotation jobs to the RAST annotation server automatically.

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use RASTserver;

my ($help, $fasta, $name, $taxid, $outdir, $username, $password); #$domain
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'n|name=s'     => \$name,
    't|taxid=s'      => \$taxid,
#    'd|domain=s'       => \$domain,
    'u|username=s'     => \$username,
    'p|password=s'     => \$password,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $name, $taxid, $username, $password){ #$domain);
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# create the rast obk
my $rastobject = RASTserver->new($username, $password);

# submit the RAST job
my $rastjob = $rastobject->submit_RAST_job({
    -taxonomyID => $taxid,
    -filetype   => 'Fasta',
    -organismName => $name,
    -keepGeneCalls => 1,
    -geneCaller => 'RAST',
    -file => $fasta
    }
);

