#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use File::Copy; 

#use Bio::SearchIO;
use Bio::SeqIO;
#use Bio::FeatureIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use IO::File;

my( 
    $thresh,
    $evalueThresh,
    $help,
    $outfile,

    );



&GetOptions( 
	     "t=s"           => \$thresh,
	     "e=s"           => \$evalueThresh,
	     
	     "o=s"           => \$outfile,
	     );
my $cmfile = shift;
my $fafile = shift;

if( $help or not $fafile ) {
    &help();
    exit(1);
}

sub help {
    print STDERR <<EOF;

$0: search a DNA fasta file against Rfam

Usage: $0 <options> cm_file fasta_file
    Options
        -h              : show this help
	-o <file>       : write the output to <file>

    Expert options
	-t <bits>       : specify cutoff in bits
	--bt <bits>     : specify blast evalue cutoff
	--local         : perform local mode search
	--global        : perform global mode search
	--nobig         : skip the large ribosomal RNAs
	--exclude [acc] : exlude family [acc] from the search
	--filter [ncbi|wu] : use wublast/ncbiblast (default ncbi)
        
EOF
}






msg("Scanning for ncRNAs... please be patient.");
my $num_ncrna = 0;
my $tool = "Infernal";
my $cmd = "cmscan --cpu 64 -E $evalueThresh --tblout /dev/stdout -o /dev/null --noali $cmfile $fafile";
msg("Running: $cmd");
open INFERNAL, '-|', $cmd;
my @f;
while (<INFERNAL>) {
  next if m/^#/;       # ignore comments
  my @x = split ' ';   # magic Perl whitespace splitter
#    msg("DEBUG: ", join("~~~", @x) );
  next unless @x > 9;  # avoid incorrect lines
  next unless defined $x[1] and $x[1] =~ m/^RF\d/;
  my $sid = $x[2];
  next unless exists $seq{$sid};
  push @{$seq{$sid}{FEATURE}}, Bio::SeqFeature::Generic->new( 
	-primary    => 'misc_RNA',
	-seq_id     => $sid,
	-source     => $tool,
	-start      => min($x[7], $x[8]),
	-end        => max($x[7], $x[8]),
	-strand     => ($x[9] eq '-' ? -1 : +1),
	-score      => undef,  # possibly x[16] but had problems here with '!'
	-frame      => 0,
	-tag        => {
	  'product' => $x[0],
	  'inference' => "COORDINATES:profile:$tool",
	}
  );
  $num_ncrna++;    
  msg("$num_ncrna ncRNA $x[0] $sid $x[7]..$x[8]");
} 
msg("Found $num_ncrna ncRNAs.");


