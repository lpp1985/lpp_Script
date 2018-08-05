#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger;
use Bio::SFF::Reader::Sequential;

# process command line arguments
my ( $basename, $chunksize, $extension, $numchunks, $verbosity );
GetOptions(
    'basename=s'  => \$basename,
    'chunksize=i' => \$chunksize,
    'extension=s' => \$extension,
    'numchunks=i' => \$numchunks,    
    'verbose+'    => \$verbosity,
);

# create logging object
my $log = Bio::Phylo::Util::Logger->new(
    '-class' => 'main',
    '-level' => $verbosity,
);

# we can open the standard input stream as a handle to read from.
# this is useful because we can then expand a compressed archive
# on the fly and pipe it into our script
my $reader = Bio::SFF::Reader::Sequential->new( 'file' => \*STDIN );

# if specified how many chunks need to be produced, compute chunk
# size based on that
if ( $numchunks ) {
    my $numreads = $reader->header->number_of_reads;
    $chunksize = int( $numreads / $numchunks ) + 1;
}

# iterate over entries, print out as FASTA
my ( $readcounter, $chunkcounter, $filename ) = ( 0, 0 );
while ( my $entry = $reader->next_entry ) {
    
    # by default, STDOUT is the output handle
    my $fh = \*STDOUT;
    
    # if a base name and chunk size have been provided, create a new
    # file handle for each chuck
    if ( $basename and $chunksize and not $readcounter % $chunksize ) {
        $filename = sprintf "%s.%i.%s", $basename, ++$chunkcounter, $extension;
        close $fh;
        open $fh, "| gzip -9 > $filename" or die $!;
    }
    
    # print table of entries in files
    my $name = $entry->name;
    $log->info( $filename, "\t", $name, "\n" );
    
    # print FASTA output
    print $fh '>', $name, "\n", $entry->bases, "\n";
    $readcounter++;
}
