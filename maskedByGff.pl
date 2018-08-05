#!/usr/bin/perl
use strict;
my $usage = "USAGE:\n\tperl maskedByGff.pl file.gff genome.fasta\n";
unless (@ARGV) { print $usage; exit }

open GFF, '<', $ARGV[0];
open FASTA, '<', $ARGV[1];

my $fasta_head;
my %fasta;
while (<FASTA>) {
	chomp;
	if (/^>/) {
		$fasta_head = $_;
		$fasta_head =~ s/>//;
	}else {
		$fasta{$fasta_head} .= $_;
	}
}

while (<GFF>) {
	if (/(.+)\t.+\t.+\t(.+)\t(.+)\t.+\t.+\t.+\t.+/) {
		my $id = $1;
		my $start = $2 - 1;
		my $length = $3 - $start;
		substr($fasta{$id},$start,$length) =~ s/\w/N/g;
	}
}

my @fasta_head = keys %fasta;
my %fasta_head;
foreach (@fasta_head) {
	if (/(\d+)/) { $fasta_head{$_} = $1 }
}

@fasta_head = sort { $fasta_head{$a} <=> $fasta_head{$b} } @fasta_head;

foreach (@fasta_head) {
	print ">$_\n$fasta{$_}\n";
}			
