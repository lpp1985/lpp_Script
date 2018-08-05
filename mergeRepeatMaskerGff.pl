#!/usr/bin/perl
use strict;
my $usage = "USAGE:\n\tperl mergeRepeatMaskerGff.pl repeatmasker1.gff repeatmasker2.gff \nIt would be better to place the repeatMasker result behind...\n\n";

unless (@ARGV) { print $usage; exit }
my %gff;
my %start;
my %geneid;
while (<>) {
	chomp;
	next if /^#/;

	if (/^(.*)\t.*\t.*\t(.*\t.*)\t.*\t.*\t.*\t.*/) {
		my $sig = "$1\t$2";
		$gff{$sig} = $_;
	}
}

my @se = keys %gff;

foreach (@se) {
	if (/(.*)\t(\d+)\t\d+/) {
		$geneid{$_} = $1;
		$start{$_} = $2;
	}
}

@se = sort { $geneid{$a} cmp $geneid{$b} or $start{$a} <=> $start{$b} } @se;

print "##gff-version 2\n";
foreach (@se) {
	print "$gff{$_}\n";
}
