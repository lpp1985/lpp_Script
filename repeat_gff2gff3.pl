#!/usr/bin/perl
use strict;

$_ =<>;
s/2/3/;
print;
while (<>) {
	s/Target (\"(\S+)\".*)/Name=\"$2\";Target=$1/;
	print 
}