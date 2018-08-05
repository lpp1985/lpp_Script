#!/usr/bin/perl -w
###_______________________________________________________________________________
###
###Program name: est_trimmer.pl
###Author:       Thomas Thiel
###Release date: 25/01/01
###Version:      04/09/02
###
###DESCRIPTION: tool for trimming EST (DNA) sequences concerning
###
###1.) Ambiguous sequences (everything apart from A,C,G,T).
###    All distal located bases containing equal or more than the predefined
###    number of ambiguous bases within a window of a defined size are deleted.
###
###2.) Removing of distal oligoN series (N stands for ONE of the characters ACGT)
###    from either the 5' or 3' end.
###    Originally thought to remove oligoT (oligoA) stretches from the 5' (3')
###    end corresponding to the bidirectional sequenced mRNA polyA tail, it is
###    also possible to remove other bases. Is done in two steps: First, all bases
###    of the set type (e.g. "T") are removed from one end (e.g. 5' end).
###    Secondly, because of potentially accumulating sequencing failures in the
###    case of the first or last positions (for very large runs) distal sequences
###    are removed containing a minimum of x repeats of the set type in a window
###    of predefined size. Setting the window size to 0 disables the second step.
###
###3.) Size cutoff.
###    All sequences being smaller than the min. size value are discarded.
###    Sequences being larger than the max. value are shortened from the 3' end.
###_______________________________________________________________________________
###
## _______________________________________________________________________________
## 
## DESCRIPTION: Tool for trimming EST (DNA) sequences
## 
## SYNTAX:   est_trimmer.pl <FASTAfile> [-amb=n,win] [-tr5=(A|C|G|T),n,win]
##                          [-tr3=(A|C|G|T),n,win] [-cut=min,max] [-id=name]
##                          [-help]
## 
##    <FASTAfile>    Single file in FASTA format containing the sequence(s).
##    [-amb=n,win]   Removes distal stretches containing "n" ambiguous bases in a
##                   "win" bp sized window.
##    [-tr5=N,n,win] Removes stretches of the given type N={A,C,G,T} from the 5'
##                   end. Value "n" defines the min. accepted repeat number of "N"
##                   in a 5' window of the size "win".
##    [-tr3=N,n,win] according to [-tr5] for the 3' end.
##    [-cut=min,max] Sets min. value for cutoff and max. sequence size.
##    [-id=name]     Optional. Final results are stored in "name".results, whereas
##                   processing steps are listed in "name".log. If not used,
##                   extensions are appended to <FASTAfile>.
##    [-help]        Further descriptions. Use "EST_trimmer.pl -help".
## 
##    Arguments can be used plurally and are processed according to their order.
## 
## EXAMPLE:  est_trimmer.pl ESTs -amb=2,50 -tr5=T,5,50 -tr3=A,5,50 -cut=100,700
## _______________________________________________________________________________
## 


# Check for arguments. If none display syntax #

if (@ARGV == 0)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Check if help is required #

if ($ARGV[0] =~ /-help/i)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\#\#(.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Open FASTA file #

open (IN,"<$ARGV[0]") || die ("\nError: File doesn't exist !\n\n");


# Checking arguments #

$arg = @ARGV;
$arg > 1 || die ("\nError: No arguments determined !\n\n");
for ($i = 1; $i < $arg; $i++)
    {
    if (($amb_n,$amb_win) = ($ARGV[$i] =~ /-amb=(\d+),(\d+)/i))
      {
      $message .= "$i. Check for ambiguous bases (search for $amb_n ambiguous bases in a $amb_win bp window).\n";
      }
    elsif (($tr3_b,$tr3_n,$tr3_win) = ($ARGV[$i] =~ /-tr3=([ACGT]),(\d+),(\d+)/i))
      {
      $message .= "$i. Trim 3' end: Remove \"$tr3_b\" tails.";
      if ($tr3_win > 0) {$message .= " Check for $tr3_n x $tr3_b in a $tr3_win bp window.\n"} else {$message .= "\n"};
      }
    elsif (($tr5_b,$tr5_n,$tr5_win) = ($ARGV[$i] =~ /-tr5=([ACGT]),(\d+),(\d+)/i))
      {
      $message .= "$i. Trim 5' end: Remove \"$tr5_b\" tails.";
      if ($tr5_win > 0) {$message .= " Check for $tr5_n x $tr5_b in a $tr5_win bp window.\n"} else {$message .= "\n"};
      }
    elsif (($cut_min,$cut_max) = ($ARGV[$i] =~ /-cut=(\d+),(\d+)/i))
      {
      $message .= "$i. Size restrictions: Size cutoff is $cut_min bp. Restrict sequence size to $cut_max bp.\n";
      }
    elsif (($file) = ($ARGV[$i] =~ /-id=(.*)/i))
      {
      }
    else {die ("\nError: Argument nr. ",++$i," is invalid !\n\n")}
    };

# Open results & log file#

if ($file) {open (OUT,">$file.results") || die ("\nError: Output file name not valid !\n\n");open (LOG,">$file.log")}
else {open (OUT,">$ARGV[0].results");open (LOG,">$ARGV[0].log")};

print "\nEST TRIMMER\n===========\n\nPerforming following steps:\n\n";
print LOG "\nLOG FILE OF ALL PERFORMED MODIFICATIONS\n";
print LOG "=======================================\n\n";
print LOG "\nPerforming following steps:\n---------------------------\n\n";
print $message;
print LOG $message;
print LOG "\n\nProcessing steps for each sequence:\n-----------------------------------\n";

# core #

$/ = ">";
while (<IN>)
  {
  next unless (($seqname,$seq) = /(.*?)\n(.*)/s);
  $seq =~ s/[\d\s>]//g;
  $seq =~ s/[^ACGTN]/N/gi; # only "acgtn" characters allowed
  $message = '';
  $discard = 'nö';
  for ($i = 1; $i < $arg; $i++)
    {
    if (($amb_n,$amb_win) = ($ARGV[$i] =~ /-amb=(\d+),(\d+)/i))
      {
      &Trim_5_n;
      &Trim_3_n;
      }
    elsif (($tr3_b,$tr3_n,$tr3_win) = ($ARGV[$i] =~ /-tr3=([ACGT]),(\d+),(\d+)/i))
      {
      &Trim_3_oligoN
      }
    elsif (($tr5_b,$tr5_n,$tr5_win) = ($ARGV[$i] =~ /-tr5=([ACGT]),(\d+),(\d+)/i))
      {
      &Trim_5_oligoN
      }
    elsif (($cut_min,$cut_max) = ($ARGV[$i] =~ /-cut=(\d+),(\d+)/i))
      {
      if ($cut_max < length $seq)
        {
        $seq = substr $seq,0,$cut_max;
        $message .= "Restrict sequence size to $cut_max bp.\n"
        };
      if ($cut_min and (length $seq < $cut_min))
        {
        $message .= "Discard sequence (below $cut_min bp cutoff).\n";
        $discard = "yo";
        }
      }
    };
  if ($discard eq 'nö')
    {
    $seq = join ("\n",grep($_,split(/(.{70})/,$seq)));
    print OUT ">$seqname\n$seq\n"
    };
  if ($message ne '') {print LOG "\n>$seqname\n$message"};
  };

print "\nDONE\n";

close (IN);
close (OUT);
close (LOG);


# subroutines #

sub Trim_5_n
  # trim 5' end until specified window contains less than n ambiguous bases
  {
  $check = '0';
  while ((length $seq > $amb_win) and !($check eq '1'))
    {
    $window = substr $seq,0,$amb_win;
    if ($window =~ /^(([^N]*N){$amb_n})/i)
      {
      $seq =~ s/^($1N*)//i;
      $message .= "Ambiguous sequence at 5' side: $1.\n"
      }
    else {$check = '1'};
    };
  };

sub Trim_3_n
  # trim 3' end until specified window contains less than n ambiguous bases
  {
  $check = '0';
  while ((length $seq > $amb_win) and !($check eq '1'))
    {
    $window = substr $seq,-$amb_win;
    if ($window =~ /.*((N[^N]*){$amb_n})/i)
      {
      $seq =~ s/(N*$1)$//i;
      $message .= "Ambiguous sequence at 3' side: $1.\n"
      }
    else {$check = '1'};
    };
  };
  
sub Trim_5_oligoN
  # remove oligoN stretches from 5' side
  {
  $check = '0';
  if ($seq =~ s/^($tr5_b+)//i) {$message .= "Remove \"$tr5_b\" tail at 5' side: $1.\n"};
  while (!($check eq '1'))
    {
    if (length $seq > $tr5_win) {$window = substr $seq,0,$tr5_win} else {$window = $seq};
    if ($window =~ /^(.*?($tr5_b){$tr5_n})/i)
      {
      $seq =~ s/^($1$tr5_b*)//i;
      $message .= "Remove \"$tr5_b\" stretch at 5' side: $1.\n"
      }
    else {$check = '1'};
    };
  };

sub Trim_3_oligoN
  # remove oligoN stretches from 3'side
  {
  $check = '0';
  if ($seq =~ s/($tr3_b+)$//i) {$message .= "Remove \"$tr3_b\" tail at 3' side: $1.\n"};
  while (!($check eq '1'))
    {
    if (length $seq > $tr3_win) {$window = substr $seq,-$tr3_win} else {$window = $seq};
    if ($window =~ /.*(($tr3_b){$tr3_n}.*)/i)
      {
      $seq =~ s/($tr3_b*$1)$//i;
      $message .= "Remove \"$tr3_b\" stretch at 3' side: $1.\n"
      }
    else {$check = '1'};
    }
  };
