#!/usr/bin/perl -w
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

use strict;
use warnings;
use Benchmark;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

$Term::ANSIColor::AUTORESET = 1;
my $time1=new Benchmark;
my $PROGRAM_NAME=`basename $0`;
chomp $PROGRAM_NAME;
my $VERSION=0.7;

sub usage {
	print CYAN << "EOF";

*************************************************************

Usage:$PROGRAM_NAME [[-s <config>][-f <forward> and -r <reverse>][-i <fragment>] [-n] [-q] [-l] [-p] [-z] [-h]

Version: $VERSION
 
Options:

        -s solexa pair-end or fragment fastq library config file
	-f solexa pair-end fastq forward file
	-r solexa piar-end fastq reverse file
	-i solexa fragment fastq file
	-o prefix of the output file prefix
	-b output the duplicate-remove reads,default F
	-n number of processors to use, default 3
	-q number of base in filename,default 2 
	-p plot the duplication data using R package
	-l set minimal duplication length,default 15
	-e set the beging base for duplication,default 1
	-z compress the outcome file
	-d delete the temp file 
	-h print this help information

	Pair-end Config: 
	../s_1_1.txt  ../s_1_2.txt
	../s_2_1.txt  ../s_2_2.txt
	.....

	Fragement Config:
	../fastq/fragment_1.txt
	../fastq/fragment_2.txt
	.....

Bug_report: Tongwu Zhang <zhangtw\@big.ac.cn>

*************************************************************

EOF
exit(1);
}

my %opts;

GetOptions(\%opts,"s:s","f:s","r:s","i:s","o:s","n:i","l:i","e:i","d","z","b","q:i","h","p");

&usage if($opts{h});
if(!$opts{o}){
	print "\nError, please set the option with -o \n";
	&usage;
}


my $type;
my @files;

my $config_file;
my $fragment_file;

if($opts{s}){
	$config_file=$opts{s};
	if($opts{i}||$opts{f}||$opts{r}){
		print "Error, please use one of input format -s or -i or (-f and -r)!!\n";
		&usage;
	}
	open(C,"$config_file")||die"Cann't open the config file!\n$!";
	@files=<C>;
	chomp @files;
	if($files[0]=~/^\S+\s+\S+/){
		$type="paired-end";
	}elsif(	$files[0]=~/^\S+/)
	{
		$type="fragment";
	}else{
		print "Error, please check the config file:\n";
		&usage;
	}

}

if($opts{i}){
	$fragment_file=$opts{i};
	if($opts{s}||$opts{f}||$opts{r}){
		print "Error, please use one of input format -s or -i or (-f and -r)!!\n";
			&usage;
	}
	$type="fragment";
	push(@files,$fragment_file);
}

if($opts{f}||$opts{r}){
	if(!$opts{f} || !$opts{r}) 
	{
		print "Error, please input correct pair-end file with -f and -r\n";
		&usage;
	}
	push(@files,$opts{f}."\t".$opts{r});
	$type="paired-end"
}

if(!$config_file && !$fragment_file && !$opts{f})
{
	print "Error, please Check the input file format!!\n";
	&usage;
}


my $output=$opts{o};
my $read_length;
my $min_length=10000;
my $parallel_num=$opts{n}||5;
my $strong=$opts{q}||2;
my $min_read_length=$opts{l}||15;
my $compress=$opts{z}||0;
my $plot=$opts{p}||0;
my $delete=defined $opts{d} ?1:0;
my $outputboth=$opts{b}||0;
my $begin=$opts{e}||1;
my @filepad= &Tandembase($strong);
my $lane=0;
my %lanelen;
my %lanepair;
my @acgt=qw(A C G T);
print "*"x20,"BEGIN","*"x20,"\n";

############# main ############

if($type eq "paired-end"){

	&Paired_end_process;

}elsif( $type eq "fragment"){

	&Fragement_process;

}
############# over main ###############

################## if it is paired-end file ##########


sub Paired_end_process{

####################### Pair end file ################################
my @filehand1;
my @filehand2;
my @filehand3;
system("rm -rf ${output}_TMPDONEFILE");
mkdir("${output}_TMPDONEFILE",0744)||die"error: Cann't creat dir$!";
#my $pad;
@filehand1=map("${output}_TMPDONEFILE/".$_.".a", @filepad);
@filehand2=map("${output}_TMPDONEFILE/".$_.".b", @filepad);
@filehand3=map("${output}_TMPDONEFILE/".$_.".c", @filepad);

foreach my $pairfile (@files){
	my ($pairfile1,$pairfile2)=$pairfile=~/(\S+)\s+(\S+)/;
	$lane++;

open(F,"$pairfile1")||die "Cann't open forward file: $pairfile1\n$!";
open(R,"$pairfile2")||die "Cann't open reverase file: $pairfile2\n$!";
	
	my $pairfiletmp1=$pairfile1;
	my $pairfiletmp2=$pairfile2;
	$pairfiletmp1=~s/^.*\///;
	$pairfiletmp2=~s/^.*\///;
	$lanepair{$lane}->[0]=$pairfiletmp1;
	$lanepair{$lane}->[1]=$pairfiletmp2;

## Get the read length of the paired-end file

<F>;<R>;
my $read1;
my $read2;
$read1 = <F>;
$read2 = <R>;
if(length(chomp($read1)) != length(chomp($read2)))
{
	print "Warning: Sequences in file $pairfile1 is not equal to sequences in file $pairfile2\n";
#	exit(1);
}
$lanelen{$lane}->[0]=length($read1)-$begin+1;
$lanelen{$lane}->[1]=length($read2)-$begin+1;

$read_length = length($read1) <length($read2)? length($read1)-$begin+1:length($read2)-$begin+1;
$min_length= $read_length < $min_length ?  $read_length : $min_length;


close F;
close R;


############## SYSTEM COMMAND #####################
`awk 'BEGIN{ORS="";OFS="\\t";}{getline var1<"$pairfile1";if(NR%4==1){name=substr(\$1,1,length(\$1)-1);}else if(NR%4==2){seq1=substr(var1,$begin);seq2=substr(\$1,$begin);gsub(/N/,"A",seq1);gsub(/N/,"A",seq2);}else if(NR%4==0){qu1=substr(var1,$begin);qu2=substr(\$1,$begin);print seq1" $lane\\n">>"${output}_TMPDONEFILE/tmp.a"; print seq2" $lane\\n">>"${output}_TMPDONEFILE/tmp.b";max=length(seq1)>length(seq2)?length(seq1):length(seq2);for(i=1;i<=max;i++){print substr(seq1,i,1)substr(seq2,i,1);}print " $lane\\t"name,qu1,qu2"\\n"}}END{close("$pairfile1")}' $pairfile2 >>${output}_TMPDONEFILE/tmp.c`;

}
############## Over each input file in @files ################



########### Close *.a *.b *.c ###################
#### sort with parallel #####
my @awkpad=(@filepad,@filepad,@filepad);
my @dir=(@filehand1,@filehand2,@filehand3);
my @ndir;
@ndir=map($_.".sort",@dir);
map(s/\w+(\.[abc])$/tmp$1/,@dir);
for(my $i=0;$i<=$#dir;){
	my $parallel=0;
	my @pida=();
	my $pid;
	while($parallel<$parallel_num && $i<=$#dir){
		defined($pida[$parallel]=fork) or die "Can not fork:$!\n";
		unless($pida[$parallel]){
			exec("awk '/^$awkpad[$i]/{print \$0 | \"sort >$ndir[$i]\"}' $dir[$i] ");
		}
		$parallel++;
		$i++;
	}
	foreach my $pid(@pida){ waitpid($pid,0);}

}

## sort -m with parallel ####

defined(my $pid1=fork) or die "Can not fork:$!\n";
unless($pid1){
	exec("sort -m ${output}_TMPDONEFILE/*.a.sort -o ${output}_TMPDONEFILE/single1.se");
}

defined(my $pid2=fork) or die "Can not fork:$!\n";
unless($pid2){
	exec("sort -m ${output}_TMPDONEFILE/*.b.sort -o ${output}_TMPDONEFILE/single2.se");
}

defined(my $pid3=fork) or die "Can not fork:$!\n";
unless($pid3){
	exec("sort -m ${output}_TMPDONEFILE/*.c.sort -o ${output}_TMPDONEFILE/pair.se");
}

waitpid($pid1,0);
waitpid($pid2,0);
waitpid($pid3,0);

#####  Delete file *.a *.b *.c ######

unlink("tmp.a");
unlink("tmp.b");
unlink("tmp.c");
map(unlink($_),@dir);
map(unlink($_),@ndir);
########################
#### Do duplication Analysis #####

#my $name1 = $pairfile1;
#my $name2 = $pairfile2;
#my $name3 = $output;
pipe(RD1,WT1)||die "Can not pipe$!\n";
pipe(RD2,WT2)||die "Can not pipe$!\n";
pipe(RD3,WT3)||die "Can not pipe$!\n";

### single1.se ###
defined($pid1=fork) or die "Can not fork:$!\n";
unless($pid1){
	close RD1;
	open(SING1,"${output}_TMPDONEFILE/single1.se")||die "Cann't not open file single1.se$!\n";
	my %string_top=();
	my %sam=();
	my $string_new;
	my %dupS1HASH;
	my %total;
	until (eof(SING1))
	{
		chomp (my $query=<SING1>);
		my ($string,$lane)=$query=~/^(\S+)\s+(\S+)/;
		$total{$lane}=0 unless(exists $total{$lane});
		$total{$lane}++;
		my $lesslen=$lanelen{$lane}[0]<$lanelen{$lane}[1]?$lanelen{$lane}[0]:$lanelen{$lane}[1];

		for(my $base_length=$min_read_length;$base_length<=$lesslen;$base_length+=5)
		{
			$string_new=substr($string,0,$base_length);
			if(!exists  $sam{$lane}->{$base_length}){ $sam{$lane}->{$base_length}=1;}
			if(!exists $string_top{$lane}->{$base_length}){ $string_top{$lane}->{$base_length}=$string_new;next;}
			if($string_new ne $string_top{$lane}->{$base_length})
			{	
				$dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
				$dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				$sam{$lane}->{$base_length}=1;
				$string_top{$lane}->{$base_length}=$string_new;
			}else{
				$sam{$lane}->{$base_length}++;

			}
			if(eof(SING1)){
				foreach $lane (keys %lanelen){
			$dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
			$dupS1HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				}
			}


		}
		
	}

	close SING1;
	open(OS1,">${output}.s1")||die"Cannot open file$!\n";
	print OS1 "lane\tdup";
	for(my $i=$min_read_length;$i<=$min_length;$i+=5){print OS1 "\tL$i";}
	print OS1 "\n";
foreach my $lanet(sort {$a<=>$b} keys %lanelen){

	my $name1= $lanepair{$lanet}[0];
	my %tmp_hash=();
	my $lesslen=$lanelen{$lanet}[0]<$lanelen{$lanet}[1]?$lanelen{$lanet}[0]:$lanelen{$lanet}[1];
	for(my $i=$min_read_length;$i<=$lesslen;$i+=5)
	{	
		my $unitype=0;
		my $uniread = exists $dupS1HASH{$lanet}->{$i}->{1} ? $dupS1HASH{$lanet}->{$i}->{1} : 0;
		map($unitype+=$_, values %{$dupS1HASH{$lanet}->{$i}});
		my $usage_ratio = sprintf("%.2f",100*$unitype / $total{$lanet});
		my $start = sprintf("%.2f",$total{$lanet}/$unitype);
		my $dup_ratio = sprintf("%.2f",100*($total{$lanet}-$uniread)/$total{$lanet});
		print WT1 "$lanet\t$name1\t$i\t$total{$lanet}\t$unitype\t$usage_ratio\t$start\t$uniread\t$dup_ratio\n";
		foreach my $t (keys %{$dupS1HASH{$lanet}->{$i}}){ $tmp_hash{$t}=1;}
	}
	
	foreach(sort {$a<=>$b} keys %tmp_hash)
	{
		print OS1 "$lanet\t$_\t";
		for(my $y=$min_read_length;$y<=$min_length;$y+=5)
		{
			if(exists $dupS1HASH{$lanet}->{$y}->{$_}){
				print OS1 $dupS1HASH{$lanet}->{$y}->{$_};
			}else{
				print OS1 "0";
			}
			print OS1 "\t";
		}
		print OS1 "\n";
	}
}
	
	close OS1;
	exit 0;
}

close WT1;

### single2.se ###
defined($pid2=fork) or die "Can not for: $!\n";
unless ($pid2){
	close RD2;
	open(SING2,"${output}_TMPDONEFILE/single2.se")||die "Cann't not open file single2.se$!\n";
	my %dupS2HASH;
	my %string_top=();
	my %sam=();
	my $string_new;
	my %total;
	until( eof(SING2))
	{
		chomp (my $query=<SING2>);
		my ($string, $lane)=$query=~/^(\S+)\s+(\S+)/;
		$total{$lane}=0 unless(exists $total{$lane});
		$total{$lane}++;
		my $lesslen=$lanelen{$lane}[0]<$lanelen{$lane}[1]?$lanelen{$lane}[0]:$lanelen{$lane}[1];
		for(my $base_length=$min_read_length;$base_length<=$lesslen;$base_length+=5)
		{
			$string_new=substr($string,0,$base_length);
			if(!exists $sam{$lane}->{$base_length}){ $sam{$lane}->{$base_length}=1;}
			if(!exists $string_top{$lane}->{$base_length}){$string_top{$lane}->{$base_length}=$string_new;next;}
			if($string_new ne $string_top{$lane}->{$base_length})
			{
				$dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
				$dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				$sam{$lane}->{$base_length}=1;
				$string_top{$lane}->{$base_length}=$string_new;
			}else{
				$sam{$lane}->{$base_length}++;
			}
			if(eof(SING2)){
				foreach $lane (keys %lanelen){
				$dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
				$dupS2HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				}
			}
		}
	}
	close SING2;
	open(OS2,">${output}.s2")||die"Can not open file$!\n";
	print OS2 "lane\tdup";
	for(my $i=$min_read_length;$i<=$min_length;$i+=5){print OS2 "\tL$i";}
	print OS2 "\n";
foreach my $lanet (sort {$a<=>$b} keys %lanelen){
	my $name2=$lanepair{$lanet}[1];
	my %tmp_hash=();
	my $lesslen=$lanelen{$lanet}[0]<$lanelen{$lanet}[1]?$lanelen{$lanet}[0]:$lanelen{$lanet}[1];
	for(my $i=$min_read_length;$i<=$lesslen;$i+=5)
	{
		
		my $unitype=0;
		my $uniread = exists $dupS2HASH{$lanet}->{$i}->{1} ? $dupS2HASH{$lanet}->{$i}->{1} : 0;
		map($unitype += $_, values %{$dupS2HASH{$lanet}->{$i}});
		my $usage_ratio = sprintf("%.2f",100*$unitype / $total{$lanet});
		my $start = sprintf("%.2f",$total{$lanet}/$unitype);
		my $dup_ratio = sprintf("%.2f",100*($total{$lanet}-$uniread)/$total{$lanet});
		print WT2 "$lanet\t$name2\t$i\t$total{$lanet}\t$unitype\t$usage_ratio\t$start\t$uniread\t$dup_ratio\n";
		foreach my $t (keys %{$dupS2HASH{$lanet}->{$i}}){ $tmp_hash{$t}=1;}
	}
	foreach(sort {$a<=>$b} keys %tmp_hash)
	{
		print OS2 "$lanet\t$_\t";
		for(my $y=$min_read_length;$y<=$min_length;$y+=5)
		{
			if(exists $dupS2HASH{$lanet}->{$y}->{$_}){
				print OS2 $dupS2HASH{$lanet}->{$y}->{$_};
			}else{
				print OS2 "0";
			}
			print OS2 "\t";
		}
		print OS2 "\n";
	}
}
	
	close OS2;
	exit 0;
}
close WT2;
############ pair.se ############
#print "lane\t$lane\n";

defined($pid3=fork) or die "Can not for: $!\n";
unless ($pid3){
	close RD3;
	open(PAIR,"${output}_TMPDONEFILE/pair.se")||die "Cann't not open file pair.se$!\n";
	my %dupS3HASH;
	my %string_top=();
	my %cstring_top=();
	my %sam=();
	my %csam=();
	my $string_new;
	my %total;
	my %ctotal;
	my %CdupSHASH;

	my $bseq_top="";
	my $bseq;
	my $base_name;
	my $qual1_top;
	my $qual2_top;
	my $qual1;
	my $qual2;
	my $seq1;
	my $seq2;
	my $qual1_bak;
	my $qual2_bak;

	if($outputboth){
	open(HQ1,">${output}_1.rdp")||die"Cann't open the file to write:$!\n";
	open(HQ2,">${output}_2.rdp")||die"Cann't open the file to write:$!\n";
	open(CP,">${output}.dup")||die"Cann't open the file to write:$!\n";
	}
	until( eof(PAIR))
	{
		chomp (my $string=<PAIR>);
		($bseq,$lane,$base_name,$qual1,$qual2)=split(/\s+/,$string);
		$total{$lane}=0 unless(exists $total{$lane});
		$total{$lane}++;
	
if($outputboth){
		unless($bseq_top)
		{
			$bseq_top=$bseq;
			$qual1_top=$qual1;
			$qual2_top=$qual2;
			$qual1_bak=$qual1;
			$qual2_bak=$qual2;
		}else{
			($seq1,$seq2)= &SplitString($bseq_top,length($qual1_top),length($qual2_top));
			if($bseq!~/^$bseq_top/i && $bseq_top!~/^$bseq/i)
			{
				print HQ1 "${base_name}1\n$seq1\n+\n$qual1_top\n";
				print HQ2 "${base_name}2\n$seq2\n+\n$qual2_top\n";
				$bseq_top=$bseq;
				$qual1_top=$qual1;
				$qual2_top=$qual2;
				$qual1_bak=$qual1;
				$qual2_bak=$qual2;
			}else{
				print CP "${base_name}1\n$seq1\n+\n$qual1_bak\n";
				print CP "${base_name}2\n$seq2\n+\n$qual2_bak\n";
				$qual1_top = &CompareQuality($qual1_top,$qual1);
				$qual2_top = &CompareQuality($qual2_top,$qual1);
				$bseq_top =length($bseq_top)>length($bseq)?$bseq_top:$bseq; 
				$qual1_bak=$qual1;
				$qual2_bak=$qual2;
				
			}
		}
		if(eof(PAIR)){
			($seq1,$seq2)= &SplitString($bseq_top,length($qual1_top),length($qual2_top));
			print HQ1 "${base_name}1\n$seq1\n+\n$qual1_top\n";
			print HQ2 "${base_name}2\n$seq2\n+\n$qual2_top\n";
		}
}
		my $lesslen=$lanelen{$lane}[0]<$lanelen{$lane}[1]?$lanelen{$lane}[0]:$lanelen{$lane}[1];
		for(my $base_length=2*$min_read_length;$base_length<=2*$lesslen;$base_length+=10)
		{
			$string_new=substr($bseq,0,$base_length);
			if(!exists $sam{$lane}->{$base_length}){ $sam{$lane}->{$base_length}=1;}
			if(!exists $string_top{$lane}->{$base_length}){$string_top{$lane}->{$base_length}=$string_new;next;}
			if($string_new ne $string_top{$lane}->{$base_length})
			{
				$dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
				$dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				$sam{$lane}->{$base_length}=1;
				$string_top{$lane}->{$base_length}=$string_new;
			}else{
				$sam{$lane}->{$base_length}++;
			}
			if(eof(PAIR)){
				foreach my $lanetmp (keys %lanelen){
				next unless(exists $dupS3HASH{$lanetmp}->{$base_length});
				$dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}}=0 if(!exists $dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}});
				$dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}}++;
				}
			}
		}
##########################################################################
#$string_new;

	for(my $clane=$lane;$clane<=scalar keys %lanelen;$clane++)
	{
		 $ctotal{$clane}=0 unless(exists $ctotal{$clane});
		 $ctotal{$clane}++;
		for(my $base_length=2*$min_read_length;$base_length<=2*$lesslen;$base_length+=10)
		{
			$string_new=substr($bseq,0,$base_length);
			if(!exists $csam{$clane}->{$base_length}){ $csam{$clane}->{$base_length}=1;}
			if(!exists $cstring_top{$clane}->{$base_length}){$cstring_top{$clane}->{$base_length}=$string_new;next;}
			if($string_new!~/$cstring_top{$clane}->{$base_length}/ && $cstring_top{$clane}->{$base_length}!~/$string_new/)
			{
				$CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}}=0 if(!exists $CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}});
				$CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}}++;
				$csam{$clane}->{$base_length}=1;
				$cstring_top{$clane}->{$base_length}=$string_new;
			}else{
				$csam{$clane}->{$base_length}++;
			}
			if(eof(PAIR) && $clane==scalar keys %lanelen){
				foreach my $clanetmp(keys %lanelen) {
				next unless(exists $CdupSHASH{$clanetmp}->{$base_length});
				$CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}}=0 if(!exists $CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}});
				$CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}}++;
				}
			}
		}

	}
	
######################################################3
	}
	close PAIR;

	if($outputboth){close HQ1;close HQ2;close CP;}
	open(OP,">${output}.p")||die"Can not open file$!\n";

	print OP "lane\tdup";
	for(my $i=$min_read_length;$i<=$min_length;$i+=5){print OP "\tL$i";}
	print OP "\n";

	open(OUT2,">${output}.du2")||die "Can not open file$!\n";
	print OUT2 "Lane\tName\tLength\tTotal\tUnitype\tUsage_ratio\tStart\tUniread\tDup_ratio\tNew\tNewPer\n";

foreach my $lanet (sort {$a<=>$b} keys %lanelen){

	my %tmp_hash=();
## Change this name to output ####
	my $name3;
	if($lanepair{$lanet}[0]=~/^(\S+)[_\.][1a].*/){
		$name3=$1;
	}else{
		$name3=$output;  }
#	print "$lanet\t$lanepair{$lanet}[0]\t$name3\n";
	my $lesslen=$lanelen{$lanet}[0]<$lanelen{$lanet}[1]?$lanelen{$lanet}[0]:$lanelen{$lanet}[1];
	for(my $i=2*$min_read_length;$i<=$lesslen*2;$i+=10)
	{
		my $tmlen=$i/2;
		my $unitype=0;
		my $uniread = exists $dupS3HASH{$lanet}->{$i}->{1} ? $dupS3HASH{$lanet}->{$i}->{1} : 0;
		map($unitype +=$_, values %{$dupS3HASH{$lanet}->{$i}});
		my $usage_ratio = sprintf("%.2f",100*$unitype / $total{$lanet});
		my $start = sprintf("%.2f",$total{$lanet}/$unitype);
		my $dup_ratio = sprintf("%.2f",100*($total{$lanet}-$uniread)/$total{$lanet});
		print WT3 "$lanet\t$name3\t$tmlen\t$total{$lanet}\t$unitype\t$usage_ratio\t$start\t$uniread\t$dup_ratio\n";
		foreach my $t (keys %{$dupS3HASH{$lanet}->{$i}}){ $tmp_hash{$t}=1;}
###################
		my $cunitype=0;
		my $cuniread = exists $CdupSHASH{$lanet}->{$i}->{1} ? $CdupSHASH{$lanet}->{$i}->{1} : 0;
		map($cunitype += $_, values %{$CdupSHASH{$lanet}->{$i}});
		my $cusage_ratio = sprintf("%.2f",100*$cunitype / $ctotal{$lanet});
		my $cstart = sprintf("%.2f",$ctotal{$lanet}/$cunitype);
		my $cdup_ratio = sprintf("%.2f",100*($ctotal{$lanet}-$cuniread)/$ctotal{$lanet});
		my $old_cunitype=0;
		if($lanet>1){
			map($old_cunitype += $_,values %{$CdupSHASH{$lanet-1}->{$i}});
		}
		my $cnew = $cunitype - $old_cunitype;
		my $cnewper = sprintf("%.2f",100*$cnew/$ctotal{$lanet});
		print OUT2 "$lanet\t$name3\t$tmlen\t$ctotal{$lanet}\t$cunitype\t$cusage_ratio\t$cstart\t$cuniread\t$cdup_ratio\t$cnew\t$cnewper\n";	
	
###################
	}
	foreach(sort {$a<=>$b} keys %tmp_hash)
	{
		print OP "$lanet\t$_\t";
		for(my $y=2*$min_read_length;$y<=2*$min_length;$y+=10)
		{
			if(exists $dupS3HASH{$lanet}->{$y}->{$_}){
				print OP $dupS3HASH{$lanet}->{$y}->{$_};
			}else{print OP "0";
			}
			print OP "\t";
		}
		print OP "\n";
	}
}
	close OP;
	close OUT2;
	exit 0;
}
close WT3;
###################### over the parallel ################
waitpid($pid1,0);
waitpid($pid2,0);
waitpid($pid3,0);
######### Print the out of duplication ###########

open(OUT,">${output}.du1")||die "Can not open file$!\n";
print OUT "Lane\tName\tLength\tTotal\tUnitype\tUsage_ratio\tStart\tUniread\tDup_ratio\n";

until(eof(RD1)||eof(RD2)||eof(RD3)){
	my $line1=<RD1>;
	my $line2=<RD2>;
	my $line3=<RD3>;
	print OUT "$line1$line2$line3";
}

close RD1;
close RD2;
close RD3;



}
################## over paired_end_process ###################################################################


sub Fragement_process{

####################### Fragement file ################################
my @filehand1;
system("rm -rf ${output}_TMPDONEFILE");
mkdir("${output}_TMPDONEFILE",0744)||die"error: Cann't creat dir. $!\n";
#my $pad;
@filehand1=map("${output}_TMPDONEFILE/".$_.".c",@filepad);

foreach my $pairfile (@files){
	if(!($pairfile=~s/^\s*(\S+)\s*/$1/)){
		print "Error with the solexa filename\n";
		exit 1;
	}
	$lane++;
	my $pairfiletmp=$pairfile;
	$pairfiletmp=~s/^.*\///;
	$lanepair{$lane}=$pairfiletmp;

open(F,"$pairfile")||die "Cann't open forward file: $pairfile $!";
## Get the read length of the paired-end file

<F>;
my $read1;
$read1 = <F>;
$lanelen{$lane}=length($read1)-$begin+1;
$min_length= $min_length < length($read1) ? $min_length : length($read1)-$begin+1;

close F;

`awk 'BEGIN{ORS="";OFS="\\t";}{if(NR%4==1){name=\$1;}else if(NR%4==2){seq1=substr(\$1,$begin);gsub(/N/,"A",seq1);}else if(NR%4==0){qu1=substr(\$1,$begin);print seq1" $lane",name,qu1"\\n">>"${output}_TMPDONEFILE/tmp.c";}}' $pairfile` ;
## Fragement file seek to begin

}
############## Over each input file in @files ################

########### Close *.c ###################

#### sort with parallel #####
my @awkpad=@filepad;
my @dir=(@filehand1);
my @ndir;
@ndir=map($_.".sort",@dir);
map(s/\w+(\.c)$/tmp$1/,@dir);
for(my $i=0;$i<=$#dir;){
	my $parallel=0;
	my @pida=();
	my $pid;
	while($parallel<$parallel_num && $i<=$#dir){
		defined($pida[$parallel]=fork) or die "Can not fork:$!\n";
		unless($pida[$parallel]){
			exec("awk '/^$awkpad[$i]/{print \$0 | \"sort >$ndir[$i]\"}' $dir[$i] ");
		}
		$parallel++;
		$i++;
	}
	foreach my $pid(@pida){ waitpid($pid,0);}

}

## sort -m with parallel ####
system("sort -m ${output}_TMPDONEFILE/*.c.sort -o ${output}_TMPDONEFILE/fragment.se");


#####  Delete file *.c ######
unlink("${output}_TMPDONEFILE/tmp.c");
map(unlink($_),@dir);
map(unlink($_),@ndir);
########################
#### Do duplication Analysis #####

############ fragment.se ############
open(OUT,">${output}.du1")||die "Can not open file$!\n";
print OUT "Lane\tName\tLength\tTotal\tUnitype\tUsage_ratio\tStart\tUniread\tDup_ratio\n";
	open(PAIR,"${output}_TMPDONEFILE/fragment.se")||die "Cann't not open file pair.se\n$!";
	my %dupS3HASH;
	my %string_top=();
	my %cstring_top=();
	my %sam=();
	my %csam=();
	my $string_new;
	my %total;
	my %ctotal;
	my %CdupSHASH;

	my $bseq_top="";
	my $bseq;
	my $base_name;
	my $qual1_top;
	my $qual1;
	my $seq1;
	my $seq2;
	my $qual1_bak;

if($outputboth){
	open(HQ1,">${output}.rdp")||die"Cann't open the file to write:$!\n";
	open(CP,">${output}.dup")||die"Cann't open the file to write:$!\n";
}
	until( eof(PAIR))
	{
		chomp (my $string=<PAIR>);
		($bseq,$lane,$base_name,$qual1)=split(/\s+/,$string);
		$total{$lane}=0 unless(exists $total{$lane});
		$total{$lane}++;

if($outputboth){
		unless($bseq_top)
		{
			$bseq_top=$bseq;
			$qual1_top=$qual1;
			$qual1_bak=$qual1;
		}else{
			if($bseq!~/^$bseq_top/i && $bseq_top!~/^$bseq/i)
			{
				print HQ1 "$base_name\n$bseq_top\n+\n$qual1_top\n";
				$bseq_top=$bseq;
				$qual1_top=$qual1;
				$qual1_bak=$qual1;
			}else{
				print CP "$base_name\n$bseq\n+\n$qual1_bak\n";
				$qual1_top = &CompareQuality($qual1_top,$qual1);
				$bseq_top =length($bseq_top)>length($bseq)?$bseq_top:$bseq; 
				$qual1_bak=$qual1;
				
			}
		}
		if(eof(PAIR)){
			print HQ1 "$base_name\n$bseq_top\n+\n$qual1_top\n";
		}
}

		my $lesslen=$lanelen{$lane};
		for(my $base_length=$min_read_length;$base_length<=$lesslen;$base_length+=5)
		{
			$string_new=substr($bseq,0,$base_length);
			if(!exists $sam{$lane}->{$base_length}){ $sam{$lane}->{$base_length}=1;}
			if(!exists $string_top{$lane}->{$base_length}){$string_top{$lane}->{$base_length}=$string_new;next;}
			if($string_new ne $string_top{$lane}->{$base_length})
			{
				$dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}=0 if(!exists $dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}});
				$dupS3HASH{$lane}->{$base_length}->{$sam{$lane}->{$base_length}}++;
				$sam{$lane}->{$base_length}=1;
				$string_top{$lane}->{$base_length}=$string_new;
			}else{
				$sam{$lane}->{$base_length}++;
			}
			if(eof(PAIR)){
				foreach my $lanetmp (keys %lanelen){
				next unless( $dupS3HASH{$lanetmp}->{$base_length});
				$dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}}=0 if(!exists $dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}});
				$dupS3HASH{$lanetmp}->{$base_length}->{$sam{$lanetmp}->{$base_length}}++;
				}
			}
		}
##########################################################################
	for(my $clane=$lane;$clane<=scalar keys %lanelen;$clane++)
	{
		 $ctotal{$clane}=0 unless(exists $ctotal{$clane});
		 $ctotal{$clane}++;
		for(my $base_length=$min_read_length;$base_length<=$lesslen;$base_length+=5)
		{
			$string_new=substr($bseq,0,$base_length);
			if(!exists $csam{$clane}->{$base_length}){ $csam{$clane}->{$base_length}=1;}
			if(!exists $cstring_top{$clane}->{$base_length}){$cstring_top{$clane}->{$base_length}=$string_new;next;}
			if($string_new!~/$cstring_top{$clane}->{$base_length}/ && $cstring_top{$clane}->{$base_length}!~/$string_new/)
			{
				$CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}}=0 if(!exists $CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}});
				$CdupSHASH{$clane}->{$base_length}->{$csam{$clane}->{$base_length}}++;
				$csam{$clane}->{$base_length}=1;
				$cstring_top{$clane}->{$base_length}=$string_new;
			}else{
				$csam{$clane}->{$base_length}++;
			}
			if(eof(PAIR) && $clane==scalar keys %lanelen){
				foreach my $clanetmp (keys %lanelen){
				next unless(exists $CdupSHASH{$clanetmp}->{$base_length} );
				$CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}}=0 if(!exists $CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}});
				$CdupSHASH{$clanetmp}->{$base_length}->{$csam{$clanetmp}->{$base_length}}++;
				}
			}
		}

	}
	
######################################################3
	}
	close PAIR;
if($outputboth){ close HQ1; close CP; }
	open(OP,">${output}.p")||die"Can not open file\n$!";
	print OP "lane\tdup";
	for(my $i=$min_read_length;$i<=$min_length;$i+=5){ print OP "\tL$i";}
	print OP "\n";
	open(OUT2,">${output}.du2")||die "Can not open file\n$!";
	print OUT2 "Lane\tName\tLength\tTotal\tUnitype\tUsage_ratio\tStart\tUniread\tDup_ratio\tNew\tNewPer\n";

foreach my $lanet (sort {$a<=>$b} keys %lanelen){
	my %tmp_hash=();
## Change this name to output ####
	my $name3;
	if($lanepair{$lanet}=~/^(\S+)[_\.][1a].*/){
		$name3=$1;
	}else{
		$name3=$output;  }
	my $lesslen=$lanelen{$lanet};
	for(my $i=$min_read_length;$i<=$lesslen;$i+=5)
	{
		my $tmlen=$i;
		my $unitype=0;
		my $uniread = exists $dupS3HASH{$lanet}->{$i}->{1} ? $dupS3HASH{$lanet}->{$i}->{1} : 0;
		map($unitype += $_, values %{$dupS3HASH{$lanet}->{$i}}) ;
		my $usage_ratio = sprintf("%.2f",100*$unitype / $total{$lanet});
		my $start = sprintf("%.2f",$total{$lanet}/$unitype);
		my $dup_ratio = sprintf("%.2f",100*($total{$lanet}-$uniread)/$total{$lanet});
		print OUT "$lanet\t$name3\t$tmlen\t$total{$lanet}\t$unitype\t$usage_ratio\t$start\t$uniread\t$dup_ratio\n";
		foreach my $t (keys %{$dupS3HASH{$lanet}->{$i}}){ $tmp_hash{$t}=1;}
###################
		my $cunitype=0;
		my $cuniread = exists $CdupSHASH{$lanet}->{$i}->{1} ? $CdupSHASH{$lanet}->{$i}->{1} : 0;
		map($cunitype += $_, values %{$CdupSHASH{$lanet}->{$i}});
		my $cusage_ratio = sprintf("%.2f",100*$cunitype / $ctotal{$lanet});
		my $cstart = sprintf("%.2f",$ctotal{$lanet}/$cunitype);
		my $cdup_ratio = sprintf("%.2f",100*($ctotal{$lanet}-$cuniread)/$ctotal{$lanet});
		my $old_cunitype=0;
		if($lanet>1){
			map($old_cunitype += $_, values %{$CdupSHASH{$lanet-1}->{$i}});
		}
		my $cnew = $cunitype - $old_cunitype;
		my $cnewper = sprintf("%.2f",100*$cnew/$ctotal{$lanet});
		print OUT2 "$lanet\t$name3\t$tmlen\t$ctotal{$lanet}\t$cunitype\t$cusage_ratio\t$cstart\t$cuniread\t$cdup_ratio\t$cnew\t$cnewper\n";	
	
###################
	}
	foreach(sort {$a<=>$b} keys %tmp_hash)
	{
		print OP "$lanet\t$_\t";
		for(my $y=$min_read_length;$y<=$min_length;$y+=5)
		{
			if(exists $dupS3HASH{$lanet}->{$y}->{$_}){
				print OP $dupS3HASH{$lanet}->{$y}->{$_};
			}else{
				print OP "0";
			}
			print OP "\t";
		}
		print OP "\n";
	}
}
	close OP;
	close OUT2;
	close OUT;
}
################## over Fragement_process ###########################################################################


sub CombinedString{
	my ($string1,$string2,)=@_;
	my @chars1=();
	my @chars2=();
	@chars1=split(//,$string1);
	@chars2=split(//,$string2);
	my $maxlen=scalar @chars1 > scalar @chars2? scalar @chars1: scalar @chars2;
	my @chars=();
	for(my $i=0;$i<$maxlen;$i++)
	{	
		unless(exists $chars1[$i]){$chars[$i]=$chars2[$i];next;}
		unless(exists $chars2[$i]){$chars[$i]=$chars1[$i];next;}
		$chars[$i]=$chars1[$i].$chars2[$i];
	}
	my $join=join('',@chars);
	return $join;

}

sub SplitString{
	my ($string,$len1,$len2) = @_;
	my @chars=();
	@chars=split(//,$string);
	my $string1="";
	my $string2="";
	my @seq=();
	my $len=$len1>$len2? $len2:$len1;
	for(my $i=0;$i<2*$len;$i+=2)
	{
		$string1.=$chars[$i];
		$string2.=$chars[$i+1];
	}
	if($len1>$len2){
		$string1.=substr($string,-1,($len1-$len2));
	}elsif($len1<$len2){
		$string2.=substr($string,-1,($len2-$len1));
	}
	push(@seq,$string1);
	push(@seq,$string2);
	return @seq;
}

sub CompareQuality{
	my ($qual1,$qual2)=@_;
	my @qu1=();
	my @qu2=();
	my @qu=();
	my $qbest="";
	my $length= length($qual1) >length($qual2)?  length($qual1):length($qual2);
	@qu1=split(//,$qual1);
	@qu2=split(//,$qual2);
	for(my $i=0;$i<$length;$i++)
	{
		unless(exists $qu1[$i]){$qu[$i]=$qu2[$i];next;}
		unless(exists $qu2[$i]){$qu[$i]=$qu1[$i];next;}
		$qu[$i]=$qu1[$i] gt $qu2[$i] ? $qu1[$i]:$qu2[$i];

	}
	$qbest=join('',@qu);
	return $qbest;
}



sub Tandembase {
	(my $num)=@_;
	@acgt=qw(A C G T);
	if($num>1) {
		my @last=Tandembase ($num-1);
		my @tmp;
		$#tmp=-1;
		for(my$i=0;$i<=$#last;$i++){
			for(my $j=0;$j<=3;$j++){
				push( @tmp,$last[$i].$acgt[$j]);
			}
		}
		return @tmp;
	}
	if($num==1) {
		return @acgt;
	}
}



################## print int to $output.R file #################
my @outputlist;

for(my $lanetmp=1;$lanetmp<=$lane;$lanetmp++)
{
	push(@outputlist,$output."_".$lanetmp);
}

#############################################################################
open(R,">${output}.R")||die "cann't open the output file:${output}.R\n$!";
#############################################################################
print R "read.table(\"${output}.du1\",header=T)->du1\n";
print R "attach(du1)\n";
 
my $lanetmp=1;

	print R "pdf(\"${output}.pdf\")\n";	
foreach my $pairfile(@outputlist){
#	print R "jpeg(\"${pairfile}_a.jpeg\",width=8000,height=6000,res=1200)\n";
#	print R "par(cex=0.7)\n";
	print R "plot(range(Length),range(0,100),type=\"n\",xlab=\"Length (bp)\",ylab=\"Scaffolding Usage (%)\",main=\"PCR duplication in solexa library:${pairfile}\")\n";
	print R "abline( h = seq(0,by=10,to=100), col = \"gray\", lwd = .3 )\n";
	print R "co=1\n";
	print R "for (i in as.character(unique(Name[Lane==$lanetmp]))){ \n";
	print R "lines(Length[Name==i & Lane==$lanetmp ],Usage_ratio[Name==i & Lane==$lanetmp],col=co,pch=20,type=\"b\")\n";
	print R "co=co+1 }\n";
	print R "legend(\"bottomleft\",as.character(unique(Name[Lane==$lanetmp])),inset=c(.05,.02), col=c(1,2,3),lty=2,pch=20,bty=\"n\")\n";
#	print R "dev.off()\n";

	print R "read.table(\"${output}.p\",header=T)->da\n";
	print R "attach(da)\n";
	print R "as.matrix(da)->ma\n";
#	print R "jpeg(\"${pairfile}_b.jpeg\",width=8000,height=6000,res=1200)\n";
#	print R "par(cex=0.7)\n";
	print R "xmax=max(ma[,2][lane==$lanetmp])\n";
	print R "ymax=max(ma[lane==$lanetmp])\n";
#	print R "plot(c(0,xmax),c(0,ymax),type=\"n\",xlab=\"Duplicate time\",ylab=\"Frequence\",main=\"Different distribution of duplicate time with different length set\")\n";
	print R "plot(c(0,30),c(0,log2(ymax)),type=\"n\",xlab=\"Duplicate depth\",ylab=\"Frequency (log2)\",main=\"Different distribution of duplicate depth with different length set\")\n";
	print R "c=1/(length(ma[1,])-2);\n";
	print R "for( i in 3:(length(ma[1,])-2))\n";
	print R "{\n";
	print R "lines(ma[,2][lane==$lanetmp],log2(ma[,i][lane==$lanetmp]),col=rgb(0,0,1,i*c))\n";
	print R "}\n";
	print R "legend(\"topright\",paste(colnames(ma)[3],\"->\",colnames(ma)[length(ma[1,])]),inset=c(.05,.02),lty=1,col=\"blue\",bty=\"n\")\n";
#	print R "dev.off()\n";

	$lanetmp++;
}

if($lane >=2){

#print R "jpeg(\"${output}_c_1.jpeg\",width=8000,height=6000,res=1200)\n";
print R "read.table(\"${output}.du2\",header=T)->du2\n";
#print R"par(cex=0.7)\n";
print R "plot(range(unique(du2\$Lane)),range(0,100),type=\"n\",xlab=\"Lane (#)\",ylab=\"Cumulative Scaffolding Usage (%)\",main=\"PCR duplication in solexa library:test\")\n";
print R "abline(h = seq(0,by=10,to=100), col = \"gray\", lwd = .3)\n";
print R "co=1\n";
print R "for (i in seq(20,by=10,to=max(du2\$Length))){ \n";
print R "lines(du2\$Lane[du2\$Length==i],du2\$Usage_ratio[du2\$Length==i],col=co,pch=20,type=\"b\")\n";
print R "co=co+1\n";
print R "}\n";
print R "legend(\"bottomleft\",paste(seq(20,by=10,to=max(du2\$Length)),\"bp\"),inset=c(.05,.02), col=1:length(unique(du2\$Length)),lty=2,pch=20,bty=\"n\")\n";
#print R "dev.off()\n";

#print R "jpeg(\"${output}_c_2.jpeg\",width=8000,height=6000,res=1200)\n";
#print R "par(cex=0.7)\n";
print R "plot(range(unique(du2\$Lane)),range(0,100),type=\"n\",xlab=\"Lane (#)\",ylab=\"Add New Lane Usage (%)\",main=\"PCR duplication in solexa library:test\")\n";
print R "abline(h = seq(0,by=10,to=100), col = \"gray\", lwd = .3)\n";
print R "co=1\n";
print R "for (i in seq(20,by=10,to=max(du2\$Length))){\n";
print R "lines(du2\$Lane[du2\$Length==i],du2\$NewPer[du2\$Length==i],col=co,pch=20,type=\"b\")\n";
print R "co=co+1\n";
print R "}\n";
print R "legend(\"bottomleft\",paste(seq(20,by=10,to=max(du2\$Length)),\"bp\"),inset=c(.05,.02), col=1:length(unique(du2\$Length)),lty=2,pch=20,bty=\"n\")\n";
#print R "dev.off()\n";

}

print R "dev.off()\n";

if($plot ||$delete){ system( "R CMD BATCH ${output}.R");}

#################################################################################

system("rm -rf ${output}_TMPDONEFILE");
my $compressname=$output.".tar.gz" if($compress);

if($delete){
#system( "R CMD BATCH ${output}.R");
#	system(" tar -zcvf $compressname ${output}.pdf ${output}.du* ${output}.s* ${output}.p ${output}.R ${output}.Rout ") if($compress);
	system(" rm -f ${output}.du* ${output}.s* ${output}.p ${output}.Rout ${output}.R");
}else{
	system(" tar -zcvf $compressname ${output}.pdf ${output}.du* ${output}.s* ${output}.p ${output}.R ${output}.Rout ") if($compress && $plot);
	system(" tar -zcvf $compressname  ${output}.du* ${output}.s* ${output}.p ${output}.R ") if($compress && !$plot);
	
}

if($compress) {

	system(" rm -rf ${output}.pdf ${output}.du* ${output}.s* ${output}.p ${output}.R ${output}.Rout ");
}


################### Over all ###################
my $time2=new Benchmark;
my $timevar1=timediff($time2,$time1);

print "*"x20,"END","*"x20,"\n";
print timestr($timevar1),"\n";

