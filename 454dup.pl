#!/usr/bin/perl -w
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use Benchmark;
my $time1=new Benchmark;
my $program = `basename $0`;
chomp $program;
my $version = "0.1";
my %opts;
GetOptions(\%opts,"i:s","o:s","q:s","s:s","t:i","g:i","h");
sub usage{
	print CYAN << "EOF";

*************************************************************************

USAGE:
	$program [[-i <fasta_in>][-q <qual_in>]or[-s <sff_file>] ] [-o <out_prefix>] [-t threads] [-g strong] [-h]

OPTIONS:
	-i	Set the FASTA file name. [use sffino -s [-n] ]
	-q	Set the QUAL file name. [use sffinfo -q [-n] ]
	-s	Set the 454_sff file.
	-o	Set the output file prefix.
	-t	Set the multi-thread in perl. default [4]
	-g	Set the base length in temp file name.default [2]
	-h	Print this help information

VERSION: $version

COMPILE: 2010.12.5

DEBUG:	Tongwu Zhang <zhangtw\@big.ac.cn>

*************************************************************************

EOF

exit(1);
}

&usage if($opts{h});
&usage if(!((($opts{i} && $opts{q} )|| $opts{s}) && $opts{o}));

if(($opts{i} || $opts{q}) &&  $opts{s}){
	print STDERR "Plase use one type data: sff file or FASTA+QUAL\n";
	&usage;
}

print "*"x20,"BEGIN","*"x20,"\n";
my $input_fasta = $opts{i};
my $input_qual = $opts{q};
my $output = $opts{o};
my $sff_file = $opts{s};
my $parallel_num=$opts{t}||4;
my $strong=$opts{g}||2;
#my $maxmis=$opts{m}||1;
my @acgt=qw(A C G T);
my @filepad = &Tandembase($strong);

if( $sff_file ) 
{
	open(FASTA,"sffinfo -s $sff_file|")||die"Cann't open the file $sff_file or sffinfo error\n$!";
	open(QUAL,"sffinfo -q $sff_file|")||die"Cann't open the file $sff_file or sffinfo error\n$!";
}else{
	open(FASTA,"$input_fasta")||die"Cann't open the file $input_fasta\n$!";
	open(QUAL,"$input_qual")||die"Cann't open the file $input_fasta\n$!"; 
}

my @filehands;

system("rm -rf ${output}_TMPDONEFILE");
mkdir("${output}_TMPDONEFILE",0744)||die"error: Cann't creat dir$!";
my $pad;
foreach $pad(@filepad){
	my $filetmp1 = "${output}_TMPDONEFILE/".$pad.".a";
	my $filehand;
	open($filehand,">$filetmp1")||die "Cann't creat file $filetmp1 to write!$!\n";
	push(@filehands,$filehand);
}


my $sequence_str = "";
my $qual_str = "" ;
my $sequence_buf = <FASTA>;
my $qual_buf = <QUAL>;
my $sequence_name;
my $qual_name;
my $sequence_ref;
my $qual_ref;
my $read_num = 0;
my $seq_len;
my $total=0;

## Read a phrap Q file ########################
while((!eof(FASTA)) && (!eof(QUAL))){

	$sequence_name = $sequence_buf;
	chomp $sequence_name;
	$qual_name = $qual_buf;

	$read_num++;

SEQ:	$sequence_buf = <FASTA>;
	if(!eof(FASTA) && $sequence_buf!~m/^>/)
	{
		chomp $sequence_buf;
		$sequence_buf=~m/^\s*(.*)\s*$/;
		$sequence_str.=$&;
		goto SEQ;
	}else{
#		print PHRAP "$sequence_str\n";
		$sequence_ref = $sequence_str;
		$sequence_str = "";
	}

QUAL:	$qual_buf = <QUAL>;
	if(!eof(QUAL) && $qual_buf!~m/^>/)
	{
		chomp $qual_buf;
		$qual_buf=~m/^\s*(.*)\s*$/;
		
		$qual_str.=" ".$&;
		goto QUAL;
	}else {
		$qual_str=~s/^\s*(.*)\s*/$1/;
		$qual_ref = $qual_str;
		$qual_str = "";

	}
	$sequence_name=~s/>(\S+).*/$1/;
	for(my $i=0;$i<=$#filepad;$i++){
		$pad=$filepad[$i];
		if($sequence_ref=~/^$pad/){
			my $seqchar= &SeqTranstoChar($sequence_ref);
			my $tmphand=$filehands[$i];
			print $tmphand "$seqchar\t$sequence_name\t$qual_ref\n";
#	print "QUAL:$seqchar\t$sequence_name\t$qual_ref\n";
		}
	}

	$total++;

}

close FASTA;
close QUAL;

########## Close *.a #############
close $_ foreach(@filehands);
#### sort with parallel #####
my @dir=glob("${output}_TMPDONEFILE/*.*");
my @ndir=map($_.".sort",@dir);
for(my $i=0;$i<=$#dir;){
	my $parallel=0;
	my @pida=();
	my $pid;
	while($parallel<$parallel_num && $i<=$#dir){
	defined($pida[$parallel]=fork) or die "Can not fork:$!\n";
	unless($pida[$parallel]){
		exec("sort $dir[$i] >$ndir[$i]");
		}
	$parallel++;
	$i++;
	}
	foreach my $pid(@pida){ waitpid($pid,0);}
}

## sort -m with parallel ####

system("sort -m ${output}_TMPDONEFILE/*.a.sort -o ${output}_TMPDONEFILE/${output}.se");

map(unlink($_), @dir);
map(unlink($_),@ndir);

open(SE,"${output}_TMPDONEFILE/${output}.se")||die"Cann't open the sort file: ${output}.se\n$!";
open(DG,">${output}.dg")||die"Cann't open the output file: ${output}.rdp\n$!";
open(DS1,">${output}.ds1")||die"Cann't open the output file: ${output}.dup\n$!";
open(DS2,">${output}.ds2")||die"Cann't open the output file: ${output}.dup\n$!";

my $Query;
my $seqbase;
my $numseq;
my $numseq_top;
my $qual;
my $name;
my $seqbase_top;

my @seqchars;
my @names;
my @quals;
my $dgroup=0;

my %duplicate;
my %misshash;


until(eof(SE)){
	$Query=<SE>;
	chomp $Query;	
	($numseq,$name,$qual)=split(/\t/,$Query);
	$seqbase=$numseq;
	$seqbase=~s/\d+//g;
	unless($seqbase_top)
	{
		push(@seqchars,$numseq);
		push(@names,$name);
		push(@quals,$qual);
		$seqbase_top=$seqbase;
		$numseq_top=$numseq;
		
	}else{
		
		if(($seqbase!~/^$seqbase_top/i && $seqbase_top!~/^$seqbase/i))
		{
			&ProcessDuplicate(\@seqchars,\@names,\@quals);
			@seqchars=();
			@names=();
			@quals=();
			push(@seqchars,$numseq);
			push(@names,$name);
			push(@quals,$qual);
			$seqbase_top=$seqbase;
		}else{
			push(@seqchars,$numseq);
			push(@names,$name);
			push(@quals,$qual);
			$seqbase_top=$seqbase if(length($seqbase_top)<length($seqbase));
		}
	}

	&ProcessDuplicate(\@seqchars,\@names,\@quals) if(eof(SE));
}

my $uniqtype=0;
map($uniqtype+=$_,values %duplicate);
my $uniqread=$duplicate{1};
my $scaf_ratio=sprintf("%.2f",$uniqtype/$total);
my $start=sprintf("%.2f",$total/$uniqtype);
my $dup_ratio=sprintf("%.2f",($total-$uniqread)/$total);

print DS1 "#Total_reads=$total\n";
print DS1 "#Unique_types=$uniqtype\n";
print DS1 "#Unique_reads=$uniqread\n";
print DS1 "#Scaffolding_ratio=$scaf_ratio\n";
print DS1 "#Duplitcate_ratio=$dup_ratio\n";
print DS1 "#start_point=$start\n";
print DS1 "Duplicate\tFreq\n";

foreach(sort {$a<=>$b}  keys %duplicate){
	print DS1 "$_\t$duplicate{$_}\n";
}

print DS2 "Miss\tFreq\n";

foreach(sort {$a<=>$b} keys %misshash){
	print DS2 "$_\t$misshash{$_}\n";
}


open(R,">${output}.R")||die"Cann't open the output file: ${output}.dup\n$!";

print R "read.table(\"${output}.ds1\",header=T)->du\n";
print R "attach(du)\n";
print R "pdf(\"${output}.pdf\")\n";
#print R "jpeg(\"${output}_1.jpeg\",width=8000,height=6500,res=1200)\n";
print R "plot(Duplicate[Duplicate!=1],Freq[Duplicate!=1],xlim=c(0,max(Duplicate)),type=\"h\",lwd=1,col=\"green4\",main=\"Frequency of duplicate in 454 sequencing\", xlab=\"Read depth\",ylab=\"Frequency\")\n";
#print R "dev.off()\n";
print R "detach(du)\n";
print R "read.table(\"${output}.ds2\",header=T)->du2\n";
print R "attach(du2)\n";
print R "lar=sum(Freq[Miss>10])\n";
#print R "jpeg(\"${output}_2.jpeg\",width=8000,height=6500,res=1200)\n";
print R "plot(Miss[Miss!=0&Miss<=10],Freq[Miss!=0&Miss<=10],xlim=c(0,11),ylim=c(0,max(Freq[Miss!=0],lar)),type=\"h\",lwd=3,col=\"red3\",main=\"Frequency of miss base in duplicate reads of 454 sequencing\", xlab=\"Miss depth\",ylab=\"Frequency\")\n";
print R "lines(11,11,type=\"h\",col=\"blue3\",lwd=3)\n";
print R "axis(1,at=11,\">10\")\n";
print R "legend(\"topright\",inset=c(.05,.02),paste(\"Miss>=10: \",lar),bty=\"n\")\n";

print R "dev.off()\n";

close R;

system("R CMD BATCH ${output}.R");
system("rm -rf ${output}_TMPDONEFILE");

my $time2=new Benchmark;
my $timevar1=timediff($time2,$time1);

print "*"x20,"END","*"x20,"\n";
print timestr($timevar1),"\n";


sub ProcessDuplicate{
	my ($ref1,$ref2,$ref3) = @_;
	my @arry1=@$ref1;
	my @arry2=@$ref2;
	my @arry3=@$ref3;
	my @arry2_sort=();
	my $group =scalar @arry1;
	my $miss=0;
	
	$dgroup++;
	$duplicate{$group}=0 unless(exists $duplicate{$group});
	$duplicate{$group}++;
	$misshash{$miss}=0 unless(exists $misshash{$miss});
	if($group==1){
		$misshash{$miss}++;
		print DG ">dgroup_id: $dgroup\tmiss: $miss\n";
		print DG "$arry2[0] ";
		print DG "\n";
		return;
	}
	my %hashbase=();
	my %hashnum=();
	my $maxb=0;
	for(my $i=0;$i<$group;$i++)
	{
		@{$hashnum{$arry2[$i]}}=split(/[a-zA-Z]+/,$arry1[$i]);
		shift(@{$hashnum{$arry2[$i]}});
		$maxb=$maxb>scalar @{$hashnum{$arry2[$i]}}? $maxb : scalar @{$hashnum{$arry2[$i]}}; 
	}
	
	for(my $i=0;$i<$maxb;$i++)
	{
		my @nums=();
		foreach(@arry2){
			push(@nums,$hashnum{$_}->[$i]) if(exists $hashnum{$_}->[$i]);
		}
		$miss+=&MissCall(\@nums) if(scalar @nums >1);
	}
	
	$misshash{$miss}++;
	
	@arry2_sort=@arry2;
	for(my $i=0;$i<$group;$i++)
	{
		my @tmparry1=();
		@tmparry1=split(//,$arry2[$i]);
		my $len=scalar @tmparry1;
		for(my $j=$i+1;$j<$group;$j++)
		{
			my @tmparr2;
			@tmparr2=split(/\s+/,$arry3[$j]);
			my $len_top=scalar @tmparr2;
			if($len<$len_top)
			{
				my $tmp=$arry2_sort[$i];
				$arry2_sort[$i]=$arry2_sort[$j];
				$arry2_sort[$j]=$tmp;
			}elsif($len=$len_top){
				my $qua=0;
				my $qub=0;
				map($qua+=$_,split(/\s+/,$arry3[$i]));
				map($qub+=$_,split(/\s+/,$arry3[$j]));
				if($qua<$qub){
					my $tmp=$arry2_sort[$i];
					$arry2_sort[$i]=$arry2_sort[$j];
					$arry2_sort[$j]=$tmp;
				}
			}
		}
	}
########## print the dgroups, sort by length and quality ###########

	print DG ">dgroup_id: $dgroup\tmiss: $miss\n";
	foreach(@arry2_sort){ print DG "$_ ";}
	print DG "\n";
}



sub MissCall{
	my $ref=shift;
	my @list=sort @$ref;
	my $mid=int((scalar @list)/2);
	my $tot=0;
	map($tot+=abs($list[$mid] - $_),@list);
	my $tmp=scalar @list;
#print "$dgroup\t$tmp\t$list[0]\t$list[1]\t$tot\n";
	return $tot;
}


	
sub CharTranstoSeq{
	my $chars=shift;
	my @chararry=$chars=~/([ATCG])/g;
	my @num=$chars=~/(\d+)/g;
	my $seq="";
	for(my $i=0;$i<=$#chararry;$i++)
	{
		$seq.=($chararry[$i] x $num[$i]);
	}
	return $seq;
}

sub SeqTranstoChar{
	my $seq = uc(shift);
	$seq =~s/[^ATCG]/A/g;
	my @chars=split(//,$seq);
	my $seqchar="";
	for(my $i = 0;$i<=$#chars;$i++)
	{
		my $num = 1;
		while($i<$#chars && ($chars[$i] eq $chars[$i+1]))
		{
			$num++;
			$i++;
			last if($i==$#chars);
		}
		$seqchar.=($chars[$i].$num);
	}
	return $seqchar;
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

	
