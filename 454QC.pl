#!/usr/bin/perl -w
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

use strict;
use warnings;
use Getopt::Long;
#use Dta::Dumper;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use Benchmark;
my $time1=new Benchmark;

my $program = `basename $0`;
chomp $program;
my $version = "1.0";
my %opts;
GetOptions(\%opts,"i:s","o:s","q:s","s:s","n:i","p","d","z","h");
sub usage{
	print CYAN << "EOF";

*************************************************************************

USAGE:
	$program [ [-i <fasta_in>] [-q <qual_in>] or [-s <sff_file>] ] [-o <out_prefix>] [-n number] [-p][-d][-z][-h]

OPTIONS:
	-i	Set the FASTA file name. [use sffino -s [-n] ]
	-o	Set the QUAL file name. [use sffinfo -q [-n] ]
	-s	Set the 454_sff file.
	-n	Set the nubmer of multiply bases(>=2). default [2]
	-p	Plot the QC, default [TURE]
	-d	Delete the QC file, just save the pictures(-p) [FALSE]
	-z	Compress the QC file and Picture default [FASLSE]
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

my $input_fasta = $opts{i};
my $input_qual = $opts{q};
my $output = $opts{o};
my $sff_file = $opts{s};
my $multiply = $opts{n} ||2;
my $plot = $opts{p} || 0;
my $compress = $opts{z} || 0;
my $delete = $opts{d} || 0; 

if( $sff_file ) 
{
	open(FASTA,"sffinfo -s $sff_file|")||die"Cann't open the file $sff_file or sffinfo error\n$!";
	open(QUAL,"sffinfo -q $sff_file|")||die"Cann't open the file $sff_file or sffinfo error\n$!";
}else{
	open(FASTA,"$input_fasta")||die"Cann't open the file $input_fasta\n$!";
	open(QUAL,"$input_qual")||die"Cann't open the file $input_fasta\n$!"; 
}

my $outfile1 = $output.".Q1";
my $outfile2 = $output.".Q2";
my $outfile3 = $output.".Q3";
my $outfile4 = $output.".Q4";

open(Q1, ">$outfile1")||die"Cann't open the file ($outfile1) to write\n$!";
open(Q2, ">$outfile2")||die"Cann't open the file ($outfile2) to write\n$!";
open(Q3, ">$outfile3")||die"Cann't open the file ($outfile3) to write\n$!";
open(Q4, ">$outfile4")||die"Cann't open the file ($outfile4) to write\n$!";

print Q1 "Readnum\tQmean\tLength\tGC\n";
print Q2 "Position\tQmean\n";
print Q3 "Position\tNumber\tQNmean\n";
#print Q4 "Readnum\tPosition\tbases\tbaselen\tQmean\tQmean1\tQlast\n";
print Q4 "Position baselen Qmean Qmean1 Qlast\n";

#open(OUT,"$output")||die"Cann't write the file $output\n$!";

### System sffino command ###



### cat fna and qual to phrap fastq like  file ###
### Note that the order must be the same ### 

#print("Convert the files to phrap fastq like format\n");
print "******** Start ***********\n";
#my $phrap_file=$output.".phrap_file";
#open(PHRAP,">$phrap_file")||die"Cann't write to the file $phrap_file\n$!";

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

my $Qmean;
my @Qpos;
my %QNpos;
my %QNNpos;

my $GC;
my $aGC;
my $agc;
my $alen;

# split it in Q4 
#my %pos_Q;
#my %num_Q;
#my %last_last;
#my %Mnum;
########################33

## Read a phrap Q file ########################
while((!eof(FASTA)) && (!eof(QUAL))){

	$sequence_name = $sequence_buf;
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
	}
	else {
		$qual_str=~s/^\s*//;
#	print PHRAP "$qual_str\n";
		$qual_ref = $qual_str;
		$qual_str = "";

	}
	%QNNpos=();
	&process_fastq($sequence_ref, $qual_ref);
############    print Q1  ####################
	print Q1 "$read_num\t$Qmean\t$seq_len\t$GC\n";
###########   print Q4 #######################
	my $bases;
	my $pos;
	my $Nlen;
	my $NNmeanlast;
	my $NNmean;
	my $NNmeanforward;
	foreach my $num(sort {$a <=> $b} keys %QNNpos)
	{
		$pos = $num + 1 ;
		my $NNallQ = 0;
		($bases)=keys %{$QNNpos{$num}};
		$Nlen = length($bases);
#		if(exists $Mnum{$Nlen}){
#			$Mnum{$Nlen}++;
#		}else{
#			$Mnum{$Nulen}=1;
#		}
		$NNmeanlast= ${${$QNNpos{$num}}{$bases}}[-1];

		map($NNallQ+=$_, @{$QNNpos{$num}{$bases}});

		$NNmean = sprintf("%.2f", $NNallQ / $Nlen);
		$NNmeanforward = sprintf("%.2f",($NNallQ - $NNmeanlast)/ ($Nlen - 1));
		if($Nlen<$multiply){
			print Q4 "#$pos $Nlen $NNmean $NNmeanforward $NNmeanlast\n";
		}else{
			print Q4 "$pos $Nlen $NNmean $NNmeanforward $NNmeanlast\n";
		}

	}

}

############  print  Q2   ######################################

my @Qpos_mean = map(sprintf("%.2f",$_/$read_num),@Qpos);

print Q2 "$_\t$Qpos_mean[-1+$_]\n" foreach(1..scalar @Qpos_mean) ;

###########  print Q3    ######################################

foreach my $num (sort {$a<=> $b} keys %QNpos)
{	
	my $pos = $num +1;
	my $QNmean=sprintf("%.2f",$QNpos{$num}[1]/$QNpos{$num}[0]);
	print  Q3 "$pos\t$QNpos{$num}[0]\t$QNmean\n";
}


########### function to process the phrap fastq like file #####################3


sub process_fastq {

	my ($sequence, $qual) = @_;
	my @seqlist = split(//,$sequence);
	my @qualist = split(/\s+/, $qual);
	$seq_len = scalar @seqlist;
	my $qua_len = scalar @qualist;
	
	if( $seq_len != $qua_len )
	{
		print STDERR "Number of bases($seq_len) is not equal to the number of quals($qua_len).\n
			Sequence: $sequence\nQual: $qual\n" ;
		exit(1);
	}

	my $gc = $sequence=~s/([gc])/$1/ig;
	$GC=sprintf("%.2f",$gc/$seq_len);

	$agc += $gc;
	$alen += $seq_len;
		
	my $Qallpos = 0;
	my $key;
	my $bac_key;
	my @value;
	my $NNpos;
	my $old_base;
	my $first  = 0;
	my $second =1;

# all pos is start 0;
	for(my $i=0; $i<=$#qualist; $i++){

	$Qallpos += $qualist[$i];
	$Qpos[$i] += $qualist[$i];
	
#	$QNpos[$i] += $qualist[$i];
	if($seqlist[$i]=~m/n/i){
		if(exists $QNpos{$i}[$first]){
			$QNpos{$i}[$first]++;
			$QNpos{$i}[$second] += $qualist[$i];
			
		}else{
			$QNpos{$i}[$first] = 1;
			$QNpos{$i}[$second] = $qualist[$i];
		}

	}

	if((uc($old_base) ne 'N') && (uc($old_base) eq $seqlist[$i]))
	{
		if(!exists $QNNpos{$NNpos})
		{
			$key=$old_base.$seqlist[$i];
			@value=();
			push(@value,($qualist[$i-1],$qualist[$i]));
			@{$QNNpos{$NNpos}{$key}}=@value;
			$bac_key = $key;

		}else{  
			push(@value,$qualist[$i]);
			($key)=keys %{$QNNpos{$NNpos}};
			$key = $bac_key.$seqlist[$i];
			delete $QNNpos{$NNpos}{$bac_key};
			@{$QNNpos{$NNpos}{$key}} = @value;
			$bac_key =$key;
		}

	}else{
		$NNpos = $i;
	}

	$old_base = $seqlist[$i] ;
	}

	$Qmean= sprintf("%.2f",$Qallpos / $qua_len);
}


close FASTA;
close QUAL;
close Q1;
close Q2;
close Q3;
close Q4;

$aGC=sprintf("%.2f",$agc/$alen);




##################### Plot QC with R software ################

my $Rfile = $output.".R";
open(R,">$Rfile")||die "Can't open the R file to plot: $Rfile\n";

##################### Print to R file ########################

print R <<"EOF";
read.table(\"$outfile1\", comment.char=\"#\", header=T,sep=\"\\t\",nrows=5)->tab5rows
sapply(tab5rows,class)->classes
read.table(\"$outfile1\", comment.char=\"#\", header=T,sep=\"\\t\",nrows=($read_num+1),colClasses=classes)->data1
read.table(\"$outfile2\", comment.char=\"#\", header=T,sep=\"\\t\")->data2
par()->old.par
pdf(\"${output}.pdf\")
#jpeg(\"${output}_1.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),cex=0.65,las=1)
density(data1\$Qmean)->d1
density(data1\$Length)->d2
plot(data2\$Qmean,ylim=c(0,40),col=\"blue\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Base quality\",main=\"Base quality of 454 GS FLX sequencing\")
abline(h=20,col=\"grey\",lty=4)
plot(range(d1\$x),range(d1\$y),xlim=c(0,40),type=\"n\",ylab=\"Density\",xlab=\"Reads mean quality\",main=\"Base quality of 454 GS FLX sequencing\")
lines(d1,col=\"blue\",lwd=0.5)
polygon(d1,col=rgb(1,0,0,0.4))
abline(v=d1\$x[which.max(d1\$y)],col=\"green3\",lty=4)
legend(\"topleft\", inset=c(.05,.05),c(\"N= $read_num\",paste(\"peak=\",round(d1\$x[which.max(d1\$y)],2))),col=c(\"black\",\"green3\"),bty=\"n\",lty=c(1,4))
#mtext(\"Base quality of 454 GS FLX sequencing\",outer=T,font=2)
#dev.off()
#jpeg(\"${output}_2.jpeg\",width=9000,height=8000,res=1200)
par(fig=c(0,1,0,0.88),oma=c(0,0,2,0),mar=c(5,4,1,1))
plot(range(d2\$x),range(d2\$y),type=\"n\",ylab=\"Density\",xlab=\"Read length (bp)\")
lines(d2,col=\"blue\",lwd=0.5)
polygon(d2,col=rgb(0,0,1,0.4))
abline(v=mean(data1\$Length),col=\"green3\",lty=4)
legend(\"topright\",inset=c(.05,.05),paste(\"Mean=\",round(mean(data1\$Length),1),\" bp\"),bty=\"n\",lty=4,col=\"green3\")
par(fig=c(0,1,0.7,1), new=TRUE)
boxplot(data1\$Length, horizontal=TRUE, axes=FALSE,col=\"green3\")
mtext(\"Read Length in 454 GS FLX sequencing\",outer=T,font=2)
#dev.off()
par(old.par)

#jpeg(\"${output}_3.jpeg\",width=8000,height=6500,res=1200)
#par(cex=0.8)
h <- hist( data1\$GC, plot = FALSE )
plot( h , border = NA, freq = FALSE, xlab = \"GC-Content\", ylab = \"Probability\",main=\"Distribution of GC-Content\")
usr <- par( \"usr\" )
ncolors <- 100
dy <- ( usr[4] - usr[3] ) / ncolors
colors <- colorRampPalette( c(\"yellow\",\"orange\",\"red\") )(ncolors)
abline( h = axTicks(2) , col = \"gray\", lwd = .5 )
for( i in 1:ncolors){clip( usr[1], usr[2], usr[3] + (i-1) * dy, usr[3] + i*dy )
plot( h, add = TRUE, axes = FALSE, ylab = \"\", xlab = \"\",col = colors[i], border = NA, freq = FALSE)}
do.call( clip, as.list( usr) )
plot( h, add = TRUE, lwd = .5 , freq = FALSE, xlab = \"\", ylab = \"\", axes = FALSE )
#rug( data1\$GC, col = \"#00000088\" )
box()
legend(\"topright\",inset=c(.05,.05),paste(\"GC-Content =\",$aGC),bty=\"n\")
#dev.off()

read.table(\"$outfile3\",comment.char=\"#\",header=T,sep=\"\\t\")->data3
#jpeg(\"${output}_4.jpeg\",width=9000,height=8000,res=1200)
#par(mar=c(5,4,4,4)+0.1)
var=max(data3\$Number)/40
var2=round(max(data3\$Number)/5)
plot(data3\$Position,data3\$Number,type='h',col=\"blue\",yaxt=\"n\",lwd=2,xlab=\"\",ylab=\"\")
lines(data3\$Position,data3\$QNmean*var,type=\"l\",col=\"red\",lty=1,lwd=2)
axis(2,at=seq(0,max(data3\$Number),var2),labels=seq(0,max(data3\$Number),var2),col.axis=\"blue\",las=2)
axis(4,at=var*seq(0,40,5),labels=seq(0,40,5),col.axis=\"red\",las=2,tck=-.01)
mtext(\"Qmean\",4,line=2.2,col=\"red\")
mtext(\"Number\",2,line=2,col=\"blue\")
title(\"Base N in GS FLX sequencing\",xlab=\"Base site\",ylab=\"\")
legend(x=\"topright\",inset=c(.05,.02),c(paste(\"N number =\",sum(data3\$Number)),paste(\" Qmean =\",round(mean(data3\$QNmean),2))),col=\"blue\",bty=\"n\")
#dev.off()
read.table(\"$outfile4\", comment.char=\"#\", header=T,sep=\" \",nrows=5)->tab5rows
sapply(tab5rows,class)->classes
read.table(\"$outfile4\", comment.char=\"#\", header=T,sep=\" \",nrows=10001,colClasses=classes)->data4
attach(data4);
#jpeg(\"${output}_5.jpeg\",width=9000,height=8000,res=1200)
#plot.new()
#par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5,4,2,2),cex=0.7)
x=sort(unique(Position))
y1=y2=y3=1:length(x)
for( i in 1:length(x))
{
	y1[i]=sum(Qmean[Position==x[i]]*baselen[Position==x[i]])/sum(baselen[Position==x[i]])
	y2[i]=sum(Qmean1[Position==x[i]]*(baselen[Position==x[i]]-1))/sum(baselen[Position==x[i]]-1)
	y3[i]=sum(Qlast[Position==x[i]])/length(Position[Position==x[i]])
}
plot(range(x),range(0:40),type=\"n\",xlab=\"First site in poly(N)\",ylab=\"Quality\",main=\"Poly(N) base quality in GS FLX 454 sequencing\")
lines(x,y1,lwd=1,col=2)
lines(x,y2,lwd=1,col=3)
lines(x,y3,lwd=1,col=4)
legend(x=\"bottomleft\",inset=c(.01,.01),c(\"All Qmean\",\"Forward Qmean\",\"Last Qmean\"),col=c(2,3,4),lty=c(1,1,1),bty=\"n\")
x=sort(unique(baselen))
y1=y2=y3=1:length(x)
for(i in 1:length(x))
{
	y1[i]=sum(Qmean[baselen==x[i]]*baselen[baselen==x[i]])/sum(baselen[baselen==x[i]])
	y2[i]=sum(Qmean1[baselen==x[i]]*(baselen[baselen==x[i]]-1))/sum(baselen[baselen==x[i]]-1)
	y3[i]=sum(Qlast[baselen==x[i]])/length(baselen[baselen==x[i]])
}
plot(range(x),range(0:40),type=\"n\",xlab=\"length of poly(N) (bp)\",ylab=\"Quality\",main=\"Poly(N) base quality in GS FLX 454 sequencing\")
lines(x,y1,lwd=1,col=2)
lines(x,y2,lwd=1,col=3)
lines(x,y3,lwd=1,col=4)
legend(x=\"topright\",inset=c(.01,.01),c(\"All Qmean\",\"Forward Qmean\",\"Last Qmean\"),col=c(2,3,4),lty=c(1,1,1),bty=\"n\")
for(i in 1:length(x))
{
	y1[i]=length(baselen[baselen==x[i]])
}
plot(x,y1,type=\"h\",xlab=\"Length of poly(N) (bp)\", ylab=\"Number\", lwd=2, col=\"green3\",main=\"Poly(N) base quality in GS FLX 454 sequencing\")
t=(Position+baselen-1)
x=sort(unique(t))
y1=1:length(x)
for(i in 1:length(x))
{
	y1[i]=sum(Qlast[x==x[i]])/length(t[x==x[i]])
}
plot(x,y1,ylim=c(0,40),type=\"h\",xlab=\"Last base site in poly(N)\",ylab=\"Last Qmean\",lwd=0.8,col=\"blue\",main=\"Poly(N) base quality in GS FLX 454 sequencing\")
#mtext(\"Poly(N) base quality in GS FLX 454 sequencing\",out=T,font=2)
dev.off()

EOF

close R;

if($plot || $delete){
	system ("R --slave --min-vsize=100M --max-vsize=10000M --min-nsize=5M --max-nsize=1000M -f $Rfile");
	system(" rm -f Rplots.pdf")
}
my $compressname=$output.".tar.gz" if($compress);
if($delete){
	system(" rm -f $outfile1 $outfile2 $outfile3 $outfile4 $Rfile");
	system(" tar -zcvf $compressname ${output}.pdf") if($compress);
}else{
	system(" tar -zcvf  $compressname $outfile1 $outfile2 $outfile3 $outfile4 $Rfile ${output}.pdf ") if($compress && $plot);
	system("tar -zcvf  $compressname $outfile1 $outfile2 $outfile3 $outfile4 $Rfile") if($compress && !$plot);
}

if($compress){ system("rm -f $outfile1 $outfile2 $outfile3 $outfile4 $Rfile ${output}.pdf");}

my $time2=new Benchmark;
my $timevar1=timediff($time2,$time1);

print "******** Finish ***********\n";
print timestr($timevar1),"\n";

