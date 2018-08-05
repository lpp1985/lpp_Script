#!/usr/bin/perl -w
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

use strict;
use warnings;
use Getopt::Long;
use Benchmark;
use Term::ANSIColor qw(:constants);

$Term::ANSIColor::AUTORESET = 1;
my $time1=new Benchmark;
my $CURRENT_VERSION="0.5";
my $PROGRAM_NAME=`basename $0`;
chomp $PROGRAM_NAME;

sub usage{
	print CYAN << "EOF";

**************************************************************************************************

 Program: $PROGRAM_NAME (preprocess of solexa fastq data)

 Version: $CURRENT_VERSION

 Contact: Tongwu Zhang <zhangtw\@big.ac.cn>

   Usage: $PROGRAM_NAME [ [-f fileF -r fileR] or [ -i file] ] [-n Number][-s Number][-p][-z][-d][-h]

 Options: 
          -f     solexa pair-end forward fastq file
	  -r     solexa pair-end reverse fastq file
	  -i     solexa fragment fastq file 
	  -o     prefix of the output jpeg file
	  -n     set length of seed for pcr duplication[50]
	  -s	 sample 1/N data to plot, default[1]
	  -q	 solexa fastq format,sanger fastq need set to 33,default [64]
	  -p 	 plot the outcome with R package [F]
          -d     delete the temp files of pic process  [F] 
	  -z 	 compress the QC file and picture [F]
	  -h     help
 Example:

          $PROGRAM_NAME -f ../s_1_1_sequence.txt -r ../s_1_2_sequence.txt -o s_1 -n 60 -p -z -s 1000

   
  Note:
  	this program is designed to analysis the quality of solexa sequncing data,
	and using R to plot. please check the R program is ok!!

**************************************************************************************************

EOF
exit(1);

}


my ($forwardFQ,$reverseFQ,$fragmentFQ,$outputP,$format,$du_length,$help,$skip,$tmpdel,$compress,$plot);
GetOptions('f|forward_fq=s'      =>  \$forwardFQ,
	   'r|reverse_fq=s'      =>  \$reverseFQ,
	   'i|fragment=s'       =>  \$fragmentFQ,
	   'o|output_prefix=s'   =>  \$outputP,
	   'n|len=i'             =>  \$du_length,
	   'q|format=i' 	 =>  \$format,
	   's|skip=i'	         =>  \$skip,
	   'z|compress'		 =>  \$compress,
	   'p|plot'		 =>  \$plot,
	   'd|del_tmp'         =>  \$tmpdel,
	   'h|help'	         =>  \$help
	   );

&usage if(!(defined $forwardFQ && defined $reverseFQ) && ! defined $fragmentFQ);
&usage if($help);
&usage if(defined $fragmentFQ && (defined $forwardFQ || defined $reverseFQ));
if(!defined $outputP){
	print "\nError, please set the option with -o\n";
	&usage;
}



my $duplen=$du_length||50;
$skip = $skip || 1;
$plot = $plot || 0;
$compress = $compress || 0;
$tmpdel = defined $tmpdel?1:0;
$format= $format || 64;

if($format != 33 &&  $format !=64){
	print "\nError, please set correct fastq format, solexa or illumina is 64, sanger is 33\n";
}


my $file1=$outputP.".qa";
my $file2=$outputP.".qp";
my $filen=$outputP.".n";
my $filec=$outputP.".fc";
my $filed=$outputP.".dg";
my $fileq=$outputP.".q20";

my $select;
my $pikn=0;
my $all_line=0;
my $pstart=0;


############################## Function to Pair-end #############################3


sub process_pair_end{

open F1 , "$forwardFQ"||die"$! cann't open solexa fowrard file:$forwardFQ\n";

open F2 , "$reverseFQ"||die"$! cann't open solexa reverse file:$reverseFQ\n";


open O1, ">$file1"||die "$!";
open O2, ">$file2"||die "$!";
open N,  ">$filen"||die "$!";
open C,  ">$filec"||die "$!";
open D,  ">$filed"||die "$!";
open QQ, ">$fileq"||die "$!";

my $line=0;
my (@posbaseA1,@posbaseC1,@posbaseG1,@posbaseT1,@posbaseN1);
my (@posbaseA2,@posbaseC2,@posbaseG2,@posbaseT2,@posbaseN2);

@posbaseA1=@posbaseC1=@posbaseG1=@posbaseT1=@posbaseN1=();
@posbaseA2=@posbaseC2=@posbaseG2=@posbaseT2=@posbaseN2=();

my (@plista, @plistb);
@plista=@plistb=();

my ($A1,$T1,$C1,$G1,$N1,$A2,$T2,$C2,$G2,$N2);

$A1=$T1=$C1=$G1=$N1=0;
$A2=$T2=$C2=$G2=$N2=0;

my ($warning,$error,$i);
$warning=$error=0;

my (%nq1,%nq2);

for( $i=0;$i<=40;$i++){

	$nq1{$i}=0;

	$nq2{$i}=0;

}

print O1 "lane flowcell x_cor y_cor qmean1 qmean2 gccont1 gccont2\n";

print D "duptime subgc1 subgc2\n";

my ($name1,$name2,$lane,$flowcell,$x_cor,$y_cor,$s1,$s2,$q1,$q2,$subs1,$subs2,$subs);
my ($len1,$len2,$p1,$p2,$qmean1,$qmean2,$gc1,$gc2,$gccont1,$gccont2);
my (@Qua1,@Qua2);
my (@seq1,@seq2);
my (%flc,%flh,%flhh,%subhash,%gchash1,%gchash2);
my ($eman,$data1,$data2);
my (@QF20,@CQF20,@QR20,@CQR20);

$/="@";
<F1>;<F2>;
$select=$skip*$pikn + int(rand($skip))+1;
while (!eof(F1) && !eof(F2)){

	$data1=<F1>;
	$data2=<F2>;
	$all_line++;
	if($all_line > $skip*($pikn+1)) {$pikn++;$pstart=0;$select=$skip*$pikn + int(rand($skip))+1;}
#$select=$skip*$pikn + int(rand($skip))+1;
	next if($select != $all_line);
	$pstart++;
	next if($pstart != 1);
#	print "$all_line\n";
	chomp $data1;
	chomp $data2;

	($name1,$s1,$eman,$q1)=split("\n",$data1);
        ($name2,$s2,$eman,$q2)=split("\n",$data2);

	$name1=~s/#.*$//;
	$name2=~s/#.*$//;

	($eman,$lane,$flowcell,$x_cor=$3,$y_cor)=split(/:/,$name1);

	if($name1 ne $name2){
		
		print "Warning: Name of FR is different. F is $name1 while R is $name2\n";

		$warning++;
	}

	$len1=length($q1);

	$len2=length($q2);

	if($len1 != $len2){

		print "Error: Length of FR is different. F is $len1 while R is $len2\n";

		$error++;

		next;

	}

	$line++;

	$p1=$p2=0;

	@Qua1=map(ord($_)-$format,split(//,$q1));

	map($p1+=$_,@Qua1);

	$qmean1=sprintf "%.2f", $p1/$len1;

	@Qua2=map(ord($_)-$format,split(//,$q2));
	
	map($p2+=$_,@Qua2);

	$qmean2=sprintf "%.2f", $p2/$len2;

	my $sig1=1;
	my $sig2=1;

	for(my $i=0; $i<$len1; $i++){
		if($Qua1[$i]>=20){ $QF20[$i]++;}else{ $sig1=0;}
		if($Qua2[$i]>=20){ $QR20[$i]++;}else{ $sig2=0;}
		$CQF20[$i]++ if($sig1);
		$CQR20[$i]++ if($sig1);
	}

	@seq1=split(//,$s1);

	@seq2=split(//,$s2);

	$gc1=$gc2=0;

	for( $i=0;$i<$len1;$i++){

		$plista[$i]+=$Qua1[$i];

		$plistb[$i]+=$Qua2[$i];

		if($seq1[$i]=~m/a/i) {
			$posbaseA1[$i]++;$A1++;
		}elsif($seq1[$i]=~m/c/i) {
			$posbaseC1[$i]++;$C1++;$gc1++;
		}elsif($seq1[$i]=~m/t/i) {
			$posbaseT1[$i]++;$T1++;
		}elsif($seq1[$i]=~m/g/i) {
			$posbaseG1[$i]++;$G1++;$gc1++;
		}elsif($seq1[$i]=~m/n/i) {$posbaseN1[$i]++;$N1++;$nq1{$Qua1[$i]}++;
		}

		if($seq2[$i]=~m/a/i) {
			$posbaseA2[$i]++;$A2++;
		}elsif($seq2[$i]=~m/c/i) {
			$posbaseC2[$i]++;$C2++; $gc2++;
		}elsif($seq2[$i]=~m/t/i) {
			$posbaseT2[$i]++;$T2++;
		}elsif($seq2[$i]=~m/g/i) {
			$posbaseG2[$i]++;$G2++; $gc2++;
		}elsif($seq2[$i]=~m/n/i) {
			$posbaseN2[$i]++;$N2++;$nq2{$Qua2[$i]}++;
		}

	}

	$gccont1=sprintf "%.2f",$gc1/$len1;
	$gccont2=sprintf "%.2f",$gc2/$len2;


	print O1 "$lane $flowcell $x_cor $y_cor $qmean1 $qmean2 $gccont1 $gccont2\n";


	if(exists $flc{$flowcell}){

		$flc{$flowcell}++;

	}else{

		$flc{$flowcell}=1;

	}

	if($qmean1>20 && $qmean2>20){

		if(exists $flh{$flowcell}){

			$flh{$flowcell}++;

			}else{
				$flh{$flowcell}=1;
			}
	}

	
	if($qmean1>30 && $qmean2>30){

		if(exists $flhh{$flowcell}){

			$flhh{$flowcell}++;

		}else{

			$flhh{$flowcell}=1;

		}
	}


	$subs1=substr($s1,0,$duplen);

	$subs2=substr($s2,0,$duplen);

	$subs=$subs1.$subs2;

	if(! exists $subhash{$subs}){

		$subhash{$subs}=1;

		$gchash1{$subs}=$gccont1;

		$gchash2{$subs}=$gccont2;

	}else{

		$subhash{$subs}++;

		$gchash1{$subs}+=$gccont1;

		$gchash2{$subs}+=$gccont2;

	}



}

$/="\n";

my ($key,$max,$min);

foreach $key (sort {$a<=>$b} keys %flc){

	if(! defined $max or !defined $min){

		$max=$min=$key;

	}

	$max=$flc{$max}>$flc{$key}?$max:$key;

	$min=$flc{$min}<$flc{$key}?$min:$key;

	$flh{$key}=0 if(! exists $flh{$key});

	$flhh{$key}=0 if(! exists $flhh{$key});

	print C "$key\t$flc{$key}\t$flh{$key}\t$flhh{$key}\n";

}


my $GC1=sprintf "%.2f",($G1+$C1)/($A1+$T1+$G1+$C1+$N1);

my $base1=$A1+$T1+$G1+$C1;

my $GC2=sprintf "%.2f",($G2+$C2)/($A2+$T2+$G2+$C2+$N2);

my $base2=$A2+$T2+$G2+$C2;

print O2 "#Forward: base=$base1; A=$A1; T=$T1; C=$C1; G=$G1; N=$N1; GC=$GC1\n";

print O2 "#Reverse: base=$base2; A=$A2; T=$T2; C=$C2; G=$G2; N=$N2; GC=$GC2\n";


print O2  "Position\tQual1\tratioA1\tratioT1\tratioC1\tratioG1\tratioN1\tQual2\tratioA2\tratioT2\tratioC2\tratioG2\tratioN2\n";

my ($ratioA1,$ratioC1,$ratioT1,$ratioG1,$ratioN1);
my ($ratioA2,$ratioC2,$ratioT2,$ratioG2,$ratioN2);
my $pos;
my (@qpa,@qpb);

for($i=0;$i<=$#Qua1;$i++){

	$qpa[$i]=sprintf "%.2f", $plista[$i]/$line;

	$qpb[$i]=sprintf "%.2f", $plistb[$i]/$line;

	$posbaseN1[$i]=0 if(!defined $posbaseN1[$i]);
	$posbaseA1[$i]=0 if(!defined $posbaseA1[$i]);
	$posbaseT1[$i]=0 if(!defined $posbaseT1[$i]);
	$posbaseC1[$i]=0 if(!defined $posbaseC1[$i]);
	$posbaseG1[$i]=0 if(!defined $posbaseG1[$i]);

	$ratioA1=sprintf "%.2f", $posbaseA1[$i]/$line;
	$ratioC1=sprintf "%.2f", $posbaseC1[$i]/$line;
	$ratioT1=sprintf "%.2f", $posbaseT1[$i]/$line;
	$ratioG1=sprintf "%.2f", $posbaseG1[$i]/$line;
	$ratioN1=sprintf "%.2f", $posbaseN1[$i]/$line;

	$posbaseN2[$i]=0 if(!defined $posbaseN2[$i]);
	$posbaseA2[$i]=0 if(!defined $posbaseA2[$i]);
	$posbaseT2[$i]=0 if(!defined $posbaseT2[$i]);
	$posbaseC2[$i]=0 if(!defined $posbaseC2[$i]);
	$posbaseG2[$i]=0 if(!defined $posbaseG2[$i]);


	$ratioA2=sprintf "%.2f", $posbaseA2[$i]/$line;
	$ratioC2=sprintf "%.2f", $posbaseC2[$i]/$line;
	$ratioT2=sprintf "%.2f", $posbaseT2[$i]/$line;
	$ratioG2=sprintf "%.2f", $posbaseG2[$i]/$line;
	$ratioN2=sprintf "%.2f", $posbaseN2[$i]/$line;

	$pos=$i+1;

	print O2 "$pos\t$qpa[$i]\t$ratioA1\t$ratioT1\t$ratioC1\t$ratioG1\t$ratioN1\t$qpb[$i]\t$ratioA2\t$ratioT2\t$ratioC2\t$ratioG2\t$ratioN2\n";

}

print N "#N1=$N1\tN2=$N2\n";

print N "Q\tN1\tN2\n";

for($i=0;$i<=40;$i++){

	print N "$i\t$nq1{$i}\t$nq2{$i}\n";
}


my ($subgc1,$subgc2);

foreach $key (keys %subhash){

	$subgc1=sprintf "%.2f", $gchash1{$key}/$subhash{$key};

	$subgc2=sprintf "%.2f", $gchash2{$key}/$subhash{$key};

	print D "$subhash{$key} $subgc1 $subgc2\n";

}

print QQ "Site\tQF20\tQR20\tCQF20\tCQR20\n";

for(my $i=0; $i<$len1; $i++)
{
	my $pos;
	$pos=$i+1;
	printf QQ  "$pos\t%.2f\t%.2f\t%.2f\t%.2f\n",$QF20[$i]/$line,$QR20[$i]/$line,$CQF20[$i]/$line,$CQR20[$i]/$line;
}	


close F1;
close F2;
close O1;
close O2;
close N;
close C;
close D;
close QQ;

print "Note: warning $warning and error $error\n";

undef %subhash;
undef %gchash1;
undef %gchash2;


################## Generate .R file to plot ###################

################## plot1: Base quality distribution of Solexa Sequencing library ##################

`awk '\$2==${max}{print \$3"\\t"\$4"\\t"\$5"\\t"\$6}' $file1 >${outputP}.max`;
`awk '\$2==${min}{print \$3"\\t"\$4"\\t"\$5"\\t"\$6}' $file1  >${outputP}.min`;
open(R,">${outputP}.R")||die"$!";
print R <<EOF;
read.table(\"$file1\",comment.char=\"#\",header=T,sep=\" \",nrows=5)->tab5rows
sapply(tab5rows,class)->classes
read.table(\"$file1\",comment.char=\"#\",header=T,sep=\" \",nrows=$line,colClasses=classes)->data1
read.table(\"$file2\",comment.char=\"#\",header=T,sep=\"\\t\")->data2
read.table(\"$filen\",comment.char=\"#\",header=T,sep=\"\\t\")->datan
read.table(\"$filec\",comment.char=\"#\",header=F,sep=\"\\t\")->slane
read.table(\"$fileq\",comment.char=\"#\",header=T,sep=\"\\t\")->Q20
attach(data1)
attach(data2)
attach(datan)
attach(Q20)

par()->old.par
pdf(\"${outputP}.pdf\")
#jpeg(\"${outputP}_1.jpeg\",width=9000,height=8000,res=1200)
density(qmean1)->d1
density(qmean2)->d2
#plot.new()
par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5,4,2,2),cex=0.7)
plot(Qual1,ylim=c(0,40),col=\"red\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Base quality\")
abline(h=20,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->F\",col=\"red\",,bty=\"n\",lty=1)
plot(range(d1\$x),range(d1\$y),xlim=c(0,40),type=\"n\",ylab=\"Density\",xlab=\"Reads mean quality\")
lines(d1,col=\"red\")
legend(x=\"topleft\",inset=c(.05,.02),\"N=$line\",col=\"red\",,bty=\"n\",lty=1)
plot(Qual2,ylim=c(0,40),col=\"blue\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Base quality\")
abline(h=20,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->R\",col=\"blue\",,bty=\"n\",lty=1)
plot(range(d2\$x),range(d2\$y),xlim=c(0,40),type=\"n\",ylab=\"Density\",xlab=\"Reads mean quality\")
lines(d2,col=\"blue\")
legend(x=\"topleft\",inset=c(.05,.02),\"N=$line\",col=\"blue\",,bty=\"n\",lty=1)
mtext(\"Base quality distribution of Solexa sequencing library\",outer=T,font=2)
#dev.off()
par(old.par)

#jpeg(\"${outputP}_2.jpeg\",width=9000,height=8000,res=1200)
#plot.new()
par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(5,4,2,2),cex=0.7)
plot(QF20,ylim=c(0,1),col=\"red\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Single Q20 Ratio\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->F\",col=\"red\",,bty=\"n\",lty=1)
plot(CQF20,ylim=c(0,1),col=\"red\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Cumulative Q20 Ratio\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->F\",col=\"red\",,bty=\"n\",lty=1)
plot(QR20,ylim=c(0,1),col=\"blue\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Single Q20 Ratio\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->R\",col=\"red\",,bty=\"n\",lty=1)
plot(CQR20,ylim=c(0,1),col=\"blue\",type=\"h\",lwd=1,xlab=\"Base site\",ylab=\"Cumulative Q20 Ratio\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->R\",col=\"blue\",,bty=\"n\",lty=1)
mtext(\"Q20 read ratio in Solexa sequencing library\",outer=T,font=2)
#dev.off()
par(old.par)


#jpeg(\"${outputP}_4.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(5,4,2,1),cex=0.6,las=1)
h <- hist( gccont1, plot = FALSE )
plot( h , border = NA, freq = FALSE, xlab = \"F GC-Content\", ylab =\"Probability\",main=\"Distribution of F GC-Content\")
usr <- par( \"usr\" )
ncolors <- 100
dy <- ( usr[4] - usr[3] ) / ncolors
colors <- colorRampPalette( c(\"yellow\",\"orange\",\"red\") )(ncolors)
abline( h = axTicks(2) , col = \"gray\", lwd = .5 )
for( i in 1:ncolors){clip( usr[1], usr[2], usr[3] + (i-1) * dy, usr[3] + i*dy )
plot( h, add = TRUE, axes = FALSE, ylab = \"\", xlab = \"\",col = colors[i], border = NA, freq = FALSE)}
do.call( clip, as.list( usr) )
plot( h, add = TRUE, lwd = .5 , freq = FALSE, xlab = \"\", ylab = \"\", axes = FALSE )
#rug( gccont1, col = \"#00000088\" )
box()
h <- hist( gccont2, plot = FALSE )
plot( h , border = NA, freq = FALSE, xlab = \"R GC-Content\", ylab = \"Probability\",main=\"Distribution of R GC-Content\")
usr <- par( \"usr\" )
ncolors <- 100
dy <- ( usr[4] - usr[3] ) / ncolors
colors <- colorRampPalette( c(\"yellow\",\"orange\",\"red\") )(ncolors)
abline( h = axTicks(2) , col = \"gray\", lwd = .5 )
for( i in 1:ncolors){clip( usr[1], usr[2], usr[3] + (i-1) * dy, usr[3] + i*dy )
plot( h, add = TRUE, axes = FALSE, ylab = \"\", xlab = \"\",col = colors[i], border = NA, freq = FALSE)}
do.call( clip, as.list( usr) )
plot( h, add = TRUE, lwd = .5 , freq = FALSE, xlab = \"\", ylab = \"\", axes = FALSE )
#rug( gccont2, col = \"#00000088\" )
box()
#dev.off()
par(old.par)

#jpeg(\"${outputP}_5.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),las=1,cex=0.48)
d1<-density(qmean1[gccont1<=0.2])
d2<-density(qmean1[gccont1>0.2 & gccont1<=0.3])
d3<-density(qmean1[gccont1>0.3 & gccont1<=0.4])
d4<-density(qmean1[gccont1>0.4 & gccont1<=0.5])
d5<-density(qmean1[gccont1>0.5 & gccont1<=0.6])
d6<-density(qmean1[gccont1>0.6])
plot(range(c(d1\$x,d2\$x,d3\$x,d4\$x,d5\$x,d6\$x)),range(c(d1\$y,d2\$y,d3\$y,d4\$y,d5\$y,d6\$y)),type=\"n\",xlab=\"F-Qmean\",ylab=\"Density\",main=\"Distribution of Qmean with different GC-Content\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
lines(d4,lwd=1,col=5)
lines(d5,lwd=1,col=6)
lines(d6,lwd=1,col=7)
legend(x=\"topleft\",inset=c(.01,.01),c(\"   0<GC<=0.2\",\"0.2<GC<=0.3\",\"0.3<GC<=0.4\",\"0.4<GC<=0.5\",\"0.5<GC<=0.6\",\"0.6<GC<1\"),col=c(2,3,4,5,6,7),lty=c(1,1,1,1,1,1),bty=\"n\")
d1<-density(qmean2[gccont2<=0.2])
d2<-density(qmean2[gccont2>0.2 & gccont1<=0.3])
d3<-density(qmean2[gccont2>0.3 & gccont1<=0.4])
d4<-density(qmean2[gccont2>0.4 & gccont1<=0.5])
d5<-density(qmean2[gccont2>0.5 & gccont1<=0.6])
d6<-density(qmean2[gccont2>0.6])
plot(range(c(d1\$x,d2\$x,d3\$x,d4\$x,d5\$x,d6\$x)),range(c(d1\$y,d2\$y,d3\$y,d4\$y,d5\$y,d6\$y)),type=\"n\",xlab=\"R-Qmean\",ylab=\"Density\",main=\"Distribution of Qmean with different GC-Content\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
lines(d4,lwd=1,col=5)
lines(d5,lwd=1,col=6)
lines(d6,lwd=1,col=7)
legend(x=\"topleft\",inset=c(.01,.01),c(\"   0<GC<=0.2\",\"0.2<GC<=0.3\",\"0.3<GC<=0.4\",\"0.4<GC<=0.5\",\"0.5<GC<=0.6\",\"0.6<GC<1\"),col=c(2,3,4,5,6,7),lty=c(1,1,1,1,1,1),bty=\"n\")
#mtext(\"Distribution of Qmean with different GC-Content\",outer=T,font=2,cex=0.6)
#dev.off()
par(old.par)

#detach (data1)
#rm(data1)
#jpeg(\"${outputP}_6.jpeg\",width=10000,height=5000,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(5,4,2,1),cex=0.6)
plot(ratioA1~Position,cex=1.2,ylim=c(0,0.5),xlab=\"Base site\",ylab=\"Ratio\",main=\"F:GC-Content=$GC1\",type=\"n\")
abline(h=0.1, col = \"gray\", lty=3)
abline(h=0.2, col = \"gray\", lty=3)
abline(h=0.3, col = \"gray\", lty=3)
abline(h=0.4, col = \"gray\", lty=3)
abline(h=0.5, col = \"gray\", lty=3)
lines(x=Position,y=ratioA1,col=\"red\")
lines(x=Position,y=ratioC1,col=\"green\")
lines(x=Position,y=ratioG1,col=\"black\")
lines(x=Position,y=ratioT1,col=\"blue\")
lines(x=Position,y=ratioN1,col=\"purple\")
legend(x=\"top\",inset=c(.05,.02),c(\"A\",\"T\",\"C\",\"G\",\"N\"),col=c(\"red\",\"blue\",\"green\",\"black\",\"purple\"),lty=c(1,1,1,1,1),horiz=T,bty=\"n\")
plot(ratioA2~Position,cex=1.2,ylim=c(0,0.5),xlab=\"Base site\",ylab=\"Ratio\",main=\"R:GC-Content=$GC2\",type=\"n\")
abline(h=0.1, col = \"gray\", lty=3)
abline(h=0.2, col = \"gray\", lty=3)
abline(h=0.3, col = \"gray\", lty=3)
abline(h=0.4, col = \"gray\", lty=3)
abline(h=0.5, col = \"gray\", lty=3)
lines(x=Position,y=ratioA2,col=\"red\")
lines(x=Position,y=ratioC2,col=\"green\")
lines(x=Position,y=ratioG2,col=\"black\")
lines(x=Position,y=ratioT2,col=\"blue\")
lines(x=Position,y=ratioN2,col=\"purple\")
legend(x=\"top\",inset=c(.05,.02),c(\"A\",\"T\",\"C\",\"G\",\"N\"),col=c(\"red\",\"blue\",\"green\",\"black\",\"purple\"),lty=c(1,1,1,1,1),horiz=T,bty=\"n\")
#dev.off()
par(old.par)

#jpeg(\"${outputP}_7.jpeg\",width=9000,height=6000,res=1200)
par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(5,4,1,2),cex=0.6)
barplot(N1,names=Q,col=\"red\",axis.lty=1,xlab=\"Q-F\",main=\"Quality of Base N in Foward Reads\")
legend(x=\"topright\",inset=c(.05,.02),\"N=$N1\",col=\"red\",pch=15,bty=\"n\")
barplot(N2,names=Q,col=\"blue\",axis.lty=1,xlab=\"Q-R\",main=\"Quality of Base N in Reversed Reads\")
legend(x=\"topright\",inset=c(.05,.02),\"N=$N2\",col=\"blue\",pch=15,bty=\"n\")
#dev.off()
par(old.par)




read.table(\"$filed\",header=T,sep=\" \")->data4
attach(data4)
#jpeg(\"${outputP}_11.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),las=1,cex=0.48)
quantile(duptime[duptime>1],probs=c(0,0.5,1))->qu
qu1<-as.numeric(qu[2])
qu2<-as.numeric(qu[3])
d1<-density(subgc1[duptime==1])
d2<-density(subgc1[duptime>1 & duptime<=qu1])
d3<-density(subgc1[duptime>qu1 & duptime<=qu2])
plot(range(c(d1\$x,d2\$x,d3\$x)),range(c(d1\$y,d2\$y,d3\$y)),type=\"n\",xlab=\"F GC-Content\",ylab=\"Density\",main=\"Distribution of GC-Content with different PCR duplication time\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
legend(x=\"topleft\",inset=c(.01,.01),c(\"Time=1\",paste(\"1<Time<=\",qu1),paste(qu1,\"<Time<=\",qu2)),col=c(2,3,4),lty=c(1,1,1),bty=\"n\")
d1<-density(subgc2[duptime==1])
d2<-density(subgc2[duptime>1 & duptime<=qu1])
d3<-density(subgc2[duptime>qu1 & duptime<=qu2])
plot(range(c(d1\$x,d2\$x,d3\$x)),range(c(d1\$y,d2\$y,d3\$y)),type=\"n\",xlab=\"R GC-Content\",ylab=\"Density\",main=\"Distribution of GC-Content with different PCR duplication time\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
legend(x=\"topleft\",inset=c(.01,.01),c(\"Time=1\",paste(\"1<Time<=\",qu1),paste(qu1,\"<Time<=\",qu2)),col=c(2,3,4),lty=c(1,1,1),bty=\"n\")
#mtext(\"Distribution of GC-Content with different PCR duplication time\",outer=T,font=2,cex=0.6)
dev.off()


jpeg(\"${outputP}_1.jpeg\",width=9000,height=8000,res=1200)
cor.test(qmean1,qmean2)->qcor
if(qcor\$p.value<0.005)  pvalue=\"p.value<0.005\"
if(qcor\$p.value>0.005)  pvalue=\"p.value>0.005\"
par(fig=c(0,0.8,0,0.8))
qqplot(qmean1,qmean2,col=\"red\",cex=0.2,xlab=\"Q-mean F\",ylab=\"Q-mean R\")
curve(x^1,col=\"blue\",add=T)
legend(\"topleft\",c(pvalue,paste(\"Cor=\",round(qcor\$estimate,digits=3))),inset=c(.05,.05))
par(fig=c(0,0.8,0.55,1), new=TRUE)
boxplot(qmean1, horizontal=TRUE, axes=FALSE,col=\"green3\",cex=0.5,pch=20)
par(fig=c(0.65,1,0,0.8),new=TRUE)
boxplot(qmean2, axes=FALSE,col=\"green3\",cex=0.5,pch=20)
mtext(\"qq-plot between Q-mean F and Q-mean R\", side=3, outer=TRUE, line=-3,cex=1.3,font=2)
dev.off()

jpeg(\"${outputP}_2.jpeg\",width=7000,height=3000,res=1200)
par(xpd=NA,cex=0.4)
le1=1.1*max(slane\$V2)
barplot(slane\$V2,names=slane\$V1,col=\"red\",axis.lty=1,xlab=\"Flowcell-tile\",main=\"Reads Number in different tiles in lane $lane\",space=0,lwd=0.8)
barplot(slane\$V3,names=slane\$V1,col=\"blue\",space=0,lwd=0.8,add=T)
barplot(slane\$V4,names=slane\$V1,col=\"green3\",space=0,lwd=0.8,add=T)
legend(1,le1,c(\"Q<=20\",\"20<Q<=30\",\"Q>30\"),fill=c(\"red\",\"blue\",\"green3\"),horiz=T,bty=\"n\")
dev.off()



read.table(\"${outputP}.max\",comment.char=\"#\",header=F,sep=\"\\t\")->maxcell
read.table(\"${outputP}.min\",comment.char=\"#\",header=F,sep=\"\\t\")->mincell
co1<-maxcell\$V1
co2<-maxcell\$V1
co3<-mincell\$V1
co4<-mincell\$V1
for(i in 1:length(maxcell\$V3)) if (maxcell\$V3[i]>20) co1[i]=\"green\" else co1[i]=\"red\"
for(i in 1:length(maxcell\$V4)) if (maxcell\$V4[i]>20) co2[i]=\"green\" else co2[i]=\"red\"
for(i in 1:length(mincell\$V3)) if (mincell\$V3[i]>20) co3[i]=\"green\" else co3[i]=\"red\"
for(i in 1:length(mincell\$V4)) if (mincell\$V4[i]>20) co4[i]=\"green\" else co4[i]=\"red\"
jpeg(\"${outputP}_3.jpeg\",width=8000,height=5000,res=1200)
par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(5,4,2,1),cex=0.6)
plot(maxcell\$V1,maxcell\$V2,cex=0.1,col=co1,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $max (MAX-F)\")
plot(maxcell\$V1,maxcell\$V2,cex=0.1,col=co2,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $max (MAX-R)\")
dev.off()
#par(old.par)

jpeg(\"${outputP}_4.jpeg\",width=8000,height=5000,res=1200)
par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(5,4,2,1),cex=0.6)
plot(mincell\$V1,mincell\$V2,cex=0.1,col=co3,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $min (MIN-F)\")
plot(mincell\$V1,mincell\$V2,cex=0.1,col=co4,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $min (MIN-R)\")
dev.off()
#par(old.par)








EOF

close R;

}
######################## process pair-end over ############################

######################## process fragment ############################

sub process_fragment{

open F1, "$fragmentFQ"||die"Cann't open solexa forward file:$forwardFQ\n$!";

open O1, ">$file1"||die "$!";
open O2, ">$file2"||die "$!";
open N,  ">$filen"||die "$!";
open C,  ">$filec"||die "$!";
open D,  ">$filed"||die "$!";
open QQ, ">$fileq"||die "$!";

my $line=0;
my (@posbaseA1,@posbaseC1,@posbaseG1,@posbaseT1,@posbaseN1);
@posbaseA1=@posbaseC1=@posbaseG1=@posbaseT1=@posbaseN1=();

my @plista;
@plista=();
my ($A1,$T1,$C1,$G1,$N1);
$A1=$T1=$C1=$G1=$N1=0;
my ($warning,$error,$i);
$warning=$error=0;
my %nq1;

for( $i=0;$i<=40;$i++){$nq1{$i}=0;}

print O1 "lane flowcell x_cor y_cor qmean1 gccont1\n";
print D "duptime subgc1\n";

my ($name1,$lane,$flowcell,$x_cor,$y_cor,$s1,$q1,$subs1,$subs);
my ($len1,$p1,$qmean1,$gc1,$gccont1);
my @Qua1;
my @seq1;
my (%flc,%flh,%flhh,%subhash,%gchash1);
my ($eman,$data1);
my (@QF20,@CQF20);

$/="@";
<F1>;

$select=$skip*$pikn + int(rand($skip))+1;
while (!eof(F1)){
	$data1=<F1>;
	$all_line++;
	if($all_line > $skip*($pikn+1)){$pikn++;$pstart=0; $select=$skip*$pikn + int(rand($skip))+1;}
#        if(!$pstart) {$select=$skip*$pikn + int(rand($skip))+1;}
	next if($select != $all_line);
	$pstart++;
	next if($pstart != 1);
#	print "$all_line\n";
	chomp $data1;

	($name1,$s1,$eman,$q1)=split("\n",$data1);
	$name1=~s/#.*$//;
	($eman,$lane,$flowcell,$x_cor,$y_cor)=split(/:/,$name1);
	
	$len1=length($q1);
	$line++;
	$p1=0;
	@Qua1=map(ord($_)-$format,split(//,$q1));
	map($p1+=$_,@Qua1);
	$qmean1=sprintf "%.2f", $p1/$len1;

	my $sig1=1;
	for(my $i=0; $i<$len1; $i++){
		if($Qua1[$i]>=20){ $QF20[$i]++;}else{ $sig1=0;}
		$CQF20[$i]++ if($sig1);
	}

	@seq1=split(//,$s1);
	$gc1=0;
	for( $i=0;$i<$len1;$i++){
		$plista[$i]+=$Qua1[$i];
		if($seq1[$i]=~m/a/i)
		{
			$posbaseA1[$i]++;$A1++;
		}elsif($seq1[$i]=~m/c/i)
		{
			$posbaseC1[$i]++;$C1++;$gc1++;
		}elsif($seq1[$i]=~m/t/i) 
		{
			$posbaseT1[$i]++;$T1++;
		}elsif($seq1[$i]=~m/g/i) 
		{
			$posbaseG1[$i]++;$G1++;$gc1++;
		}elsif($seq1[$i]=~m/n/i) 
		{
			$posbaseN1[$i]++;$N1++;$nq1{$Qua1[$i]}++;
		}
	}

	$gccont1=sprintf "%.2f",$gc1/$len1;
	print O1 "$lane $flowcell $x_cor $y_cor $qmean1 $gccont1\n";

	if(exists $flc{$flowcell}){
		$flc{$flowcell}++;
	}else{
		$flc{$flowcell}=1;
	}
	
	if($qmean1>20){

		if(exists $flh{$flowcell}){
			$flh{$flowcell}++;
		}else{
			$flh{$flowcell}=1;
		}

	}

	if($qmean1>30){
		if(exists $flhh{$flowcell}){
			$flhh{$flowcell}++;
		}else{
			$flhh{$flowcell}=1;
		}
	}

	$subs1=substr($s1,0,$duplen);
	$subs=$subs1;
	
	if(! exists $subhash{$subs}){
		$subhash{$subs}=1;
		$gchash1{$subs}=$gccont1;
	}else{
		$subhash{$subs}++;
		$gchash1{$subs}+=$gccont1;
	}

}

##################################

$/="\n";
my ($key,$max,$min);

foreach $key (sort {$a<=>$b} keys %flc){
	if(! defined $max or !defined $min){
		$max=$min=$key;
	}
	
	$max=$flc{$max}>$flc{$key}?$max:$key;
	$min=$flc{$min}<$flc{$key}?$min:$key;
	$flh{$key}=0 if(! exists $flh{$key});
	$flhh{$key}=0 if(! exists $flhh{$key});
	print C "$key\t$flc{$key}\t$flh{$key}\t$flhh{$key}\n";
}

my $GC1=sprintf "%.2f",($G1+$C1)/($A1+$T1+$G1+$C1+$N1);
my $base1=$A1+$T1+$G1+$C1;
print O2 "#fragment: base=$base1; A=$A1; T=$T1; C=$C1; G=$G1; N=$N1; GC=$GC1\n";
print O2  "Position\tQual1\tratioA1\tratioT1\tratioC1\tratioG1\tratioN1\n";

my ($ratioA1,$ratioC1,$ratioT1,$ratioG1,$ratioN1);
my $pos;
my @qpa;

for($i=0;$i<=$#Qua1;$i++){
	
	$qpa[$i]=sprintf "%.2f", $plista[$i]/$line;
	$posbaseN1[$i]=0 if(!defined $posbaseN1[$i]);
	$posbaseA1[$i]=0 if(!defined $posbaseA1[$i]);
	$posbaseT1[$i]=0 if(!defined $posbaseT1[$i]);
	$posbaseC1[$i]=0 if(!defined $posbaseC1[$i]);
	$posbaseG1[$i]=0 if(!defined $posbaseG1[$i]);
	
	$ratioA1=sprintf "%.2f", $posbaseA1[$i]/$line;
	$ratioC1=sprintf "%.2f", $posbaseC1[$i]/$line;
	$ratioT1=sprintf "%.2f", $posbaseT1[$i]/$line;
	$ratioG1=sprintf "%.2f", $posbaseG1[$i]/$line;
	$ratioN1=sprintf "%.2f", $posbaseN1[$i]/$line;
	$pos=$i+1;
	print O2 "$pos\t$qpa[$i]\t$ratioA1\t$ratioT1\t$ratioC1\t$ratioG1\t$ratioN1\n";
}

print N "#N1=$N1\n";
print N "Q\tN1\n";

for($i=0;$i<=40;$i++){
	print N "$i\t$nq1{$i}\n";
}

my $subgc1;

foreach $key (keys %subhash){
	$subgc1=sprintf "%.2f", $gchash1{$key}/$subhash{$key};
	print D "$subhash{$key} $subgc1\n";
}


print QQ "Site\tQF20\tCQF20\n";

for(my $i=0; $i<$len1; $i++)
{
	my $pos;
	$pos=$i+1;
	printf QQ  "$pos\t%.2f\t%.2f\n",$QF20[$i]/$line,$CQF20[$i]/$line;
}	





close F1;
close O1;
close O2;
close N;
close C;
close D;
close QQ;

undef %subhash;
undef %gchash1;

#################################### fragment ####################################################
################## plot1: Base quality distribution of Solexa Sequencing library ##################

`awk '\$2==${max}{print \$3"\\t"\$4"\\t"\$5}' $file1 >${outputP}.max`;
`awk '\$2==${min}{print \$3"\\t"\$4"\\t"\$5}' $file1  >${outputP}.min`;
open(R,">${outputP}.R")||die"$!";
print R <<EOF;
read.table(\"$file1\",comment.char=\"#\",header=T,sep=\" \",nrows=5)->tab5rows
sapply(tab5rows,class)->classes
read.table(\"$file1\",comment.char=\"#\",header=T,sep=\" \",nrows=$line,colClasses=classes)->data1
read.table(\"$file2\",comment.char=\"#\",header=T,sep=\"\\t\")->data2
read.table(\"$filen\",comment.char=\"#\",header=T,sep=\"\\t\")->datan
read.table(\"$filec\",comment.char=\"#\",header=F,sep=\"\\t\")->slane
read.table(\"$fileq\",comment.char=\"#\",header=T,sep=\"\\t\")->Q20
attach(data1)
attach(data2)
attach(datan)
attach(Q20)
par()->old.par
pdf(\"${outputP}.pdf\")
#jpeg(\"${outputP}_1.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),cex=0.6,las=1)
density(qmean1)->d1
plot(Qual1,ylim=c(0,40),col=\"blue\",type=\"h\",lwd=2,xlab=\"Base site\",ylab=\"Base quality\",main=\"Base quality distribution of Solexa sequencing library I\")
abline(h=20,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->S\",col=\"blue\",bty=\"n\",lty=1)
plot(range(d1\$x),range(d1\$y),xlim=c(0,40),type=\"n\",ylab=\"Density\",xlab=\"Reads mean quality\",main=\"Base quality distribution of Solexa sequencing library II\")
lines(d1,col=\"blue\")
legend(x=\"topleft\",inset=c(.05,.02),\"N=$line\",col=\"red\",,bty=\"n\",lty=1)
#mtext(\"Base quality distribution of Solexa sequencing library\",outer=T,font=2)
#dev.off()

par(old.par)
#jpeg(\"${outputP}_2.jpeg\",width=10000,height=4500,res=1200)
#par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),cex=0.6,las=1)
plot(QF20,ylim=c(0,1),col=\"blue\",type=\"h\",lwd=2,xlab=\"Base site\",ylab=\"Single Q20 Ratio\",main=\"Q20 read ratio in Solexa sequencing library\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->S\",col=\"blue\",,bty=\"n\",lty=1)
plot(CQF20,ylim=c(0,1),col=\"blue\",type=\"h\",lwd=2,xlab=\"Base site\",ylab=\"Cumulative Q20 Ratio\",main=\"Q20 read ratio in Solexa sequencing library\")
abline(h=0.5,col=\"black\",lty=4)
legend(x=\"topright\",inset=c(.05,.02),\"${outputP}->S\",col=\"blue\",,bty=\"n\",lty=1)
#mtext(\"Q20 read ratio in Solexa sequencing library\",outer=T,font=2)
#dev.off()

#par(old.par)
#jpeg(\"${outputP}_3.jpeg\",width=8000,height=6500,res=1200)
#par(cex=0.8)
h <- hist( gccont1, plot = FALSE )
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
#rug( gccont1, col = \"#00000088\" )
box()
#dev.off()
par(old.par)

#jpeg(\"${outputP}_4.jpeg\",width=8000,height=6500,res=1200)
par(cex=0.8)
d1<-density(qmean1[gccont1<=0.2])
d2<-density(qmean1[gccont1>0.2 & gccont1<=0.3])
d3<-density(qmean1[gccont1>0.3 & gccont1<=0.4])
d4<-density(qmean1[gccont1>0.4 & gccont1<=0.5])
d5<-density(qmean1[gccont1>0.5 & gccont1<=0.6])
d6<-density(qmean1[gccont1>0.6])
plot(range(c(d1\$x,d2\$x,d3\$x,d4\$x,d5\$x,d6\$x)),range(c(d1\$y,d2\$y,d3\$y,d4\$y,d5\$y,d6\$y)),type=\"n\",xlab=\"Qmean\",ylab=\"Density\",main=\"Distribution of Qmean with different GC-Content\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
lines(d4,lwd=1,col=5)
lines(d5,lwd=1,col=6)
lines(d6,lwd=1,col=7)
legend(x=\"topleft\",inset=c(.01,.01),c(\"   0<GC<=0.2\",\"0.2<GC<=0.3\",\"0.3<GC<=0.4\",\"0.4<GC<=0.5\",\"0.5<GC<=0.6\",\"0.6<GC<1\"),col=c(2,3,4,5,6,7),lty=c(1,1,1,1,1,1),bty=\"n\")
#dev.off()

#detach (data1)
#rm(data1)
par(old.par)

#jpeg(\"${outputP}_5.jpeg\",width=8000,height=6500,res=1200)
par(cex=0.8)
plot(ratioA1~Position,cex=1.2,ylim=c(0,0.5),xlab=\"Base site\",ylab=\"Ratio\",main=\"GC-Content=$GC1\",type=\"n\")
abline(h=0.1, col = \"gray\", lty=3)
abline(h=0.2, col = \"gray\", lty=3)
abline(h=0.3, col = \"gray\", lty=3)
abline(h=0.4, col = \"gray\", lty=3)
abline(h=0.5, col = \"gray\", lty=3)
lines(x=Position,y=ratioA1,col=\"red\")
lines(x=Position,y=ratioC1,col=\"green\")
lines(x=Position,y=ratioG1,col=\"black\")
lines(x=Position,y=ratioT1,col=\"blue\")
lines(x=Position,y=ratioN1,col=\"purple\")
legend(x=\"top\",inset=c(.05,.02),c(\"A\",\"T\",\"C\",\"G\",\"N\"),col=c(\"red\",\"blue\",\"green\",\"black\",\"purple\"),lty=c(1,1,1,1,1),horiz=T,bty=\"n\")
#dev.off()
par(old.par)

#jpeg(\"${outputP}_6.jpeg\",width=9000,height=5000,res=1200)
par(cex=0.8)
barplot(N1,names=Q,col=\"blue\",axis.lty=1,xlab=\"Qmean\")
legend(x=\"topright\",inset=c(.05,.02),\"N=$N1\",col=\"blue\",pch=15,bty=\"n\")
#dev.off()
par(old.par)



read.table(\"$filed\",header=T,sep=\" \")->data4
attach(data4)

#jpeg(\"${outputP}_9.jpeg\",width=8000,height=6500,res=1200)
#par(cex=0.8)
quantile(duptime[duptime>1],probs=c(0,0.5,1))->qu
qu1<-as.numeric(qu[2])
qu2<-as.numeric(qu[3])
d1<-density(subgc1[duptime==1])
d2<-density(subgc1[duptime>1 & duptime<=qu1])
d3<-density(subgc1[duptime>qu1 & duptime<=qu2])
plot(range(c(d1\$x,d2\$x,d3\$x)),range(c(d1\$y,d2\$y,d3\$y)),type=\"n\",xlab=\"GC-Content\",ylab=\"Density\",main=\"Distribution of GC-Content with different PCR duplication time\")
lines(d1,lwd=1,col=2)
lines(d2,lwd=1,col=3)
lines(d3,lwd=1,col=4)
legend(x=\"topleft\",inset=c(.01,.01),c(\"Time=1\",paste(\"1<Time<=\",qu1),paste(qu1,\"<Time<=\",qu2)),col=c(2,3,4),lty=c(1,1,1),bty=\"n\")
dev.off()




jpeg(\"${outputP}_1.jpeg\",width=7000,height=3000,res=1200)
par(xpd=NA,cex=0.4)
le1=1.1*max(slane\$V2)
barplot(slane\$V2,names=slane\$V1,col=\"red\",axis.lty=1,xlab=\"Flowcell-tile\",main=\"Reads Number in different tiles in lane $lane\",space=0,lwd=0.8)
barplot(slane\$V3,names=slane\$V1,col=\"blue\",space=0,lwd=0.8,add=T)
barplot(slane\$V4,names=slane\$V1,col=\"green3\",space=0,lwd=0.8,add=T)
legend(1,le1,c(\"Q<=20\",\"20<Q<=30\",\"Q>30\"),fill=c(\"red\",\"blue\",\"green3\"),horiz=T,bty=\"n\")
dev.off()


read.table(\"${outputP}.max\",comment.char=\"#\",header=F,sep=\"\\t\")->maxcell
read.table(\"${outputP}.min\",comment.char=\"#\",header=F,sep=\"\\t\")->mincell
co1<-maxcell\$V1
co3<-mincell\$V1
for(i in 1:length(maxcell\$V3)) if (maxcell\$V3[i]>20) co1[i]=\"green\" else co1[i]=\"red\"
for(i in 1:length(mincell\$V3)) if (mincell\$V3[i]>20) co3[i]=\"green\" else co3[i]=\"red\"
jpeg(\"${outputP}_2.jpeg\",width=8000,height=5000,res=1200)
par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(5,4,2,1),cex=0.6)
plot(maxcell\$V1,maxcell\$V2,cex=0.1,col=co1,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $max (MAX)\")
plot(mincell\$V1,mincell\$V2,cex=0.1,col=co3,xlab=\"X-axis\",ylab=\"Y-axis\",main=\"XY coordinate in Solexa lane $lane tile $min (MIN)\")
dev.off()

EOF

close R;

}

######################## fragment over ########################



if(defined $fragmentFQ){

	&process_fragment();

}else{
	&process_pair_end();

}


########### Run the R script ###########
if($plot||$tmpdel){
	
	system ("R  --slave --min-vsize=100M --max-vsize=10000M --min-nsize=5M --max-nsize=1000M -f ${outputP}.R");
	system("rm -f Rplots.pdf")
}

my $compressname=$outputP.".tar.gz" if($compress);

if($tmpdel)
{
	system(" rm -f $file1 $file2 $filen $filec $fileq $filed ${outputP}.R ${outputP}.max ${outputP}.min");
	system(" tar -zcvf $compressname ${outputP}.pdf ${outputP}*.jpeg") if($compress);
}else{
	system(" tar -zcvf $compressname ${outputP}.pdf ${outputP}*.jpeg $file1 $file2 $filen $fileq $filec $filed ${outputP}.R ${outputP}.max ${outputP}.min ") if($compress && $plot);
	system(" tar -zcvf $compressname $file1 $file2 $filen $fileq $filec $filed ${outputP}.R ${outputP}.max ${outputP}.min ") if($compress && !$plot);
}

if($compress){ system("rm -f $file1 $file2 $fileq $filen $filec $filed ${outputP}.R ${outputP}.pdf");}

############################################


my $time2=new Benchmark;

my $timevar1=timediff($time2,$time1);

print "******** Finish ***********\n";

print timestr($timevar1),"\n";

############################################
