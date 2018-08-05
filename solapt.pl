#!/usr/bin/perl   -w
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

use warnings;
use Getopt::Std;
use Getopt::Long;

my $allusage=qq(
Usage: solapt.pl 

Command:

solsize      => Estimate the solexa sequencing library insert size
soljoin      => Join the forward and reverse read to single longer read for solexa small sequencing library.
solclq       => Statistic and filter the continue low base quality reads in solexa sequencing
solfilter1   => Trim and filter low base quality reads of solexa fragment reads
solfilter2   => Trim and filter low base quality reads of solexa paired-end reads.
soldtrf      => Find and filter the tandem nucleotide repeat  reads in solexa sequencing reads
solin1      => Remove internal adapter for forward reads
solin2      => Remove internal adapter for reverse reads.

\n);

die($allusage) if(@ARGV<1);

my $cmd =$ARGV[0];
my %cmd_hash =(solsize=>\&solsize,soljoin=>\&soljoin,solclq=>\&solclq,solfilter1=>\&solfilter1,solfilter2=>\&solfilter2,soldtrf=>\&soldtrf,solin1=>\&solin1,solin2=>\&solin2);

if(defined($cmd_hash{$cmd})){
	&{$cmd_hash{$cmd}};
}else {
	die ("** Unrecongized command $cmd");
}

################################################################

############# SOLSIZE ##################
sub solsize {

#print $#ARGV;

if($#ARGV !=2) {
	print "\nUsage: solapt.pl  solsize  PE_MAP_file EXP_Insert_Size\n";

print <<"EOF";

This perl script is designed to extract information of soap out, you can access the quality of you 
solexa pair-end reads

Note: the PE file should be sort in reads name! if no, try to give a hash to \@lines!

This version is V1.0
write by Tongwu Zhang at March,23,2010
	
if there is any problem about this script, please contact the author
	
Email: zhangtw\@big.ac.cn
	
EOF


	exit (0);
}


$TMP=$ARGV[1];

$TMP=~s/\..*$//g;

@atcg=("A","T","C","G");

@tandem= Tandembase (2);


open(O,">$TMP.dat")||die"$!";

open(F,$ARGV[1])||die"$!";

%time=();

while(<F>){
 
 @tmp=split(/\t/,$_);

 if(exists $ht{$tmp[0]}) {
	 
   if($time{$tmp[0]} ==1) {

	   $query_last=$ht{$tmp[0]}."\t".$hash_p{$tmp[0]}."\t".$hash_d{$tmp[0]};
	   
	   @lines=();
	   
	   push (@lines,$query_last);
          
         }
          
	   $query=$tmp[7]."\t".$tmp[8]."\t".$tmp[6];

	   $time{$tmp[0]}++;
    
          push (@lines,$query);


	  @{$repeat{$tmp[0]}}=@lines;



   }else {

         $time{$tmp[0]}=1;

         $ht{$tmp[0]}=$tmp[7];             ### %ht name--->bast ref
 
	  $hash_p{$tmp[0]}=$tmp[8];         ### %hash_p name--->position of ref

	  $hash_d{$tmp[0]}=$tmp[6];         ### %hash_d name--->direction of ref
   }

$name=$tmp[0];

if(! exists $dtrf{$name}){

	$sequence=$tmp[1];

 $Lpan=0;

 foreach $repeatsub (@tandem){

	 $pan=multiTRF($sequence,$repeatsub);

	 if($pan) {

		 if($pan>$Lpan) {

			 $Lpan=$pan;

			 $Lrepeat=$repeatsub;
		 }
	 }
 }

$dtrf{$name}=$Lrepeat."\t".$Lpan;

}
$len{$name}=$tmp[5];

       $name=~s/#.*$//;                  ### %list name of mated-pair reads

       $list_name{$name}=1;

       $list{$tmp[0]}=$time{$tmp[0]};      
}
close F;
foreach  $col ( sort keys %list_name){

	$key1=$col."#0/1";

	$key2=$col."#0/2";

if($list{$key1}==1 and $list{$key2}==1 ) {

	if(exists $ht{$key1} and exists $ht{$key2}) {
 
		if($ht{$key1} eq $ht{$key2}) {

			if($hash_p{$key1}>$hash_p{$key2}){

				 $abs=$hash_p{$key1}-$hash_p{$key2}+$len{$key1};

			}elsif($hash_p{$key1}<$hash_p{$key2}){

				$abs=$hash_p{$key2}-$hash_p{$key1}+$len{$key2};

			}else{

				$abs=$len{$key1}>=$len{$key2}?$len{$key1}:$len{$key2}
			}

#		$abs=abs($hash_p{$key1}-$hash_p{$key2});

			if($hash_d{$key1} eq $hash_d{$key2}){

				$header="#";
				
			}else {
				$header="" ;
			}

			$col1=$col."#0/1";

			$col2=$col."#0/2";

			$dtrfv1=$dtrf{$col1};

			$dtrfv2=$dtrf{$col2};

			$dtrfv1=~s/^.*\t//;

			$dtrfv2=~s/^.*\t//;

			$header="#" if($dtrfv1>0.9 || $dtrfv1>0.9) ;

	print  O "$header$col\t$ht{$key1}\t$abs\t$hash_d{$key1}/$hash_d{$key2}\t$dtrf{$col1}\t$dtrf{$col2}\n";
           }
         }
}

else {

if($list{$key1}==1) { @lines=();push(@lines,$ht{$key1}."\t".$hash_p{$key1}."\t".$hash_d{$key1} ) ; @{$repeat{$key1}}=@lines;}

if($list{$key2}==1) { @lines=();push(@lines,$ht{$key2}."\t".$hash_p{$key2}."\t".$hash_d{$key2} ) ; @{$repeat{$key2}}=@lines;}

@arryref1=();

%ref1_pos=%ref1_all=%ref1_name=();

foreach $element ( @{$repeat{$key1}} ) {


    
    @region=split(/\t+/,$element);
     
    $ref1_name{$region[0]}=$region[0];

    push(@arryref1,$region[1]);
    
   @{$ref1_pos{$region[0]}}=@arryref1;
    
    $ref1_all{$region[1]}=$element;
}

@arryref2=();

%ref2_pos=%ref2_all=%ref2_name=();

foreach $element ( @{$repeat{$key1}} ) {

    @region=split(/\t/,$element);
     
    $ref2_name{$region[0]}=$region[0];

    push(@arryref2,$region[1]);
        
  @{$ref2_pos{$region[0]}}=@arryref2;

    $ref2_all{$region[1]}=$element;

 }

foreach $refname ( keys %ref1_name) {
 
	next if(! exists $ref2_name{$refname} );

        @ele1=@{$ref1_pos{$refname}};

        @ele2=@{$ref2_pos{$refname}};

	for($i=0;$i<=$#ele1;$i++) {

		for($j=0;$j<=$#ele2;$j++) {

			if($ele1[$i]>$ele2[$j]){

				$tmpvalue=$ele1[$i]-$ele2[$j]+$len{$key1};

			}elsif($ele1[$i]<$ele2[$j]){

				$tmpvalue=$ele2[$i]-$ele1[$j]+$len{$key2};

			}else{

				$tmpvalue=$len{$key1}>=$len{$key2}?$len{$key1}:$len{$key2};
			}
#			$value=abs($ARGV[1]-abs($ele1[$i]-$ele2[$j]));

			$value=abs($ARGV[1]-$tmpvalue);
	
			if($i==0 and $j==0) { $min=$value;$imin=$i;$jmin=$j;next;}

			if($min > $value) {
				
				$min=$value;

				$imin=$i;

				$jmin=$j;
			}
		}
	}

@last1=split(/\t+/, $ref1_all{$ele1[$imin]});

@last2=split(/\t+/, $ref2_all{$ele2[$jmin]});

if($last1[1]>$last2[1]){

	$abs=$last1[1]-$last2[1]+$len{$key1};

}elsif($last1[1]<$last2[1]){

	$abs=$last2[1]-$last1[1]+$len{$key2}
}else{

	$abs=$len{$key1}>=$len{$key2}?$len{$key1}:$len{$key2};
}
#$abs=abs($last1[1]-$last2[1]);

# $gam=$abs/$ARGV[1];

#############################Changed by tongwu zhang at 27/3/2010###########

##### if($gam<1.5 && $gam >0.5)  {

#	print  O "$col\t$last1[0]\t$abs\t$last1[2]/$last2[2]\n";
#  }
#################################

$col1=$col."#0/1";

$col2=$col."#0/2";

print O "#$col\t$last1[0]\t$abs\t$last1[2]/$last2[2]\t$dtrf{$col1}\t$dtrf{$col2}\n";
		
}


}

}
close O;
sub multiTRF {
	
	($seq,$mutibase)=@_;

	$seq=uc($seq);

	$length=length($seq);

	$num=$seq=~s/$mutibase/$mutibase/g;

	$value=$num*2/$length;

	return  $value;
}
sub Tandembase {

	($num)=@_;

	if($num>1) {

		@last=Tandembase ($num-1);

		$#tmpsub=-1;

		for($i=0;$i<=$#last;$i++){

			for($j=0;$j<=3;$j++){

				push( @tmpsub,$last[$i].$atcg[$j]);

				}
		}
		
		return @tmpsub;
	}

	if($num==1) {

		return @atcg;

		}
}
################# Using R to plot the Distrubtion of Solexa libarary Insert size ########

####### Write the R script file ############


##system(" grep -v \"^#\" $TMP.dat >$TMP.dat.uniq");

open(R,">$TMP.R")||die "$!";

#print R "library(lattice)\n";

print R "pdf(\"$TMP.pdf\")\n";
#print R "jpeg(\"$TMP.jpeg\",width=9000,height=5000,res=1200)\n";

print R "read.table(\"$TMP.dat\",comment.char=\"\")->x\n";

print R "read.table(\"$TMP.dat\")->u\n";

print R "d1<-density(x\$V3)\n";

print R "d2<-density(u\$V3)\n";

$m1=`wc -l  $TMP.dat|cut -d " " -f 1`;

$m2=`grep -v "^#" $TMP.dat|wc -l |cut -d " " -f 1`;

chomp $m1;

chomp $m2;

#print R "plot.new()\n";

#print R "par(mfcol=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,2),cex=0.65)\n";

print R "plot(range(d1\$x,d2\$x),range(d1\$y,d2\$y),type=\"n\",xlab=\"Insert Size (bp)\",ylab=\"Density\",main=\"Distribution of library insert size in Solexa sequencing\")\n";

print R "lines(d1,col=\"blue\")\n";

print R "segments(d1\$x[which.max(d1\$y)],0,d1\$x[which.max(d1\$y)],max(d1\$y),col=\"green3\",lty=2)\n";

print R "text(d1\$x[which.max(d1\$y)],max(d1\$y),labels=round(d1\$x[which.max(d1\$y)],1),pos=4)\n";

print R "legend(x=\"topright\",inset=c(.05,.02),c(\"${TMP}_all\",\"N=$m1\"),col=c(\"blue\",\"blue\"),bty=\"n\",lty=c(1,0))\n";

print R "rug(x\$V3,col=\"blue\",lwd=0.5)\n";

print R "plot(range(d1\$x,d2\$x),range(d1\$y,d2\$y),type=\"n\",xlab=\"Insert Size (bp)\",ylab=\"Density\",main=\"Distribution of library insert size in Solexa sequencing\")\n";

print R "lines(d2,col=\"red\")\n";

print R "segments(d2\$x[which.max(d2\$y)],0,d2\$x[which.max(d2\$y)],max(d2\$y),col=\"green3\",lty=2)\n";

print R "text(d2\$x[which.max(d2\$y)],max(d2\$y),labels=round(d2\$x[which.max(d2\$y)],1),pos=4)\n";

print R "legend(x=\"topright\",inset=c(.05,.02),c(\"${TMP}_uniq\",\"N=$m2\"),col=c(\"red\",\"red\"),bty=\"n\",lty=c(1,0))\n";

print R "rug(u\$V3,col=\"red\",lwd=0.5)\n";

#print R "mtext(\"Distribution of library insert size in Solexa sequencing\",outer=T,font=4)\n";

print R "dev.off()\n";

close R;

######## run the R script file ##########
 system ("R  --slave -f $TMP.R");


 }

################## OVER SOLSIZE ########################

#################  SOLJOIN #############################

sub soljoin {


$CURRENT_VERSION="0.1";

$PROGRAM_NAME=$0;

$PROGRAM_NAME=~s|.*/||;

GetOptions('f|forward_fq=s'      => \$forwardFQ,
	   'r|reverse_fq=s'      => \$reverseFQ,
	   'o|output_prefix=s'   => \$outputP, 
	   'm|min_overlap:i'     => \$min_overlap,
	   'x|max_overlap:i'     => \$max_overlap,
	   'n|mismatch:i'        => \$misMatch,
	   'q|format:i'		 => \$format,
	   'h|help'              => \$help,
	   'v|version'           => \$version,
	  );

if($version){

	print STDERR "$PROGRAM_NAME, Version $CURRENT_VERSION\n";

	exit;
}

if($help || !defined $forwardFQ || !defined $reverseFQ){

	&PrintUsage();

	exit;
}

$format ||=64;
$misMatch||= 1;
$min_overlap ||= 10;
$maxstar=0;
$maxstar=1 unless($max_overlap);

sub PrintUsage {

	print STDERR <<'END';

Usage: solapt.pl soljoin [options] -f <FORWARD_FILE> -r <REVERSE_FILE> -o <OUTPUT_PREFIX>

Options: 
          -f or --forward_fqs      <str>     : Forward fastq file 
	  
	  -r or --reverse_fq       <str>     : Reverse fastq file

	  -o or --output           <str>     : Output prefix

	  -m or --min_overlap      <int>     : Min overlap between forward and reverse reads, default 10

	  -x or --max_ovrlap       <int>     : Max overlap between forward and reverse reads,default minmum read length.

	  -n or --mismatch         <int>     : Mismatch between overlap,default 1

	  -q or --format	   <int>     : fastq format: set 64 for illumina or solexa and 33 for sanger fastq,default 64

	  -v or --version                    : Prints version of the program

	  -h or --help                       : Prints this usage summary

Note that default values are not stringent and will result in very littel filtering
the max_overlap is set to larger than length of FR reads when you want to see more
overlap between FR reads. this program use R package to plot the density distribution plot
of joined length and overlap length.if find any question and bugs, please sent it to zhangtw
Email: zhangtw@big.ac.cn

END
	return;

}

if($format !=64 && !$format !=33){
	print "\nError, please set the correct fastq, 64 for solexa or illumina and 33 for sange fastq\n";
	exit 0;
}


open(F1,$forwardFQ)||die"cann't open the input fastq forward file '$forwardFQ':$!\n";

open(F2,$reverseFQ)||die"cann't open the input fastq reverse file '$reverseFQ': $!\n";

$OUT1=$outputP.".fqa";

$OUT2=$outputP.".fqb";

$OUT=$outputP.".fqc";

$OUTT=$outputP.".inf";

open(O,">$OUT")||die"cann't open the output file '$OUT':$!\n";

open(O1,">$OUT1")||die"cann't open the output file '$OUT1':$!\n";

open(O2,">$OUT2")||die"cann't open the output file '$OUT2':$!\n";

open(T,">$OUTT")||die"cann't open the output file '$OUTT':$!\n";

print T "length_combine\tlength_overlap\tmismatch\n";

$comn=$pair=0;

while(!eof(F1) && !eof(F2)){

	$name1=<F1>;

	$name2=<F2>;

	$name=$name1;

	chomp $name;

	$name=~s/#.*$//;

	$seq1=<F1>;

	$seq2=<F2>;

	chomp $seq1;

	chomp $seq2;

	if($maxstar){

		$max_overlap=length($seq1) > length($seq2) ? length($seq2) : length($seq1)
	}

	<F1>;<F2>;

	$qual1=<F1>;

	$qual2=<F2>;

	chomp $qual1;

	chomp $qual2;

	$revcom_seq2=reverse $seq2;

	$revcom_seq2=~tr/ATCG/TAGC/;

	$rev_qual2=reverse $qual2;

	$len=length($seq1);

#$len2=length($seq2);

	@seq_arry1=split(//,$seq1);

	@seq_arry2=split(//,$revcom_seq2);

	@qu_arry1=map(ord($_)-$format,split(//,$qual1));

	@qu_arry2=map(ord($_)-$format,split(//,$rev_qual2));

	$combine=0;

	for($i=($len-$max_overlap);$i<=($len-$min_overlap);$i++){

		$i=0 if($i<0);

		$j=$i-1;$k=-1;$mis=0;

		while($mis<=$misMatch && $j<$len-1){

			$j++;$k++;

			$mis++ if($seq_arry1[$j] ne $seq_arry2[$k]);

		}


		
		if($j==$len-1 && $mis<=$misMatch ){

			
			$combine=1;

			$mseq="";

			$mqu="";

			$k=0;

			for($s=$i;$s<=$len-1;$s++){

				if($seq_arry1[$s] ne $seq_arry2[$k]){

					$mseq.=$seq_arry1[$s];

					$maxq=$qu_arry1[$s]>=$qu_arry2[$k]?$qu_arry1[$s]:$qu_arry2[$k];

					$minq=$qu_arry1[$s]>=$qu_arry2[$k]?$qu_arry2[$k]:$qu_arry1[$s];

					if($minq>=30){
						$w=0.4;
					}elsif($minq>=20){
						$w=0.3;
					}elsif($minq>=10){
						$w=0.2;
					}else{  $w=0.1;
					}

					$qvalue=$maxq+$minq*$w;

					$qvalue=40 if( $qvalue >40);

					$mqu.=chr($qvalue+$format);

				}else{
					$mseq.=$seq_arry1[$s] if($qu_arry1[$s]>=$qu_arry2[$k]);

					$mseq.=$seq_arry2[$k] if($qu_arry1[$s]<$qu_arry2[$k]);

					$qvalue=abs($qu_arry1[$s]-$qu_arry2[$k]);

					$qvalue=2 if($qvalue<2);

					$mqu.=chr($qvalue+$format);

				}

				$k++;
			}

			$seq=substr($seq1,0,$i).$mseq.substr($revcom_seq2,$k);

			$qual=substr($qual1,0,$i).$mqu.substr($rev_qual2,$k);

			last;

		}

	}


	if($combine){

		$comn++;

		$lencom=length($seq);

		$lenoverlap=length($mseq);

		$name.="#".$lencom.":".$lenoverlap;

		print O "$name\n$seq\n+\n$qual\n";

		print T "$lencom\t$lenoverlap\t$mis\n";

		

	}else{ 
		$pair++;

		print O1 "$name1$seq1\n+\n$qual1\n";

		print O2 "$name2$seq2\n+\n$qual2\n";

	}


}

$total=$comn+$pair;

$percent=sprintf '%.2f', 100*$comn/$total;

print T "#Total\tcombined\tuncombined\tpercent\n";

print T "#$total\t$comn\t$pair\t$percent\n";

close F1;

close F2;

close O;

close O1;

close O2;

close T;

$OUTR=$outputP.".R";

open(R,">$OUTR")||die"$!";

print R "pdf(\"${outputP}.pdf\")\n";
#print R "jpeg(\"${outputP}.jpeg\",width=10000,height=4500,res=1200)\n";

#print R "par(mfrow=c(1,2),oma=c(0,0,2,0),mar=c(5,4,2,1),cex=0.6)\n";

print R "read.table(\"$OUTT\",comment.char=\"#\",header=T,sep=\"\\t\")->dat\n";

print R "d1<-density(dat\$length_combine)\n";

print R "d2<-density(dat\$length_overlap)\n";

print R "plot(d1\$x,d1\$y,type=\"l\",col=\"blue\",cex=0.6,xlab=\"Length of joined FR reads (bp)\",ylab=\"Density\",main=\"Length distribution of joined FR reads\")\n";

print R "plot(d2\$x,d2\$y,type=\"l\",col=\"red\",cex=0.6,xlab=\"Length of overlap (bp)\",ylab=\"Density\",main=\"Length distribution of overlap\")\n";

#print R "mtext(\"Length distribution of joined FR reads and overlap\",outer=T,font=2,cex=0.6)\n";

print R "dev.off()\n";

system ("R --min-vsize=100M --max-vsize=10000M --min-nsize=5M --max-nsize=1000M -f $OUTR");

#`rm $OUTR`;


}
##################### OVER SOLJOIN #######################
#####################  SOLCLQ  ########################

sub solclq{

if(@ARGV<2){

	print "Usage: solapt.pl solclq <ILLUMINA_FASTQ> <OUTPUTFILE> <PHRAP_QUALITY> <HOLD_LENGTH>\n";
	exit;

}

$phrap_q=$ARGV[3]||20;

$read=0;

$read_last=0;

$holdlength=$ARGV[4]||20;

open(O,">$ARGV[2]")||die"$!";

open(F1,"$ARGV[1]")||die"$!";

while(!eof(F1)){

	$name=<F1>;

	chomp $name;

	$seq=<F1>;

	chomp $seq;

	<F1>;

	$qual=<F1>;

	chomp $qual;

	@q=map((ord($_)-64),split(//,$qual));

	$read++;

	$count=0;

	$find=0;

	for($i=0;$i<=$#q;$i++){

		$num=0;

		while($q[$i] < $phrap_q){

			$num++;

			$flq=$i if($num==1 && $find== 0);

			$i++;
			
			last if($i>$#q);
		}

		 if($num>=2){

			 $find=1;
			 
			 $count++;

		 }

	}

	if($find){

		if($flq>=$holdlength){

			$seq1=substr($seq,0,$flq);

			$qual1=substr($qual,0,$flq);

			print O "$name\n$seq1\n+\n$qual1\n";

			$read_last++;

		}
	}else{
		print O "$name\n$seq\n+\n$qual\n";

		$read_last++;
	}



	if(!exists $time{$count}){
		
		$time{$count}=1;

	}else{

		$time{$count}++;

	}
}

$per=100*$time{0}/$read;

$per1=100*$read_last/$read;

print "$ARGV[0]\n";

print "All Reads:$read\tGood Reads:$time{0}\tPercent:$per\n";

print "All Reads:$read\tLast Reads:\t$read_last\tPercent:$per1\n";

foreach $key (sort {$a<=>$b} keys %time){

	print "$key\t$time{$key}\n";

}

close F1;

}
################## OVER SOLCLQ #################
#################  SOLFILTER1 #################

sub solfilter1{

my %opts;
GetOptions(\%opts,"i=s","o=s","n=i","q=i","l=i","h!");
my $Ver="0.1";
my $usage=<<"USAGE";
Program:solapt.pl solfilter1
Version:$Ver
Contact:Tongwu Zhang <zhangtw\@big.ac.cn>

  Usage:solapt.pl solfilter1 
  	   -i <input> -o <output_prefix> 
           -q <low_quality,default:20> 
           -n <number of low_qualtiy,default:3>
           -l <length_to_trimmed,default:35>
	   -h Display this usage information

USAGE

die $usage if ( $opts{h});
die "Error of the parameter!!\n\n$usage\n" if(!defined $opts{i}||!defined $opts{o});
my $input=$opts{i};
my $lowq=$opts{q}||20;
my $mism=$opts{n}||3;
my $output=$opts{o};
my $data1=0;
my $data2=0;
my $read1=0;
my $read2=0;
my $length=$opts{l}||35;

open(OUT,">${output}.flq")||die"cann't open ${output}.clean\n";
open(IN,"$input")||die"$!";
while(!eof(IN)){

	$id=<IN>;
	$seq=<IN>;
	<IN>;
	$qual=<IN>;
	chomp $qual;
	chomp $seq;

	$read1++;
	$data1+=length($seq);

	$seq1=substr($seq,0,$length);
	$qual1=substr($qual,0,$length);

	@Q=map(ord($_)-64,split(//,$qual1));
	@S=split(//,$seq1);
	$low=0;
	$posbak=-1;
	$subfind=0;

	for($i=0;$i<$length;$i++){

		if($Q[$i]<$lowq || $S[$i]=~m/N/i){

			$low++;

			$sub=$i-$posbak;

			if($sub==1){ $subfind=1;last}else{$posbak=$i;}
		}


			
	}
	
	next if($low>$mism|| $subfind);
	
	print OUT "$id$seq1\n+\n$qual1\n";
	$read2++;
	$data2+=length($seq1);

	
}
close IN;
close OUT;

print "$output\tInput_reads:$read1\tOutput_reads:$read2\tInput_data:$data1\tOutput_data:$data2\n";

}

##################### OVER SOLFILTER ######################
####################  SOLFILTER2 #########################

sub solfilter2 {

my %opts;
GetOptions(\%opts,"f=s","r=s","o=s","n=i","q=i","l=i","h!");
my $Ver="0.1";
my $usage=<<"USAGE";
Program:solapt.pl solfilter2
Version:$Ver
Contact:Tongwu Zhang <zhangtw\@big.ac.cn>

  Usage:solapt.pl solfilter2
           -f <forward_input> 
           -r <reverse_input> 
	   -o <output_prefix>
           -q <low_quality,default:20> 
           -n <number of low_qualtiy,default:3>
           -l <min_length,default:35>
	   -h Display this usage information

USAGE

die $usage if ( $opts{h});
die "Error of the parameter!!\n\n$usage\n" if(!defined $opts{f}||!defined $opts{r}||!defined $opts{o});
my $finput=$opts{f};
my $rinput=$opts{r};
local $lowq2=$opts{q}||20;
local $mism2=$opts{n}||3;
my $output=$opts{o};
my $fdata1=0;
my $fdata2=0;
my $fread1=0;
my $fread2=0;
my $rdata1=0;
my $rdata2=0;
my $rread1=0;
my $rread2=0;
my $length1;
my $length2;
my $mlength=$opts{l}||35;

open(P,">${output}.flq.pair")||die"cann't open ${output}.flq.pair\n";
open(S,">${output}.flq.single")||die"cann't open ${output}.flq.single\n";

open(INF,"$finput")||die"$!";
open(INR,"$rinput")||die"$!";
while(!eof(INF) && !eof(INR)){
	$id=<INF>;
	$seq=<INF>;
	<INF>;
	$qual=<INF>;
	chomp $qual;
	chomp $seq;
	$fread1++;
	$fdata1+=length($seq);
	$length1= &getlength($seq,$qual);
	if($length1>=$mlength){
		$pr1=$id.substr($seq,0,$length1)."\n+\n".substr($qual,0,$length1)."\n";
		$fread2++;
		$fdata2+=$length1;
	}
	
	$id=<INR>;
	$seq=<INR>;
	<INR>;
	$qual=<INR>;
        chomp $qual;
	chomp $seq;
	$rread1++;
	$rdata1+=length($seq);
	$length2= &getlength($seq,$qual);
	if($length2>=$mlength){
		$pr2=$id.substr($seq,0,$length2)."\n+\n".substr($qual,0,$length2)."\n";
		$rread2++;
		$rdata2+=$length2;
	}
	
	if($length1>=$mlength && $length2>=$mlength){
		print P "$pr1$pr2";
	}elsif($length1>=$mlength){
		print S "$pr1";
	}elsif($length2>=$mlength){ 
		print S "$pr2";}
}


sub getlength {

	my $seqt;
	my $qualt;
	my @Q=();
	my @S=();
	
	($seqt,$qualt)=@_;
	@Q=map(ord($_)-64,split(//,$qualt));
	@S=split(//,$seqt);
	my $low=0;
	my $posbak=-1;
	my $subfind=0;
	my $i;

	for( $i=0;$i<length($seqt);$i++){

		if($Q[$i]<$lowq2 || $S[$i]!~m/[atcg]/i){

			$low++;

			last if($low>$mism2);

			$sub=$i-$posbak;

			if($sub==1){ $subfind=1;last}else{$posbak=$i;}
		}


			
	}

	return $i;

}

close INF;
close INR;
close P;
close S;


print "$finput\tInput_reads $fread1\tOutput_reads $fread2\tInput_data $fdata1\tOutput_data\t$fdata2\n";
print "$rinput\tInput_reads $rread1\tOutput_reads $rread2\tInput_data $rdata1\tOutput_data\t$rdata2\n";


}

################### OVER SOLFILER2 ###################

###################  SOLDTRF ######################

sub soldtrf {

if($#ARGV < 3) {
	
	print "Usage: solapt.pl  Dtrf Solexa_fastq_file Repeat_base_Number cut_off output_prefix\n";

	exit 1;
}

if($#ARGV == 4){

	$Num=$ARGV[2];

	open($filehand,"$ARGV[1]")||die"$!";

	$output=$ARGV[4];

	$cutoff=$ARGV[3];

}elsif($#ARGV == 3){

	$filehand="STDIN";

	$Num=$ARGV[1];

	$output=$ARGV[3];

	$cutoff=$ARGV[2];
}else{

	print "Error with the options, please use the help\n";

	exit;

}

$repeatfile=$output."_re";

$uniqfile=$output."_un";

$report=$output.".rereport";

open(O1,">$repeatfile")||die "$!";

open(O2,">$uniqfile")||die "$!";

open(O3,">$report")||die "$!";

@atcg=("A","T","C","G");

@tandem= Tandembase1 ($Num);

while(<$filehand>){

	chomp;

	s/^\s+//g;

	s/\s+$//g;

	$seqname=$_;

	$sequence=<$filehand>;

	chomp $sequence;

	<$filehand>;

	$qual=<$filehand>;

	chomp $qual;

        $Lpan=0;

	foreach $repeat (@tandem){

			$pan=multiTRF1($sequence,$repeat);

			if($pan) {
                                 
                                 if($pan>$Lpan) { 
                                                  
                                                  $Lpan=$pan;
                                                  $Lrepeat=$repeat;
                                      }




			}

	
             	}

	print O1 "$seqname\n$sequence\n+\n$qual\n"  if($Lpan >$cutoff);

	print O2 "$seqname\n$sequence\n+\n$qual\n"  if($Lpan <=$cutoff);

	print O3 "$seqname\t$Lrepeat\t$Lpan\n";


}

sub multiTRF1 {

    ($seq,$mutibase)=@_;

    $seq=uc($seq);
    
    $length=length($seq);

    $num=$seq=~s/$mutibase/$mutibase/g;

    $value=$num*$Num/$length;

    return  $value;



}



sub Tandembase1 {

		($num)=@_;

		if($num>1) {

			@last=Tandembase ($num-1);

			$#tmp=-1;

			for($i=0;$i<=$#last;$i++){

					for($j=0;$j<=3;$j++){

								push( @tmp,$last[$i].$atcg[$j]);

					}}
			return @tmp;
		}

		if($num==1) {

			return @atcg;
		}
	}

if($#ARGV == 3) { close $filehand;}

close O1;
close O2;
close O3;

}

######################## OVER SOLDTRF ###################
####################### SOLIN1 ########################

sub solin1{


if( $#ARGV <3){

	print "Usage: <STDIN>|solapt.pl solin1 <adapter_sequence> <STDOUT> option: <0 or NUM>\n";
	
	print "Usage: solapt.pl  solin1  <INPUT_FILE> <adapter_sequence> <STDOUT> option: <0 or NUM>\n";

	exit;
}

if($#ARGV == 3){

$adapter=$ARGV[1];

$filehand="STDIN";

$outfile=$ARGV[2];

$filter=$ARGV[3];

}elsif($#ARGV == 4){

	$adapter=$ARGV[2];

	$outfile=$ARGV[3];

	$filter=$ARGV[4];

	open($filehand,"$ARGV[1]")||die"$!";

}else{

	print "Error with the option, please check the option with help!\n";

}

open(O,">$outfile")||die"$!";

open(O1,">${outfile}.add1")||die"$!";

open(O2,">${outfile}.add2")||die"$!";

$inputnum=0;

$outputnum=0;

$addpair=0;

$addsingle=0;

while(<$filehand>){
	
	chomp;

	$name=$_;

	$seq=<$filehand>;

	chomp $seq;

	$seq_bc=$seq;

	$seq=uc $seq;

	<$filehand>;

	$qual=<$filehand>;

	chomp $qual;

	$inputnum++;

	next if($seq=~s/$adapter/$adapter/g >1);

	$pos=index($seq,$adapter);

	if($pos == -1){

		 print O "$name\n$seq_bc\n+\n$qual\n";

		 $outputnum++;

		 next;

	}

	next if($filter == 0);

	$len1=$pos;

	$len2=length($seq)-$len1-length($adapter);

	if($len1>$filter){

		$name1=$name."/1";

		$seq1=substr($seq_bc,0,$len1);

		$qual1=substr($qual,0,$len1);

		print O1 "$name1\n$seq1\n+\n$qual1\n";

	}
	
	if($len2>$filter){

		$name2=$name."/2";

		$seq2=substr($seq_bc,$pos+length($adapter),$len2);

		$seq2=reverse $seq2;

		$seq2=~tr/ATCGNatcgn/TAGCNtagcn/;

		$qual2=substr($qual,$pos+length($adapter),$len2);

		$qual2=reverse $qual2;
		
		print O2 "$name2\n$seq2\n+\n$qual2\n";

	}

	if($len1>$filter && $len2>$filter){

		        $addpair++;

	}elsif($len1>$filter || $len2>$filter){

		$addsingle++;
	}

  }

close O1;

close O2;

close O;

print "Input Reads:$inputnum\n";
print "Output Reads:$outputnum\n";
print "Add paired reads:$addpair\n";
print "Add single reads:$addsingle\n";

}
###################### OVER SOLIN1 #####################
######################  SOLIN2 ########################

sub solin2{

if( $#ARGV <3){

	print "Usage: <STDIN>|solapt.pl solin2  <adapter_sequence> <STDOUT> option: <0 or NUM>\n";
	
	print "Usage: solapt.pl solin2 <INPUT_FILE> <adapter_sequence> <STDOUT> option: <0 or NUM>\n";

	exit;
}

if($#ARGV == 3){

$adapter=$ARGV[1];

$filehand="STDIN";

$outfile=$ARGV[2];

$filter=$ARGV[3];

}elsif($#ARGV == 4){

	$adapter=$ARGV[2];

	$outfile=$ARGV[3];

	$filter=$ARGV[4];

	open($filehand,"$ARGV[1]")||die"$!";

}else{

	print "Error with the option, please check the option with help!\n";

}

open(O,">$outfile")||die"$!";

open(O1,">${outfile}.add1")||die"$!";

open(O2,">${outfile}.add2")||die"$!";

$inputnum=0;

$outputnum=0;

$addpair=0;

$addsingle=0;

while(<$filehand>){
	
	chomp;

	$name=$_;

	$seq=<$filehand>;

	chomp $seq;

	$seq_bc=$seq;

	$seq=uc $seq;

	<$filehand>;

	$qual=<$filehand>;

	chomp $qual;

	$inputnum++;

	next if($seq=~s/$adapter/$adapter/g >1);

	$pos=index($seq,$adapter);

	if($pos == -1){

		 print O "$name\n$seq_bc\n+\n$qual\n";

		 $outputnum++;

		 next;

	}

	next if($filter == 0);

	$len1=$pos;

	$len2=length($seq)-$len1-length($adapter);

	if($len1>$filter){

		$name1=$name."/1";

		$seq1=substr($seq_bc,0,$len1);

		$seq1=reverse $seq1;

		$seq1=~tr/ATCGNatcgn/TAGCNtagcn/;

		$qual1=substr($qual,0,$len1);

		$qual1=reverse $qual1;


		print O1 "$name1\n$seq1\n+\n$qual1\n";

	}
	
	if($len2>$filter){

		$name2=$name."/2";

		$seq2=substr($seq_bc,$pos+length($adapter),$len2);

		$qual2=substr($qual,$pos+length($adapter),$len2);

		print O2 "$name2\n$seq2\n+\n$qual2\n";

	}

	if($len1>$filter && $len2>$filter){

		        $addpair++;

	}elsif($len1>$filter || $len2>$filter){

		$addsingle++;
	}

  }

close O1;

close O2;

close O;

print "Input Reads:$inputnum\n";
print "Output Reads:$outputnum\n";
print "Add paired reads:$addpair\n";
print "Add single reads:$addsingle\n";

}

############# OVER SOLIN2 ############################

################# OVER ALL ####################





	




