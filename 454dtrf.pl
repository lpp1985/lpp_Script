#!/usr/bin/perl  -w
### write this perl script for 454 seqeuncing dinucletide tandem repeat finde in fastq format,and you can 
### change to find three,four,or more hight dinucletide TRF!!
### Tongwu Zhang
### Key Laboratory of Genome Sciences and Information, Beijing Institute of Genomics,
### Email contact <zhangtw@big.ac.cn>
### January 2011

# Released under GNU General Public License version 3

if($#ARGV != 3) {
	
	print "Usage: perl 454dtrf.pl <sff2fna_file> <outprefix> <Repeat_Number> <cut_off>\n";
	exit 1;
}

$Num=$ARGV[2];

open(F,"$ARGV[0]")||die"$!";
open(O,">$ARGV[1]_di.dat")||die"$!";
print O "Read_name\tRepeat_base\tdi\n";

$mincutoff = $ARGV[3] || 0;
@atcg=("A","T","C","G");
$finder=0;
$first=1;

$sequence="";

 @tandem= Tandembase ($Num);

while(<F>){

	chomp;

	s/^\s+//g;

	s/\s+$//g;

	if(/^>/) {

		m/^(>\S+)/;
		
		if($first){$seqname=$1;$first=0; next;}

	        $last=$seqname;

		$seqname=$1;
		
		$finder=1;

		}
         
	if(! $finder){ $sequence.=$_; next;

	}else{
		
                $Lpan=0;

		foreach $repeat (@tandem){

			$pan=multiTRF($sequence,$repeat);

			if($pan) {
                                
                                 if($pan>$Lpan) { 
                                                  
                                                  $Lpan=$pan;
                                                  $Lrepeat=$repeat;
                                      }
			}
             	}
#				print O "$last\t$Lrepeat\t$Lpan\n$sequence\n" if($Lpan);
				print O "$last\t$Lrepeat\t$Lpan\n" if($Lpan);
	}

	$sequence="";

	$finder=0;

}

 $Lpan=0;

  foreach $repeat (@tandem){

	  $pan=multiTRF($sequence,$repeat);

	  if($pan) {
	                  
		  if($pan>$Lpan) {
			  
			  $Lpan=$pan;
		
			  $Lrepeat=$repeat;
		  }
       }
}

#print O "$seqname\t$Lrepeat\t$Lpan\n$sequence\n" if($Lpan);
print O "$last\t$Lrepeat\t$Lpan\n" if($Lpan);

close F;
close O;

sub multiTRF {

    ($seq,$mutibase)=@_;

    $seq=uc($seq);
    
    $length=length($seq);

    @bases=split(//,$seq);

    @match=split(//,$mutibase);

    $num=0;

    for( $i=0;$i<=$#bases-$Num+1;$i++){
      
	$find=1;

         for($j=0;$j<=$#match;$j++){

		 $ni=$i+$j;

		 if($bases[$ni] ne $match[$j]) { $find=0;last;}
	 }

	 if($find) { $num++;$i=$ni;}
    }
    $value=$num*$Num/$length;

    if($value >= $mincutoff ) {

         return  $value;
    }else{
        
	 return  0;
    }
}


sub Tandembase {
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


$NAME=$ARGV[1];
$m1=`wc -l  ${NAME}_di.dat|cut -d " " -f 1`;
chomp $m1;
open(R, ">${NAME}_di.R")||die"$!";
print R "pdf(\"${NAME}.pdf\")\n";
#print R  "jpeg(\"${NAME}_di.jpeg\",width=8000,heigh=6500,res=1200)\n";
#print R  "plot.new()\n";
#print O  "library(Rcmdr)\n";
#print R "par(cex=0.8)\n";
print R  "read.table(\"${NAME}_di.dat\",head=T)->f\n";
print R  "hist(f\$di,freq=T,xlim=c(0,1), breaks=100, col=\"#0000FF\",xlab=\"Ratio of dinucleotide sequence\",ylab=\"Frequency\",main=\"Distribution of dinucleotide ratio in reads of GS FLX/454 sequencing\")\n";
print R "legend(x=\"topright\",inset=c(.05,.02),c(\"${NAME}\",\"N=$m1\"),col=c(\"blue\",\"blue\"),bty=\"n\",lty=c(0,0))\n";
print R "box()\n";
print R  "dev.off()\n";
close R;
close R;
system(qq(R -f ${NAME}_di.R));
