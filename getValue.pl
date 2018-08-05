open(IN,"@ARGV[0]");
open(OUT,">@ARGV[2]");
$cout=0;
$countReads=0;
@valueArray=(0)x@ARGV[1];

while(<IN>){
        $cout++;
        if($cout==4){
                chomp;
                @tmp=split(//,$_);
                for ($i=0;$i<=$#valueArray;$i++){
             $qValue = ord($tmp[$i]) - 64;
             $valueArray[$i]=$valueArray[$i]+$qValue;
            }
          $countReads++;
          $cout=0;
                }
        }

        print OUT "Bp Qvalue\n";

for ($i=0;$i<=$#valueArray;$i++){
        $qValueAverage=$valueArray[$i]/$countReads;
        $qValueAverage=sprintf "%d",$qValueAverage;
        $j=$i+1;
        print OUT "$j $qValueAverage\n";
        }
close IN;
close OUT;
