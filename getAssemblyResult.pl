opendir(NOW,".");
@filelist=readdir(NOW);
open(OUT,">assemblyResult.sta");
foreach $element (@filelist){
if($element=~/run(.+)\.log$/){
   $kmer=$1;
   $count=1;
   open(IN,"$element");
   while(<IN>){
    if(/scaffolds from/){
       $isprint=1;
        print OUT "$kmer\n$_";
        next;
      }
    if($isprint==1){
      $count++;
       print OUT "$_";
      }
    if($count>4){
      $isprint=0;
      last;
      }
   }
   close IN;
 }
}
close OUT;
closedir NOW;
