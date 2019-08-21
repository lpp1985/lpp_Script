java -XX:-UseGCOverheadLimit -Xmx307200m -jar /home/nfs/SOFTWARE/Other/pilon-1.16.jar    --frags *.bam  --genome $1  --output  $2    --threads 64     --changes   --minmq 40
