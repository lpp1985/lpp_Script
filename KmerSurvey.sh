Data=$1
/home/lpp/pub/SOFTWARE/Pacbio/canu/Linux-amd64/bin/gatekeeperCreate -o ./$Data.gkpStore ./$Data.gkp
/pub/SOFTWARE/Pacbio/canu/Linux-amd64/bin/meryl   -B -C -L 2 -v -m 17 -threads 64 -memory 65536   -s ./$Data.gkpStore   -o ./$Data.ms17
/pub/SOFTWARE/Pacbio/canu/Linux-amd64/bin/meryl      -Dh -s /$Data.ms17  >  $Data.ms17.histogram
