#!/bin/csh
#$ -j y
#$ -o /xraid/habitat/mcipriano/seq_align/rio01.log
#$ -N rio
#$ -cwd

alias rm=rm;
setenv BLASTMAT /opt/BioBrew/NCBI/6.1.0/data;
setenv PATH ${PATH}:/opt/BioBrew/Fasta/3.4/bin:/opt/BioBrew/ClustalW/1.8.3/bin/clustalw:/opt/BioBrew/NCBI/6.1.0/bin;
/xraid/habitat/mcipriano/seq_align/riostart.pl;
