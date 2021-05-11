#!/bin/csh
#$ -j y
#$ -o /xraid/habitat/mcipriano/seq_align/caligb01.log
#$ -N alignments
#$ -cwd

alias rm=rm;
setenv BLASTMAT /opt/BioBrew/NCBI/6.1.0/data;
setenv PATH ${PATH}:/opt/BioBrew/Fasta/3.4/bin:/opt/BioBrew/ClustalW/1.8.3/bin/:/opt/BioBrew/NCBI/6.1.0/bin;
hostname;
/xraid/habitat/mcipriano/seq_align/cmsa.pl;
