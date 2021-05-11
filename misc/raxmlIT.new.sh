#!/bin/sh

DIR=$1;
FILE=$2;

mkdir -p $DIR/trees
mkdir -p $DIR/trees/raxml

 mkdir -p $DIR/trees/raxml/full
 cd $DIR/trees/raxml/full
 /data/projects/kinetochore/data2/scripts/changedups.pl $FILE full.fasta
 /data/projects/kinetochore/data2/scripts/alignment_format.pl fasta phylip < full.fasta > full.phy
 # Trim the alignment
 trimal -in full.fasta -out full.trim.fasta -automated1
 /data/projects/kinetochore/data2/scripts/alignment_format.pl fasta phylip < full.trim.fasta > full.trim.phy
 mkdir modeltest
 cd modeltest
 cp ../full.trim.phy .
 /data/projects/kinetochore/data2/scripts/ProteinModelSelection.pl full.trim.phy > bestmodel
 cd ..
 raxml -f a -x 12345 -p 12345 -N 100 -m `cat modeltest/bestmodel` -s full.trim.phy  -n full.phy-RaxML -T 7


