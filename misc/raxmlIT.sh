#!/bin/sh

DIR=$1;
NAME=$2;

mkdir -p $DIR/trees
mkdir -p $DIR/trees/raxml

 mkdir $DIR/trees/raxml/full
 /data/projects/kinetochore/data2/scripts/alignment_format.pl fasta phylip < $DIR/$NAME > $DIR/trees/raxml/full/full.phy
 cd $DIR/trees/raxml/full
 raxml -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAJTT -s full.phy  -n full.phy-RaxML -T 7


