mira --project=giardiawbc6 --fasta --job=denovo,genome,accurate,solexa,sanger -GE:not=15 SANGER_SETTINGS -CL:msvs=yes -LR:mxti=yes


cd bchoc_assembly1
cp /path/where/the/data/resides/sequences.fasta bchoc_in.sanger.fasta
cp /path/where/the/data/resides/qualities.someextension bchoc_in.sanger.fasta.qual
ssaha2 -output ssaha2 -kmer 8 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 6 /path/where/the/vector/data/resides/vector.fasta bchoc_in.sanger.fasta > bchoc_ssaha2vectorscreen_in.txt
mira -project=bchoc -fasta -job=denovo,genome,normal,sanger SANGER_SETTINGS -CL:msvs=yes
