#! /bin/bash 

rdir='/home/eric/R/aiv'


sdir=/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924
    

echo $sdir/


fa=$(ls ${sdir}/*_sampling.fa| sed s/_sampling.fa//)

for a in $fa
    do
    echo "Processing file: ${a}\_sampling.fa"
    echo "store in ${a}\_sampling_ft.tree"
    # echo ${fadir}/gisaid_IRD_merged_${a}\1000.fasta
    # echo ${sdir}/${a}_gtr.tree
    fasttree -gtr -nt $a\_sampling.fa > $a\_sampling_ft.tree
    
    done   