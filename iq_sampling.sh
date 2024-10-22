#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge/0307/sampling_subset/0922


fa=("PB1" 'PA' 'NP' 'MP')

for a in "${fa[@]}"
do
    echo ~/Analysis/aiv/merge/0307/distance/iq/$a\_5p_sampling.fa
    # echo $a\_aligned_fill_trimed_outlier.fa
    # iqtree -s $a\_aligned_fill_trimed.fa -st DNA -m MFP -nt 10 -bb 1000 
    # iqtree -s  ~/Analysis/aiv/merge/0307/distance/iq/$a\_5p_sampling.fa -st DNA -m GTR+I+G -nt AUTO -bb 1000 
  
done


# fa=$(ls ${sdir}/*_aligned_fill_trimed.fa| sed s/_aligned_fill_trimed.fa//)
#  for a in $fa 
#         do
        
#         # echo $a
#         # echo $a\_aligned_fill_trimed.fa
#         # iqtree -s $a\_aligned_fill_trimed.fa -st DNA -m MFP -nt 10 -bb 1000 
#         iqtree -s $a\_aligned_fill_trimed.fa -st DNA -m GTR+G4+F -nt 10

#         done   




# for i in 

# 

# $a\_aligned.fa