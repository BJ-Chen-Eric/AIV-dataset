#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge
    


# elements=("MP" "NS" "PB1" "PA" "NP")
# elements=("NS" "PB2")
elements=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9")



# Loop through the elements
for i in "${elements[@]}"
do
  echo "Processing element: $i"
  
  Rscript ${rdir}/agreement_grouping_0924.R \
        -seg $i \
        -iqt ~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_$i\_aligned_iqtree.nexus \
        -iqg ~/Analysis/aiv/merge/0307/distance/iq/$i\_iq_patristic_group.RData \
        -ftg ~/Analysis/aiv/merge/0307/distance/ft/$i\_ft_patristic_group.RData \
        -tree ~/Analysis/aiv/merge/0307/iq_tree/ \
        -tp _aligned_iqtree.nexus \
        -o ~/Analysis/aiv/merge/0307/distance/agree_pd/

  # Rscript ${rdir}/agreement_grouping_0924.R \
  #       -seg $i \
  #       -dis ~/Analysis/aiv/merge/0307/iq_raw/gisaid_IRD_merged_$i\_aligned_fill_trimed_remove_outlier.fa.mldist \
  #       -iqt ~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_$i\_aligned_iqtree.nexus \
  #       -iqg ~/Analysis/aiv/merge/0307/distance/iq/$i\_iq_patristic_group.RData \
  #       -ftg ~/Analysis/aiv/merge/0307/distance/ft/$i\_ft_patristic_group.RData \
  #       -tree ~/Analysis/aiv/merge/0307/iq_tree/ \
  #       -tp _aligned_iqtree.nexus \
  #       -o ~/Analysis/aiv/merge/0307/distance/agree_pd/
done



