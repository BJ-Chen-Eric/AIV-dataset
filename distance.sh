#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge/0307/sampling_subset/iq
    


# elements=("MP" "NS" "PB2" "PB1" "PA" "NP")  #
# elements=("PB2")
# elements=("H5" "N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9" "MP" "NS" "PB2" "PB1" "PA" "NP")
# elements=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9")

tree=$(ls ${sdir}/*_iqtree.nexus | xargs -n  1 basename | sed s/_iqtree.nexus//)


# Loop through the elements
for i in $tree
do
  echo "Processing element: $i"
  # echo "Performing iq distance"
  # Rscript ${rdir}/distance_shell.R \
  #       -seg $i \
  #       -tree /home/eric/Analysis/aiv/merge/0307/sampling_subset/iq/ \
  #       -p _iqtree.nexus \
  #       -o /home/eric/Analysis/aiv/merge/0307/distance/ns_iq_subset/ \
  #       -of _ns_iq_subset_patristic_dist_matrix

  echo "Performing iq partition"
  Rscript ${rdir}/partistic_partition.R \
        -seg $i \
        -pd ~/Analysis/aiv/merge/0307/distance/ns_iq_subset/$i\_ns_iq_subset_patristic_dist_matrix.RData \
        -o ~/Analysis/aiv/merge/0307/distance/ns_iq_subset/ \
        -of _ns_iq_subset_patristic_group_bympd_percentile 

#   echo "Performing ft distance"        
#   Rscript ${rdir}/distance_shell.R \
#         -seg $i \
#         -tree /home/eric/Analysis/aiv/merge/0307/sampling_subset/NS_sam/ \
#         -p _fasttree.nexus \
#         -o /home/eric/Analysis/aiv/merge/0307/distance/ns_ft_subset/ \
#         -of _ns_ft_subset_patristic_dist_matrix

  echo "Performing ft partition"
  Rscript ${rdir}/partistic_partition_quantile.R \
        -seg $i \
        -pd ~/Analysis/aiv/merge/0307/distance/ns_ft_subset/$i\_ns_ft_subset_patristic_dist_matrix.RData \
        -o ~/Analysis/aiv/merge/0307/distance/ns_ft_subset/ \
        -of _ns_ft_subset_patristic_group_bympd_percentile 

  echo "Performing agreement"
  Rscript ${rdir}/agreement_grouping_0924.R \
        -seg $i \
        -iqt ~/Analysis/aiv/merge/0307/sampling_subset/iq/$i\_iqtree.nexus \
        -iqg ~/Analysis/aiv/merge/0307/distance/ns_iq_subset/$i\_ns_iq_subset_patristic_group_bympd_percentile.RData \
        -ftg ~/Analysis/aiv/merge/0307/distance/ns_ft_subset/$i\_ns_ft_subset_patristic_group_bympd_percentile.RData \
        -op _ns_sampling_bympd_percentile \
        -o ~/Analysis/aiv/merge/0307/distance/ns_sampling/
done
