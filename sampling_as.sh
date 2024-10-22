#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge/0307/sampling_subset/iq
    


elements=("MP" "NS" "PB2" "PB1" "PA" "NP")
# elements=("N9")
# elements=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9" "MP" "NS" "PB2" "PB1" "PA" "NP")
# elements=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9")

# tree=$(ls ${sdir}/*_iqtree.nexus | xargs -n  1 basename | sed s/_iqtree.nexus//)

for i in "${elements[@]}"
do
  # Rscript ${rdir}/patristic_distance_analysis.R \
  #   --segment $i \
  #   --tree_dir "/home/eric/Analysis/aiv/merge/0307/iq_tree/" \
  #   --prefix _aligned_iqtree.nexus \
  #   --out_fix "_iq" \
  #   --out_dir "/home/eric/Analysis/aiv/merge/0307/distance/NA_1009/"

  # Rscript ${rdir}/patristic_distance_analysis.R \
  #   --segment $i \
  #   --tree_dir "/home/eric/Analysis/aiv/merge/0307/ft_tree/" \
  #   --prefix _aligned_fasttree.tree \
  #   --out_fix _ft \
  #   --out_dir "/home/eric/Analysis/aiv/merge/0307/distance/NA_1009/"

  # Rscript ${rdir}/PB2_NS_sampling_result.R  \
  #   --segment $i \
  #   --iqtree /home/eric/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_$i\_aligned_iqtree.nexus \
  #   -iqg ~/Analysis/aiv/merge/0307/distance/NA_1009/$i\_iq_MPD_groups.RData \
  #   -ftg ~/Analysis/aiv/merge/0307/distance/NA_1009/$i\_ft_MPD_groups.RData \
  #   -op '' \
  #   -o ~/Analysis/aiv/merge/0307/distance/NA_1009/

  Rscript ${rdir}/PB2_NS_sampling_result.R  \
    --segment $i \
    --iqtree /home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/$i\_5p_sampling.nexus \
    -iqg ~/Analysis/aiv/merge/0307/distance/sampling_0924/$i\_iq_MPD_groups.RData \
    -ftg ~/Analysis/aiv/merge/0307/distance/sampling_0924/$i\_ft_MPD_groups.RData \
    -op _as_sampling \
    -o ~/Analysis/aiv/merge/0307/distance/sampling_0924/
done
