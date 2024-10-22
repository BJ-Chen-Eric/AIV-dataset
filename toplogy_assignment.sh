#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge
    


elements=('NS' "PB2" "MP" "PB1" "PA" "NP")
# elements=("NS")
# elements=()
# elements=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9")



# Loop through the elements
for i in "${elements[@]}"
do
  echo "Processing element: $i"
  
#   Rscript ${rdir}/topology_assignment.R \
#         -seg $i \
#         -iqt ~/Analysis/aiv/merge/0307/sampling_subset/0924/iq/$i\_5p_sampling.nexus \
#         -iqg ~/Analysis/aiv/merge/0307/distance/sampling_0924/$i\_iq_MPD_groups.RData \
#         -op '_as_sampling' \
#         -o ~/Analysis/aiv/merge/0307/distance/sampling_0924/
  Rscript ${rdir}/sampling_aa_Identicle_reflection.R \
        -seg $i \
        -iqt ~/Analysis/aiv/merge/0307/sampling_subset/0924/iq/$i\_5p_sampling.nexus \
        -iqg ~/Analysis/aiv/merge/0307/distance/sampling_0924/$i\_lineage_last.RData \
        -op '_as_sampling' \
        -o ~/Analysis/aiv/merge/0307/distance/sampling_0924/

done



