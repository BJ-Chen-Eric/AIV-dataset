#! /bin/bash 

rdir='/home/eric/R/aiv'

sdir=/home/eric/Analysis/aiv/merge/0307
    
cd ${sdir}/

echo $sdir/

# Rscript ${rdir}/Command_meta_clean.R -m /home/eric/Analysis/aiv/all/meta/ \
#         -gis /home/eric/Analysis/aiv/all/all.fasta \
#         -ird /home/eric/Analysis/aiv/ird/IRD_Sequence.fa \
#         -bc 10 \
#         -p gisaid_IRD_merged_ \
#         -o ${sdir}/



fa=$(ls ${sdir}/*.fa| sed s/.fa//)

    for a in $fa 
        do
        
        mafft --auto --anysymbol --thread -8  $a\.fa > $a\_aligned.fa

        done   
    


# sdir='/home/eric/Analysis/aiv/all'

# cd ${sdir}/

# echo $sdir/

# Rscript ${rdir}/Command_meta_clean.R -m ${sdir}/meta/ \
#         -s ${sdir}/all.fasta \
#         -p all_clean_ \
#         -o ${sdir}/



# fa=$(ls ${sdir}/*.fa| sed s/.fa//)

# for a in $fa 
#     do

#     mafft --auto --anysymbol  $a\.fa > $a\_aligned.fa

#     done 


# Rscript ${rdir}/Command_seq_clean_fillgap.R -i ${sdir}/ \
#         -o ${sdir}/


# hs=$(ls ~/Analysis/aiv/Non_H5_data/)
# for h in $hs
#     do

#     sdir='/home/eric/Analysis/aiv/Non_H5_data'

#     sdir=${sdir}/$h

#     cd ${sdir}/

#     echo $sdir/

#     Rscript ${rdir}/Command_meta_clean.R -m ${sdir}/gisaid_epiflu_isolates.csv \
#             -s ${sdir}/gisaid_epiflu_sequence.fasta \
#             -p $h\_clean_ \
#             -o ${sdir}/



#     fa=$(ls ${sdir}/*.fa| sed s/.fa//)

#     for a in $fa 
#         do

#         mafft --auto --anysymbol  $a\.fa > $a\_aligned.fa

#         done 


#     Rscript ${rdir}/Command_seq_clean_fillgap.R -i ${sdir}/ \
#             -o ${sdir}/

#     done


