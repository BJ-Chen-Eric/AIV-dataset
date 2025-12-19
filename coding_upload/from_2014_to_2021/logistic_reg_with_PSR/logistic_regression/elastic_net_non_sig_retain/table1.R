library(data.table)
library(readr)
library(dplyr)
library(pROC)
library(ggplot2)
library(vegan)




read_csv_as_list <- function(file_OwO){
  all_file_name <- list.files(file_OwO)
  OwO <- c()
  for (i in all_file_name) {
    file_name <- paste(file_OwO, i, sep = "")
    OwO[[i]] <- read.csv(file_name)
  }
  for (j in seq_along(OwO)) {
    OwO1 <- strsplit(names(OwO)[j], split = "\\.")
    OwO2 <- unlist(OwO1)
    names(OwO)[j] <- OwO2[[1]]
  }
  return(OwO)
}


results <- read_csv_as_list(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/mul_log_csv/")

results_OwO <- results[(grepl(pattern = 'DrWu', x = names(results)))&(grepl(pattern = 'mix', x = names(results)))]


df_colname <- c('Model', 'Variable', 'AOR', '95% CI', 'P-value', 'VIF')
df_final <- matrix(data = NA, nrow = 0, ncol = length(df_colname)) %>% as.data.frame()
colnames(df_final) <- df_colname



# test <- results_OwO[[1]]
# 
# i <- 1
for (i in seq_along(results_OwO)) {
  target_df <- results_OwO[[i]]
  target_name <- names(results_OwO)[[i]]
  colnames(target_df) <- c('Variable', 'AOR', '95% CI', 'P-value', 'VIF')
  target_df[, 'Model'] <- target_name
  df_final <- rbind(target_df, df_final)
}


df_final$Model %>% table()

df_final[, 'VIF'] <- df_final$VIF %>% round(digits = 4)


df_final[df_final$Model %in% 'clog_AS_group3_mix_DrWu', 'Model'] <- 'Asia group 2'
df_final[df_final$Model %in% 'clog_AS_group2_mix_DrWu', 'Model'] <- 'Asia group 1'
df_final[df_final$Model %in% 'clog_EU_group2_mix_DrWu', 'Model'] <- 'Europe group 1'
df_final[df_final$Model %in% 'clog_EU_group3_mix_DrWu', 'Model'] <- 'Europe group 2'



write.csv(df_final, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/table1.csv', row.names = F)








# Wed May 21 15:36:16 2025 ------------------------------
#  haven't run below yet







# rename overall table for paper sup. --------------------------------------
library(data.table)

four_index_compete_overall_table <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/Over_all_table_20250402.csv') %>% as.data.frame()
four_index_compete_PN_VAR <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/PN_VAR.csv') %>% as.data.frame()



for (i in colnames(four_index_compete_overall_table)) {
  
  if (class(four_index_compete_overall_table[, i]) %in% 'numeric') {
    four_index_compete_overall_table[, i] <- four_index_compete_overall_table[, i] %>% round(digits = 4)
  } else {
    next
  }
  
  
}








# four_index_compete_overall_table
four_index_compete_overall_table1 <- four_index_compete_overall_table
four_index_compete_overall_table1$Model_Index <- c('Asia group 1', 'Asia group 2', 'Europe group 1', 'Europe group 2')
colnames(four_index_compete_overall_table1)[colnames(four_index_compete_overall_table1) %in% 'Variable_num'] <- 'Variable number'
four_index_compete_overall_table1[nrow(four_index_compete_overall_table1)+1, ] <- 0
four_index_compete_overall_table1[2:5, ] <- four_index_compete_overall_table1[1:4, ] 
four_index_compete_overall_table1[1, ] <- ''
four_index_compete_overall_table1[1, 1] <- 'Model'
colnames(four_index_compete_overall_table1)[colnames(four_index_compete_overall_table1) %in% 'Model_Index'] <- ''
write.csv(four_index_compete_overall_table1,
          file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/for_paper_sup/over_all_table.csv',
          row.names = F)


# four_index_compete_PN_VAR
four_index_compete_PN_VAR1 <- four_index_compete_PN_VAR

four_index_compete_PN_VAR2 <- four_index_compete_PN_VAR1[c(3,4,1,2), ]
four_index_compete_PN_VAR2$Model <- c('Asia group 1', 'Asia group 2', 'Europe group 1', 'Europe group 2')

four_index_compete_PN_VAR2[, ncol(four_index_compete_PN_VAR2)+1] <- four_index_compete_PN_VAR2[, 15]
four_index_compete_PN_VAR3 <- four_index_compete_PN_VAR2[, -15]

colnames(four_index_compete_PN_VAR3)[19] <- 'Logit'

four_index_compete_PN_VAR3[nrow(four_index_compete_PN_VAR3)+1, ] <- 0
four_index_compete_PN_VAR3[2:5, ] <- four_index_compete_PN_VAR3[1:4, ] 
four_index_compete_PN_VAR3[1, ] <- ''
four_index_compete_PN_VAR3[1, 1] <- 'Model'
colnames(four_index_compete_PN_VAR3)[colnames(four_index_compete_PN_VAR3) %in% 'Model'] <- ''
write.csv(four_index_compete_PN_VAR3,
          file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/for_paper_sup/PN_VAR.csv',
          row.names = F)







# ----------------------------------------------------------------------------




single_PRS_overall_table <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/Over_all_table_20250402.csv') %>% as.data.frame()
single_PRS_PN_VAR <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/PN_VAR.csv') %>% as.data.frame()



i <- 2
for (i in colnames(single_PRS_overall_table)) {
  
  if (class(single_PRS_overall_table[, i]) %in% 'numeric') {
    single_PRS_overall_table[, i] <- single_PRS_overall_table[, i] %>% round(digits = 4)
  } else {
    next
  }
  
  
}




single_PRS_overall_table1 <- single_PRS_overall_table

i <- 1
for (i in 1:nrow(single_PRS_overall_table1)) {
  target_c <- single_PRS_overall_table1[i, 'Model_Index']
  
  if (grepl(x = target_c, pattern = 'AS')) {
    region <- 'Asia'
  }
  if (grepl(x = target_c, pattern = 'EU')) {
    region <- 'Europe'
  }
  
  
  
  if (grepl(x = target_c, pattern = 'group2')) {
    group <- 'group 1'
  }
  if (grepl(x = target_c, pattern = 'group3')) {
    group <- 'group 2'
  }
  
  
  if (grepl(x = target_c, pattern = 'DrWu')) {
    index <- 'logit'
  }
  if (grepl(x = target_c, pattern = 'richness')) {
    index <- 'richness'
  }
  if (grepl(x = target_c, pattern = 'shannon')) {
    index <- 'shannon'
  }
  if (grepl(x = target_c, pattern = 'simpson')) {
    index <- 'simpson'
  }
  
  single_PRS_overall_table1[i, 'Model_Index'] <- paste(region, group, index, sep = ' ')
  
  
  
}


single_PRS_overall_table1[is.na(single_PRS_overall_table1$Youden_Index), 'Youden_Index'] <- ''

single_PRS_overall_table1[nrow(single_PRS_overall_table1)+1, ] <- 0
single_PRS_overall_table1[2:nrow(single_PRS_overall_table1), ] <- single_PRS_overall_table1[1:(nrow(single_PRS_overall_table1)-1), ]
single_PRS_overall_table1[1, ] <- ''
colnames(single_PRS_overall_table1)[1] <- ''
single_PRS_overall_table1[1, 1] <- 'Model and PRS'

colnames(single_PRS_overall_table1) <- c('',
                                         'AUC',
                                         "AIC", 
                                         "Variable number", 
                                         'Include bird number', 
                                         "P value threshold",
                                         "Youden index",
                                         "Alpha", 
                                         "Lambda")

write.csv(single_PRS_overall_table1,
          file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/for_paper_sup/over_all_table.csv',
          row.names = F)
















# single_PRS_PN_VAR

single_PRS_PN_VAR1 <- single_PRS_PN_VAR

i <- 1
for (i in 1:nrow(single_PRS_PN_VAR1)) {
  target_c <- single_PRS_PN_VAR1[i, 'Model']
  
  if (grepl(x = target_c, pattern = 'AS')) {
    region <- 'Asia'
  }
  if (grepl(x = target_c, pattern = 'EU')) {
    region <- 'Europe'
  }
  
  
  
  if (grepl(x = target_c, pattern = 'group2')) {
    group <- 'group 1'
  }
  if (grepl(x = target_c, pattern = 'group3')) {
    group <- 'group 2'
  }
  
  
  if (grepl(x = target_c, pattern = 'DrWu')) {
    index <- 'logit'
  }
  if (grepl(x = target_c, pattern = 'richness')) {
    index <- 'richness'
  }
  if (grepl(x = target_c, pattern = 'shannon')) {
    index <- 'shannon'
  }
  if (grepl(x = target_c, pattern = 'simpson')) {
    index <- 'simpson'
  }
  
  single_PRS_PN_VAR1[i, 'Model'] <- paste(region, group, index, sep = ' ')
  
  
  
}



single_PRS_PN_VAR1[nrow(single_PRS_PN_VAR1)+1, ] <- 0
single_PRS_PN_VAR1[2:nrow(single_PRS_PN_VAR1), ] <- single_PRS_PN_VAR1[1:(nrow(single_PRS_PN_VAR1)-1), ]
single_PRS_PN_VAR1[1, ] <- ''
colnames(single_PRS_PN_VAR1)[1] <- ''
single_PRS_PN_VAR1[1, 1] <- 'Model and PRS'


write.csv(single_PRS_PN_VAR1,
          file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/for_paper_sup/PN_VAR.csv',
          row.names = F)



