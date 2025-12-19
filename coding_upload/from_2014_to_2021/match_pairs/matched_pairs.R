library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(survival)
library(metafor)

setwd("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/")

complete_grid_data <- read.csv("母資料/complete_grid_data_20250519.csv")
# grid_data_M1 <- complete_grid_data[complete_grid_data$FAO_or_GSGD_outbreak_2016_to_2021==1,]
# grid_data_M2 <- complete_grid_data
Asia <- complete_grid_data[complete_grid_data$Asia ==1,]
Europe <- complete_grid_data[complete_grid_data$Europe ==1,]





Asia_PS_env <- read.csv("母資料/Asia_PS_env.csv")
Europe_PS_env <- read.csv("母資料/Europe_PS_env.csv")


Asia_truth <- read.csv("母資料/Asia_truth.csv")
Europe_truth <- read.csv("母資料/Europe_truth.csv")



Europe_group2 <- Europe[Europe$Group_2 %in% 1, ]
Europe_group3 <- Europe[Europe$Group_3 %in% 1, ]

table(Europe_group2$EU_HH)
table(Europe_group3$EU_HH)


Asia_group2 <- Asia[Asia$Group_2 %in% 1, ]
Asia_group3 <- Asia[Asia$Group_3 %in% 1, ]
table(Asia_group2$AS_HH)
table(Asia_group3$AS_HH)

# test-----------------
analysis_df <-  Europe
location_PS <-  Europe_PS_env
location_truth <-  Europe_truth
case_col <- 'EU_HH'
birdnum <- 1

clog_with_risk_set_sampling <- function(analysis_df, location_PS, location_truth, case_col){
  birdlist <- colnames(location_PS)[-1]
  # out <- matrix(0,length(birdlist),7) %>% as.data.frame()
  # row.names(out) <- birdlist
  # colnames(out) <- c("Npa", "Nsp(in Npa)", "Nng", "Nsn(in Nng)", "Nnon", "ORisNA", "Max_cscn_diff")
  sub_ID <- analysis_df$Id
  case <- analysis_df[, c('Id', case_col)]
  PS <- location_PS[location_PS$Id %in% sub_ID,]
  truth <- location_truth[location_truth$Id %in% sub_ID,]
  ORtable <- matrix(0,length(birdlist),5) %>% as.data.frame()
  colnames(ORtable) <- c("Bird_name", "OR", "P-value", "95CI_lb", "95CI_ub")
  options(scipen=999)
  pb <- txtProgressBar(min = 0, max = length(birdlist), style = 3)
  for (birdnum in seq_along(birdlist)) {
    birdname <- birdlist[birdnum]
    PS1 <- PS[,c(1,which(colnames(PS) %in% birdname))]
    truth1 <- truth[,c(1,which(colnames(truth) %in% birdname))]
    AAA <- merge(x=case, y=PS1, by = "Id")
    BBB <- merge(x=AAA, y=truth1, by = "Id")
    colnames(BBB)[3] <- "PS"
    colnames(BBB)[4] <- "Truth"
    outbreak_1 <- BBB[BBB[, case_col] %in% 1,]
    outbreak_0 <- BBB[BBB[, case_col] %in% 0,]
    # Npa <- c()
    # Nsp <- c()
    # Nng <- c()
    # Nsn <- c()
    # Nnon <- c()
    # ORisNA <- c()
    diff_v <- c()
    PSdiffs <- c()
    for (j in 1:nrow(outbreak_1)) {
      # PSdiffs[[j]] <- min(abs(outbreak_1[j, 3] - outbreak_0[ ,3]))
      for (k in 1:nrow(outbreak_0)) {
        diffs <- outbreak_1[j, 3] - outbreak_0[k ,3] ##爆發禽場傾向分數-沒有爆發禽場傾向分數
        diff_v[k] <- diffs
      }
      PSdiffs[[j]] <- min(abs(diff_v))
    }
    PSdiffs <- unlist(PSdiffs)
    
    dataforclog <- matrix(0,0,4)
    colnames(dataforclog) <- c('Match_group',case_col,'PS','Truth')
    for (j in 1:nrow(outbreak_1)) {
      case_PS <- outbreak_1[j,3] #抓出 case 的 ps
      DL<- round((case_PS-max(PSdiffs)), digits = 12)
      UL<- round((case_PS+max(PSdiffs)), digits = 12)
      outbreak_0$PS <- round((outbreak_0$PS), digits = 12)
      control_risk_set <- outbreak_0[outbreak_0$PS >= DL & outbreak_0$PS <= UL,]
      ran_row <- sample(nrow(control_risk_set),1)
      dat <- rbind.data.frame(outbreak_1[j,],control_risk_set[ran_row,]) 
      colnames(dat)[1] <- 'Match_group'
      dat$Match_group <- j
      dataforclog <- rbind.data.frame(dataforclog, dat)
      }
    dataforclog[, case_col] <- ifelse(dataforclog[, case_col] %in% 1,TRUE,FALSE)
    formula_OwO <- paste(case_col, ' ~ Truth + strata(Match_group)', sep = '') %>% as.formula() 
    
    clogmodel <- clogit(formula_OwO , data = dataforclog)
    
    OR <- summary(clogmodel)$coefficients[,2]
    if (OR %in% NA) {
      ORtable[birdnum,"Bird_name"] <- birdname
      ORtable[birdnum,"OR"] <- NA
      ORtable[birdnum,"P-value"] <- NA
      ORtable[birdnum,"95CI_lb"] <- NA
      ORtable[birdnum,"95CI_ub"] <- NA
      next
    }
    
    ORtable[birdnum,"Bird_name"] <- birdname
    ORtable[birdnum,"OR"] <- summary(clogmodel)$coefficients[,2]
    ORtable[birdnum,"P-value"] <- summary(clogmodel)$coefficients[,5]
    ORtable[birdnum,"95CI_lb"] <- exp(summary(clogmodel)$coefficients[,1] - 1.96*summary(clogmodel)$coefficients[,3])
    ORtable[birdnum,"95CI_ub"] <- exp(summary(clogmodel)$coefficients[,1] + 1.96*summary(clogmodel)$coefficients[,3])
    
    setTxtProgressBar(pb, birdnum)
  }
  close(pb)
  return(ORtable)
}

clog_EU_all_grid <- clog_with_risk_set_sampling(analysis_df = Europe, 
                                                location_PS = Europe_PS_env, 
                                                location_truth = Europe_truth, 
                                                case_col = 'EU_HH')


clog_EU_group2 <- clog_with_risk_set_sampling(analysis_df = Europe_group2, 
                                                location_PS = Europe_PS_env, 
                                                location_truth = Europe_truth, 
                                                case_col = 'EU_HH')



clog_EU_group3 <- clog_with_risk_set_sampling(analysis_df = Europe_group3, 
                                              location_PS = Europe_PS_env, 
                                              location_truth = Europe_truth, 
                                              case_col = 'EU_HH')


clog_AS_all_grid <- clog_with_risk_set_sampling(analysis_df = Asia, 
                                                location_PS = Asia_PS_env, 
                                                location_truth = Asia_truth, 
                                                case_col = 'AS_HH')


clog_AS_group2 <- clog_with_risk_set_sampling(analysis_df = Asia_group2, 
                                                location_PS = Asia_PS_env, 
                                                location_truth = Asia_truth, 
                                                case_col = 'AS_HH')


clog_AS_group3 <- clog_with_risk_set_sampling(analysis_df = Asia_group3, 
                                              location_PS = Asia_PS_env, 
                                              location_truth = Asia_truth, 
                                              case_col = 'AS_HH')





write.csv(clog_EU_all_grid, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_EU_all_grid.csv', row.names = F)
write.csv(clog_EU_group2, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_EU_group2.csv', row.names = F)
write.csv(clog_EU_group3, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_EU_group3.csv', row.names = F)


write.csv(clog_AS_all_grid, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_AS_all_grid.csv', row.names = F)
write.csv(clog_AS_group2, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_AS_group2.csv', row.names = F)
write.csv(clog_AS_group3, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/clog_AS_group3.csv', row.names = F)


