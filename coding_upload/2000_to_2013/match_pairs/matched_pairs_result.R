library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggVennDiagram)
library(ggplotify)
library(aplot)
library(viridis)



file_OwO <- "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/result/"

match_pair_result <- function(file_OwO){
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

results <- match_pair_result(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/result/")

# 
# clog_AS_all_grid <- results[[1]]
# clog_AS_group2 <- results[[2]]
# clog_AS_group3 <- results[[3]]
# clog_EU_all_grid <- results[[4]]
# clog_EU_group2 <- results[[5]]
# clog_EU_group3 <- results[[6]]




pos_sig_results <- results

for (i in seq_along(pos_sig_results)) {
  ORtable <- pos_sig_results[[i]]
  ORtable1 <- ORtable[!(is.na(ORtable$OR)), ]
  ORtable2 <- ORtable1[ORtable1$OR > 1 & ORtable1$P.value < 0.05, ]
  ORtable3 <- ORtable2[order(-ORtable2$OR),]
  pos_sig_results[[i]] <- ORtable3
}



pos_sig_results_AS <- pos_sig_results[1:3]
pos_sig_results_EU <- pos_sig_results[4:6]








# ---------------------------------------------------------------------------
# AS

AS_pos_sig_bird <- union(pos_sig_results_AS[[1]]$Bird_name, pos_sig_results_AS[[2]]$Bird_name)
AS_pos_sig_bird <- union(AS_pos_sig_bird, pos_sig_results_AS[[3]]$Bird_name)


AS_plot_df <- matrix(NA, nrow = length(AS_pos_sig_bird)*3, ncol = 3) %>% as.data.frame()
colnames(AS_plot_df) <- c('Bird', 'Control_group', 'OR')

AS_plot_df$Bird <- AS_pos_sig_bird
# table(AS_plot_df$Bird)

# i <- 1
for (i in 1:3) {
  ii <- i %>% as.character()
  len <- length(AS_pos_sig_bird)
  AS_plot_df[(1+len*(i-1)):(i*len), 'Control_group'] <- paste('group_', ii, sep = '')
}


# i <- 1
# k <- 1
for (i in seq_along(pos_sig_results_AS)) {
  result_OwO <- pos_sig_results_AS[[i]]
  df_name <- names(pos_sig_results_AS)[[i]]
  if (df_name %in% 'clog_AS_all_grid') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      AS_plot_df[(AS_plot_df$Bird %in% the_bird) & (AS_plot_df$`Control_group` %in% 'group_1'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
  if (df_name %in% 'clog_AS_group2') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      AS_plot_df[(AS_plot_df$Bird %in% the_bird) & (AS_plot_df$`Control_group` %in% 'group_2'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
  if (df_name %in% 'clog_AS_group3') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      AS_plot_df[(AS_plot_df$Bird %in% the_bird) & (AS_plot_df$`Control_group` %in% 'group_3'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
}

colnames(AS_plot_df)

AS_plot_df <- AS_plot_df[order(AS_plot_df$Bird),]




# -------------------------------------------------------------------------------------
# group bar plot, but no a good idea, so use heat map below
# 將 category 轉為數值並計算邊界
AS_plot_df <- AS_plot_df %>%
  mutate(category_numeric = as.numeric(as.factor(Bird)))

# 計算每個 group 的邊界 (取最大值 + 0.5)
group_breaks <- AS_plot_df %>%
  group_by(Bird) %>%
  summarise(break_position = max(category_numeric) + 0.5) %>%
  pull(break_position)


# group_breaks <- head(group_breaks, -1)  # 移除最後一組的邊界


ggplot(AS_plot_df, aes(fill= Control_group, y=Bird, x=OR)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_hline(yintercept = group_breaks, linetype = "solid", color = "black")





# ----------------------------------------------------------
# Heat map



AAA <- ggplot(AS_plot_df, aes(fill= OR, y=Bird, x=Control_group)) + 
  geom_tile() +
  labs(title = '', x = '', y = 'Bird') +
  scale_fill_gradient(low = "white", high = "red")+
  theme(
    panel.background = element_rect(fill = "white", color = NA), # 繪圖區背景
    plot.background = element_rect(fill = "white", color = NA),  # 整體背景
    panel.border = element_blank(),                              # 移除面板邊框
    axis.text.x = element_text(size = 18),       # x 軸標籤
    axis.text.y = element_text(size = 20),       # y 軸標籤
    axis.title.x = element_text(size = 14),      # x 軸標題
    axis.title.y = element_text(size = 14),      # y 軸標題
  ) +
  scale_x_discrete(labels = c("group_1" = "Group 1", "group_2" = "Group 2", "group_3" = "Group 3"))
  
AAA

# for (i in seq_along(pos_sig_results)) {
#   print(pos_sig_results[[i]])
# }


# ------------------------------------------------------------------------------
# combine with bird information

# bird information
asia_bird_information <- read_excel("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/bird/asia_bird_information.xlsx")
europe_bird_information <- read_excel("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/bird/europe_bird_information.xlsx")

asia_bird_information1 <- asia_bird_information[,c(2,3,4,6,12,13,14)]
europe_bird_information1 <- europe_bird_information[,c(2,3,4,6,14,15,16)]


asia_bird_information2 <- asia_bird_information1[!(asia_bird_information1$`order name` %in% "NA"),]

sciname <- asia_bird_information2$SCIENTIFIC_NAME
sciname1 <- sub(pattern = " ", replacement = "\\.", sciname)
asia_bird_information2$SCIENTIFIC_NAME <- sciname1
# unique(asia_bird_information2$`Family name`)

# 檢查大家 Family name 都可用 " (" 做分割
OwO <- asia_bird_information2$`Family name`
TF1 <- str_detect(string = asia_bird_information2$`Family name`, pattern = " \\(") 
OwO1 <- OwO[!(TF1)]

OwO2 <- str_split(string = OwO, pattern = " \\(")
for (i in seq_along(OwO2)) {
  asia_bird_information2[i,"Family name"] <- OwO2[[i]][[1]]
}

setdiff(AS_pos_sig_bird, asia_bird_information2$SCIENTIFIC_NAME)
# all birds name be correct

# birdname_fix <- read_excel("D://梁雋承/研究所/集大成HH_as_case_2014to2023/result_visualization/bird/鳥種名稱修正.xlsx")
# for (i in birdname_fix$incorrect) {
#   target_AS[target_AS %in% i] <- birdname_fix[birdname_fix$incorrect %in% i, "correct"]
#   clog_AS_result_out[clog_AS_result_out$X %in% i,"X"] <- birdname_fix[birdname_fix$incorrect %in% i, "correct"]
#   MH_AS_result_out[MH_AS_result_out$X %in% i,"X"] <- birdname_fix[birdname_fix$incorrect %in% i, "correct"]
# }
# target_AS <- unlist(target_AS)

HM_df <- matrix(NA, nrow = length(AS_pos_sig_bird), ncol = 1) %>% as.data.frame()
colnames(HM_df) <- 'Bird'
HM_df$Bird <- AS_pos_sig_bird


HM_df$B_Order <- "Order"
HM_df$B_Family <- "Family"
HM_df$B_Migratory <- "Migratory"
HM_df$B_AIV <- "AIV"

HM_df$Order <- NA
HM_df$Family <- NA
HM_df$Migratory <- NA
HM_df$AIV <- "Yes"


for (i in AS_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Order'] <- asia_bird_information2[asia_bird_information2$SCIENTIFIC_NAME %in% i, 'order name']
}

for (i in AS_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Family'] <- asia_bird_information2[asia_bird_information2$SCIENTIFIC_NAME %in% i, 'Family name']
}

for (i in AS_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Migratory'] <- asia_bird_information2[asia_bird_information2$SCIENTIFIC_NAME %in% i, 'Migratory']
}

HM_df[HM_df$Migratory %in% "T", 'Migratory'] <- 'Yes'
HM_df[HM_df$Migratory %in% "F", 'Migratory'] <- 'No'


colnames(asia_bird_information2)[5] <- "Aiv"
AS_infor <- asia_bird_information2[,c("SCIENTIFIC_NAME", "Aiv", "science", "FAO")]
HM_df1 <- merge(x = HM_df, y = AS_infor, by.x = "Bird", by.y = "SCIENTIFIC_NAME", all.x = T)
HM_df1[(HM_df1$Aiv %in% NA)&
                 (HM_df1$science %in% NA)&
                 (HM_df1$FAO %in% NA),"AIV"] <- "No"


# factor
# null_table_HM1$Bird <- factor(null_table_HM1$Bird, levels = null_table_mix1$Bird)

# O_and_F <- null_table_HM1[,c("Order", "Family")]
# O_and_F1 <- O_and_F[order(O_and_F$Order),]
# null_table_HM1$Family <- factor(null_table_HM1$Family, levels = unique(O_and_F1$Family))

HM_df1$Migratory <- factor(HM_df1$Migratory, levels = c("Yes", "No"))
HM_df1$AIV <- factor(HM_df1$AIV, levels = c("Yes", "No"))


table(HM_df1$Order)

HM11 <- ggplot(HM_df1, aes(B_Order, Bird, fill= Order)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Anseriformes" ="#01A01C",
                               "Charadriiformes" ="black",
                               "Passeriformes" = "#8900DD",
                               "Coraciiformes" = "#0CF5E0",
                               "Pelecaniformes" = "red"))


table(HM_df1$Family)

HM12 <- ggplot(HM_df1, aes(B_Family, Bird, fill= Family)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Anatidae" ="#01A01C",
                               "Laridae" = "black",
                               "Scolopacidae" = "#666666",
                               "Sturnidae" = "#244BCF",
                               "Coraciidae" = "#0CF5E0",
                               "Hirundinidae" = "#7331C6",
                               "Threskiornithidae" = "#F323C7",
                               "Ardeidae" = "red"
                               ))


HM13 <- ggplot(HM_df1, aes(B_AIV, Bird, fill= AIV)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Yes" = "#FF8075", "No" = "Black"))




HM14 <- ggplot(HM_df1, aes(B_Migratory, Bird, fill= Migratory)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Yes" = "#66FF82", "No" = "Black"))


AS_final_plot <- AAA %>% insert_right(HM11, width = 0.5) %>% insert_right(HM12, width = 0.5) %>% insert_right(HM13, width = 0.5) %>% insert_right(HM14, width = 0.5) %>% as.ggplot()
AS_final_plot + labs(title = "Asia match-pairs result") 





# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# EU

EU_pos_sig_bird <- union(pos_sig_results_EU[[1]]$Bird_name, pos_sig_results_EU[[2]]$Bird_name)
EU_pos_sig_bird <- union(EU_pos_sig_bird, pos_sig_results_EU[[3]]$Bird_name)


EU_plot_df <- matrix(NA, nrow = length(EU_pos_sig_bird)*3, ncol = 3) %>% as.data.frame()
colnames(EU_plot_df) <- c('Bird', 'Control_group', 'OR')

EU_plot_df$Bird <- EU_pos_sig_bird
# table(AS_plot_df$Bird)

# i <- 1
for (i in 1:3) {
  ii <- i %>% as.character()
  len <- length(EU_pos_sig_bird)
  EU_plot_df[(1+len*(i-1)):(i*len), 'Control_group'] <- paste('group_', ii, sep = '')
}


# i <- 1
# k <- 1
for (i in seq_along(pos_sig_results_EU)) {
  result_OwO <- pos_sig_results_EU[[i]]
  df_name <- names(pos_sig_results_EU)[[i]]
  if (df_name %in% 'clog_EU_all_grid') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      EU_plot_df[(EU_plot_df$Bird %in% the_bird) & (EU_plot_df$`Control_group` %in% 'group_1'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
  if (df_name %in% 'clog_EU_group2') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      EU_plot_df[(EU_plot_df$Bird %in% the_bird) & (EU_plot_df$`Control_group` %in% 'group_2'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
  if (df_name %in% 'clog_EU_group3') {
    for (k in 1:nrow(result_OwO)) {
      the_bird <- result_OwO[k, 'Bird_name']
      EU_plot_df[(EU_plot_df$Bird %in% the_bird) & (EU_plot_df$`Control_group` %in% 'group_3'), 'OR'] <- result_OwO[k, 'OR']
    }
  }
}

colnames(EU_plot_df)

EU_plot_df <- EU_plot_df[order(EU_plot_df$Bird),]





# ------------------------------------------------------------------------------
# combine with bird information

# bird information
# asia_bird_information <- read_excel("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/bird/asia_bird_information.xlsx")
europe_bird_information <- read_excel("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/match_pairs/bird/europe_bird_information.xlsx")

# asia_bird_information1 <- asia_bird_information[,c(2,3,4,6,12,13,14)]
europe_bird_information1 <- europe_bird_information[,c(2,3,4,6,14,15,16)]


europe_bird_information2 <- europe_bird_information1[!(europe_bird_information1$`order name` %in% "NA"),]

sciname <- europe_bird_information2$SCIENTIFIC_NAME
sciname1 <- sub(pattern = " ", replacement = "\\.", sciname)
europe_bird_information2$SCIENTIFIC_NAME <- sciname1
# unique(asia_bird_information2$`Family name`)

# 檢查大家 Family name 都可用 " (" 做分割
OwO <- europe_bird_information2$`Family name`
TF1 <- str_detect(string = europe_bird_information2$`Family name`, pattern = " \\(") 
OwO1 <- OwO[!(TF1)]

OwO2 <- str_split(string = OwO, pattern = " \\(")
for (i in seq_along(OwO2)) {
  europe_bird_information2[i,"Family name"] <- OwO2[[i]][[1]]
}

name_noncorrect <- setdiff(EU_pos_sig_bird, europe_bird_information2$SCIENTIFIC_NAME)

# correct bird name from europe_bird_information2
for (i in name_noncorrect) {
  target_row <- grep(pattern = i, x = europe_bird_information2$SCIENTIFIC_NAME)
  EU_plot_df[EU_plot_df$Bird %in% i, 'Bird'] <- europe_bird_information2[target_row, 'SCIENTIFIC_NAME']
}

# Heat map
AAA <- ggplot(EU_plot_df, aes(fill= OR, y=Bird, x=Control_group)) + 
  geom_tile() +
  labs(title = '', x = '', y = 'Bird') +
  scale_fill_gradient(low = "white", high = "red")+
  theme(
    panel.background = element_rect(fill = "white", color = NA), # 繪圖區背景
    plot.background = element_rect(fill = "white", color = NA),  # 整體背景
    panel.border = element_blank(),                              # 移除面板邊框
    axis.text.x = element_text(size = 18),       # x 軸標籤
    axis.text.y = element_text(size = 12),       # y 軸標籤
    axis.title.x = element_text(size = 14),      # x 軸標題
    axis.title.y = element_text(size = 14),      # y 軸標題
  ) +
  scale_x_discrete(labels = c("group_1" = "Group 1", "group_2" = "Group 2", "group_3" = "Group 3"))

AAA


# correct bird name from europe_bird_information2
for (i in name_noncorrect) {
  target_row <- grep(pattern = i, x = europe_bird_information2$SCIENTIFIC_NAME)
  EU_pos_sig_bird[EU_pos_sig_bird %in% i] <- europe_bird_information2[target_row, 'SCIENTIFIC_NAME']
}

EU_pos_sig_bird <- unlist(EU_pos_sig_bird)





HM_df <- matrix(NA, nrow = length(EU_pos_sig_bird), ncol = 1) %>% as.data.frame()
colnames(HM_df) <- 'Bird'
HM_df$Bird <- EU_pos_sig_bird


HM_df$B_Order <- "Order"
HM_df$B_Family <- "Family"
HM_df$B_Migratory <- "Migratory"
HM_df$B_AIV <- "AIV"

HM_df$Order <- NA
HM_df$Family <- NA
HM_df$Migratory <- NA
HM_df$AIV <- "Yes"


for (i in EU_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Order'] <- europe_bird_information2[europe_bird_information2$SCIENTIFIC_NAME %in% i, 'order name']
}

for (i in EU_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Family'] <- europe_bird_information2[europe_bird_information2$SCIENTIFIC_NAME %in% i, 'Family name']
}

for (i in EU_pos_sig_bird) {
  HM_df[HM_df$Bird %in% i, 'Migratory'] <- europe_bird_information2[europe_bird_information2$SCIENTIFIC_NAME %in% i, 'Migratory']
}

HM_df[HM_df$Migratory %in% "T", 'Migratory'] <- 'Yes'
HM_df[HM_df$Migratory %in% "F", 'Migratory'] <- 'No'


colnames(europe_bird_information2)[5] <- "Aiv"
EU_infor <- europe_bird_information2[,c("SCIENTIFIC_NAME", "Aiv", "science", "FAO")]
HM_df1 <- merge(x = HM_df, y = EU_infor, by.x = "Bird", by.y = "SCIENTIFIC_NAME", all.x = T)
HM_df1[(HM_df1$Aiv %in% NA)&
         (HM_df1$science %in% NA)&
         (HM_df1$FAO %in% NA),"AIV"] <- "No"


# factor
# null_table_HM1$Bird <- factor(null_table_HM1$Bird, levels = null_table_mix1$Bird)

# O_and_F <- null_table_HM1[,c("Order", "Family")]
# O_and_F1 <- O_and_F[order(O_and_F$Order),]
# null_table_HM1$Family <- factor(null_table_HM1$Family, levels = unique(O_and_F1$Family))

HM_df1$Migratory <- factor(HM_df1$Migratory, levels = c("Yes", "No"))
HM_df1$AIV <- factor(HM_df1$AIV, levels = c("Yes", "No"))


table(HM_df1$Order)

HM11 <- ggplot(HM_df1, aes(B_Order, Bird, fill= Order)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Anseriformes" ="#01A01C",
                               "Charadriiformes" ="black",
                               "Columbiformes" = "#fbff00",
                               "Coraciiformes" = "#0CF5E0",
                               "Galliformes" = "#0032ff",
                               "Gruiformes" = "#ffa200",
                               "Passeriformes" = "#8900DD",
                               "Pelecaniformes" = "red",
                               "Strigiformes" = "#895700"))


table(HM_df1$Family)
major_family <- c('Anatidae', 'Laridae', 'Scolopacidae')
HM_df1[!(HM_df1$Family %in% major_family), 'Family'] <- "others"
HM_df1$Family <- factor(HM_df1$Family, levels = c('Anatidae', 'Laridae', 'Scolopacidae', 'others'))



HM12 <- ggplot(HM_df1, aes(B_Family, Bird, fill= Family)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Anatidae" ="#01A01C",
                               "Laridae" = "black",
                               "Scolopacidae" = "#666666",
                               "others" = "#E6E6E6"
  ))


HM13 <- ggplot(HM_df1, aes(B_AIV, Bird, fill= AIV)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Yes" = "#FF8075", "No" = "Black"))




HM14 <- ggplot(HM_df1, aes(B_Migratory, Bird, fill= Migratory)) + 
  geom_tile()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 18),
        panel.background = element_rect(fill = "white"),
        legend.key.height = unit(1, "pt"),
        legend.key.width  = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))+
  scale_fill_manual(values = c("Yes" = "#66FF82", "No" = "Black"))


EU_final_plot <- AAA %>% insert_right(HM11, width = 0.5) %>% insert_right(HM12, width = 0.5) %>% insert_right(HM13, width = 0.5) %>% insert_right(HM14, width = 0.5) %>% as.ggplot()
EU_final_plot + labs(title = "Europe match-pairs result") 




























