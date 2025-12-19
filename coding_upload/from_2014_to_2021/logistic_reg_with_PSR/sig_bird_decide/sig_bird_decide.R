library(data.table)
library(readr)
library(dplyr)
library(pROC)
library(ggplot2)
library(vegan)

file_OwO <- "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/"

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
results <- match_pair_result(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/match_pairs/result/")
results_OwO <- results[c(2, 3, 5, 6)]



complete_grid_data <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/FAO_and_GISAID_record_grid_data/complete_grid_data_20250519.csv') %>% as.data.frame()



bird_YN_asia <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_YN_asia.csv') %>% as.data.frame()
bird_YN_europe <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_YN_europe.csv') %>% as.data.frame()


bird_num_asia <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_num_asia.csv')
bird_num_europe <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_num_europe.csv')




ID_table <- bird_num_asia[, 'Id'] %>% as.data.frame()
AS_tran <- apply(bird_num_asia[,-1], MARGIN = 2, FUN = function(X){ceiling(log((X+1)))}) %>% as.data.frame()
bird_num_asia1 <- cbind(ID_table, AS_tran)



ID_table <- bird_num_europe[, 'Id'] %>% as.data.frame()
EU_tran <- apply(bird_num_europe[,-1], MARGIN = 2, FUN = function(X){ceiling(log((X+1)))}) %>% as.data.frame()
bird_num_europe1 <- cbind(ID_table, EU_tran)




# Tue May 20 17:32:33 2025 ------------------------------

# use the integrate index: Dr_Wu_with_AIC_method, shannon index, richness, gini-simpson
# use different subset from conditional logistic regression:
# point estimation > 1, point estimation < 1, mix 
# go to see the plot

# before this, please run the code from 1-48


model_subset_list <- c()
# subset_list <- c('pos', 'neg', 'mix')
index_list <- c('DrWu', 'shannon', 'richness', 'simpson')

# i <- 1
# j <- 'pos'
# k <- 'DrWu'


# create all table
for (i in seq_along(results_OwO)) {
  target_df <- results_OwO[[i]]
  target_df_name <- names(results_OwO)[[i]]
  # pos subset
  target_df_pos <- target_df[target_df$OR > 1, ]
  target_df_pos1 <- target_df_pos[!(is.na(target_df_pos$OR)), ]
  target_df_pos2 <- target_df_pos1[order(target_df_pos1$P.value), ]
  plot_table <- matrix(data = NA, nrow = nrow(target_df_pos2), ncol = 7) %>% as.data.frame()
  colnames(plot_table) <- c('Bird_name', 'Cumulate_bird_num', 'AUC', 'AUC_UP', 'AUC_DW', 'p_value', 'model_AIC')
  plot_table[, 'Bird_name'] <- target_df_pos2[, 'Bird_name']
  plot_table[, 'p_value'] <- target_df_pos2[, 'P.value']
  plot_table[, 'Cumulate_bird_num'] <- 1:nrow(target_df_pos2)
  
  for (k in index_list) {
    tablename <- paste(target_df_name, '_pos_', k, sep = '')
    model_subset_list[[tablename]] <- plot_table
  }
  
  # neg subset
  target_df_neg <- target_df[target_df$OR < 1, ]
  target_df_neg1 <- target_df_neg[!(is.na(target_df_neg$OR)), ]
  target_df_neg2 <- target_df_neg1[order(target_df_neg1$P.value), ]
  plot_table <- matrix(data = NA, nrow = nrow(target_df_neg2), ncol = 7) %>% as.data.frame()
  colnames(plot_table) <- c('Bird_name', 'Cumulate_bird_num', 'AUC', 'AUC_UP', 'AUC_DW', 'p_value', 'model_AIC')
  plot_table[, 'Bird_name'] <- target_df_neg2[, 'Bird_name']
  plot_table[, 'p_value'] <- target_df_neg2[, 'P.value']
  plot_table[, 'Cumulate_bird_num'] <- 1:nrow(target_df_neg2)
  
  for (k in index_list) {
    tablename <- paste(target_df_name, '_neg_', k, sep = '')
    model_subset_list[[tablename]] <- plot_table
  }
  
  
  # mix 
  target_df_all <- target_df[!(is.na(target_df$OR)), ]
  target_df_all1 <- target_df_all[order(target_df_all$P.value), ]
  plot_table <- matrix(data = NA, nrow = nrow(target_df_all1), ncol = 7) %>% as.data.frame()
  colnames(plot_table) <- c('Bird_name', 'Cumulate_bird_num', 'AUC', 'AUC_UP', 'AUC_DW', 'p_value', 'model_AIC')
  plot_table[, 'Bird_name'] <- target_df_all1[, 'Bird_name']
  plot_table[, 'p_value'] <- target_df_all1[, 'P.value']
  plot_table[, 'Cumulate_bird_num'] <- 1:nrow(target_df_all1)
  
  for (k in index_list) {
    tablename <- paste(target_df_name, '_mix_', k, sep = '')
    model_subset_list[[tablename]] <- plot_table
  }
}




# calculation four indices


model_subset_index_list <- c()

# i <- 1


total_steps <- length(model_subset_list)
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)

all_start_time <- Sys.time()

for (i in seq_along(model_subset_list)) {
  
  start_time <- Sys.time()
  
  target_df <- model_subset_list[[i]]
  target_df_name <- names(model_subset_list)[[i]]
  
  # filter the continent and study group grid
  if (grepl(x = target_df_name, pattern = 'AS')) {
    continent_id <- complete_grid_data[complete_grid_data$Asia == 1, ]
    colnames(continent_id)[colnames(continent_id) %in% 'AS_HH'] <- 'case'
    bird_YN <- bird_YN_asia
    bird_num <- bird_num_asia1
  } 
  if (grepl(x = target_df_name, pattern = 'EU')) {
    continent_id <- complete_grid_data[complete_grid_data$Europe == 1, ]
    colnames(continent_id)[colnames(continent_id) %in% 'EU_HH'] <- 'case'
    bird_YN <- bird_YN_europe
    bird_num <- bird_num_europe1
  } 
  if (grepl(x = target_df_name, pattern = 'group2')) {
    continent_group_id <- continent_id[continent_id$Group_2 == 1, ]
  } 
  if (grepl(x = target_df_name, pattern = 'group3')) {
    continent_group_id <- continent_id[continent_id$Group_3 == 1, ]
  }
  
  continent_group_id1 <- continent_group_id[, c('Id', 'case')]
  
  
  if (grepl(x = target_df_name, pattern = 'DrWu')) {
    Dr_Wu_plot <- target_df
    for (j in 1:nrow(Dr_Wu_plot)) {
      include_bird_list <- Dr_Wu_plot[1:j, 'Bird_name']
      target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% include_bird_list)|(colnames(bird_YN)) %in% 'Id')]
      log_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      log_table1 <- log_table[, !(colnames(log_table) %in% 'Id')] %>% as.data.frame()
      model <- glm(formula = case ~ ., data = log_table1, family = 'binomial')
      pred_prob <- predict(model, type = "response") 
      for_roc <- continent_group_id1
      roc_obj <- roc(for_roc$case, pred_prob, levels = c(0, 1), direction = "<")
      ci_auc <- ci.auc(roc_obj)
      Dr_Wu_plot[j, 'AUC_UP'] <- ci_auc[3]
      Dr_Wu_plot[j, 'AUC_DW'] <- ci_auc[1]
      Dr_Wu_plot[j, 'AUC'] <- auc(roc_obj)
      Dr_Wu_plot[j, 'model_AIC'] <- AIC(model)
    }
    out_put_table <- Dr_Wu_plot
  }
  
  if (grepl(x = target_df_name, pattern = 'shannon')) {
    shannon_plot <- target_df
    for (k in 1:nrow(shannon_plot)) {
      include_bird_list <- shannon_plot[1:k, 'Bird_name']
      target_bird_df <- bird_num[, ((colnames(bird_num) %in% include_bird_list)|(colnames(bird_num)) %in% 'Id')]
      shannon_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      shannon_table[, 'shannon_index'] <- NA
      for (l in 1:nrow(shannon_table)) {
        shannon_table[l,"shannon_index"] <- diversity(shannon_table[l, which(!(colnames(shannon_table) %in% c('Id', 'case', 'shannon_index')))], index = "shannon")
      }
      roc_obj <- roc(shannon_table$case, shannon_table$shannon_index, levels = c(0, 1), direction = "<")
      ci_auc <- ci.auc(roc_obj)
      shannon_plot[k, 'AUC_UP'] <- ci_auc[3]
      shannon_plot[k, 'AUC_DW'] <- ci_auc[1]
      shannon_plot[k, 'AUC'] <- auc(roc_obj)
    }
    out_put_table <- shannon_plot
  }
  
  if (grepl(x = target_df_name, pattern = 'richness')) {
    richness_plot <- target_df
    for (k in 1:nrow(richness_plot)) {
      include_bird_list <- richness_plot[1:k, 'Bird_name']
      target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% include_bird_list)|(colnames(bird_YN)) %in% 'Id')]
      richness_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      richness_table[, 'richness'] <- NA
      for (l in 1:nrow(richness_table)) {
        # sum(1, 2, 3)
        richness_table[l,"richness"] <- sum(richness_table[l, which(!(colnames(richness_table) %in% c('Id', 'case', 'richness')))])
      }
      roc_obj <- roc(richness_table$case, richness_table$richness, levels = c(0, 1), direction = "<")
      ci_auc <- ci.auc(roc_obj)
      richness_plot[k, 'AUC_UP'] <- ci_auc[3]
      richness_plot[k, 'AUC_DW'] <- ci_auc[1]
      richness_plot[k, 'AUC'] <- auc(roc_obj)
    }
    out_put_table <- richness_plot
  }
  
  if (grepl(x = target_df_name, pattern = 'simpson')) {
    gini_simpson_plot <- target_df
    for (k in 1:nrow(gini_simpson_plot)) {
      include_bird_list <- gini_simpson_plot[1:k, 'Bird_name']
      target_bird_df <- bird_num[, ((colnames(bird_num) %in% include_bird_list)|(colnames(bird_num)) %in% 'Id')]
      gini_simpson_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      gini_simpson_table[, 'gini_simpson_index'] <- NA
      for (l in 1:nrow(gini_simpson_table)) {
        gini_simpson_table[l,"gini_simpson_index"] <- diversity(gini_simpson_table[l, which(!(colnames(gini_simpson_table) %in% c('Id', 'case', 'gini_simpson_index')))], index = "simpson")
      }
      roc_obj <- roc(gini_simpson_table$case, gini_simpson_table$gini_simpson_index, levels = c(0, 1), direction = "<")
      ci_auc <- ci.auc(roc_obj)
      gini_simpson_plot[k, 'AUC_UP'] <- ci_auc[3]
      gini_simpson_plot[k, 'AUC_DW'] <- ci_auc[1]
      gini_simpson_plot[k, 'AUC'] <- auc(roc_obj)
    }
    out_put_table <- gini_simpson_plot
  }
  
  model_subset_index_list[[target_df_name]] <- out_put_table
  
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  setTxtProgressBar(pb, i)
  cat("         Iteration", i, "took", elapsed_time, "seconds")
  
}
close(pb)

all_end_time <- Sys.time()

all_elapsed_time <- all_end_time - all_start_time









# write out csv
 for (k in seq_along(model_subset_index_list)) {
   df_OwO <- model_subset_index_list[[k]]
   df_name <- names(model_subset_index_list)[[k]]
   route_df <- paste('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/sig_bird_decide/AUC_df/3_subset_4_model_4_index/', df_name, '.csv', sep = '')
   write.csv(df_OwO, file = route_df, row.names = F)
 }




# draw plot

k <- 1

library(tidyr)
library(scales)

for (k in seq_along(model_subset_index_list)) {
  df_OwO <- model_subset_index_list[[k]]
  df_name <- names(model_subset_index_list)[[k]]
  
  if (grepl(x = df_name, pattern = 'DrWu')) {
    df_OwO1 <- df_OwO
    df_OwO1[df_OwO1$model_AIC > 1000, 'model_AIC'] <- NA
    # plot_table_OwO <- df_OwO1 %>%
    #   pivot_longer(cols = c(AUC, p_value, model_AIC), names_to = "Line", values_to = "Value") %>% as.data.frame()
    list_OwO1 <- df_OwO1[, 'model_AIC'][!(is.na(df_OwO1[, 'model_AIC']))]
    df_OwO1[is.na(df_OwO1$model_AIC), 'model_AIC'] <- (max(list_OwO1)+1)
    
    plot_table_OwO <- df_OwO1 %>%
      mutate(
        model_AIC_scaled = scales::rescale(model_AIC, to = c(0, 1))
      ) %>% 
      pivot_longer(cols = c(AUC, p_value, model_AIC_scaled), names_to = "Line", values_to = "Value") %>% 
      as.data.frame()
    
    for_point_plot <- df_OwO1[df_OwO1$model_AIC %in% max(df_OwO1$model_AIC), ]
    for_point_plot[, 'on_plot'] <- 1
    
    
    plot_need_to_save <- ggplot()+
      geom_line(data = plot_table_OwO, aes(x = Cumulate_bird_num, y = Value, color = Line), size = 0.5)+
      scale_y_continuous(
        sec.axis = sec_axis(
          transform = ~ . * (max(df_OwO1$model_AIC) - min(df_OwO1$model_AIC)) + min(df_OwO1$model_AIC),
          name = "Value of model AIC"
        )) +
      geom_point(data = for_point_plot, aes(x = Cumulate_bird_num, y = on_plot, shape = 'Point'), 
                 color = "#138b50", size = 3) +
      scale_color_manual(values = c("AUC" = "black", "p_value" = "red", 'model_AIC_scaled'= '#138b00'),
                         labels = c("AUC" = "AUC", "p_value" = "P-Value", 'model_AIC_scaled'='AIC'))+
      scale_shape_manual(values = c("Point" = 10), # for shape
                         labels = c("Point" = "AIC > 1000")) +
      
      geom_ribbon(data = df_OwO, aes(x = Cumulate_bird_num, ymin=AUC_DW, ymax=AUC_UP), fill = 'blue', alpha = 0.2)+
      
      labs(x = 'Include bird number', y = 'Value of AUC and P-value', title = df_name, shape = "Point")+
      theme(
        plot.title = element_text(size = 30, face = "bold"),  # Title text size
        axis.title.x = element_text(size = 25),  # X-axis title text size
        axis.title.y = element_text(size = 25),  # Y-axis title text size
        axis.text.x = element_text(size = 20),  # X-axis tick labels size
        axis.text.y = element_text(size = 20),  # Y-axis tick labels size
        legend.title = element_text(size = 22),  # Legend title size
        legend.text = element_text(size = 15),  # Legend text size
        panel.background = element_blank(),  # Remove background
        # plot.background = element_blank(),  # Remove plot background
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 1),  # Show x and y axes
        axis.ticks = element_line(color = "black", linewidth = 1)  # Show axis ticks
      )
  } else {
    plot_table_OwO <- df_OwO %>%
      pivot_longer(cols = c(AUC, p_value), names_to = "Line", values_to = "Value") %>% as.data.frame()
    plot_need_to_save <- ggplot()+
      geom_line(data = plot_table_OwO, aes(x = Cumulate_bird_num, y = Value, color = Line), size = 0.5)+
      scale_color_manual(values = c("AUC" = "black", "p_value" = "red"),
                         labels = c("AUC" = "AUC", "p_value" = "P-Value"))+
      geom_ribbon(data = df_OwO, aes(x = Cumulate_bird_num, ymin=AUC_DW, ymax=AUC_UP), fill = 'blue', alpha = 0.2)+
      labs(x = 'Include bird number', y = 'Value', title = df_name)+
      theme(
        plot.title = element_text(size = 30, face = "bold"),  # Title text size
        axis.title.x = element_text(size = 25),  # X-axis title text size
        axis.title.y = element_text(size = 25),  # Y-axis title text size
        axis.text.x = element_text(size = 20),  # X-axis tick labels size
        axis.text.y = element_text(size = 20),  # Y-axis tick labels size
        legend.title = element_text(size = 22),  # Legend title size
        legend.text = element_text(size = 15),  # Legend text size
        panel.background = element_blank(),  # Remove background
        # plot.background = element_blank(),  # Remove plot background
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 1),  # Show x and y axes
        axis.ticks = element_line(color = "black", linewidth = 1)  # Show axis ticks
      )
  }
  
  save_route <- paste('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/sig_bird_decide/AUC_plot/3_subset_4_model_4_index/', df_name, '.png', sep = '')
  ggsave(plot = plot_need_to_save, filename = save_route, width = 12, height = 9, dpi = 300)
}




















