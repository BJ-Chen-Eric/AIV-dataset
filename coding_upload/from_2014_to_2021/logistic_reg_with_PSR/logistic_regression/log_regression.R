# setwd('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/')
library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(dendextend)
library(car) #vif()
library(formattable)
library(pROC)
library(vegan)
library(data.table)

complete_grid_data <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/FAO_and_GISAID_record_grid_data/complete_grid_data_20250519.csv') %>% as.data.frame()
bird_YN_asia <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_YN_asia.csv') %>% as.data.frame()
bird_YN_europe <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_YN_europe.csv') %>% as.data.frame()


bird_num_asia <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_num_asia.csv')
bird_num_europe <- fread('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/鳥種母資料/bird_num_europe.csv')




ID_table <- bird_num_asia[, 'Id'] %>% as.data.frame()
AS_tran <- apply(bird_num_asia[,-1], MARGIN = 2, FUN = function(X){ceiling(log((X+1)))}) %>% as.data.frame()
bird_num_asia1 <- cbind(ID_table, AS_tran) # num transform



ID_table <- bird_num_europe[, 'Id'] %>% as.data.frame()
EU_tran <- apply(bird_num_europe[,-1], MARGIN = 2, FUN = function(X){ceiling(log((X+1)))}) %>% as.data.frame()
bird_num_europe1 <- cbind(ID_table, EU_tran) # num transform




# other ENV factors
H5_HPAI_farm <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/整理好的/H5_HPAI_farm.csv")
landcovertype_threshold_10_percent <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/整理好的/landcovertype_threshold_10_percent.csv")
livestock <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/整理好的/livestock.csv")




file_OwO <- "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/sig_bird_decide/AUC_df/3_subset_4_model_4_index/"


# j <- 1
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
results <- match_pair_result(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/sig_bird_decide/AUC_df/3_subset_4_model_4_index/")





# --------------------------------------------------------------------------------------
# Use elastic net for model selection, don't get rid of unsig. variable from model
# PRS must be in the model

# One model for four index selection

# criteria on Dr.Wu method: 
# bird num < 50
# AIC min

# criteria on other method: 
# bird num < 50
# AUC max

# PRS must less then 50 birds, since too much bird will lead to other variable couldn't get into final model



# before this, run the code from 1-64








DrWu_mix <- intersect(grep(x = names(results), pattern = 'DrWu'), grep(x = names(results), pattern = 'mix'))
richness_pos <- intersect(grep(x = names(results), pattern = 'richness'), grep(x = names(results), pattern = 'pos'))
shannon_pos <- intersect(grep(x = names(results), pattern = 'shannon'), grep(x = names(results), pattern = 'pos'))
simpson_pos <- intersect(grep(x = names(results), pattern = 'simpson'), grep(x = names(results), pattern = 'pos'))

target_num <- c(DrWu_mix, richness_pos, shannon_pos, simpson_pos)

target_result <- results[target_num]

# names(target_result)


# data <- HHH
# y_col <- 'case'


simple_logistic_regression <- function(data, y_col){
  all_variable <- colnames(data)[-1]
  simple_result <- data.frame(matrix(0, length(all_variable), 4))
  colnames(simple_result) <- c("Variable", "OR", "95% CI", "P-value")
  simple_result[,"Variable"] <- all_variable
  
  for (i in all_variable) {
    y <- y_col
    x <- i
    simple_form <- paste(y, " ~ ", x, sep = "") %>% as.formula()
    simple_log_model <- glm(simple_form, data = data, family = 'binomial')
    simple_result[simple_result$Variable %in% i,"OR"] <- summary(simple_log_model)$coefficients[i,"Estimate"] %>% exp() %>% round(digits = 4)
    if (is.na(confint(simple_log_model)[i,'2.5 %'])) {
      CI_Low <- 0
      CI_Up <- confint(simple_log_model)[i,'97.5 %']  %>% exp() %>% round(digits = 4)
      simple_result[simple_result$Variable %in% i,"95% CI"] <- paste(CI_Low, '-', CI_Up, sep = '')
    } else {
      CI_Low <- confint(simple_log_model)[i,'2.5 %']  %>% exp() %>% round(digits = 4)
      CI_Up <- confint(simple_log_model)[i,'97.5 %']  %>% exp() %>% round(digits = 4)
      simple_result[simple_result$Variable %in% i,"95% CI"] <- paste(CI_Low, '-', CI_Up, sep = '')
    }
    simple_result[simple_result$Variable %in% i,"P-value"] <- summary(simple_log_model)$coefficients[i,"Pr(>|z|)"] %>% round(digits = 4)
  }
  
  for (i in varname$original) {
    simple_result[simple_result$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  
  simple_result[simple_result$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  
  library('gt')
  library('webshot2')
  simple_table <- simple_result %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Odds ratio"),
                 locations = cells_column_labels(columns = 'OR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI')) 
  output <- c()
  output[[1]] <- simple_result
  output[[2]] <- simple_table
  return(output)
}



# model 外存

# model <- model_1
# summary(model)

model_save_OR_OwO <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  VIF_model <- car::vif(model) %>% as.data.frame()
  colnames(VIF_model) <- "VIF"
  VIF_model$Variable <- row.names(VIF_model)
  OvO <- merge(x = coefficient, y = VIF_model, by = "Variable", all.x = TRUE)
  for (i in varname$original) {
    OvO[OvO$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  bird <- OvO[grep(pattern = '\\.', x = OvO$Variable),'Variable']
  OvO$Variable <- factor(OvO$Variable, levels = c(bird, OvO$Variable[!(OvO$Variable %in% c(bird, "(Intercept)"))],  '(Intercept)'))
  OvO <- OvO[order(OvO$Variable),]
  OvO$Variable <- OvO$Variable %>% as.character()
  OvO[OvO$Variable %in% '(Intercept)', 'Variable'] <- 'Intercept'
  
  OvO[OvO$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  colnames(OvO)[2] <- 'AOR'
  OvO <- OvO[!(OvO$Variable %in% 'Intercept'), ]
  
  library('gt')
  library('webshot2')
  OvO_table <- OvO %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Adjusted odds ratio"),
                 locations = cells_column_labels(columns = 'AOR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI')) %>% 
    tab_footnote(footnote = md("Variance inflation factor"),
                 locations = cells_column_labels(columns = 'VIF'))
  output <- c()
  output[[1]] <- OvO
  output[[2]] <- OvO_table
  return(output)
}



# model 外存 for only one var.
model_save_OR_OwO_for_onevar <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  
  for (i in varname$original) {
    coefficient[coefficient$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  bird <- coefficient[grep(pattern = '\\.', x = coefficient$Variable),'Variable']
  coefficient$Variable <- factor(coefficient$Variable, levels = c(bird, coefficient$Variable[!(coefficient$Variable %in% c(bird, "(Intercept)"))],  '(Intercept)'))
  coefficient <- coefficient[order(coefficient$Variable),]
  coefficient$Variable <- coefficient$Variable %>% as.character()
  coefficient[coefficient$Variable %in% '(Intercept)', 'Variable'] <- 'Intercept'
  
  coefficient[coefficient$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  colnames(coefficient)[2] <- 'OR'
  coefficient <- coefficient[!(coefficient$Variable %in% 'Intercept'), ]
  
  library('gt')
  library('webshot2')
  OvO_table <- coefficient %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Odds ratio"),
                 locations = cells_column_labels(columns = 'OR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI'))
  output <- c()
  output[[1]] <- coefficient
  output[[2]] <- OvO_table
  return(output)
}



# model 外存 become a table
# model <- forward_log
# summary(model)
# confint.default(model)[,'2.5 %']


model_save_OwO <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  VIF_model <- car::vif(model) %>% as.data.frame()
  colnames(VIF_model) <- "VIF"
  VIF_model$Variable <- row.names(VIF_model)
  OvO <- merge(x = coefficient, y = VIF_model, by = "Variable", all.x = TRUE)
  return(OvO)
}

# model 把不顯著者踢掉

# model <- model_AAA
# case_col <- 'case'

model_final <- function(model, case_col){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient_sig <- coefficient[coefficient$`Pr(>|z|)` < 0.05,]
  coefficient_sig1 <- coefficient_sig[-1,]
  variableAA <- row.names(coefficient_sig1)
  
  
  if ((length(variableAA) == 1) & ('PRS' %in% variableAA)) {
    formA <- paste(case_col, 'PRS', sep = " ~ ") %>% as.formula() # c
    model_simplify <- glm(formA, data = GGG, family = 'binomial')
    return(model_simplify)
  } else {
    variableAA <- variableAA[!(variableAA %in% 'PRS')] # c
    AAAA <- paste(variableAA, collapse = " + ")
    formA <- paste(case_col, AAAA, sep = " ~ PRS+") %>% as.formula() # c
    model_simplify <- glm(formA, data = GGG, family = 'binomial')
    return(model_simplify)
  }
  
}




# elastic net
# data <- HHH
# y_col <- 'case'

elastic_net_variable_selection_log <- function(data, y_col){
  library(glmnet)
  y <- data[, y_col]
  X <- data[, -which(colnames(data) %in% y_col)] %>% as.matrix()
  
  
  # Define Alpha Values to Test
  alpha_values <- seq(0, 1, length.out= 200)  # Test alpha from 0 to 1
  cv_errors <- numeric(length(alpha_values))  # Store cross-validation errors
  lambda_values <- numeric(length(alpha_values))
  
  
  # Loop over each alpha value
  for (i in seq_along(alpha_values)) {
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = alpha_values[i], lambda = seq(0, 1, length.out=200))
    cv_errors[i] <- min(cv_fit$cvm)  # Store minimum cross-validation error
    lambda_values[i] <- cv_fit$lambda[cv_fit$cvm %in% min(cv_fit$cvm)]
  }
  
  # Find the Best Alpha
  best_alpha <- alpha_values[which.min(cv_errors)]
  best_lambda <- lambda_values[which.min(cv_errors)]
  
  cv.fit <- cv.glmnet(X, y, family = "binomial", alpha = best_alpha, lambda = seq(0, 1, length.out=200))
  
  
  
  
  elastic_result <- 
    data.frame(log_lambda = log(cv.fit$lambda),mse = cv.fit$cvm, 
               upper = cv.fit$cvm + cv.fit$cvsd, lower = cv.fit$cvm - cv.fit$cvsd, 
               n_features = cv.fit$nzero ) 
  
  
  # visualize for best alpha with different lambda
  
  
  best_alpha_txt <- best_alpha %>% round(digits = 4) %>% as.character()
  best_lambda_txt <- best_lambda %>% round(digits = 4) %>% as.character()
  
  plot_name <- paste('Cross-Validation Error vs Log Lambda   alpha = ', best_alpha_txt, ' lambda = ', best_lambda_txt, sep = '')
  
  best_alpha_plot <- elastic_result %>% mutate(label=ifelse(seq_len(200) %% 10 == 0, yes=n_features, no="")) %>% 
    ggplot(aes(x = log_lambda, y = mse)) +
    geom_point(color = "black") +  # MSE points
    geom_line(color = "black") +   # Line connecting points
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "gray") +  # Error bars
    geom_vline(xintercept = log(cv.fit$lambda.min), color = "green", linetype = "dashed") +  # Lambda.min
    geom_vline(xintercept =  log(cv.fit$lambda.1se), color = "blue", linetype = "dashed") +  # Lambda.1se
    # geom_vline(xintercept = log(min_lambda), color = "red", linetype = "dashed") +  # Custom lambda
    geom_text(aes(label = label, y = max(mse) + 0.0001), color = "black", size = 5, vjust = 0) +  # Feature count
    labs(title = plot_name,
         x = "Log(Lambda)", y = "Mean Squared Error") +
    theme_bw()
  
  
  best.model <- glmnet(X, y, alpha = best_alpha, lambda = best_lambda, family = "binomial") # model selection
  coefficients <- coef(best.model) %>% as.matrix() %>% as.data.frame()
  
  coefficients$var <- row.names(coefficients)
  
  coefficients1 <- coefficients[!(coefficients$s0 %in% 0), ] %>% as.data.frame()
  coefficients2 <- coefficients1[-1, ]
  
  
  
  output <- c()
  output[['alpha']] <- best_alpha
  output[['lambda']] <- best_lambda
  output[['variable']] <- coefficients2$var
  output[['best_alpha_with_lambda_plot']] <- best_alpha_plot
  
  return(output)
}


# test <- elastic_net_variable_selection_log(data = HHH, y_col = 'case')






list_OwO <- target_result
# names(list_OwO)
# i <- 3

# for output route
sim_log_csv <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/sim_log_csv/'
sim_log_table <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/sim_log_table/'

mul_log_csv <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/mul_log_csv/'
mul_log_table <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/mul_log_table/'



final_output_df_colname <- c('Model_Index', 'AUC', 'AIC', 'Variable_num', 'Include_bird_num', 'p_value_thre', 'Youden_Index', 'Alpha', 'Lambda')
final_output_df <- matrix(data = NA, nrow = length(list_OwO), ncol = length(final_output_df_colname)) %>% as.data.frame()
colnames(final_output_df) <- final_output_df_colname


# i <- 11




total_steps <- length(list_OwO)
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)


for (i in seq_along(list_OwO)) {
  target_df <- list_OwO[[i]]
  target_df_name <- names(list_OwO)[[i]]
  
  # select the bird subgroup by method
  if (grepl(x = target_df_name, pattern = 'DrWu')) {
    target_df1 <- target_df[target_df$Cumulate_bird_num <=50, ]
    target_df1 <- target_df1[!(target_df1$AUC %in% 1),]
    include_bird_num <- which(target_df1$model_AIC %in% min(target_df1$model_AIC))
    include_pvalue <- target_df1[include_bird_num, 'p_value']
    target_df2 <- target_df1[1:include_bird_num, ]
  } else {
    target_df1 <- target_df[target_df$Cumulate_bird_num <=50, ]
    include_bird_num <- which(target_df1$AUC %in% max(target_df1$AUC))
    include_pvalue <- target_df1[include_bird_num, 'p_value']
    target_df2 <- target_df1[1:include_bird_num, ]
  }
  # target_df2
  
  
  
  
  
  
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
  # continent_group_id1
  
  
  
  
  # calculate the indecies
  
  youden_index <- NA
  if (grepl(x = target_df_name, pattern = 'DrWu')) {
    # target_df2
    target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% target_df2$Bird_name)|(colnames(bird_YN)) %in% 'Id')]
    log_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
    log_table1 <- log_table[, !(colnames(log_table) %in% 'Id')] %>% as.data.frame()
    model <- glm(formula = case ~ ., data = log_table1, family = 'binomial')
    pred_prob <- predict(model, type = "response") 
    for_roc <- continent_group_id1
    roc_obj <- roc(for_roc$case, pred_prob, levels = c(0, 1), direction = "<")
    youden_index <- roc_obj$thresholds[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)]
    renew_index <- ifelse(pred_prob >= youden_index, 1, 0)
    continent_group_id2 <- continent_group_id1
    continent_group_id2[, 'PRS'] <- renew_index
  } 
  if (grepl(x = target_df_name, pattern = 'richness')) {
    target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% target_df2$Bird_name)|(colnames(bird_YN)) %in% 'Id')]
    richness_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
    richness_table[, 'richness'] <- NA
    for (l in 1:nrow(richness_table)) {
      # sum(1, 2, 3)
      richness_table[l,"richness"] <- sum(richness_table[l, which(!(colnames(richness_table) %in% c('Id', 'case', 'richness')))])
    }
    continent_group_id2 <- richness_table[, colnames(richness_table) %in% c('Id', 'case', 'richness')]
    colnames(continent_group_id2)[3] <- 'PRS'
  } 
  if (grepl(x = target_df_name, pattern = 'shannon')) {
    target_bird_df <- bird_num[, ((colnames(bird_num) %in% target_df2$Bird_name)|(colnames(bird_num)) %in% 'Id')]
    shannon_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
    shannon_table[, 'shannon_index'] <- NA
    for (l in 1:nrow(shannon_table)) {
      shannon_table[l,"shannon_index"] <- diversity(shannon_table[l, which(!(colnames(shannon_table) %in% c('Id', 'case', 'shannon_index')))], index = "shannon")
    }
    continent_group_id2 <- shannon_table[, colnames(shannon_table) %in% c('Id', 'case', 'shannon_index')]
    colnames(continent_group_id2)[3] <- 'PRS'
  } 
  if (grepl(x = target_df_name, pattern = 'simpson')) {
    target_bird_df <- bird_num[, ((colnames(bird_num) %in% target_df2$Bird_name)|(colnames(bird_num)) %in% 'Id')]
    gini_simpson_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
    gini_simpson_table[, 'gini_simpson_index'] <- NA
    for (l in 1:nrow(gini_simpson_table)) {
      gini_simpson_table[l,"gini_simpson_index"] <- diversity(gini_simpson_table[l, which(!(colnames(gini_simpson_table) %in% c('Id', 'case', 'gini_simpson_index')))], index = "simpson")
    }
    continent_group_id2 <- gini_simpson_table[, colnames(gini_simpson_table) %in% c('Id', 'case', 'gini_simpson_index')]
    colnames(continent_group_id2)[3] <- 'PRS'
  } 
  # continent_group_id2
  
  # data combine for further simple and multiple log. reg.
  CCC <- merge(x = continent_group_id2, y = H5_HPAI_farm, by = "Id", all.x = TRUE)
  DDD <- merge(x = CCC, y = landcovertype_threshold_10_percent, by = "Id", all.x = TRUE)
  EEE <- merge(x = DDD, y = livestock[,c(1,5,6,7)], by = "Id", all.x = TRUE)
  GGG <- EEE[,-1]
  
  
  # remove all 0 variables
  # j <- 1
  all_var_0_TF <- c()
  for (j in 1:ncol(GGG)) {
    TF_all <- ifelse(GGG[, j] %in% 0, T, F)
    if (FALSE %in% TF_all) {
      all_var_0_TF[[j]] <- TRUE
    } else {
      all_var_0_TF[[j]] <- FALSE
    }
  }
  all_var_0_TF <- all_var_0_TF %>% unlist()
  
  HHH <- GGG[, all_var_0_TF]
  
  
  
  # 參數對照表
  varname <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_regression/data/參數對照.csv")
  
  # simple_logistic_regression
  simple_log_result <- simple_logistic_regression(data = HHH, y_col = 'case')
  
  route <- paste(sim_log_csv, target_df_name, '.csv', sep = '')
  write.csv(simple_log_result[[1]], route, row.names = F)
  
  route1 <- paste(sim_log_table, target_df_name, '.png', sep = '')
  gtsave(simple_log_result[[2]], route1)
  
  
  # multiple regression
  # null = glm(case ~ 1, data = HHH, family = 'binomial')
  # full = glm(case ~ ., data = HHH, family = 'binomial')
  # forward_log <- step(null, scope=list(lower=null, upper=full), direction="forward")
  # summary(forward_log)
  
  
  
  # elastic net for model selection
  ES_result <- elastic_net_variable_selection_log(data = HHH, y_col = 'case')
  # names(ES_result)
  # ES_result[[4]]
  
  if (!('PRS' %in% ES_result[['variable']])) {
    all_var <- c('PRS', ES_result[['variable']])
  } else {
    all_var <- ES_result[['variable']]
  }
  
  
  target_variable <- all_var %>% paste(collapse = '+') 
  form <- paste('case ~ ', target_variable) %>% as.formula()
  model_1 <- glm(form, data = HHH, family = 'binomial')
  # summary(model_1)
  
  # forward_log <- model_1
  
  
  
  
  
  
  # some model only include PSR
  if (length(model_1$coefficients) == 2) {
    
    model_save_result <- model_save_OR_OwO_for_onevar(model_1)
    
    route2 <- paste(mul_log_csv, target_df_name, '.csv', sep = '')
    write.csv(model_save_result[[1]], route2, row.names = F)
    
    route3 <- paste(mul_log_table, target_df_name, '.png', sep = '')
    gtsave(model_save_result[[2]], route3)
    
    # final auc
    pred_prob_final <- predict(model_1, type = "response") 
    roc_obj_final <- roc(HHH$case, pred_prob_final, levels = c(0, 1), direction = "<")
    
    final_output_df[i, 'Model_Index'] <- target_df_name
    final_output_df[i, 'AUC'] <- auc(roc_obj_final)
    final_output_df[i, 'AIC'] <- model_1$aic
    final_output_df[i, 'Variable_num'] <- 1
    # final_output_df[i, 'Exclude_num'] <- 0
    final_output_df[i, 'Include_bird_num'] <- include_bird_num
    final_output_df[i, 'p_value_thre'] <- include_pvalue
    final_output_df[i, 'Youden_Index'] <- youden_index
    final_output_df[i, 'Alpha'] <- ES_result[['alpha']]
    final_output_df[i, 'Lambda'] <- ES_result[['lambda']]
  } else {
    
    model_save_result <- model_save_OR_OwO(model = model_1)
    
    route2 <- paste(mul_log_csv, target_df_name, '.csv', sep = '')
    write.csv(model_save_result[[1]], route2, row.names = F)
    
    route3 <- paste(mul_log_table, target_df_name, '.png', sep = '')
    gtsave(model_save_result[[2]], route3)
    
    
    
    
    final_model_df <- model_save_result[[1]]
    
    
    # final auc
    pred_prob_final <- predict(model_1, type = "response") 
    roc_obj_final <- roc(HHH$case, pred_prob_final, levels = c(0, 1), direction = "<")
    
    
    
    
    final_output_df[i, 'Model_Index'] <- target_df_name
    final_output_df[i, 'AUC'] <- auc(roc_obj_final)
    final_output_df[i, 'AIC'] <- model_1$aic
    final_output_df[i, 'Variable_num'] <- nrow(final_model_df)
    # final_output_df[i, 'Exclude_num'] <- exclude_num
    final_output_df[i, 'Include_bird_num'] <- include_bird_num
    final_output_df[i, 'p_value_thre'] <- include_pvalue
    final_output_df[i, 'Youden_Index'] <- youden_index
    final_output_df[i, 'Alpha'] <- ES_result[['alpha']]
    final_output_df[i, 'Lambda'] <- ES_result[['lambda']]
    
    
    
    
  }
  setTxtProgressBar(pb, i)
}
close(pb)







write.csv(final_output_df, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/Over_all_table_20250521.csv', row.names = F)




final_output_df_simple <- final_output_df





# Coefficient P/N table in mul. reg. simple PRS--------------------------------------------
mul_log_results <- match_pair_result(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/mul_log_csv/")


OwO_AAAA <- c()

for (i in seq_along(mul_log_results)) {
  OwO_aa <- mul_log_results[[i]]
  OwO_aaa <- OwO_aa$Variable
  OwO_AAAA <- union(OwO_AAAA, OwO_aaa)
}



over_all_tt <- matrix(data = 'x', nrow = length(mul_log_results), ncol = length(OwO_AAAA)) %>% as.data.frame()
colnames(over_all_tt) <- OwO_AAAA
row.names(over_all_tt) <- names(mul_log_results)


# j <- 5
for (j in seq_along(mul_log_results)) {
  OwO_aa <- mul_log_results[[j]]
  OwO_name <- names(mul_log_results)[[j]]
  
  for (k in 1:nrow(OwO_aa)) {
    if (OwO_aa[k, 'AOR'] > 1) {
      
      if (OwO_aa[k, 'P.value'] < 0.05) {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'P/S'
      } else {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'P/US'
      }
      
    } else {
      
      if (OwO_aa[k, 'P.value'] < 0.05) {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'N/S'
      } else {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'N/US'
      }
      
    }
  }
  
  
  
  
} 





write.csv(over_all_tt, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/PN_VAR.csv')








































# four index compete in same model-------------------------------------------------------------







simple_logistic_regression <- function(data, y_col){
  all_variable <- colnames(data)[-1]
  simple_result <- data.frame(matrix(0, length(all_variable), 4))
  colnames(simple_result) <- c("Variable", "OR", "95% CI", "P-value")
  simple_result[,"Variable"] <- all_variable
  
  for (i in all_variable) {
    y <- y_col
    x <- i
    simple_form <- paste(y, " ~ ", x, sep = "") %>% as.formula()
    simple_log_model <- glm(simple_form, data = data, family = 'binomial')
    simple_result[simple_result$Variable %in% i,"OR"] <- summary(simple_log_model)$coefficients[i,"Estimate"] %>% exp() %>% round(digits = 4)
    if (is.na(confint(simple_log_model)[i,'2.5 %'])) {
      CI_Low <- 0
      CI_Up <- confint(simple_log_model)[i,'97.5 %']  %>% exp() %>% round(digits = 4)
      simple_result[simple_result$Variable %in% i,"95% CI"] <- paste(CI_Low, '-', CI_Up, sep = '')
    } else {
      CI_Low <- confint(simple_log_model)[i,'2.5 %']  %>% exp() %>% round(digits = 4)
      CI_Up <- confint(simple_log_model)[i,'97.5 %']  %>% exp() %>% round(digits = 4)
      simple_result[simple_result$Variable %in% i,"95% CI"] <- paste(CI_Low, '-', CI_Up, sep = '')
    }
    simple_result[simple_result$Variable %in% i,"P-value"] <- summary(simple_log_model)$coefficients[i,"Pr(>|z|)"] %>% round(digits = 4)
  }
  
  for (i in varname$original) {
    simple_result[simple_result$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  
  simple_result[simple_result$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  
  library('gt')
  library('webshot2')
  simple_table <- simple_result %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Odds ratio"),
                 locations = cells_column_labels(columns = 'OR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI')) 
  output <- c()
  output[[1]] <- simple_result
  output[[2]] <- simple_table
  return(output)
}



# model 外存
model_save_OR_OwO <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  VIF_model <- car::vif(model) %>% as.data.frame()
  colnames(VIF_model) <- "VIF"
  VIF_model$Variable <- row.names(VIF_model)
  OvO <- merge(x = coefficient, y = VIF_model, by = "Variable", all.x = TRUE)
  for (i in varname$original) {
    OvO[OvO$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  bird <- OvO[grep(pattern = '\\.', x = OvO$Variable),'Variable']
  OvO$Variable <- factor(OvO$Variable, levels = c(bird, OvO$Variable[!(OvO$Variable %in% c(bird, "(Intercept)"))],  '(Intercept)'))
  OvO <- OvO[order(OvO$Variable),]
  OvO$Variable <- OvO$Variable %>% as.character()
  OvO[OvO$Variable %in% '(Intercept)', 'Variable'] <- 'Intercept'
  
  OvO[OvO$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  colnames(OvO)[2] <- 'AOR'
  OvO <- OvO[!(OvO$Variable %in% 'Intercept'), ]
  
  library('gt')
  library('webshot2')
  OvO_table <- OvO %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Adjusted odds ratio"),
                 locations = cells_column_labels(columns = 'AOR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI')) %>% 
    tab_footnote(footnote = md("Variance inflation factor"),
                 locations = cells_column_labels(columns = 'VIF'))
  output <- c()
  output[[1]] <- OvO
  output[[2]] <- OvO_table
  return(output)
}



# model 外存 for only one var.
model_save_OR_OwO_for_onevar <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  
  for (i in varname$original) {
    coefficient[coefficient$Variable %in% i,"Variable"] <- varname[varname$original %in% i,"full_name"]
  }
  bird <- coefficient[grep(pattern = '\\.', x = coefficient$Variable),'Variable']
  coefficient$Variable <- factor(coefficient$Variable, levels = c(bird, coefficient$Variable[!(coefficient$Variable %in% c(bird, "(Intercept)"))],  '(Intercept)'))
  coefficient <- coefficient[order(coefficient$Variable),]
  coefficient$Variable <- coefficient$Variable %>% as.character()
  coefficient[coefficient$Variable %in% '(Intercept)', 'Variable'] <- 'Intercept'
  
  coefficient[coefficient$`P-value` %in% 0, 'P-value'] <- '<0.0001'
  colnames(coefficient)[2] <- 'OR'
  coefficient <- coefficient[!(coefficient$Variable %in% 'Intercept'), ]
  
  library('gt')
  library('webshot2')
  OvO_table <- coefficient %>%
    gt() %>%
    cols_align(
      align = "center", # Set alignment to "center"
      columns = everything() # Apply to all columns
    ) %>% 
    tab_options(
      table.border.top.width = px(3),           # Set top border width
      table.border.top.color = "black",         # Change top border color
      table.border.bottom.width = px(3),        # Set bottom border width
      table.border.bottom.color = "black",      # Change bottom border color
      column_labels.border.bottom.width = px(2), # Set border width under column labels
      column_labels.border.bottom.color = "black", # Change border color under column labels
      table_body.hlines.width = 0,               # Remove horizontal lines in table body
      table_body.border.bottom.color = 'black'
    ) %>% 
    tab_footnote(footnote = md("Odds ratio"),
                 locations = cells_column_labels(columns = 'OR')) %>% 
    tab_footnote(footnote = md("95% confidence interval"),
                 locations = cells_column_labels(columns = '95% CI'))
  output <- c()
  output[[1]] <- coefficient
  output[[2]] <- OvO_table
  return(output)
}



# model 外存 become a table
model_save_OwO <- function(model){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient$Variable <- row.names(coefficient)
  coefficient$OR <- coefficient$Estimate %>% exp() %>% round(digits = 4)
  CI_Low <- confint.default(model)[,'2.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  CI_Up <- confint.default(model)[,'97.5 %']  %>% exp() %>% round(digits = 4) %>% format(scientific = FALSE)
  coefficient$'95% CI' <- paste(CI_Low, '-', CI_Up, sep = '')
  coefficient$'P-value' <- coefficient$`Pr(>|z|)` %>% round(digits = 4)
  coefficient <- coefficient[,5:8]
  VIF_model <- car::vif(model) %>% as.data.frame()
  colnames(VIF_model) <- "VIF"
  VIF_model$Variable <- row.names(VIF_model)
  OvO <- merge(x = coefficient, y = VIF_model, by = "Variable", all.x = TRUE)
  return(OvO)
}

# model 把不顯著者踢掉
model_final <- function(model, case_col){
  coefficient <- summary(model)$coefficient %>% as.data.frame()
  coefficient_sig <- coefficient[coefficient$`Pr(>|z|)` < 0.05,]
  coefficient_sig1 <- coefficient_sig[-1,]
  variableAA <- row.names(coefficient_sig1)
  AAAA <- paste(variableAA, collapse = " + ")
  formA <- paste(case_col, AAAA, sep = " ~ ") %>% as.formula()
  model_simplify <- glm(formA, data = GGG, family = 'binomial')
  return(model_simplify)
}





# elastic net
# data <- HHH
# y_col <- 'case'

elastic_net_variable_selection_log <- function(data, y_col){
  library(glmnet)
  y <- data[, y_col]
  X <- data[, -which(colnames(data) %in% y_col)] %>% as.matrix()
  
  
  # Define Alpha Values to Test
  alpha_values <- seq(0, 1, length.out= 200)  # Test alpha from 0 to 1
  cv_errors <- numeric(length(alpha_values))  # Store cross-validation errors
  lambda_values <- numeric(length(alpha_values))
  
  
  # Loop over each alpha value
  for (i in seq_along(alpha_values)) {
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = alpha_values[i], lambda = seq(0, 1, length.out=200))
    cv_errors[i] <- min(cv_fit$cvm)  # Store minimum cross-validation error
    lambda_values[i] <- cv_fit$lambda[cv_fit$cvm %in% min(cv_fit$cvm)]
  }
  
  # Find the Best Alpha
  best_alpha <- alpha_values[which.min(cv_errors)]
  best_lambda <- lambda_values[which.min(cv_errors)]
  
  cv.fit <- cv.glmnet(X, y, family = "binomial", alpha = best_alpha, lambda = seq(0, 1, length.out=200))
  
  
  
  
  elastic_result <- 
    data.frame(log_lambda = log(cv.fit$lambda),mse = cv.fit$cvm, 
               upper = cv.fit$cvm + cv.fit$cvsd, lower = cv.fit$cvm - cv.fit$cvsd, 
               n_features = cv.fit$nzero ) 
  
  
  # visualize for best alpha with different lambda
  
  
  best_alpha_txt <- best_alpha %>% round(digits = 4) %>% as.character()
  best_lambda_txt <- best_lambda %>% round(digits = 4) %>% as.character()
  
  plot_name <- paste('Cross-Validation Error vs Log Lambda   alpha = ', best_alpha_txt, ' lambda = ', best_lambda_txt, sep = '')
  
  best_alpha_plot <- elastic_result %>% mutate(label=ifelse(seq_len(200) %% 10 == 0, yes=n_features, no="")) %>% 
    ggplot(aes(x = log_lambda, y = mse)) +
    geom_point(color = "black") +  # MSE points
    geom_line(color = "black") +   # Line connecting points
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "gray") +  # Error bars
    geom_vline(xintercept = log(cv.fit$lambda.min), color = "green", linetype = "dashed") +  # Lambda.min
    geom_vline(xintercept =  log(cv.fit$lambda.1se), color = "blue", linetype = "dashed") +  # Lambda.1se
    # geom_vline(xintercept = log(min_lambda), color = "red", linetype = "dashed") +  # Custom lambda
    geom_text(aes(label = label, y = max(mse) + 0.0001), color = "black", size = 5, vjust = 0) +  # Feature count
    labs(title = plot_name,
         x = "Log(Lambda)", y = "Mean Squared Error") +
    theme_bw()
  
  
  best.model <- glmnet(X, y, alpha = best_alpha, lambda = best_lambda, family = "binomial") # model selection
  coefficients <- coef(best.model) %>% as.matrix() %>% as.data.frame()
  
  coefficients$var <- row.names(coefficients)
  
  coefficients1 <- coefficients[!(coefficients$s0 %in% 0), ] %>% as.data.frame()
  coefficients2 <- coefficients1[-1, ]
  
  
  
  output <- c()
  output[['alpha']] <- best_alpha
  output[['lambda']] <- best_lambda
  output[['variable']] <- coefficients2$var
  output[['best_alpha_with_lambda_plot']] <- best_alpha_plot
  
  return(output)
}








list_OwO <- target_result
# names(list_OwO)
# i <- 3

# for output route
sim_log_csv <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/sim_log_csv/'
sim_log_table <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/sim_log_table/'

mul_log_csv <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/mul_log_csv/'
mul_log_table <- '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/mul_log_table/'



model_list <- c('log_AS_group2', 'log_AS_group3', 'clog_EU_group2', 'clog_EU_group3')



final_output_df_colname <- c('Model_Index', 'AUC', 'AIC', 'Variable_num', 'Alpha', 'Lambda')
final_output_df <- matrix(data = NA, nrow = length(model_list), ncol = length(final_output_df_colname)) %>% as.data.frame()
colnames(final_output_df) <- final_output_df_colname


Over_all_table_20250521 <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/PRS_separate/Over_all_table/Over_all_table_20250521.csv')

# j <- 3
# i <- 3

for (j in seq_along(model_list)) {
  model_name <- model_list[[j]]
  refer_table <- Over_all_table_20250521[grep(pattern = model_name, x = Over_all_table_20250521$Model_Index), ]
  
  
  list_OwO1 <- list_OwO[grep(pattern = model_name, x = names(list_OwO))]
  # names(list_OwO1)
  
  
  for (i in seq_along(list_OwO1)) {
    target_df <- list_OwO1[[i]]
    target_df_name <- names(list_OwO1)[[i]]
    
    # select the bird subgroup by method
    if (grepl(x = target_df_name, pattern = 'DrWu')) {
      target_df1 <- target_df[target_df$Cumulate_bird_num <=100, ]
      target_df1 <- target_df1[!(target_df1$AUC %in% 1),]
      include_bird_num <- which(target_df1$model_AIC %in% min(target_df1$model_AIC))
      include_pvalue <- target_df1[include_bird_num, 'p_value']
      target_df2 <- target_df1[1:include_bird_num, ]
    } else {
      target_df1 <- target_df[target_df$Cumulate_bird_num <=100, ]
      include_bird_num <- which(target_df1$AUC %in% max(target_df1$AUC))
      include_pvalue <- target_df1[include_bird_num, 'p_value']
      target_df2 <- target_df1[1:include_bird_num, ]
    }
    # target_df2
    
    
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
    # continent_group_id1
    
    
    
    
    # calculate the indecies
    
    youden_index <- NA
    if (grepl(x = target_df_name, pattern = 'DrWu')) {
      # target_df2
      target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% target_df2$Bird_name)|(colnames(bird_YN)) %in% 'Id')]
      log_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      log_table1 <- log_table[, !(colnames(log_table) %in% 'Id')] %>% as.data.frame()
      model <- glm(formula = case ~ ., data = log_table1, family = 'binomial')
      pred_prob <- predict(model, type = "response") 
      for_roc <- continent_group_id1
      roc_obj <- roc(for_roc$case, pred_prob, levels = c(0, 1), direction = "<")
      youden_index <- roc_obj$thresholds[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)]
      renew_index <- ifelse(pred_prob >= youden_index, 1, 0)
      continent_group_id2 <- continent_group_id1
      continent_group_id2[, 'Logit'] <- renew_index
      Logit_T <- continent_group_id2
    } 
    if (grepl(x = target_df_name, pattern = 'richness')) {
      target_bird_df <- bird_YN[, ((colnames(bird_YN) %in% target_df2$Bird_name)|(colnames(bird_YN)) %in% 'Id')]
      richness_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      richness_table[, 'richness'] <- NA
      for (l in 1:nrow(richness_table)) {
        # sum(1, 2, 3)
        richness_table[l,"richness"] <- sum(richness_table[l, which(!(colnames(richness_table) %in% c('Id', 'case', 'richness')))])
      }
      continent_group_id2 <- richness_table[, colnames(richness_table) %in% c('Id', 'case', 'richness')]
      colnames(continent_group_id2)[3] <- 'Richness'
      Richness_T <- continent_group_id2
    } 
    if (grepl(x = target_df_name, pattern = 'shannon')) {
      target_bird_df <- bird_num[, ((colnames(bird_num) %in% target_df2$Bird_name)|(colnames(bird_num)) %in% 'Id')]
      shannon_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      shannon_table[, 'shannon_index'] <- NA
      for (l in 1:nrow(shannon_table)) {
        shannon_table[l,"shannon_index"] <- diversity(shannon_table[l, which(!(colnames(shannon_table) %in% c('Id', 'case', 'shannon_index')))], index = "shannon")
      }
      continent_group_id2 <- shannon_table[, colnames(shannon_table) %in% c('Id', 'case', 'shannon_index')]
      colnames(continent_group_id2)[3] <- 'Shannon'
      Shannon_T <- continent_group_id2
    } 
    if (grepl(x = target_df_name, pattern = 'simpson')) {
      target_bird_df <- bird_num[, ((colnames(bird_num) %in% target_df2$Bird_name)|(colnames(bird_num)) %in% 'Id')]
      gini_simpson_table <- merge(x = continent_group_id1, y = target_bird_df, all.x = T, by = 'Id')
      gini_simpson_table[, 'gini_simpson_index'] <- NA
      for (l in 1:nrow(gini_simpson_table)) {
        gini_simpson_table[l,"gini_simpson_index"] <- diversity(gini_simpson_table[l, which(!(colnames(gini_simpson_table) %in% c('Id', 'case', 'gini_simpson_index')))], index = "simpson")
      }
      continent_group_id2 <- gini_simpson_table[, colnames(gini_simpson_table) %in% c('Id', 'case', 'gini_simpson_index')]
      colnames(continent_group_id2)[3] <- 'Simpson'
      Simpson_T <- continent_group_id2
    } 
    # continent_group_id2
  }
  
  
  AA <- merge(x = Logit_T, y = Richness_T[, c('Id', 'Richness')], by = 'Id')
  BB <- merge(x = AA, y = Shannon_T[, c('Id', 'Shannon')], by = 'Id')
  CC <- merge(x = BB, y = Simpson_T[, c('Id', 'Simpson')], by = 'Id')
  
  
  
  CCC <- merge(x = CC, y = H5_HPAI_farm, by = "Id", all.x = TRUE)
  DDD <- merge(x = CCC, y = landcovertype_threshold_10_percent, by = "Id", all.x = TRUE)
  EEE <- merge(x = DDD, y = livestock[,c(1,5,6,7)], by = "Id", all.x = TRUE)
  GGG <- EEE[,-1]
  
  # remove all 0 variables
  # k <- 1
  all_var_0_TF <- c()
  for (k in 1:ncol(GGG)) {
    TF_all <- ifelse(GGG[, k] %in% 0, T, F)
    if (FALSE %in% TF_all) {
      all_var_0_TF[[k]] <- TRUE
    } else {
      all_var_0_TF[[k]] <- FALSE
    }
  }
  all_var_0_TF <- all_var_0_TF %>% unlist()
  
  HHH <- GGG[, all_var_0_TF]
  
  # 參數對照表
  varname <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_regression/data/參數對照.csv")
  
  # simple_logistic_regression
  simple_log_result <- simple_logistic_regression(data = HHH, y_col = 'case')
  
  route <- paste(sim_log_csv, model_name, '.csv', sep = '')
  write.csv(simple_log_result[[1]], route, row.names = F)
  
  route1 <- paste(sim_log_table, model_name, '.png', sep = '')
  gtsave(simple_log_result[[2]], route1)
  
  
  # multiple regression
  # null = glm(case ~ 1, data = HHH, family = 'binomial')
  # full = glm(case ~ ., data = HHH, family = 'binomial')
  # forward_log <- step(null, scope=list(lower=null, upper=full), direction="forward")
  # summary(forward_log)
  
  
  
  # elastic net for model selection
  ES_result <- elastic_net_variable_selection_log(data = HHH, y_col = 'case')
  # names(ES_result)
  
  
  
  target_variable <- ES_result[['variable']] %>% paste(collapse = '+') 
  form <- paste('case ~ ', target_variable) %>% as.formula()
  model_1 <- glm(form, data = HHH, family = 'binomial')
  # summary(model_1)
  
  forward_log <- model_1
  
  
  
  
  
  
  # some model only include PSR
  if (length(forward_log$coefficients) == 2) {
    
    model_save_result <- model_save_OR_OwO_for_onevar(forward_log)
    
    route2 <- paste(mul_log_csv, model_name, '.csv', sep = '')
    write.csv(model_save_result[[1]], route2, row.names = F)
    
    route3 <- paste(mul_log_table, model_name, '.png', sep = '')
    gtsave(model_save_result[[2]], route3)
    
    # final auc
    pred_prob_final <- predict(forward_log, type = "response") 
    roc_obj_final <- roc(HHH$case, pred_prob_final, levels = c(0, 1), direction = "<")
    
    final_output_df[j, 'Model_Index'] <- model_name
    final_output_df[j, 'AUC'] <- auc(roc_obj_final)
    final_output_df[j, 'AIC'] <- forward_log$aic
    final_output_df[j, 'Variable_num'] <- 1
    # final_output_df[j, 'Exclude_num'] <- 0
    # final_output_df[j, 'Youden_Index'] <- youden_index
    final_output_df[j, 'Alpha'] <- ES_result[['alpha']]
    final_output_df[j, 'Lambda'] <- ES_result[['lambda']]
  } else {
    
    model_save_result <- model_save_OR_OwO(model = forward_log)
    
    route2 <- paste(mul_log_csv, model_name, '.csv', sep = '')
    write.csv(model_save_result[[1]], route2, row.names = F)
    
    route3 <- paste(mul_log_table, model_name, '.png', sep = '')
    gtsave(model_save_result[[2]], route3)
    
    
    
    
    final_model_df <- model_save_result[[1]]
    
    
    # final auc
    pred_prob_final <- predict(forward_log, type = "response") 
    roc_obj_final <- roc(HHH$case, pred_prob_final, levels = c(0, 1), direction = "<")
    
    
    
    
    final_output_df[j, 'Model_Index'] <- model_name
    final_output_df[j, 'AUC'] <- auc(roc_obj_final)
    final_output_df[j, 'AIC'] <- forward_log$aic
    final_output_df[j, 'Variable_num'] <- nrow(final_model_df)
    # final_output_df[j, 'Exclude_num'] <- exclude_num
    # final_output_df[j, 'Youden_Index'] <- youden_index
    final_output_df[j, 'Alpha'] <- ES_result[['alpha']]
    final_output_df[j, 'Lambda'] <- ES_result[['lambda']]
    
    
    
  }
  
  
  
}



write.csv(final_output_df, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/Over_all_table_20250521.csv', row.names = F)




# Coefficient P/N table in mul. reg.--------------------------------------------
mul_log_results <- match_pair_result(file_OwO = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/mul_log_csv/")

# Over_all_table_20250402 <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/Over_all_table_20250402.csv')


OwO_AAAA <- c()

for (i in seq_along(mul_log_results)) {
  OwO_aa <- mul_log_results[[i]]
  OwO_aaa <- OwO_aa$Variable
  OwO_AAAA <- union(OwO_AAAA, OwO_aaa)
}



over_all_tt <- matrix(data = 'x', nrow = length(mul_log_results), ncol = length(OwO_AAAA)) %>% as.data.frame()
colnames(over_all_tt) <- OwO_AAAA
row.names(over_all_tt) <- names(mul_log_results)


# j <- 5
for (j in seq_along(mul_log_results)) {
  OwO_aa <- mul_log_results[[j]]
  OwO_name <- names(mul_log_results)[[j]]
  
  
  for (k in 1:nrow(OwO_aa)) {
    if (OwO_aa[k, 'AOR'] > 1) {
      
      if (OwO_aa[k, 'P.value'] < 0.05) {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'P/S'
      } else {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'P/US'
      }
      
    } else {
      
      if (OwO_aa[k, 'P.value'] < 0.05) {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'N/S'
      } else {
        over_all_tt[row.names(over_all_tt) %in% OwO_name, colnames(over_all_tt) %in% OwO_aa[k, 'Variable']] <- 'N/US'
      }
      
    }
  }
  
  
  
  
  
  
} 





write.csv(over_all_tt, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/logistic_reg_with_PSR/logistic_regression/elastic_net_non_sig_retain/4_index_compete/Over_all_table/PN_VAR.csv')


