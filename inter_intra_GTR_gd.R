source('~/R/aiv/function.R')
suppressPackageStartupMessages(library(castor))
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("-dis", "--distance", type='character', default=TRUE,
                    help="GTR genetic distance matrix")

parser$add_argument("-seg", "--segment", type='character', default=TRUE,
                    help="Operating segment")

parser$add_argument("-iqg", "--iq_grouping", type='character', default=TRUE,
                    help="Iqtree grouping reuslt")

parser$add_argument("-op", "--output_prefix", type='character', default='agreement',
                    help="output file prefix")

parser$add_argument("-o", "--out_dir", type='character', default=TRUE,
                    help="Generated distance matrix store direction")
args <- parser$parse_args()


# args$distance <- '~/Analysis/aiv/merge/0307/sampling_subset/1015/NP_aa_sampling.fa.mldist'
# args$segment <- 'NP'
# args$iq_grouping <- '~/Analysis/aiv/merge/0307/lineages/NP_aa_ident_lineage.csv'
# args$output_prefix <- '_GTR_genetic_distance_aa_ident'
# args$out_dir <- '~/Analysis/aiv/merge/0307/distance/inter_intra_gd/'


dir.create(args$out_dir)
# data preparing -------------------------------------------------------------------
dis <- fread(args$distance) %>% as.data.frame()
colnames(dis)[1] <- 'seq'
dis <- dis %>% tibble::column_to_rownames('seq')
colnames(dis) <- rownames(dis) <- str_extract(rownames(dis), pattern = '[0-9]+\\|H|[0-9]+_H') %>% str_remove(pattern = '\\|H|_H')


iqtree_out <- fread(args$iq_grouping) %>% as.data.frame() %>% select(c('Strain_number', 'correct_out1'))
iqtree_out[colnames(iqtree_out)] <- lapply(iqtree_out[colnames(iqtree_out)], as.character)

# intra -------------------------------------------------------------------
intra_list <- list()
for(i in colnames(iqtree_out)[2:ncol(iqtree_out)]) {
  cat(paste('performing ', i, '\n', sep = ''))
  group <- iqtree_out[, c('Strain_number', i)]
  colnames(group)[2] <- 'group'
  result <- c()
  x <- 0
  for(j in unique(group$group)) {
    x <- x+1
    cat(paste('performing ', x, ',', length(unique(group$group)), '\n', sep = ''))
    clade <- group[group$group %in% j, 1] %>% as.character()
    n <- length(clade)
    clade <- which(rownames(dis) %in% clade)
    A <- dis[clade, clade]
    # A[1:5, 1:5] %>% print()
    result[j] <- (A %>% sum())/(n*n-n)*100
  }
  intra_list[[i]] <- result
}

print(intra_list)


# inter group -------------------------------------------------------------------
inter_list <- list()
for(i in colnames(iqtree_out)[2:ncol(iqtree_out)]) {
  cat(paste('performing ', i, '\n', sep = ''))
  group <- iqtree_out[, c('Strain_number', i)]
  colnames(group)[2] <- 'group'
  
  group_DF <- table_DF(group$group)
  result <- c()
  x <- 0
  for(j in 1:nrow(group_DF)) {
    x <- x+1
    cat(paste('performing ', x, ',', nrow(group_DF), '\n', sep = ''))
    for(k in 2:nrow(group_DF)) {
      if(k > j) {
        sub_dis <- dis[which(rownames(dis) %in% group[group$group %in% group_DF[j, 1], 1]), 
                       which(rownames(dis) %in% group[group$group %in% group_DF[k, 1], 1])]
        result[paste(group_DF$x[j], group_DF$x[k], sep = '_')] <- mean(sub_dis %>% as.matrix())*100
      }
    }
  }
  inter_list[[i]] <- result
}
print(inter_list)

saveRDS(list(intra=intra_list, inter=inter_list), 
        paste(args$out_dir, args$segment, args$output_prefix, '.RData', sep = ''))

 
# # intra_list <- readRDS(paste(args$out_dir, args$segment, args$output_prefix, '.RData', sep = ''))[[1]]
# dis_intra <- do.call(rbind, intra_list) %>% as.data.frame()
# for(i in seq_along(intra_list)) {
#   if(length(intra_list[[i]])==ncol(dis_intra)) {next}
#   else{dis_intra[i, (length(intra_list[[i]])+1):ncol(dis_intra)] <- ''}
# }
# dis_intra$m <- lapply(intra_list, mean) %>% unlist()
# dis_intra <- t(dis_intra) %>% as.data.frame() %>% tibble::rownames_to_column('group')
# dis_intra <- dis_intra %>% tidyr::gather(key, value, c(2:ncol(dis_intra)))
# dis_intra <- dis_intra[!dis_intra$value %in% '', ] %>% mutate(int='Intra', value=as.numeric(value))
# 
# 
# # inter_list <- readRDS(paste(args$out_dir, args$segment, args$output_prefix, '.RData', sep = ''))[[2]]
# result <- do.call(inter_list, what=rbind) %>% as.data.frame()
# for(i in seq_along(inter_list)) {
#   if(length(inter_list[[i]])==ncol(result)) {next}
#   else{result[i, (length(inter_list[[i]])+1):ncol(result)] <- ''}
# }
# result$m <- lapply(inter_list, mean) %>% unlist()
# inter <- t(result) %>% as.data.frame() %>% tibble::rownames_to_column('group')
# inter <- inter %>% tidyr::gather(key, value, c(2:ncol(inter)))
# inter <- inter[!inter$value %in% '', ] %>% mutate(int='Inter', value=as.numeric(value))
# 
# 
# A <- rbind(dis_intra, inter) 
# A$int <- factor(A$int, levels = c('Intra', 'Inter')) 
# A <- replace(A, is.na(A), 0)
# A$key <- (1-(str_remove(A$key, pattern = 'd_') %>% as.numeric())) %>% as.character()
# A$label <- paste(A$value %>% round(digits = 2), '%', sep = '')
# # A$key <- factor(A$key, levels = rev(A$key %>% unique()))
# # if(args$segment %in% 'PB2') {
# #   A$key <- (1-(str_remove(A$key, pattern = 'd_') %>% as.numeric())) %>% as.character()
# # }
# # if(args$segment != 'PB2') {
# #   A$key <- (str_remove(A$key, pattern = 'd_')) %>% as.character()
# # }
# 
# 
# png(paste('~/Analysis/aiv/merge/0307/result/', args$segment, args$output_prefix, '.png', sep = ''), 
#     width = 9600, height = 5400, res = 600)
# ggplot(A %>% filter(group != 'm'))+
#   # geom_boxplot(aes(x=key, y=value, color=key))+
#   geom_boxplot(aes(x=key, y=value, fill=int),  position = position_dodge(1))+
#   geom_point(aes(x=key, y=value, fill=int),  position = position_dodge(1), size = 2)+
#   geom_text(data = A %>% filter(group %in% 'm', int == 'Intra')
#             , mapping=aes(x=key, y=-1, label=label), color='#E68865'
#             , size=4, show.legend = F)+ #1.65*6/length(unique(A$key))
#   geom_text(data = A %>% filter(group %in% 'm', int == 'Inter')
#             , mapping=aes(x=key, y=(max((A %>% filter(int == 'Inter'))$value))*1.1
#                           , label=label), color='#3BAEDA'
#             , size=4, show.legend = F)+ #0.6*6/length(unique(A$key))
#   scale_y_continuous(limits = c(-1, (max((A %>% filter(int == 'Inter'))$value))*1.15), 
#                      breaks =  c(0, seq(2.5, (max(A$value))+2, 2.5) ))+
#   scale_fill_manual(values=c('#E68865', '#3BAEDA'), labels=c('Intra (within groups)', 'Inter (between groups)'))+
#   labs(x='Patristic Distance Threshold', y="Genetic distance (%)", color='Inter or\nIntra groups'
#        , title = paste(args$segment, ", Genetic distance (Kimura's 2 parameter)", sep = ''), fill='')+
#   theme_bw()+gg_theme+theme(axis.text.x = element_text(size = 10))
# dev.off()
# # plot_list <- list()
# # for(i in unique(A$key)) {
# #   # plot_list[[i]] <- 
# #   
# # }
# 
# 
# 



