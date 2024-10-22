source('~/R/aiv/function.R')
suppressPackageStartupMessages(library(castor))
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("-seg", "--segment", type='character', default=TRUE,
                    help="Segments compute the distance matrix")

parser$add_argument("-iqt", "--iqtree", type='character', default=TRUE,
                    help="Iqtree direction")

parser$add_argument("-iqg", "--iq_grouping", type='character', default=TRUE,
                    help="Iqtree grouping reuslt")

parser$add_argument("-ftg", "--ft_grouping", type='character', default=TRUE,
                    help="fasttree grouping reuslt")

parser$add_argument("-op", "--output_prefix", type='character', default='agreement',
                    help="output file prefix")

parser$add_argument("-o", "--out_dir", type='character', default=TRUE,
                    help="Generated distance matrix store direction")
args <- parser$parse_args()


# args$segment <- 'MP'
# args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/MP_5p_sampling.nexus'
# args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/MP_iq_MPD_groups.RData'
# args$ft_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/MP_ft_MPD_groups.RData'
# args$output_prefix <- '_as_sampling'
# args$out_dir <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/'



args$segment <- 'PB2'
args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/PB2_5p_sampling.nexus'
args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/PB2_iq_MPD_groups.RData'
args$ft_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/PB2_ft_MPD_groups.RData'
args$output_prefix <- '_as_sampling'
args$out_dir <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/'


# args$segment <- 'N8'
# args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_N8_aligned_iqtree.nexus'
# args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/NA_1009/N8_iq_MPD_groups.RData'
# args$ft_grouping <- '~/Analysis/aiv/merge/0307/distance/NA_1009/N8_ft_MPD_groups.RData'
# args$output_prefix <- ''
# args$out_dir <- '~/Analysis/aiv/merge/0307/distance/NA_1009/'


# mergeresult ------------------------------------------------------------
dir.create(args$out_dir)


tree <- treeio::read.beast(args$iqtree)

out <- readRDS(args$ft_grouping)[[2]]

rank <- colnames(out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V2=as.numeric(V2), name=paste(V1, V2, sep = '_')) %>%
  arrange(V2) %>% select(name) %>% unlist()

out <- out[, c('Strain_number', rank)]
if(all(out[, ncol(out)] == out[, (ncol(out)-1)])) {out <- out[, -ncol(out)]}

# apply(out, MARGIN = 2, function(x){x %>% is.na() %>% table()})



iqtree_out <- readRDS(args$iq_grouping)[[2]]

while(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {
  iqtree_out <- iqtree_out[, -ncol(iqtree_out)]
}
# if(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {iqtree_out <- iqtree_out[, -ncol(iqtree_out)]}
if(colnames(iqtree_out)[2] %in% 'd_0.0075') {
  iqtree_out <- iqtree_out %>% select(-c('d_0.0075'))
}


nc <- min(ncol(out), ncol(iqtree_out))
# u <- list()
# for(j in 2:nc) {
#   A <- out[, c(1, j)]
#   colnames(A)[2] <- 'ft'
#   B <- iqtree_out[, c(1, j)]
#   colnames(B)[2] <- 'iq'
# 
#   cat(paste(colnames(out)[j], ', ', j, ', ', nc, ', ', nrow(out), '\n', sep = ''))
# 
#   test <- merge(A, B, by = 'Strain_number')
#   if(all(complete.cases(test)) == FALSE) next
#   rownames(test) <- test$Strain_number
#   # test <- test[rank, ]
# 
#   # nc <- min(ncol(out), ncol(iqtree_out))
#   # u_sub <- matrix(0, nc-1, ncol = 5) %>% as.data.frame()
#   # colnames(u_sub) <- c('u', 'u-same', 'u-distinct', 'ft_groups', 'iq_groups')
#   # rownames(u_sub) <- colnames(out)[2:nc]
# 
#   result <- list()
#   for(i in 1:nrow(test)) {
#     # sub_ft <-
#     ft_gg_mat <- ifelse(test[i, 'ft']==test[i:nrow(test), 'ft'], yes = 1, no = 0)
#     iq_gg_mat <- ifelse(test[i, 'iq']==test[i:nrow(test), 'iq'], yes = 1, no = 0)
# 
#     result[[i]] <- (ft_gg_mat+iq_gg_mat) %>% table_DF() %>% mutate(group=paste(x, Freq, sep = '_')) %>%
#       select(group) %>% unlist()
#   }
# 
# 
#   result <- result %>% unlist()
#   n <- nrow(out)
#   grid <- (n*(n-1)/2)
#   agree_diff <- (str_remove(result[grepl(pattern='0_', result, )], pattern = '0_') %>% as.numeric() %>% sum())
#   disagree <- str_remove(result[grepl(pattern='1_', result, )], pattern = '1_') %>% as.numeric() %>% sum()
#   agree_same <- (str_remove(result[grepl(pattern='2_', result, )], pattern = '2_') %>% as.numeric() %>% sum())-n
# 
# 
#   bs <- tree@data %>% mutate(node_new=as.numeric(node)-length(tree@phylo$tip.label)) %>% as.data.frame() %>%
#     filter(node_new %in% unique(test$iq)) %>% mutate(dis=colnames(out)[j], group='bs') %>% rename('y'='bs') %>%
#     select(c('dis', 'group', y))
# 
# 
#   name <- colnames(out)[j]
#   agreement <- (agree_diff+agree_same)/grid
#   u_sub <- matrix(0, nrow = 5, ncol = 3) %>% as.data.frame()
#   colnames(u_sub) <- c('dis', 'group', 'y')
#   u_sub$dis <- colnames(out)[j]
#   u_sub$group <- c('u', 'u-same', 'u-distinct', 'ft_groups', 'iq_groups')
#   u_sub$y <- c(agreement, (agree_same)/grid, (agree_diff)/grid, out[, c(j)] %>% unique() %>% length(), iqtree_out[, c(j)] %>% unique() %>% length())
#   u[[(j-1)]] <- rbind(u_sub, bs)
#   # u[(j-1), 1] <- agreement # consistence amoung iqtree and fasttree
#   # u[(j-1), 2] <- (agree_same)/grid#/agreement
#   # u[(j-1), 3] <- (agree_diff)/grid#/agreement
#   # u[(j-1), 4] <- out[, c(j)] %>% unique() %>% length()
#   # u[(j-1), 5] <- iqtree_out[, c(j)] %>% unique() %>% length()
#   gc()
# }


A <- readRDS(paste(args$out_dir, args$segment, args$output_prefix, '_agreement_result.RData', sep = ''))
u <- A$agreement_table


plot_data <- u %>% do.call(what=rbind) %>% as.data.frame()
# plot_data <- u %>% tibble::rownames_to_column('dis') %>% tidyr::gather(key, value, 2:(ncol(u)+1))
# plot_data$dis <- str_remove(string = plot_data$dis, pattern = 'd_')
# g <- plot_data[plot_data$group %in% 'u', ] %>% nrow()
plot_data[plot_data$group %in% "u", 'group'] <- "Agreement"
plot_data[plot_data$group %in% "u-same", 'group'] <- "Agreement-same"
plot_data[plot_data$group %in% "u-distinct", 'group'] <- "Agreement-diff"
colnames(plot_data)[2:3] <- c('key', 'value')
plot_data[plot_data$key %in% c('Agreement', 'Agreement-same', 'Agreement-diff'), 'value'] <-
  plot_data[plot_data$key %in% c('Agreement', 'Agreement-same', 'Agreement-diff'), 'value']*100
plot_data$dis <- str_remove(plot_data$dis, pattern = 'd_') %>% as.numeric()

plot_data <- rbind(plot_data,
                   data.frame(dis=plot_data$dis %>% unique(),
                              key='Disagreement',
                              value=max(plot_data[plot_data$key %in% c('Agreement'), 'value']) -plot_data[plot_data$key %in% c('Agreement'), 'value']))
plot_data$dis <- factor(plot_data$dis, levels = unique(plot_data$dis))
# plot_data[plot_data$key %in% c('Agreement'), 'value'] <-
# plot_data[plot_data$key %in% "Agreement", 'key'] <- "Disagreement"

dis_bs <- data.frame(mean=plot_data %>% filter(key == c('bs')) %>% group_split(dis) %>%
                       lapply(FUN = function(x){median(x$value) %>% round(digits = 1)}) %>% unlist(),
                     x=(plot_data %>% filter(key == c('bs')))$dis %>% unique())
iq_groups <- plot_data %>% filter(key == c('iq_groups'))

as <-
  plot_data %>% filter(key %in% c('Agreement', 'Agreement-same', 'iq_groups')) %>%
  reshape(idvar = "dis", timevar = "key", direction = "wide")
colnames(as) <- c('dis', 'AS', 'AS_same', 'iq_groups')

if(args$segment %in% c(paste('N', 1:9, sep = ''))){
  as <- as %>% mutate(diff1=c(AS_same[-1], AS_same[nrow(as)])-AS_same, 
                      diff2=c(0, diff1[-nrow(as)])) %>% filter(AS > 80, iq_groups >= 3, iq_groups <= 20)
}else{
  as <- as %>% mutate(diff1=c(AS_same[-1], AS_same[nrow(as)])-AS_same, 
                      diff2=c(0, diff1[-nrow(as)])) %>% filter(AS > 75, iq_groups > 3, iq_groups <= 20)
}


# elbow point exist
if(any(as$diff1 >= 6 & as$diff2 <2)) {
  as_pick <- (as[as$diff1 >=6 & as$diff2 <2, ])
  as_pick <- as_pick[which.max(as_pick$AS), 1] %>% unlist() %>% as.character()
  }
# no elbow point
if(isFALSE(any(as$diff1 >= 6 & as$diff2 <2))) {as_pick <- as[which.max(as$AS_same+as$AS), 1] %>% as.character()}


rec <- grep(plot_data$dis %>% unique(), pattern=paste(as_pick, '$', sep = ''))
as_pick <- grep(colnames(iqtree_out), pattern=paste(as_pick, '$', sep = ''))-1


if(args$segment %in% c(paste('N', 1:9, sep = ''))){
  dis <- plot_data %>% filter(key %in% 'iq_groups', value >=3, value<= 20) %>% select(dis) %>% unlist %>% as.character()
}else{
  dis <- plot_data %>% filter(key %in% 'iq_groups', value >3, value<= 20) %>% select(dis) %>% unlist %>% as.character()
}

dmax <- grep(plot_data$dis %>% unique(), pattern=paste(max(dis), '$', sep = ''))
dmin <- grep(plot_data$dis %>% unique(), pattern=paste(min(dis), '$', sep = ''))

# as <- plot_data %>% filter(key %in% c('Agreement', 'Agreement-same')) %>% 
#   reshape(idvar = "dis", timevar = "key", direction = "wide") %>% filter(`value.Agreement-same` <= 50 ) %>% 
#   mutate(sum=c(`value.Agreement` +`value.Agreement-same`)) 
# as <- which.max(as$sum)
# if(identical(args$segment, 'MP')) {as <- 8}

max_g <- max(plot_data[plot_data$key == 'iq_groups', 'value'])
coeff <- max_g/100
print(coeff)
# sec_y <- max(plot_data[1:g, 'value'])
agreement <-
  ggplot()+
  # geom_violin(plot_data %>% filter(key %in% c('bs')), mapping=aes(x=dis, y=value, group=dis), fill='#958BE3', linewidth=1.2, alpha=0.4)+
  # geom_jitter(plot_data %>% filter(key %in% c('bs')), mapping=aes(x=dis, y=value, group=dis), fill='black')+
  geom_col(plot_data %>% filter(key %in% c('iq_groups')), mapping=aes(x=dis, y=value, group=dis), fill='grey75')+
  geom_hline(mapping = aes(yintercept = 80* coeff), linetype = 'dashed', linewidth = 1.5)+
  geom_point(plot_data %>% filter(key %in% c('Agreement', 'Agreement-same')), mapping=aes(x=dis, y=value* coeff, color=key), size=4)+
  geom_line(plot_data %>% filter(key %in% c('Agreement', 'Agreement-same')), mapping=aes(x=dis, y=value* coeff, group=key, color=key), linewidth=1.8)+
  geom_text(plot_data %>% filter(key %in% c('iq_groups')), mapping=aes(x=dis, y=value, label = value), vjust=-0.5, size=5)+
  annotate("rect", xmin=rec-0.47 , xmax=rec+0.5, ymin=-2 , ymax=max_g+1, alpha=0.2, color="red", fill=NA, linewidth=1)+
  annotate("rect", xmin=dmin-0.525 , xmax=dmax+0.525, ymin=-4 , ymax=max_g+5, alpha=1, color="grey25", fill=NA, linewidth=1)+
  # annotate(geom="text", x=0, y=104, size=5, label="Bootstrap value\nmedian", color="black")+
  # geom_point(plot_data %>% filter(key %in% c('ft_groups', 'iq_groups')), mapping=aes(x=dis, y=value * coeff, color=key), size=3)+
  # geom_line(plot_data %>% filter(key %in% c('ft_groups', 'iq_groups')), mapping=aes(x=dis, y=value * coeff, group=key, color=key), linewidth=1.2)+
  # geom_hline(yintercept = 100, color='red', linetype='dashed')+
  # geom_text(dis_bs, mapping=aes(x=x, y=104, label=mean), size=5)+
  # geom_text(plot_data %>% filter(key == 'iq_groups'), mapping=aes(x=dis, y=104, label=value), size=5)+
  # geom_text(iq_groups, mapping=aes(x=dis, y=-3, label=value), size=5)+
  labs(x='MPD thresholds', title = paste(args$segment, ' fasttree and iqtree partition consistence evaluation', sep = ''),
       color='')+
  scale_y_continuous(name = "Assigned lineages", 
                     sec.axis = sec_axis(~./coeff, name="Agreement Statistics", breaks = seq(0, 100, 10), ), 
                     limits = c(-5, max_g+10), labels= c(seq(0, max_g, 10)), breaks= c(seq(0, max_g, 10)))+
  scale_color_manual(values = c('#FFB161', '#66B3FF'))+
  # scale_y_continuous(name = "Agreement/Bootstrap", limits = c(-5, 106))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=18, angle = 45, vjust=1.05, hjust=1), legend.position = 'bottom')
gc()




# for(i in names(plot_list)) {
#   # png(paste('~/Analysis/aiv/merge/0307/plot/', i, "_ft_tree.png"), width = 1600, height = 900)
#   # print(ft_list[[i]])
#   # dev.off()
#   png(paste('~/Analysis/aiv/merge/0307/plot/', i, "_agree_same_dsitinct.png"), width = 1600, height = 900)
#   print(plot_list[[i]])
#   dev.off()
#   # png(paste('~/Analysis/aiv/merge/0307/plot/', i, "_iq_tree.png"), width = 1600, height = 900)
#   # print(iq_list[[i]])
#   # dev.off()
# }

png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment, '_agreement_1015.png', sep = ''), width=9600, height=5400, res = 600)
agreement
dev.off()
# 
# cat(paste('Saving agreement\n', sep = ''))
# saveRDS(list(agreement_table=u, agreement_plot=agreement),
#         file = paste(args$out_dir, args$segment, '_agreement_result.RData', sep = ''))



# iqtree tree vis ---------------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame() %>% mutate(Strain_number=as.character(Strain_number))
tree_list <- readRDS('/home/eric/Analysis/aiv/merge/0307/iq_tree/internal_NA_tree.RData')
# colnames(meta) <- c('Isolate_Id', 'Strain_number', 'Subtype', 'segment', 'Collection_Date', 'Location', 'Host'
#                     , 'seq', 'Date')

clade_table <- table_DF(meta[grepl(meta$Subtype, pattern='H5'), 'Clade'])
clade_table[, 'Clade'] <- 'GsGD-others'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4'), 'Clade'] <- '2.3.4.4'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4b'), 'Clade'] <- '2.3.4.4b'
clade_table[grep(clade_table$x, pattern = 'nonGsGD'), 'Clade'] <- 'nonGsGD'
clade_table[grep(clade_table$x, pattern = 'IRD_unlabel'), 'Clade'] <- 'IRD_unlabel'
nc <- min(ncol(out), ncol(iqtree_out))

cat(paste('Loading tree\n', sep = ''))



# iqtree_out <- readRDS(args$iq_grouping)[[2]]
# rank <- colnames(iqtree_out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V2=as.numeric(V2), name=paste(V1, V2, sep = '_')) %>%
#   arrange(V2) %>% select(name) %>% unlist()
# 
# iqtree_out <- iqtree_out[, c('Strain_number', rank)]
# if(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {iqtree_out <- iqtree_out[, -ncol(iqtree_out)]}
# if(colnames(iqtree_out)[2] %in% 'd_0.0075') {
#   iqtree_out <- iqtree_out %>% select(-c('d_0.0075'))
# }

# tree <- read_as_list(path = args$tree_dir, file_type = 'nexus',
#                      prefix = paste(args$segment, args$tree_prefix , '$', sep = ''))[[1]]


tree <- treeio::read.beast(args$iqtree)
tip_name <- tree@phylo$tip.label %>% str_remove_all(pattern = "'")
# tree <- ape::read.nexus(paste("~/Analysis/aiv/merge/0307/gisaid_IRD_merged_", seg,"_aligned_iqtree.nwk", sep = ''))

info <- data.frame(Strain_number=tip_name %>% str_extract(pattern = '_[0-9]+_H|\\|[0-9]+\\|H') %>% str_remove(pattern = '_|\\|') %>%
                     str_remove(pattern = '_H|\\|H'))

rank <- info %>% unlist()

info <- merge(info, meta[, c('Strain_number', 'new_header')],
              by = 'Strain_number', all.x = T)

info <- header_cleaning(info$new_header, pattern = '\\|')

colnames(info) <- c("Accesion_number","Strain_number","Subtype", 'Clade', "Segment", 'Collection_date'
                    ,"Location","Host", "Header_Host", "Host_type")

country <- meta[meta$Strain_number %in% info$Strain_number, c('Strain_number', 'Location')]
rownames(country) <- country$Strain_number
country <- country[info$Strain_number, ]
info$Location <- country$Location

info[info$Accesion_number %in% 'EPI_ISL_66116', 'Collection_date'] <- 1970
info$Host_type <- replace(info$Host_type, info$Host_type %in% '', 'Wild')
info$Year <- word(info$Collection_date, sep = '-', 1) %>% as.numeric()

info$continent <- word(info$Location, 1, sep = "/")

info[, 'time_group'] <- 0
info[info$Year < 1996, 'time_group'] <- "before_1995"
info[info$Year <= 2016 & info$Year >= 1996, 'time_group'] <- "1996~2013"
info[info$Year <= 2020 & info$Year >= 2016, 'time_group'] <- "2014~2020"
info[info$Year >=2021, 'time_group'] <- "after_2021"


info$time_group <- factor(info$time_group, levels = c("before_1995", "1996~2013", "2014~2020", "after_2021"))

info$h <- str_extract(info$Subtype, pattern = 'H[0-9]+')
info$h <- factor(info$h, levels = data.frame(l=str_extract(unique(info$h), pattern = '[A-Z]+')
                                             , n=str_extract(unique(info$h), pattern = '[0-9]+') %>% as.numeric()) %>% arrange(n) %>% mutate(level=paste(l, n, sep = '')) %>% select(level))
info$n <- str_extract(info$Subtype, pattern = 'N[0-9]')
info$n <- factor(info$n, levels = data.frame(l=str_extract(unique(info$n), pattern = '[A-Z]+'), n=str_extract(unique(info$n), pattern = '[0-9]+') %>% as.numeric()) %>%
                   arrange(n) %>% mutate(level=paste(l, n, sep = '')) %>% select(level))

info[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  info[info$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
info$Clade_new <- factor(info$Clade_new , levels = c("2.3.4.4","2.3.4.4b", "GsGD-others","nonGsGD", 'IRD_unlabel'))



info[, 'Subtype_new'] <- 'Other subtype'

# info[info$Subtype %in% c('MPN1', 'MPN2', 'MPN6', 'MPN8', 'H7N9', 'H9N2'), 'Subtype_new'] <-
#   info[info$Subtype %in% c('MPN1', 'MPN2', 'MPN6', 'MPN8', 'H7N9', 'H9N2'), 'Subtype']

retain <- (info$Subtype %>% table_DF() %>% arrange(desc(Freq)))[1:8, ] %>% select(x) %>% unlist
info[info$Subtype %in% retain, 'Subtype_new'] <-
  info[info$Subtype %in% retain, 'Subtype']
info[, 'Subtype_new'] <- factor(info[, 'Subtype_new'], levels = c(retain, 'Other subtype'))

## group assignment (manually) sequences groups
# g <- (min(ncol(out), ncol(iqtree_out)))
# out_sub <- out[, 1:g]
# iqtree_out_sub <- iqtree_out[, 1:g]
# colnames(out_sub) <- c('Strain_number', paste('FT', str_replace(colnames(out_sub)[-1], pattern = 'd_', replacement = '\n'), sep = ''))
# colnames(iqtree_out_sub) <- c('Strain_number', paste('IQ', str_replace(colnames(iqtree_out_sub)[-1], pattern = 'd_', replacement = '\n'), sep = ''))
# FT_IQ_group <- data.frame(Strain_number=out_sub$Strain_number)
# for(i in 2:g) {
#   FT_IQ_group <- merge(FT_IQ_group, out_sub[, c(1, i)], by = 'Strain_number')
#   FT_IQ_group <- merge(FT_IQ_group, iqtree_out_sub[, c(1, i)], by = 'Strain_number')
# }


# iqtree_out <- iqtree_out[, c(1, c(1+as))]
# colnames(iqtree_out)[2] <- 'Lineage'
# info <- merge(info, iqtree_out, by='Strain_number', all=T)
# rownames(info) <- info$Strain_number
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# group <- table_DF(info$Lineage) %>% filter(Freq > nrow(iqtree_out)*0.01)
# info[!info$Lineage %in% group$x, 'Lineage'] <- NA
#
# for(k in is.na(info[, 'Lineage']) %>% which()) {
#   if(k==1) {info[k, 'Lineage'] <- info[k+1, 'Lineage']}
#   else(info[k, 'Lineage'] <- info[k-1, 'Lineage'])
# }

# info <- merge(info, FT_IQ_group, by='Strain_number', all=T)
# rownames(info) <- info$Strain_number
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))

info <- merge(info, iqtree_out[, c(1, as_pick+1)], by = 'Strain_number')
rownames(info) <- info$Strain_number
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))

bs <- tree@data %>% mutate(node=as.numeric(node)-length(tree@phylo$tip.label), index=1:nrow(tree@data)) %>% as.data.frame()

bs1 <- merge((table_DF(info[, 18]) %>% mutate(x=as.numeric(x)-1)), bs, by.x='x', by.y='node', all.y=T) %>%
  mutate(diff= x-index)

nodes <- bs1[rownames(bs1[!is.na(bs1$Freq) & bs1$diff <1, ]), ]
A <- rep(NA, nrow(bs))
A[(nodes$index %>% as.numeric())] <- nodes$bs %>% as.numeric()
tree@data$nodes <- A


top20 <- bs1[(which(bs1$x <= 12)), ]
A <- rep(NA, nrow(bs))
A[(top20$index %>% as.numeric())] <- top20$bs %>% as.numeric()
tree@data$top20 <- A


for(j in 18:ncol(info)) {
  group <- table_DF(info[, j]) %>% arrange(desc(Freq))
  if(nrow(group) > 18) {
    group <- group[1:18, ]
  }# else(group <- group %>% filter(Freq>100))
  rownames(group) <- group$x
  info[!info[, j] %in% group$x, j] <- NA
  remain_g <- info[, j] %>% unique()
  group <- group[remain_g[!is.na(remain_g)], ]
  group[, 'edit'] <- paste('L', 1:nrow(group), sep = '')


  for(k in seq_len(nrow(group))) {
    info[info[, j] %in% group[k, 'x'], j] <- group[k, 3]
  }
  # info[!info$g %in% group$x, 'clade'] <- 'Minor'
  info[, j] <- factor(info[, j], levels = c(paste('L', 1:18, sep = ''), NA))
}

# info <- merge(info, segment_groups[[args$segment]], by = 'Strain_number')


rownames(info) <- info$Strain_number %>% as.character()
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
info$Strain_number <- factor(info$Strain_number, levels =rev(rank))
# info$data_source <- str_remove_all(info$Accesion_number, pattern = '_[0-9]+')
# info[!(info$data_source %in% 'IRD'), 'data_source'] <- 'GISAID'

#change tree header to match the data
new_header <- info$Strain_number %>% as.character()
tree@phylo$tip.label <- new_header
#
clade <- info[, c('Strain_number', 'continent')]
# tree <- tree_list[[args$segment]]
ph_tree <- ggtree(tree)
n <- length(clade$continent %>% unique())
set.seed(12)
# head(clade)
# phylo <-
#   ph_tree %<+% clade +
#   geom_tippoint(aes(color = continent), size = 3) +
#   scale_color_manual(values=rainbow(n)[sample(1:n, size = n, replace = F)],
#                      na.value = "transparent")+
#   # scale_color_brewer(palette='Accent')+# scale_colo(palette = 'Accent')+
#   theme(legend.position = "right")  # Display the legend on the right

# tree1 <- tree
# tree1@phylo$tip.label <- paste(info$Strain_number %>% as.character(), info[, 'lineages'] %>% as.character(), sep = '_')
# write.tree(tree1 %>% as.phylo(),
#            file = paste('~/Analysis/aiv/merge/0307/sampling_subset/0924/iq/', args$segment, '.nexus', sep = ''))

clade <- info[, c(1, 18)]
# clade <- info[, c('Strain_number', 'lineages')]
colnames(clade)[2] <- 'lineage'
lineage_rank <-
  clade$lineage %>% unique() %>% header_cleaning('L') %>% mutate(V1='L', V2=as.numeric(V2)) %>% arrange(V2) %>%
  mutate(p=paste(V1, V2, sep = '')) %>% select(p) %>% unlist()
clade$lineage <- factor(clade$lineage, levels = lineage_rank)

# head(clade)
phylo_bs <-
  ph_tree %<+% clade +
  geom_tippoint(aes(color = lineage), size = 3) +
  geom_text(aes(label=top20), color='blue', hjust = 1.1, vjust=-0.4) +
  geom_text(aes(label=nodes), color='red', hjust = 1.1, vjust=-0.4) +
  scale_color_manual(values=custom_color,
                     na.value = "transparent")+
  # scale_color_brewer(palette='Accent')+# scale_colo(palette = 'Accent')+
  theme(legend.position = "right")

# png(paste('/home/eric/Analysis/aiv/merge/0307/result/new_', args$segment, '_', args$output_prefix,'_ftiq_tree_phylo.png', sep = ''), width=9600, height=5400, res = 600)
# # (Heatmap_6 %>% insert_left(phylo, width = 12) %>% insert_right(Heatmap_1))
# phylo
# dev.off()


cat(paste('Plotting tree\n', sep = ''))
#Heatmap from sequence information
Heatmap_1 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=continent, color=continent), linewidth=0.05)+
  labs(x="continent", fill='continent', color='continent') + ylab(NULL)+
  scale_color_brewer(palette = 'Accent')+
  scale_fill_brewer(palette = 'Accent')+
  phylo_theme

Heatmap_2 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Clade_new, color=Clade_new), linewidth=0.05)+
  labs(x="Clade", fill="Clade", color="Clade") + ylab(NULL)+
  scale_color_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  scale_fill_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  # scale_color_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGd-others'='#5AA977','nonGsGd'='#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  # scale_fill_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGd-others'='#5AA977','nonGsGd'='#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  phylo_theme

Heatmap_3 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Subtype_new, color=Subtype_new), linewidth=0.05)+
  labs(x="Subtype", fill="Subtype", color="Subtype") + ylab(NULL)+
  # scale_color_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
  # scale_fill_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
  scale_color_brewer(palette = 'Set2')+
  scale_fill_brewer(palette = 'Set2')+
  phylo_theme


# Heatmap_4 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Host_type, color=Host_type), linewidth=0.05)+
#   labs(x="Host", fill="Host", color="Host") + ylab(NULL)+
#   scale_color_brewer(palette = 'Set3')+
#   scale_fill_brewer(palette = 'Set3')+
#   phylo_theme


# Heatmap_5 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=order_new, color=order_new), linewidth=0.05)+
#   labs(x="Order", fill="Order", color="Order") + ylab(NULL)+
#   scale_color_manual(values = c('#FAA508','#647531','#314975','#6B3175','#EBF2F7'))+
#   scale_fill_manual(values = c('#FAA508','#647531','#314975','#6B3175','#EBF2F7'))+
#   phylo_theme


Heatmap_6 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=time_group, color=time_group), linewidth=0.05)+
  labs(x="Year", fill="Year", color="Year") + ylab(NULL)+
  scale_color_brewer(palette = 'Blues')+
  scale_fill_brewer(palette = 'Blues')+
  phylo_theme


Heatmap_7 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=`FT\n0.01`, color=`FT\n0.01`), linewidth=0.05)+
  labs(x='FT\n1', fill="Grouping", color="Grouping") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(axis.title.x = element_text(size=5))

# ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=`IQ\n0.05`, color=`IQ\n0.05`), linewidth=0.05)+
#   labs(x='FT\n0.025', fill="Grouping", color="Grouping") + ylab(NULL)+
#   scale_color_brewer(palette = 'Set3')+
#   scale_fill_brewer(palette = 'Set3')+
#   phylo_theme
#
# group_Heatmaps <- list()
# for(i in colnames(info)[19:(ncol(info))]) {
#   A <- info[, c('Strain_number', i)]
#   colnames(A)[2] <- 'test'
#   group_Heatmaps[[i]] <-
#     ggplot(A)+geom_tile(aes(x="",y=Strain_number, fill=test, color=test), linewidth=0.05)+
#     labs(x=i, fill="d_0.05", color="d_0.05") + ylab(NULL)+
#     scale_color_brewer(palette = 'Set3')+
#     scale_fill_brewer(palette = 'Set3')+
#     phylo_theme+theme(legend.position = "none", axis.title.x = element_text(size=5))
# }

# (Heatmap_6 %>% insert_left(phylo, width = 12) %>% insert_right(Heatmap_1))
p <- (Heatmap_1 %>% insert_left(phylo_bs, width = 12))
# for(i in seq_along(group_Heatmaps)) {
#   p <- p %>% insert_right(group_Heatmaps[[i]])
# }

# p <- (Heatmap_1 %>% insert_left(phylo, width = 12))
p <- as.ggplot(p %>% insert_right(Heatmap_3) %>% insert_right(Heatmap_2) %>%
                 insert_right(Heatmap_6) )+ #%>% insert_right(Heatmap_7)
  ggtitle(paste(args$segment, "segment phylogenetic tree (iqtree)"))+
  theme(plot.title = element_text(face = "bold",size = 25 ,hjust = 0.5))

# cat(paste('Saving whole result', sep = ''))
# saveRDS(list(agreement_table=u, agreement_plot=agreement, tree=p, FT_IQ_group=FT_IQ_group),
#         file = paste(args$out_dir, args$segment,'_agreement_result_iqtree.RData', sep = ''))
png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment,
          '_iq_tree_1015.png', sep = ''), width=9600, height=5400, res = 600)
p
dev.off()
# 
# 
# 
# 
# tree_list <- readRDS('/home/eric/Analysis/aiv/merge/0307/iq_tree/internal_NA_tree.RData')
# segment <- args$segment
# iqtree <- paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', args$segment, '_aligned_iqtree.nexus', sep = '')
# output_prefix <- 'lineage_0919'
# out_dir <- '~/Analysis/aiv/merge/0307/result/'
# 
# 
# cat(paste('Loading tree\n', sep = ''))
# tree <- ape::read.nexus(iqtree)
# 
# tip_name <- tree$tip.label %>% str_remove_all(pattern = "'") # 286980
# 
# info <- data.frame(Strain_number=tip_name %>% str_extract(pattern = '_[0-9]+_H|\\|[0-9]+\\|H') %>% str_remove(pattern = '_|\\|') %>%
#                      str_remove(pattern = '_H|\\|H'))
# 
# rank <- info %>% unlist()
# 
# info <- merge(info, meta[, c('Strain_number', 'new_header')],
#               by = 'Strain_number', all.x = T)
# 
# info <- header_cleaning(info$new_header, pattern = '\\|')
# 
# colnames(info) <- c("Accesion_number","Strain_number","Subtype", 'Clade', "Segment", 'Collection_date'
#                     ,"Location","Host", "Header_Host", "Host_type")
# 
# country <- meta[meta$Strain_number %in% info$Strain_number, c('Strain_number', 'Location')]
# rownames(country) <- country$Strain_number
# country <- country[info$Strain_number, ]
# info$Location <- country$Location
# 
# info[info$Accesion_number %in% 'EPI_ISL_66116', 'Collection_date'] <- 1970
# info$Host_type <- replace(info$Host_type, info$Host_type %in% '', 'Wild')
# info$Year <- word(info$Collection_date, sep = '-', 1) %>% as.numeric()
# 
# info$continent <- word(info$Location, 1, sep = "/")
# 
# info[, 'time_group'] <- 0
# info[info$Year < 1996, 'time_group'] <- "before_1995"
# info[info$Year <= 2016 & info$Year >= 1996, 'time_group'] <- "1996~2013"
# info[info$Year <= 2020 & info$Year >= 2016, 'time_group'] <- "2014~2020"
# info[info$Year >=2021, 'time_group'] <- "after_2021"
# 
# 
# info$time_group <- factor(info$time_group, levels = c("before_1995", "1996~2013", "2014~2020", "after_2021"))
# 
# info$h <- str_extract(info$Subtype, pattern = 'H[0-9]+')
# info$h <- factor(info$h, levels = data.frame(l=str_extract(unique(info$h), pattern = '[A-Z]+')
#                                              , n=str_extract(unique(info$h), pattern = '[0-9]+') %>% as.numeric()) %>% arrange(n) %>% mutate(level=paste(l, n, sep = '')) %>% select(level))
# info$n <- str_extract(info$Subtype, pattern = 'N[0-9]')
# info$n <- factor(info$n, levels = data.frame(l=str_extract(unique(info$n), pattern = '[A-Z]+'), n=str_extract(unique(info$n), pattern = '[0-9]+') %>% as.numeric()) %>%
#                    arrange(n) %>% mutate(level=paste(l, n, sep = '')) %>% select(level))
# 
# info[, 'Clade_new'] <- 'nonGsGD'
# for(k in seq_len(nrow(clade_table)))  {
#   info[info$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
# }
# info$Clade_new <- factor(info$Clade_new , levels = c("2.3.4.4","2.3.4.4b", "GsGD-others","nonGsGD", 'IRD_unlabel'))
# 
# 
# 
# info[, 'Subtype_new'] <- 'Other subtype'
# 
# # info[info$Subtype %in% c('MPN1', 'MPN2', 'MPN6', 'MPN8', 'H7N9', 'H9N2'), 'Subtype_new'] <-
# #   info[info$Subtype %in% c('MPN1', 'MPN2', 'MPN6', 'MPN8', 'H7N9', 'H9N2'), 'Subtype']
# 
# retain <- (info$Subtype %>% table_DF() %>% arrange(desc(Freq)))[1:8, ] %>% select(x) %>% unlist
# info[info$Subtype %in% retain, 'Subtype_new'] <-
#   info[info$Subtype %in% retain, 'Subtype']
# info[, 'Subtype_new'] <- factor(info[, 'Subtype_new'], levels = c(retain, 'Other subtype'))
# 
# 
# ## group assignment (manually) sequences groups
# info <- merge(info, iqtree_out[, c(1, as_pick + 1)], by='Strain_number', all=T)
# 
# colnames(info)[ncol(info)] <- 'small'
# # info[!info$small %in% (info$small %>% table_DF() %>% arrange(desc(Freq)))[1:18, 'x'], 'small'] <- NA
# rownames(info) <- info$Strain_number
# # info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# # info$small %>% unique()
# 
# # group <- (table_DF(info[, 'clades']) %>% arrange(desc(Freq)))[1:12, ]
# # info[!info$clades %in% group$x, 'clades'] <- NA
# group <- table_DF(info[, 'small']) %>% arrange(x %>% as.numeric())
# if(nrow(group) > 19) {
#   group <- group[1:18, ]
# }else(group <- group)
# rownames(group) <- group$x
# info[!info[, 'small'] %in% group$x, 'small'] <- NA
# # remain_g <- info[, 'small'] %>% unique()
# # group <- group[remain_g[!is.na(remain_g)], ]
# group[, 'edit'] <- paste('s_l', 1:nrow(group), sep = '')
# 
# for(k in seq_len(nrow(group))) {
#   info[info[, 'small'] %in% group[k, 'x'], 'small'] <- group[k, 3]
# }
# # info[!info$g %in% group$x, 'clade'] <- 'Minor'
# info[, 'small'] <- factor(info[, 'small'], levels = c(paste('s_l', 1:18, sep = ''), NA))
# 
# rownames(info) <- info$Strain_number %>% as.character()
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# # info$Strain_number <- factor(info$Strain_number, levels =rev(rank))
# 
# 
# 
# rownames(info) <- info$Strain_number %>% as.character()
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# info$Strain_number <- factor(info$Strain_number, levels =rev(rank))
# 
# 
# # info$clades <- factor(info$clades, levels = c(unique(info[!is.na(info$clades), 'clades']), NA))
# # info$data_source <- str_remove_all(info$Accesion_number, pattern = '_[0-9]+')
# # info[!(info$data_source %in% 'IRD'), 'data_source'] <- 'GISAID'
# 
# #change tree header to match the data
# new_header <- info$Strain_number %>% as.character()
# tree$tip.label <- new_header
# #
# 
# # phylo_raw <- ggtree(tree)
# # tree_list[[i]] <- phylo_raw
# phylo_raw <- tree_list[[args$segment]]
# clade <- info[, c('Strain_number', 'small')]
# colnames(clade)[2] <- 'Clade'
# phylo <- 
#   phylo_raw %<+% clade +
#   geom_tippoint(aes(color = Clade), size = 2) +
#   scale_color_manual(values=custom_color,
#                      na.value = "transparent")+
#   # scale_color_brewer(palette='Set3')+
#   labs(title = paste(args$segment, 'Global tree', (iqtree_out[, c(as_pick + 1)] %>% unique() %>% length()), sep = ''), 
#        color='Sampling tree\nlineages')+ 
#   theme(legend.position = "right", legend.text = element_text(size=20-10),
#         legend.title = element_text(size=22-10), title = element_text(size=24))  
# 
# png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment,'_iq_tree_from_sampling.png',sep = ''),
#     width=9600, height=5400, res = 600)
# phylo
# dev.off()
# 
# # cat(paste('Plotting tree\n', sep = ''))
# # #Heatmap from sequence information
# # Heatmap_1 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=continent, color=continent), linewidth=0.05)+
# #   labs(x="Continent", fill='Continent', color='Continent') + ylab(NULL)+
# #   scale_color_brewer(palette = 'Accent')+
# #   scale_fill_brewer(palette = 'Accent')+
# #   phylo_theme
# # 
# # Heatmap_2 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Clade_new, color=Clade_new), linewidth=0.05)+
# #   labs(x="GISAID\nClade", fill="GISAID\nClade", color="GISAID\nClade") + ylab(NULL)+
# #   # scale_color_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
# #   # scale_fill_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
# #   scale_color_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
# #   scale_fill_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
# #   phylo_theme
# # 
# # Heatmap_3 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Subtype_new, color=Subtype_new), linewidth=0.05)+
# #   labs(x="Serotype", fill="Serotype", color="Serotype") + ylab(NULL)+
# #   # scale_color_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
# #   # scale_fill_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
# #   scale_color_brewer(palette = 'Set2')+
# #   scale_fill_brewer(palette = 'Set2')+
# #   phylo_theme
# # 
# # 
# # # Heatmap_4 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Host_type, color=Host_type), linewidth=0.05)+
# # #   labs(x="Host", fill="Host", color="Host") + ylab(NULL)+
# # #   scale_color_brewer(palette = 'Set3')+
# # #   scale_fill_brewer(palette = 'Set3')+
# # #   phylo_theme
# # 
# # 
# # # Heatmap_5 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=order_new, color=order_new), linewidth=0.05)+
# # #   labs(x="Order", fill="Order", color="Order") + ylab(NULL)+
# # #   scale_color_manual(values = c('#FAA508','#647531','#314975','#6B3175','#EBF2F7'))+
# # #   scale_fill_manual(values = c('#FAA508','#647531','#314975','#6B3175','#EBF2F7'))+
# # #   phylo_theme
# # 
# # 
# # Heatmap_6 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=time_group, color=time_group), linewidth=0.05)+
# #   labs(x="Year", fill="Year", color="Year") + ylab(NULL)+
# #   scale_color_brewer(palette = 'Blues')+
# #   scale_fill_brewer(palette = 'Blues')+
# #   phylo_theme
# # 
# # test <- info[, c(1, 19)]
# # colnames(test)[2] <- 'color'
# # Heatmap_7 <- ggplot(test)+geom_tile(aes(x="",y=Strain_number, fill=color, color=color), linewidth=0.05)+
# #   labs(x=colnames(info)[19], fill="Grouping", color="Grouping") + ylab(NULL)+
# #   scale_color_brewer(palette = 'Set3')+
# #   scale_fill_brewer(palette = 'Set3')+
# #   phylo_theme+theme(axis.title.x = element_text(size = 10))
# # 
# # 
# # group_Heatmaps <- list()
# # for(i in colnames(info)[20:ncol(info)]) {
# #   A <- info[, c('Strain_number', i)]
# #   colnames(A)[2] <- 'test'
# #   group_Heatmaps[[i]] <-
# #     ggplot(A)+geom_tile(aes(x="",y=Strain_number, fill=test, color=test), linewidth=0.05)+
# #     labs(x=i, fill="d_0.05", color="d_0.05") + ylab(NULL)+
# #     scale_color_brewer(palette = 'Set3')+
# #     scale_fill_brewer(palette = 'Set3')+
# #     phylo_theme+theme(legend.position = "none", axis.title.x = element_text(size = 10))
# # }
# # 
# # 
# # # p <- (Heatmap_1 %>% insert_left(phylo, width = 12))
# # p <- (Heatmap_7 %>% insert_left(phylo, width = 12))
# # for(i in seq_along(group_Heatmaps)) {
# #   p <- p %>% insert_right(group_Heatmaps[[i]])
# # }
# # 
# # p <-as.ggplot(p %>% insert_right(Heatmap_1) %>% insert_right(Heatmap_3) %>% insert_right(Heatmap_2) %>%
# #                 insert_right(Heatmap_6) )+ #%>% insert_right(Heatmap_7)
# #   ggtitle(paste(args$segment, "segment phylogenetic tree (iqtree, lineages from sampling tree)"))+
# #   theme(plot.title = element_text(face = "bold",size = 25 ,hjust = 0.5))
# # 
# # 
# # png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment,'_iq_tree_from_samll.png',sep = ''),
# #     width=9600, height=5400, res = 600)
# # p
# # dev.off()
# # saveRDS(tree_list, '/home/eric/Analysis/aiv/merge/0307/iq_tree/internal_NA_tree.RData')
# # 
