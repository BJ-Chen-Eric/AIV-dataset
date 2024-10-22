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

parser$add_argument("-op", "--output_prefix", type='character', default='agreement',
                    help="output file prefix")

parser$add_argument("-o", "--out_dir", type='character', default=TRUE,
                    help="Generated distance matrix store direction")
args <- parser$parse_args()


# args$segment <- 'MP'
# args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/MP_5p_sampling.nexus'
# args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/MP_iq_MPD_groups.RData'
# args$ft_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/MP_ft_MPD_groups.RData'
# args$tree_dir <- '~/Analysis/aiv/merge/0307/sampling_subset/0924/iq/'
# args$tree_prefix <- '_5p_sampling.nexus'
# args$out_dir <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/'



args$segment <- 'PA'
args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/PA_5p_sampling.nexus'
args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/PA_iq_MPD_groups.RData'
args$output_prefix <- '_as_sampling'
args$out_dir <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/'


# merge result ------------------------------------------------------------
dir.create(args$out_dir)


tree <- treeio::read.beast(args$iqtree)


iqtree_out <- readRDS(args$iq_grouping)[[2]]
rank <- colnames(iqtree_out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V2=as.numeric(V2), name=paste(V1, V2, sep = '_')) %>%
  arrange(V2) %>% select(name) %>% unlist()

iqtree_out <- iqtree_out[, c('Strain_number', rank)]
if(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {iqtree_out <- iqtree_out[, -ncol(iqtree_out)]}
if(colnames(iqtree_out)[2] %in% 'd_0.0075') {
  iqtree_out <- iqtree_out %>% select(-c('d_0.0075'))
}


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

as <- as %>% mutate(diff=c(AS_same[-1], AS_same[nrow(as)])-AS_same) %>% filter(AS > 90, iq_groups > 3)
# elbow point exist
if(any(as$diff >= 6)) {as_pick <- (as[as$diff >=6 , ][, 1] %>% as.character() %>% rev())[1]}
# no elbow point
if(isFALSE(any(as$diff >= 6))) {as_pick <- as[which.max(as$AS_same), 1] %>% as.character()}


rec <- grep(plot_data$dis %>% unique(), pattern=paste(as_pick, '$', sep = ''))
as_pick <- grep(colnames(iqtree_out), pattern=paste(as_pick, '$', sep = ''))-1




# topology assignment -----------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'
meta$p <- paste(meta$Isolate_Id, meta$Location, meta$Collection_Date, sep = '_')
meta$Strain_number <- as.character(meta$Strain_number)

iqtree <- paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', args$segment, '_aligned_iqtree.nexus', sep = '')
tree <- ape::read.nexus(iqtree)

tip_name <- tree$tip.label %>% str_remove_all(pattern = "'") # 286980

info <- data.frame(Strain_number=tip_name %>% str_extract(pattern = '_[0-9]+_H|\\|[0-9]+\\|H') %>% str_remove(pattern = '_|\\|') %>%
                     str_remove(pattern = '_H|\\|H'))

rank <- info %>% unlist()

info <- merge(info, meta[, c('Strain_number', 'new_header')],
              by = 'Strain_number', all.x = T)
# head(info) %>% print()
info <- header_cleaning(info$new_header, pattern = '\\|')


colnames(info) <- c("Accesion_number","Strain_number","Subtype", 'Clade', "Segment", 'Collection_date'
                    ,"Location","Host", "Header_Host", "Host_type")

info_sampling_lineage <- merge(info, iqtree_out[, c(1, as_pick + 1)], by='Strain_number', all=T)
colnames(info_sampling_lineage)[11] <- 'sample_lineage'


# node_leafe --------------------------------------------------------------
pd_list <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', args$segment,'_iq_patristic_dist_matrix.RData', sep = ''))

out_list <- list()


leafe <- pd_list[[2]]#[-c(1:3)]
names(leafe) <- 1:length(leafe)

leafe_all <- leafe %>% lapply(FUN = function(x){data.frame(epi=x)}) %>% do.call(what=rbind) %>%
  as.data.frame() %>% tibble::rownames_to_column('node') %>% mutate(node=str_remove(node, patter='\\.[0-9]+') %>% as.numeric())


leafe_all$Strain_number <- strain_number_extract(leafe_all$epi)


leafe_all <- merge(info_sampling_lineage[!is.na(info_sampling_lineage$sample_lineage),
                                         c("Strain_number", 'sample_lineage')], leafe_all, by = 'Strain_number', all=T)

all_leafe <- leafe_all[, 'node'] %>% table_DF()

leafe_all_list <- leafe_all %>% group_split(node) %>% as.list() %>% lapply(as.data.frame)

names(leafe_all_list) <- lapply(leafe_all_list, function(x){x[1, 3]})


leafe_all_list <- leafe_all_list[lapply(leafe_all_list, function(x){any(!x[, 2] %>% is.na())}) %>% unlist() %>% which()]

leaves <- lapply(leafe_all_list, function(x){nrow(x)}) %>% do.call(what=rbind) %>% as.data.frame() %>%
  tibble::rownames_to_column('node') %>% rename('tips'='V1')

sample_leaves <- lapply(leafe_all_list, function(x){((!(x[, 2]) %>% is.na()) %>% sum())}) %>%
  do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('node') %>% rename('Sample'='V1')


test <- lapply(leafe_all_list, function(x){(x[, 2] %>% table_DF(prop = T) %>% mutate(Prop=Prop*100))})# %>% do.call(what=rbind) %>% as.data.frame()
test <- lapply(test, function(x){x[which.max(x$Freq), ]}) %>% do.call(what=rbind) %>% as.data.frame() %>%
  tibble::rownames_to_column('node')
# test <- test %>% do.call(what=rbind) %>% as.data.frame() %>%
#   tibble::rownames_to_column('node')
test <- merge(test, leaves, by = 'node', all=T)
test <- merge(test, sample_leaves, by = 'node', all=T)



test <- test %>% group_split(x) %>% as.list() %>% lapply(as.data.frame)


node_list <- list()
for(i in seq_along(test)) {
  A <- test[[i]] %>% filter(Prop >= 50)

  take_node <- A %>% arrange(tips)
  leafe_subset <- leafe_all[leafe_all$node %in% A$node, ]

  nodes <- c()
  x <- 0
  while (nrow(leafe_subset) > 0) {
    x <- x+1
    seed <- leafe_subset[leafe_subset$node %in% take_node[1, 1], 'epi'] %>% unique() #
    seed_node <- leafe_subset[leafe_subset$epi %in% seed, 'node'] %>% unique()
    seed_node_epi <- leafe_subset[leafe_subset$node %in% seed_node, 'epi'] %>% unique()
    node <- leafe_subset[leafe_subset$epi %in% seed_node_epi, 'node'] %>% table_DF() %>% arrange(desc(Freq))

    print(node[1, ])
    take_node <- take_node[!take_node$node %in% (node$x %>% unlist), ]
    leafe_subset <- leafe_subset[!leafe_subset$node %in% (node$x %>% unlist), ]

    nodes[x] <- node$x[1]
    print(leafe_subset %>% nrow())
  }


  node_list[[i]] <- data.frame(sample=i, global=nodes)
  # result[[i]] <- A[A$node %in% max(A[A$Freq %in% max(A$Freq), ][, 'node']),  ]
}

result <- node_list %>% do.call(what=rbind)


assigned_lineage <- list()
for(i in 1:nrow(result)) {
  A <- leafe_all[leafe_all$node %in% result[i, 'global'], ]
  A$sample <- result[i, 'sample']

  assigned_lineage[[result[i, 'global']]] <- leafe_all[leafe_all$node %in% result[i, 'global'], ]
}

assigned_lineage <- assigned_lineage %>% do.call(what=rbind) %>% as.data.frame()

test1 <- assigned_lineage

duplicate <- table_DF(test1$epi) %>% filter(Freq != 1) %>% select(x) %>% unlist()
duplicate_node <- test1[test1$epi %in% duplicate, 'node'] %>% unique()

complete <- test1[!test1$node %in% duplicate_node, ] %>% unique()
complete$assigned <- complete$node

take_node <- leaves[leaves$node %in% duplicate_node, ] %>% arrange(tips)
leafe_subset <- leafe_all[leafe_all$node %in% duplicate_node, ]

nodes <- list()
x <- 0
while (nrow(leafe_subset) > 0) {
  seed <- leafe_subset[leafe_subset$node %in% take_node[1, 1], 'epi'] %>% unique() #
  seed_node <- leafe_subset[leafe_subset$epi %in% seed, 'node'] %>% unique()
  seed_node_epi <- leafe_subset[leafe_subset$node %in% seed_node, 'epi'] %>% unique()
  node <- leafe_subset[leafe_subset$epi %in% seed_node_epi, 'node'] %>% table_DF() %>% arrange(Freq)

  for(i in 1:nrow(node)) {
    x <- x+1
    A <- leafe_subset[leafe_subset$node %in% node[i, 'x'], ]
    A$assigned <- node[i, 'x']
    nodes[[x]] <- A
    leafe_subset <- leafe_subset[!leafe_subset$epi %in% A$epi, ]
  }

  # test2 <- complete[complete$epi %in% seed_node_epi, ]
  take_node <- take_node[!take_node$node %in% (node$x %>% unlist), ]
  print(leafe_subset %>% nrow())
}
nodes <- nodes %>% do.call(what=rbind) %>% as.data.frame()

all_assign <- rbind(complete, nodes)
A <- leafe_all[leafe_all$Strain_number %in% setdiff(info$Strain_number, all_assign$Strain_number), ] %>%
  arrange(desc(node)) %>% distinct(epi, .keep_all = T) %>% mutate(assigned=NA)
all_assign <- rbind(all_assign, A)
rownames(all_assign) <- all_assign$Strain_number
all_assign <- all_assign[rank, ]

unwant <- all_assign$assigned %>% table_DF() %>% filter(Freq <=14) %>% select(x) %>% unlist()
all_assign[all_assign$assigned %in% unwant, 'assigned'] <- NA

for(k in is.na(all_assign[, 'assigned']) %>% which()) {
  if(k==1) {all_assign[k, 'assigned'] <- all_assign[k+1, 'assigned']}
  else(all_assign[k, 'assigned'] <- all_assign[k-1, 'assigned'])
}


all_assign$correct <- ''
all_assign$correct[1] <- 1
x <- 1
for(k in 2:nrow(all_assign)) {
  if(all_assign[k-1, 'assigned'] == all_assign[k, 'assigned']) {
    all_assign[k, 'correct'] <- x
  }else{
    x <- x+1
    all_assign[k, 'correct'] <- x
  }
}

unwant <- all_assign$correct %>% table_DF() %>% filter(Freq <3) %>% select(x) %>% unlist()
all_assign[all_assign$correct %in% unwant, 'correct'] <- NA

for(k in is.na(all_assign[, 'correct']) %>% which()) {
  if(k==1) {all_assign[k, 'correct'] <- all_assign[k+1, 'correct']}
  else(all_assign[k, 'correct'] <- all_assign[k-1, 'correct'])
}


all_assign$correct_out <- ''
all_assign$correct_out[1] <- 1
x <- 1
for(k in 2:nrow(all_assign)) {
  if(all_assign[k-1, 'correct'] == all_assign[k, 'correct']) {
    all_assign[k, 'correct_out'] <- x
  }else{
    x <- x+1
    all_assign[k, 'correct_out'] <- x
  }
}


all_assign$correct_out <- paste('L', all_assign$correct_out, sep = '')
out <- all_assign[, c('Strain_number', 'correct_out')]
colnames(out)[2] <- 'lineages'
write.csv(all_assign[, c(1, 6)], file = paste('~/Analysis/aiv/merge/0307/lineages/', args$segment, '_lineage.csv', sep = ''), row.names = F)

# all_assign$index <- 1:nrow(all_assign)
# 
# test2 <- all_assign %>% group_split(correct) %>% as.list() %>% lapply(function(x){x$index}) %>%
#   lapply(function(x) {c(min(x), max(x))})
# names(test2) <- 1:length(test2)
# test2 <- test2 %>% do.call(what=rbind) %>% as.data.frame() %>% arrange(V1) %>% mutate(p=paste(V1, V2, sep = '_'))
# colnames(test2)[1:2] <- c('before_s', 'before_e')
# 
# 
# test3 <- all_assign %>% group_split(assigned) %>% as.list() %>% lapply(function(x){x$index}) %>%
#   lapply(function(x) {c(min(x), max(x))})
# names(test3) <- 1:length(test3)
# test3 <- test3 %>% do.call(what=rbind) %>% as.data.frame() %>% arrange(V1) %>% mutate(p=paste(V1, V2, sep = '_'))
# colnames(test3)[1:2] <- c('correct_s', 'correct_e')
# 
# test4 <- merge(test2, test3, by = 'p', all = T)


# vis ---------------------------------------------------------------------
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

segment <- args$segment
iqtree <- paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', args$segment, '_aligned_iqtree.nexus', sep = '')
output_prefix <- 'lineage_0919'
out_dir <- '~/Analysis/aiv/merge/0307/result/'


cat(paste('Loading tree\n', sep = ''))
tree <- ape::read.nexus(iqtree)

tip_name <- tree$tip.label %>% str_remove_all(pattern = "'") # 286980

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
info <- merge(info, all_assign[, c('Strain_number', 'correct')], by='Strain_number', all=T)

colnames(info)[ncol(info)] <- 'Assigned'
info[!info$Assigned %in% (info$Assigned %>% table_DF() %>% arrange(desc(Freq)))[1:18, 'x'], 'Assigned'] <- NA
rownames(info) <- info$Strain_number
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# info$Assigned %>% unique()

# group <- (table_DF(info[, 'clades']) %>% arrange(desc(Freq)))[1:12, ]
# info[!info$clades %in% group$x, 'clades'] <- NA
group <- table_DF(info[, 'Assigned']) %>% arrange(x %>% as.numeric())
if(nrow(group) > 19) {
  group <- group[1:18, ]
}else(group <- group)
rownames(group) <- group$x
info[!info[, 'Assigned'] %in% group$x, 'Assigned'] <- NA
# remain_g <- info[, 'Assigned'] %>% unique()
# group <- group[remain_g[!is.na(remain_g)], ]
group[, 'edit'] <- paste('L', 1:nrow(group), sep = '')

for(k in seq_len(nrow(group))) {
  info[info[, 'Assigned'] %in% group[k, 'x'], 'Assigned'] <- group[k, 3]
}
# info[!info$g %in% group$x, 'clade'] <- 'Minor'
info[, 'Assigned'] <- factor(info[, 'Assigned'], levels = c(paste('L', 1:18, sep = ''), NA))




# sampling ----------------------------------------------------------------
info <- merge(info, iqtree_out[, c(1, as_pick + 1)], by='Strain_number', all=T)

colnames(info)[ncol(info)] <- 'small'
info[!info$small %in% (info$small %>% table_DF() %>% arrange(desc(Freq)))[1:18, 'x'], 'small'] <- NA
rownames(info) <- info$Strain_number
# info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# info$small %>% unique()

# group <- (table_DF(info[, 'clades']) %>% arrange(desc(Freq)))[1:12, ]
# info[!info$clades %in% group$x, 'clades'] <- NA
group <- table_DF(info[, 'small']) %>% arrange(x %>% as.numeric())
if(nrow(group) > 19) {
  group <- group[1:18, ]
}else(group <- group)
rownames(group) <- group$x
info[!info[, 'small'] %in% group$x, 'small'] <- NA
# remain_g <- info[, 'small'] %>% unique()
# group <- group[remain_g[!is.na(remain_g)], ]
group[, 'edit'] <- paste('s_l', 1:nrow(group), sep = '')

for(k in seq_len(nrow(group))) {
  info[info[, 'small'] %in% group[k, 'x'], 'small'] <- group[k, 3]
}
# info[!info$g %in% group$x, 'clade'] <- 'Minor'
info[, 'small'] <- factor(info[, 'small'], levels = c(paste('s_l', 1:18, sep = ''), NA))

rownames(info) <- info$Strain_number %>% as.character()
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# info$Strain_number <- factor(info$Strain_number, levels =rev(rank))




rownames(info) <- info$Strain_number %>% as.character()
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
# info$Strain_number <- factor(info$Strain_number, levels =rev(rank))

rownames(info) <- info$Strain_number %>% as.character()
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
info$Strain_number <- factor(info$Strain_number, levels =rev(rank))


# info$clades <- factor(info$clades, levels = c(unique(info[!is.na(info$clades), 'clades']), NA))
# info$data_source <- str_remove_all(info$Accesion_number, pattern = '_[0-9]+')
# info[!(info$data_source %in% 'IRD'), 'data_source'] <- 'GISAID'

#change tree header to match the data
new_header <- info$Strain_number %>% as.character()
tree$tip.label <- new_header


# phylo_raw <- ggtree(tree)
# tree_list[[i]] <- phylo_raw
phylo_raw <- tree_list[[args$segment]]
clade <- info[, c('Strain_number', 'small')]
colnames(clade)[2] <- 'small'
n <- length(clade$small %>% unique())

phylo <- phylo_raw %<+% clade +
  geom_tippoint(aes(color = small), size = 2) +
  scale_color_manual(values=custom_color,
                     na.value = "transparent")+
  # scale_color_brewer(palette='Set3')+
  labs(color='small\nlineage')+# scale_colo(palette = 'Accent')+
  labs(color='Sampling tree\nlineages')+# scale_colo(palette = 'Accent')+
  theme(legend.position = "right", legend.text = element_text(size=20-10),
        legend.title = element_text(size=22-10))  # Display the legend on the right


cat(paste('Plotting tree\n', sep = ''))
#Heatmap from sequence information
Heatmap_1 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=continent, color=continent), linewidth=0.05)+
  labs(x="Continent", fill='Continent', color='Continent') + ylab(NULL)+
  scale_color_brewer(palette = 'Accent')+
  scale_fill_brewer(palette = 'Accent')+
  phylo_theme

Heatmap_2 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Clade_new, color=Clade_new), linewidth=0.05)+
  labs(x="GISAID\nClade", fill="GISAID\nClade", color="GISAID\nClade") + ylab(NULL)+
  # scale_color_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  # scale_fill_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  scale_color_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  scale_fill_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  phylo_theme

Heatmap_3 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Subtype_new, color=Subtype_new), linewidth=0.05)+
  labs(x="Serotype", fill="Serotype", color="Serotype") + ylab(NULL)+
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

test <- info[, c(1, 18)]
colnames(test)[2] <- 'color'
Heatmap_7 <- ggplot(test)+geom_tile(aes(x="",y=Strain_number, fill=color, color=color), linewidth=0.05)+
  labs(x=colnames(info)[18], fill="Assigned", color="Assigned") + ylab(NULL)+
  scale_color_manual(values=custom_color, na.value = "transparent")+
  scale_fill_manual(values=custom_color, na.value = "transparent")+
  phylo_theme+theme(axis.title.x = element_text(size = 10))



# p <- (Heatmap_1 %>% insert_left(phylo, width = 12))
p <- (Heatmap_7 %>% insert_left(phylo, width = 12))

p <- as.ggplot(p %>% insert_right(Heatmap_1) %>% insert_right(Heatmap_3) %>% insert_right(Heatmap_2) %>%
                insert_right(Heatmap_6) )+ #%>% insert_right(Heatmap_7)
  ggtitle(paste(args$segment, "segment phylogenetic tree (iqtree, lineages from sampling tree)"))+
  theme(plot.title = element_text(face = "bold",size = 25 ,hjust = 0.5))


png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment, '_iq_topology_assignemnt.png', sep = ''), width=1600, height=900)
p
dev.off()


