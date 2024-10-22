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



args$segment <- 'NS'
args$iqtree <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0924/iq/NS_5p_sampling.nexus'
args$iq_grouping <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/NS_lineage_last.RData'
args$output_prefix <- '_as_sampling'
args$out_dir <- '~/Analysis/aiv/merge/0307/distance/sampling_0924/'


# merge result ------------------------------------------------------------
# dir.create(args$out_dir)
# 
# 
# tree <- treeio::read.beast(args$iqtree)
# 
iqtree_out <- readRDS(args$iq_grouping)
 
# while(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {
#   iqtree_out <- iqtree_out[, -ncol(iqtree_out)]
# }
# # if(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {iqtree_out <- iqtree_out[, -ncol(iqtree_out)]}
# 
# 
# 
# A <- readRDS(paste(args$out_dir, args$segment, args$output_prefix, '_agreement_result.RData', sep = ''))
# u <- A$agreement_table
# 
# 
# plot_data <- u %>% do.call(what=rbind) %>% as.data.frame()
# # plot_data <- u %>% tibble::rownames_to_column('dis') %>% tidyr::gather(key, value, 2:(ncol(u)+1))
# # plot_data$dis <- str_remove(string = plot_data$dis, pattern = 'd_')
# # g <- plot_data[plot_data$group %in% 'u', ] %>% nrow()
# plot_data[plot_data$group %in% "u", 'group'] <- "Agreement"
# plot_data[plot_data$group %in% "u-same", 'group'] <- "Agreement-same"
# plot_data[plot_data$group %in% "u-distinct", 'group'] <- "Agreement-diff"
# colnames(plot_data)[2:3] <- c('key', 'value')
# plot_data[plot_data$key %in% c('Agreement', 'Agreement-same', 'Agreement-diff'), 'value'] <-
#   plot_data[plot_data$key %in% c('Agreement', 'Agreement-same', 'Agreement-diff'), 'value']*100
# plot_data$dis <- str_remove(plot_data$dis, pattern = 'd_') %>% as.numeric()
# 
# plot_data <- rbind(plot_data,
#                    data.frame(dis=plot_data$dis %>% unique(),
#                               key='Disagreement',
#                               value=max(plot_data[plot_data$key %in% c('Agreement'), 'value']) -plot_data[plot_data$key %in% c('Agreement'), 'value']))
# plot_data$dis <- factor(plot_data$dis, levels = unique(plot_data$dis))
# # plot_data[plot_data$key %in% c('Agreement'), 'value'] <-
# # plot_data[plot_data$key %in% "Agreement", 'key'] <- "Disagreement"
# 
# dis_bs <- data.frame(mean=plot_data %>% filter(key == c('bs')) %>% group_split(dis) %>%
#                        lapply(FUN = function(x){median(x$value) %>% round(digits = 1)}) %>% unlist(),
#                      x=(plot_data %>% filter(key == c('bs')))$dis %>% unique())
# iq_groups <- plot_data %>% filter(key == c('iq_groups'))
# 
# as <-
#   plot_data %>% filter(key %in% c('Agreement', 'Agreement-same', 'iq_groups')) %>%
#   reshape(idvar = "dis", timevar = "key", direction = "wide")
# colnames(as) <- c('dis', 'AS', 'AS_same', 'iq_groups')
# 
# if(args$segment %in% c(paste('N', 1:9, sep = ''))){
#   as <- as %>% mutate(diff1=c(AS_same[-1], AS_same[nrow(as)])-AS_same, 
#                       diff2=c(0, diff1[-nrow(as)])) %>% filter(AS > 80, iq_groups >= 3, iq_groups <= 20)
# }else{
#   as <- as %>% mutate(diff1=c(AS_same[-1], AS_same[nrow(as)])-AS_same, 
#                       diff2=c(0, diff1[-nrow(as)])) %>% filter(AS > 75, iq_groups > 3, iq_groups <= 20)
# }
# 
# 
# # elbow point exist
# # if(any(as$diff1 >= 6 & as$diff2 <2)) {as_pick <- (as[as$diff1 >=6 & as$diff2 <2, ][, 1] %>% as.character() %>% rev())[1]}
# if(any(as$diff1 >= 6 & as$diff2 <2)) {
#   as_pick <- (as[as$diff1 >=6 & as$diff2 <2, ])
#   as_pick <- as_pick[which.max(as_pick$AS), 1] %>% unlist() %>% as.character()
# }
# # no elbow point
# if(isFALSE(any(as$diff1 >= 6 & as$diff2 <2))) {as_pick <- as[which.max(as$AS_same+as$AS), 1] %>% as.character()}
# 
# 
# rec <- grep(plot_data$dis %>% unique(), pattern=paste(as_pick, '$', sep = ''))
# as_pick <- grep(colnames(iqtree_out), pattern=paste(as_pick, '$', sep = ''))-1



# topology assignment -----------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'
meta$p <- paste(meta$Isolate_Id, meta$Location, meta$Collection_Date, sep = '_')
meta$Strain_number <- as.character(meta$Strain_number)

# iqtree <- paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', args$segment, '_aligned_iqtree.nexus', sep = '')
iqtree <- paste('~/Analysis/aiv/merge/0307/sampling_subset/1015/iq/', args$segment, '_aa_sampling.nexus', sep = '')
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

# info_sampling_lineage <- merge(info, iqtree_out[, c(1, as_pick + 1)], by='Strain_number', all.x = T)
info_sampling_lineage <- merge(info, iqtree_out, by='Strain_number', all.x = T)
# setdiff(iqtree_out[, c(1)], info$Strain_number) %>% length()
colnames(info_sampling_lineage)[11] <- 'sample_lineage'


# node_leafe --------------------------------------------------------------
pd_list <- extract_leaf_nodes(tree)

# pd_list <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', args$segment,'_iq_patristic_dist_matrix.RData', sep = ''))

leafe <- pd_list
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
  A <- test[[i]] %>% filter(Prop >= 50) %>% top_n(1, tips)

  # A <- test[[i]] %>% filter(Prop >= 50)
  # take_node <- A %>% arrange(tips)
  # leafe_subset <- leafe_all[leafe_all$node %in% A$node, ]
  # nodes <- c()
  # x <- 0
  # while (nrow(leafe_subset) > 0) {
  #   x <- x+1
  #   seed <- leafe_subset[leafe_subset$node %in% take_node[1, 1], 'epi'] %>% unique() #
  #   seed_node <- leafe_subset[leafe_subset$epi %in% seed, 'node'] %>% unique()
  #   seed_node_epi <- leafe_subset[leafe_subset$node %in% seed_node, 'epi'] %>% unique()
  #   node <- leafe_subset[leafe_subset$epi %in% seed_node_epi, 'node'] %>% table_DF() %>% arrange(desc(Freq))
  #   
  #   print(node[1, ])
  #   take_node <- take_node[!take_node$node %in% (node$x %>% unlist), ]
  #   leafe_subset <- leafe_subset[!leafe_subset$node %in% (node$x %>% unlist), ]
  #   
  #   nodes[x] <- node$x[1]
  #   print(leafe_subset %>% nrow())
  # }
  
  
  node_list[[i]] <- data.frame(sample=i, global=A$node, tips=A$tips)
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



# node up and down relation -----------------------------------------------
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


# remove minor lineages ---------------------------------------------------
unwant <- all_assign$assigned %>% table_DF() %>% filter(Freq <3) %>% select(x) %>% unlist()

all_assign[all_assign$assigned %in% unwant, 'assigned'] <- NA
for(k in is.na(all_assign[, 'assigned']) %>% which()) {
  if(k==1) {all_assign[k, 'assigned'] <- all_assign[k+1, 'assigned']}
  else(all_assign[k, 'assigned'] <- all_assign[k-1, 'assigned'])
}



all_assign$correct_out <- ''
all_assign$correct_out[1] <- 1
x <- 1
for(k in 2:nrow(all_assign)) {
  if(all_assign[k-1, 'assigned'] == all_assign[k, 'assigned']) {
    all_assign[k, 'correct_out'] <- x
  }else{
    x <- x+1
    all_assign[k, 'correct_out'] <- x
  }
}


remain <- c()
for(i in all_assign$assigned %>% unique()) {
  remain[i] <- all_assign[all_assign$assigned == i, 'correct_out'] %>% table_DF %>% 
    top_n(1, Freq) %>% select(x) %>% unlist()
}


write.csv((all_assign[all_assign$correct_out %in% remain, 'correct_out'] %>% table_DF() %>% 
             select(Freq) %>% sum())/nrow(info), 
          file = paste('~/Analysis/aiv/merge/0307/lineages/', args$segment, '_unassigned_ratio.csv', sep = ''), 
          row.names = F, quote = F)


all_assign[!all_assign$correct_out %in% remain, c('assigned', 'correct_out')] <- NA

for(k in is.na(all_assign[, 'correct_out']) %>% which()) {
  if(k==1) {all_assign[k, 'correct_out'] <- all_assign[k+1, 'correct_out']}
  else(all_assign[k, 'correct_out'] <- all_assign[k-1, 'correct_out'])
}


all_assign$correct_out1 <- ''
all_assign$correct_out1[1] <- 1
x <- 1
for(k in 2:nrow(all_assign)) {
  if(all_assign[k-1, 'correct_out'] == all_assign[k, 'correct_out']) {
    all_assign[k, 'correct_out1'] <- x
  }else{
    x <- x+1
    all_assign[k, 'correct_out1'] <- x
  }
}

all_assign$correct_out1 <- paste('L', all_assign$correct_out1, sep = '')


for(k in is.na(all_assign[, 'assigned']) %>% which()) {
  if(k==1) {all_assign[k, 'assigned'] <- all_assign[k+1, 'assigned']}
  else(all_assign[k, 'assigned'] <- all_assign[k-1, 'assigned'])
}

write.csv(all_assign, file = paste('~/Analysis/aiv/merge/0307/lineages/', args$segment, '_aa_ident_lineage.csv', sep = ''), row.names = F)



# assign identical aa back ------------------------------------------------
a <- read.fasta(paste('/home/eric/Analysis/aiv/merge/0307/ORF_filter_outlier/gisaid_IRD_merged_', args$segment,'_aligned_fill_trimed_remove_outlier.fa', sep = ''), as.string = T) %>% 
  do.call(what=rbind) %>% as.data.frame() 
colnames(a) <- 'sequence'
a <- a %>% mutate(V2=header_cleaning(rownames(a), pattern = '\\|') %>% select(V2) %>% unlist(), 
                  sequence=as.character(sequence)) %>% rename('Strain_number'='V2')


rownames(a) <- 1:nrow(a)
a <- a %>% mutate(sequence=str_replace_all(sequence, pattern = 'u|p|x|e|i|o|l', replacement = 'n'))

a$aa <-
  Biostrings::DNAStringSet(a$sequence %>% str_remove_all(pattern = '-')) %>%
  Biostrings::translate(if.fuzzy.codon = 'X') %>% as.character()

seg <- merge(all_assign[, c('Strain_number', 'correct_out1')], a, by = 'Strain_number', all.y = T) %>% 
  group_split(aa) %>% as.list() %>% lapply(FUN = function(x){as.data.frame(x[, 1:2])})

multiple <- seg[which(lapply(seg, nrow) != 1)]

for(i in seq_along(multiple)) {
  multiple[[i]][, 2] <- (multiple[[i]] %>% arrange(correct_out1))[1, 2]
}

out <- rbind(seg[which(lapply(seg, nrow) == 1)] %>% do.call(what=rbind) %>% as.data.frame(), 
             multiple %>% do.call(what=rbind) %>% as.data.frame())
colnames(out)[2] <- 'lineages'
write.csv(out, file = paste('~/Analysis/aiv/merge/0307/lineages/', args$segment, '_lineage.csv', sep = ''), row.names = F)




meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame() %>% 
  mutate(Strain_number=as.character(Strain_number))

clade_table <- table_DF(meta[grepl(meta$Subtype, pattern='H5'), 'Clade'])
clade_table[, 'Clade'] <- 'GsGD-others'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4'), 'Clade'] <- '2.3.4.4'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4b'), 'Clade'] <- '2.3.4.4b'
clade_table[grep(clade_table$x, pattern = 'nonGsGD'), 'Clade'] <- 'nonGsGD'
clade_table[grep(clade_table$x, pattern = 'IRD_unlabel'), 'Clade'] <- 'IRD_unlabel'


i <- args$segment
cat(paste('Loading tree\n', sep = ''))
# tree <- ape::read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i, '_aligned_iqtree.nexus', sep = ''))
tree <- ape::read.nexus(paste('~/Analysis/aiv/merge/0307/sampling_subset/1015/iq/', i,'_aa_sampling.nexus', sep = ''))

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
info[info$Subtype %in% c('H5N1', 'H5N2', 'H5N6', 'H5N8', 'H7N9', 'H9N2'), 'Subtype_new'] <-
  info[info$Subtype %in% c('H5N1', 'H5N2', 'H5N6', 'H5N8', 'H7N9', 'H9N2'), 'Subtype']


info[, 'Subtype_new'] <- factor(info[, 'Subtype_new'], levels = c('H5N1', 'H5N2', 'H5N6', 'H5N8', 'H7N9', 'H9N2', 'Other subtype'))


## group assignment (manually) sequences groups

# info <- merge(info, out, by='Strain_number', all=T)
info <- merge(info, all_assign[, c('Strain_number', 'sample_lineage', 'correct_out1')], by='Strain_number', all=T)
rownames(info) <- info$Strain_number
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
info[, 19] <- factor(info[, 19], 
                     levels = table(info$sample_lineage, info$correct_out1) %>% as.data.frame(stringsAsFactors=F) %>% 
                       filter(Freq !=0) %>% mutate(Var1=as.numeric(Var1)) %>% arrange(Var1) %>% select(Var2) %>% 
                       unlist() %>% unique())

rownames(info) <- info$Strain_number %>% as.character()
info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))

table(info$sample_lineage, info$correct_out1)

info$Strain_number <- factor(info$Strain_number, levels =rev(rank))
# info$data_source <- str_remove_all(info$Accesion_number, pattern = '_[0-9]+')
# info[!(info$data_source %in% 'IRD'), 'data_source'] <- 'GISAID'

#change tree header to match the data
new_header <- info$Strain_number %>% as.character()
tree$tip.label <- new_header


# tree2 <- groupClade(tree, .node = c((unique(info$assigned) %>% as.numeric())+Ntip(tree)))
# Check that strain numbers match the tree tip labels


# Call groupOTU with the correct focus
# tree2 <- groupOTU(tree, .node = group)

# phylo_raw <- 
# ggtree(tree2, aes(color=group))+
# scale_color_manual(values=c(custom_color), na.value = "transparent")

phylo_raw <- ggtree(tree)
# tree_list[[i]] <- phylo_raw
# phylo_raw <- tree_list[[i]]
clade <- info[, c(1, ncol(info))]
colnames(clade)[2] <- 'lineages'
phylo <-
  phylo_raw %<+% clade +
  geom_tippoint(aes(color = lineages), size = 2) +
  scale_color_manual(values=custom_color, na.value = "transparent")+
  # scale_color_brewer(palette='Set3')+
  labs(color='Assigned\nlineage')+# scale_colo(palette = 'Accent')+
  theme(legend.position = "right", legend.text = element_text(size=20-8),
        legend.title = element_text(size=22-8))  # Display the legend on the right


cat(paste('Plotting tree', sep = ''))
#Heatmap from sequence information
Heatmap_1 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=continent, color=continent), linewidth=0.05)+
  labs(x="Continent", fill='Continent', color='Continent') + ylab(NULL)+
  scale_color_brewer(palette = 'Accent')+
  scale_fill_brewer(palette = 'Accent')+
  guides(fill = guide_legend(nrow = 1, direction = "vertical"))+
  phylo_theme+theme(legend.position = "bottom") 

Heatmap_2 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Clade_new, color=Clade_new), linewidth=0.05)+
  labs(x="GISAID\nClade", fill="GISAID\nClade", color="GISAID\nClade") + ylab(NULL)+
  # scale_color_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  # scale_fill_manual(values = c('#4FC0E8','#7EEDA7','#5AA977','#EBF2F7', '#B63B3D'))+ #'#DF9E01','#DF3A01','#7B4322'
  scale_color_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  scale_fill_manual(values = c("2.3.4.4"='#4FC0E8', "2.3.4.4b"='#7EEDA7','GsGD-others'='#5AA977','nonGsGD'='grey70', 'IRD_unlabel'= '#EBF2F7'))+ #'#DF9E01','#DF3A01','#7B4322'
  guides(fill = guide_legend(nrow = 1, direction = "vertical"))+
  phylo_theme+theme(legend.position = "bottom") 

Heatmap_3 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=Subtype_new, color=Subtype_new), linewidth=0.05)+
  labs(x="Subtype", fill="Subtype", color="Subtype") + ylab(NULL)+
  scale_color_manual(values = c('H5N1'='#FFD97D','H5N2'='#FF8F54','H5N6'='#FF5356','H5N8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
  scale_fill_manual(values = c('H5N1'='#FFD97D','H5N2'='#FF8F54','H5N6'='#FF5356','H5N8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
  # scale_color_brewer(palette = 'Set2')+
  # scale_fill_brewer(palette = 'Set2')+
  phylo_theme


Heatmap_6 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=time_group, color=time_group), linewidth=0.05)+
  labs(x="Year", fill="Year", color="Year") + ylab(NULL)+
  scale_color_brewer(palette = 'Blues')+
  scale_fill_brewer(palette = 'Blues')+
  guides(fill = guide_legend(nrow = 1, direction = "vertical"))+
  phylo_theme+theme(legend.position = "bottom") 

# Load patchwork library
library(patchwork)

combined_plot <-
  (phylo + Heatmap_1 + Heatmap_3 + Heatmap_2 + Heatmap_6) +
  plot_layout(guides = "collect", widths = c(0.6, 0.1, 0.1, 0.1, 0.1)) & 
  theme(legend.position = "bottom",
        legend.direction = 'horizontal',
        legend.box = 'vertical')

# Add title and display
png(paste('/home/eric/Analysis/aiv/merge/0307/result/', args$segment, '_iqtree_aa_identical_lineage.png', sep = ''), width=9600, height=5400, res = 600)
combined_plot + labs(title = paste(args$segment, "segment phylogenetic tree (iqtree)")) +
  theme(plot.title = element_text(face = "bold", size = 25, hjust=2.5))
dev.off()


