.libPaths(c("/home/eric/R/x86_64-pc-linux-gnu-library/4.4/"))
source('~/R/aiv/function.R')
class <- list(Internal=seg_sub_level[c(1:3, 20, 30, 31)], HA=seg_sub_level[c(4:19)], N=seg_sub_level[c(21:29)])
bird <- fread(file = '~/Analysis/aiv/merge/0307/birdlist.csv') %>% as.data.frame()
out <- fread(file = '~/Analysis/aiv/merge/0307/outlier.csv') %>% as.data.frame() %>% 
  apply(MARGIN = 2, FUN = function(x){strain_number_extract(x)}) %>% as.vector() %>% na.omit()

meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame() %>% 
  filter(!Strain_number %in% out)
meta[is.na(meta$segment), 'segment'] <- 'NA'
meta$p <- paste(meta$Isolate_Id, meta$Location, meta$Collection_Date, sep = '_')


# mpd_thres <- list(
#   # PB2=0.25,
#   # PB1=0.175,
#   # PA=0.125,
#   # NP=0.2,
#   # MP=0.075,
#   # NS=0.075,
#   N1=0.19, #0.175,
#   N2=0.15, #0.225,
#   N3=0.19,#0.175,
#   N4=0.08,#0.075,
#   N5=0.06,#0.075,
#   N6=0.13,#0.1,
#   N7=0.07,#0.075,
#   N8=0.08,#0.075,
#   N9=0.04 #0.05
# )# %>% do.call(what=rbind)  %>% as.data.frame() %>% tibble::rownames_to_column('seg') %>% rename('thres'='V1')

segment_groups <- list()
for(i in c(class$Internal)) {
  
  # # iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group.RData', sep = ''))[[2]]
  # iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/NA_1009/', i, '_iq_MPD_groups.RData', sep = ''))[[2]]
  # 
  # colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
  # iq_grouping <- iq_grouping[, colN]
  # 
  # tree <- read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i,'_aligned_iqtree.nexus', sep = ''))
  # tip_name <- tree$tip.label %>% str_remove_all(pattern = "'") # 286980
  # info <- data.frame(Strain_number=tip_name %>% str_extract(pattern = '_[0-9]+_H|\\|[0-9]+\\|H') %>% str_remove(pattern = '_|\\|') %>%
  #                      str_remove(pattern = '_H|\\|H'))
  # rank <- info %>% unlist()
  # 
  # rownames(iq_grouping) <- iq_grouping$Strain_number
  # iq_grouping <- iq_grouping[rank, ]
  # colnames(iq_grouping)[2] <- 'lineages'
  # 
  # group <- table_DF(iq_grouping[, 2]) #%>% filter(Freq > nrow(iq_grouping)*0.01)
  # # iq_grouping[!iq_grouping$clades %in% group$x, 'clades'] <- NA
  # # for(j in 2:ncol(iq_grouping)) {
  # #   for(k in is.na(iq_grouping[, j]) %>% which()) {
  # #     if(k==1) {iq_grouping[k, j] <- iq_grouping[k+1, j]}
  # #     else(iq_grouping[k, j] <- iq_grouping[k-1, j])
  # #   }
  # # }
  # # segment_groups[[i]] <- group
  # rownames(group) <- group$x
  # group <- group[iq_grouping[,2] %>% unique(), ]
  # group[, 'edit'] <- paste('L', 1:nrow(group), sep = '')
  # 
  # 
  # for(k in seq_len(nrow(group))) {
  #   iq_grouping[iq_grouping[, 2] %in% group[k, 'x'], 2] <- group[k, 3]
  # }
  # 
  # if(i %in% c(class$N)) {iq_grouping$lineages <- paste(i, iq_grouping$lineages, sep = '')}
  # write.csv(iq_grouping, file = paste('~/Analysis/aiv/merge/0307/lineages/', i, '_lineage.csv', sep = ''), quote = F, row.names = F)
  
  iq_grouping <- fread(paste('~/Analysis/aiv/merge/0307/lineages/', i, '_lineage.csv', sep = '')) %>% as.data.frame()
  segment_groups[[i]] <- iq_grouping
  # meta_geno <- merge(meta_geno, iq_grouping, by = 'Strain_number', all = T)
  # colnames(meta_geno)[ncol(meta_geno)] <- paste(i, '_group', sep = '')
}

meta$p <- paste(meta$Isolate_Id, meta$Location, meta$Collection_Date, sep = '_')
eight <- list()
# group_out
for(i in c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'))  {
  if((i %in% c('NA')))  {
    sub_meta <- meta[meta$segment %in% i, ]
    A <- data.frame(Strain_number=sub_meta$Strain_number, lineages=sub_meta$sub)
    A <- merge(A, sub_meta %>% mutate(Strain_number=as.character(Strain_number)), by = 'Strain_number', all.x = T) %>% 
      distinct(Strain_number, .keep_all = T)
  } 
  if((i %in% c('HA')))  {
    sub_meta <- meta[meta$segment %in% i, ]
    A <- data.frame(Strain_number=sub_meta$Strain_number, lineages=sub_meta$sub)
    A <- merge(A, sub_meta %>% mutate(Strain_number=as.character(Strain_number)), by = 'Strain_number', all.x = T) %>% 
      distinct(Strain_number, .keep_all = T)
  } 
  if(!(i %in% c('HA', 'NA'))) {
    A <- segment_groups[[i]]
    sub_meta <- meta[meta$segment %in% i, ]
    A <- merge(A, sub_meta, by = 'Strain_number', all.x = T) %>% 
      distinct(Strain_number, .keep_all = T) #%>% mutate(sub_clade=paste(sub, lineages, sep = '_'))
  }
  # A <- A %>% arrange(desc(Collection_Date)) %>% distinct(Isolate_Id, .keep_all = T) # remove duplicate epi
  eight[[i]] <- A %>% select(-c('seq'))
}



complete_virus <- Reduce(intersect, lapply(eight, function(x){x$Isolate_Id}))

more8_virus <- intersect(table_DF(meta$Isolate_Id) %>% filter(Freq>8) %>% select(x) %>% unlist, 
                            meta$Isolate_Id %>% unique())
full <- meta$Isolate_Id %>% unique()
genotype <- data.frame(Isolate_Id=full)
for(i in names(eight))  {
  A <- eight[[i]] 
  colnames(A)[2] <- i
  genotype <- merge(genotype, A[, c('Isolate_Id', i)], by = 'Isolate_Id', all.x=T) %>% distinct(Isolate_Id, .keep_all = T)
  genotype[, i] <- replace(genotype[, i], is.na(genotype[, i]), 'X')
}

genotype$complete  <- 0
genotype[genotype$Isolate_Id %in%  complete_virus, 'complete'] <- 1
genotype[genotype$Isolate_Id %in%  more8_virus, 'complete'] <- 2
  
  
genotype$genotype <- paste(genotype$PB2, genotype$PB1, genotype$PA, genotype$HA,
                           genotype$NP, genotype$`NA`, genotype$MP, genotype$NS, sep = '_')



A <- lapply(eight, function(x) {x[x$Isolate_Id %in% more8_virus, ]}) %>% do.call(what=rbind) %>% as.data.frame() %>% 
  group_split(Isolate_Id) %>% as.list() %>% lapply(as.data.frame)
names(A) <- lapply(A, FUN = function(x){x[1, 'Isolate_Id']})


A <- lapply(A, function(x) {table(x$lineages, x$segment) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq !=0)})  


for(i in names(A)) {
  x <- merge(A[[i]], data.frame(Var2=seg_level), by = 'Var2', all = T)
  times <- prod(table_DF(x$Var2) %>% select(Freq) %>% unlist() %>% unique())
  x <- replace(x, is.na(x), values = 'X')
  A[[i]] <- 
    data.frame(PB2=rep(x[x$Var2 == 'PB2', 'Var1'], each=times/(x[x$Var2 == 'PB2', 'Var1'] %>% length)), 
               PB1=rep(x[x$Var2 == 'PB1', 'Var1'], each=times/(x[x$Var2 == 'PB1', 'Var1'] %>% length)), 
               PA=rep(x[x$Var2 == 'PA', 'Var1'], each=times/(x[x$Var2 == 'PA', 'Var1'] %>% length)), 
               HA=rep(x[x$Var2 == 'HA', 'Var1'], each=times/(x[x$Var2 == 'HA', 'Var1'] %>% length)), 
               NP=rep(x[x$Var2 == 'NP', 'Var1'], each=times/(x[x$Var2 == 'NP', 'Var1'] %>% length)), 
               NA_gene=rep(x[x$Var2 == 'NA', 'Var1'], each=times/(x[x$Var2 == 'NA', 'Var1'] %>% length)), 
               MP=rep(x[x$Var2 == 'MP', 'Var1'], each=times/(x[x$Var2 == 'MP', 'Var1'] %>% length)), 
               NS=rep(x[x$Var2 == 'NS', 'Var1'], each=times/(x[x$Var2 == 'NS', 'Var1'] %>% length))) %>% 
    mutate(genotypes=paste(PB2, PB1, PA, HA, NP, NA_gene, MP, NS, sep = '_')) %>% 
    select(genotypes) %>% unlist() %>% paste(collapse = ', ')
}

more8_geno <- A %>% do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('Isolate_Id')
colnames(more8_geno)[2] <- 'multiples_geno'

genotype <- merge(genotype, more8_geno, by = 'Isolate_Id', all = T)
genotype$multiples_geno <- replace(genotype$multiples_geno, is.na(genotype$multiples_geno), values = '')



meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno$Host <- meta_geno$Host %>% str_replace_all(pattern = '_', replacement = ' ')

host <- meta_geno %>% group_split(Host) %>% as.list() %>% lapply(as.data.frame)
for(i in seq_along(host)) {
  host[[i]][, c('scentific_name', 'H_order', 'H_family')] <- 
    bird[bird$name %in% host[[i]][, 'Host'], c('scientific_name', 'order', 'family', 'w_d')]
}
meta_geno <- host %>% do.call(what=rbind) %>% as.data.frame()



# meta[meta$Strain_number %in% 1386870758, c('Subtype', 'sub')] <- c('H4N2', 'H4')
# write.csv(meta, file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv', quote = F, row.names = F)

# tree_list <- read_as_list(path = '~/Analysis/aiv/merge/0307/', file_type = 'nwk', prefix = 'aligned.nwk$')
# names(tree_list) <- str_extract(names(tree_list), pattern = '[A-Z]+[0-9]+_a|[a-z]+[A-Z]+_a|[A-Z]+_a|NS1true_a') %>% 
#   str_remove_all(pattern = '_a')

# group <- list()
# for(i in list.files(path = '/home/eric/Analysis/aiv/merge/0307/group/groups/'))  {
#   clade <- read_as_list(paste('/home/eric/Analysis/aiv/merge/0307/group/groups/' ,i , '/', sep = ''), file_type = 'txt')
#   for(j in names(clade))  {
#     clade[[j]] <- sub(pattern = 'tree tree_1 \\= \\[\\&R\\] ', clade[[j]][2,1], replacement = '') %>%
#       strsplit(split='\'') %>% unlist() %>% grep(pattern = '^EPI|^IRD', value = T) %>% as.data.frame()
#   }
#   names(clade) <- LETTERS[1:length(clade)]
#   clade <- do.call(rbind, clade) %>% as.data.frame() %>% tibble::rownames_to_column('group')
#   clade$group <- sub(clade$group, pattern = '\\..*', replacement = '') %>% str_remove(pattern = 'up')
#   clade[, c('epi', 'strain_number')] <- (clade$. %>% strsplit(split = '\\|') %>% do.call(what=rbind))[, c(1,2)]
#   clade$ui <- paste(clade$epi, clade$strain_number, sep = '|')
#   clade$s_g <- paste(i, clade$group, sep = '')
#   group[[i]] <- clade
# }
# genotype <- assign_geno(group_list = group, meta = meta)



clade_table <- table_DF(meta[grepl(meta$Subtype, pattern='H5'), 'Clade'])
clade_table[, 'Clade'] <- 'GsGD-others'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4'), 'Clade'] <- '2.3.4.4'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4b'), 'Clade'] <- '2.3.4.4b'
clade_table[grep(clade_table$x, pattern = 'nonGsGD'), 'Clade'] <- 'nonGsGD'
clade_table[grep(clade_table$x, pattern = 'IRD_unlabel'), 'Clade'] <- 'IRD_unlabel'


# GISAID IRD merge --------------------------------------------------------
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-m", "--meta", type='character', default=TRUE,
                    help="Requied input bed file")

parser$add_argument("-gis", "--gisaid_seq_fasta", type='character', default=TRUE,
                    help="Generated out Rdata")

parser$add_argument("-ird", "--ird_seq_fasta", type='character', default=TRUE,
                    help="Generated out Rdata")

parser$add_argument("-bc", "--box_cut", type='numeric', default=10,
                    help="Generated out Rdata")

parser$add_argument("-prop", "--non_nt_prop", type='numeric', default=5,
                    help="Generated out Rdata")

parser$add_argument("-ns", "--numbers_segments", type='numeric', default=4,
                    help="Generated out Rdata")

parser$add_argument("-p", "--out_fasta_prefix", type='character', default=TRUE,
                    help="Generated out Rdata")

parser$add_argument("-o", "--out_fasta_dir", type='character', default=TRUE,
                    help="Generated out Rdata")

# 
args <- parser$parse_args()
args$meta <- '~/Analysis/aiv/all/meta/'
args$gisaid_seq_fasta <- '~/Analysis/aiv/all/all.fasta'
args$ird_seq_fasta <- '/home/eric/Analysis/aiv/ird/IRD_Sequence.fa'
args$length_down_limit <- 0.8
args$known_outlier <- 0.8
args$out_fasta_prefix <- 'gisaid_IRD_merged_'
args$out_fasta_dir <- '~/Analysis/aiv/merge/0307/'
cat(paste('End of sequence preprocessing', sep = ''))

GISAID <- read_as_list(path = args$meta , prefix = '.csv', header = T, file_type = 'txt') %>% 
  do.call(what=rbind) %>% distinct(Isolate_Id, .keep_all = T) %>% select(matches("INSDC")) %>% unlist %>% 
  str_split(pattern = ', ') %>% unlist


# processing IRD sequences: import
ird <- read.fasta(args$ird_seq_fasta, as.string = T, whole.header = T) %>% 
  do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('header') %>% 
  mutate(header=str_remove_all(header, pattern=' ')
         , index = word(header, 1, sep = '/')) %>% select(header)

IRD <- header_cleaning(header_vector = ird$header, pattern = '\\/') %>% select(V2) %>% unlist


# lt <- list(set1 = ncbi %>% unlist(), set2 = ird$V2 %>% unlist())meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'

# ComplexHeatmap::list_to_matrix(lt)
catelog <- data.frame(cate=union(GISAID, IRD), GISAID=0, IRD=0) %>% tidyr::drop_na(cate)
catelog[catelog$cate %in% GISAID %>% unlist(), 'GISAID'] <- 1
catelog[catelog$cate %in% IRD %>% unlist(), 'IRD'] <- 1
catelog <- catelog %>% tibble::column_to_rownames('cate')

# library(ComplexUpset)
# upset(catelog ,colnames(catelog)[1:2]
#       , width_ratio = 0.4 #ratio between horizen and vertical bar
#       , name='Conditions' # names of verticle bar
#       , base_annotations=list('Intersection size'=
#                                 intersection_size(text=list(vjust=-0.1,hjust=0.5, size=5))) # text in group bar
#       , set_sizes=(upset_set_size()+ theme(axis.text.x=element_text(size=16)
#                                            , axis.title = element_text(size=20)))
#       , themes=upset_modify_themes(
#         list('Intersection size'=theme(axis.text=element_text(size=15)
#                                        , axis.title=element_text(size=30))
#              , 'overall_sizes'=theme(axis.text.x=element_text(size = 24))
#              , 'intersections_matrix'=theme(text=element_text(size=24))))
#       ,  wrap=TRUE
#       , matrix=intersection_matrix(geom=geom_point(stroke=4))
#       , stripes=upset_stripes(geom = geom_segment(size=15)))+
#   labs(title = 'Translated sgmRNA among conditions'
#        , subtitle = 'Recombination more than 1 ')+
#   theme(title = element_text(size=36))

# sequences number --------------------------------------------------------
library("VennDiagram") 

png(paste('~/Analysis/aiv/merge/0307/result/sequences_number_venn.png', sep = ''), width = 9600, height = 5400, res=600)
# move to new plotting page 
grid.newpage() 
draw.pairwise.venn(area1=320698, area2=221754, cross.area=181948, 
                   category=c("GISAID", "IRD"), fill=c("#80DFAE", "#3BAEDA"),
                   cat.cex = c(2, 2), cex= c(2), fontfamily=('NimbusMon'))
dev.off()


# genome length -----------------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'
segment <- organ_subtype(in_seq = meta)


plot_genome_l <- meta[, c('sub', 'nt_number')] %>% mutate(g='after filtering')
filter_min <- plot_genome_l %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(FUN = function(x) {min(x$nt_number)}) %>% unlist()
names(filter_min) <- plot_genome_l %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(FUN = function(x) {unique(x$sub)}) %>% unlist()

genome_median <- plot_genome_l %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(FUN = function(x) {median(x$nt_number)}) %>% unlist()
names(genome_median) <- plot_genome_l %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(FUN = function(x) {unique(x$sub)}) %>% unlist()


genome_thres <- genome_median*0.8

# meta %>% 
#   ggplot()+
#   geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number))+
#   scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
#   labs(x='Segment', y='Genome length', title='Genome length distribution after filtering')+
#   theme_bw()+gg_theme

ORF_length <- c(2280, 2274 ,2148, rep(1707, 16) ,1497 , rep(1410, 9), 759, 693)
names(ORF_length) <- seg_sub_level

meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information_raw.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'
segment <- organ_subtype(in_seq = meta)
for(i in names(segment)) {
  A <- segment[[i]][, c('sub', 'nt_number', 'segment')] %>% mutate(color='Raw')
  # A[A$nt_number < filter_min[i], 'color'] <- 'Removed'
  # A[(nrow(A)+1), ] <- c(i, ORF_length[i], i, 'ORF')
  segment[[i]] <- A
}

# plot_genome_l <- rbind(plot_genome_l, meta[, c('sub', 'nt_number')] %>% mutate(g='before filtering')) %>% 
#   mutate(g=factor(g, levels=c('before filtering', 'after filtering')))
plot_genome_l <- do.call(rbind, segment) %>% as.data.frame() %>% mutate(nt_number=as.numeric(nt_number))

genome_l <- matrix(lapply(segment, function(x){summary(x[, 'nt_number'])}) %>% unlist(), nrow = length(segment), ncol = 6, byrow = T) %>% 
  as.data.frame() %>% mutate('medain_80%'=V3*0.8)

rownames(genome_l) <- names(segment)
colnames(genome_l)[1:6] <- names(segment[[31]][, 'nt_number'] %>% summary)
write.csv(genome_l[seg_sub_level, ], file = '~/Analysis/aiv/merge/0307/result/genome_length.csv')

png('~/Analysis/aiv/merge/0307/result/length_box.png', width = 9600, height = 5400, res = 600)
plot_genome_l %>% as.data.table() %>% 
  ggplot()+
  geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=segment))+
  # geom_hline(plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color=='ORF'), 
  #            mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  # geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  # geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', title='Genome length distribution', fill='')+
  # facet_wrap(~factor(sub, level=seg_sub_level))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16))
dev.off()

png('~/Analysis/aiv/merge/0307/result/length_box_internal.png', width = 9600, height = 5400, res = 600)
plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color!='ORF') %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  geom_hline(plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color=='ORF'), 
             mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', title='Genome length distribution before/after filtering', fill='')+
  facet_wrap(~factor(sub, level=seg_sub_level))+
  theme_bw()+gg_theme
dev.off()

png('~/Analysis/aiv/merge/0307/result/length_box_HA.png', width = 9600, height = 5400, res = 600)
plot_genome_l %>% filter(sub %in% c(paste('H', seq(1, 8, 1), sep = '')), color!='ORF') %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  geom_hline(plot_genome_l %>% filter(sub %in% c(paste('H', seq(1, 8, 1), sep = '')), color=='ORF'), 
             mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', title='Genome length distribution before/after filtering', fill='')+
  facet_wrap(~factor(sub, level=seg_sub_level), nrow = 2)+
  theme_bw()+gg_theme+theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14))

plot_genome_l %>% filter(sub %in% c(paste('H', seq(9, 16, 1), sep = '')), color!='ORF') %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  geom_hline(plot_genome_l %>% filter(sub %in% c(paste('H', seq(9, 16, 1), sep = '')), color=='ORF'), 
             mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Gephylo_bsnome length', title='Genome length distribution before/after filtering', fill='')+
  facet_wrap(~factor(sub, level=seg_sub_level), nrow = 2)+
  theme_bw()+gg_theme+theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14))
dev.off()

png('~/Analysis/aiv/merge/0307/result/length_box_NA.png', width = 9600, height = 5400, res = 600)
plot_genome_l %>% filter(sub %in% c(paste('N', seq(1, 9, 1), sep = '')), color!='ORF') %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  geom_hline(plot_genome_l %>% filter(sub %in% c(paste('N', seq(1, 9, 1), sep = '')), color=='ORF'), 
             mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', title='Genome length distribution before/after filtering', fill='')+
  facet_wrap(~factor(sub, level=seg_sub_level), nrow = 3)+
  theme_bw()+gg_theme+theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14))
dev.off()


length_diff <- 
  merge(ORF_length %>% as.data.frame() %>% tibble::rownames_to_column('sub') %>% rename('ORF length'='.'), 
        genome_thres %>% as.data.frame() %>% tibble::rownames_to_column('sub') %>% rename('80% median genome length'='.'), by = 'sub') 
# tidyr::gather(key, value, 2:3) %>% mutate(sub=factor(sub, levels=seg_sub_level))

length_diff <-
  merge(genome_median %>% as.data.frame() %>% tibble::rownames_to_column('sub') %>% rename('Genome median length'='.'),
        length_diff, by = 'sub') %>% tidyr::gather(key, value, 2:4)


meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information_raw.csv') %>% as.data.frame() %>% 
  mutate(key='violin') %>% select(c('sub', 'key', 'nt_number'))
colnames(meta) <- colnames(length_diff)
length_diff <- rbind(meta, length_diff) %>% 
  mutate(sub=factor(sub, levels=seg_sub_level), 
         key=factor(key, levels=c('Genome median length', 'ORF length', '80% median genome length', 'violin')))


png('~/Analysis/aiv/merge/0307/result/ORF_length_filter_length.png', width = 9600, height = 5400, res = 600)
ggplot()+
  # geom_violin(length_diff %>% filter(key == 'violin'), mapping=aes(x=sub, y=value), linewidth=0.5)+
  geom_boxplot(length_diff %>% filter(key == 'violin'), mapping=aes(x=sub, y=value), linewidth=0.5)+
  geom_point(length_diff %>% filter(key != 'violin'), mapping=aes(x=sub, y=value, color=key, shape=key), size=5)+
  # geom_line(length_diff %>% filter(key != 'violin'), mapping=aes(x=sub, y=value, group=key, color=key), linewidth=1.2)+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  scale_color_manual(values = c('80% median genome length'='#FFA888', 'ORF length'='#034FDE', 'Genome median length'='#48CFAE'))+
  labs(x='', y='Genome length', color='', shape='', 
       title = 'Genomic Length Variation Across Genes: ORF and Median Comparisons')+
  theme_bw()+gg_theme+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=14), legend.position="bottom")
dev.off()


length_diff %>% filter(sub %in%  c(paste('H', seq(1, 16, 1), sep = ''))) %>% 
  ggplot()+
  geom_point(mapping=aes(x=sub, y=value, color=key), size=3)+
  geom_line(mapping=aes(x=sub, y=value, group=key, color=key), linewidth=1.2)+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='HA genes', y='Genome length', color='')+
  theme_bw()+gg_theme


length_diff %>%filter(sub %in%  c(paste('N', seq(1, 9, 1), sep = ''))) %>% 
  ggplot()+
  geom_point(mapping=aes(x=sub, y=value, color=key), size=3)+
  geom_line(mapping=aes(x=sub, y=value, group=key, color=key), linewidth=1.2)+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='NA genes', y='Genome length', color='')+
  theme_bw()+gg_theme


length_diff %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS')) %>% 
  ggplot()+
  geom_point(mapping=aes(x=sub, y=value, color=key), size=3)+
  geom_line(mapping=aes(x=sub, y=value, group=key, color=key), linewidth=1.2)+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', color='')+
  theme_bw()+gg_theme


# ambiguous characters ----------------------------------------------------
meta_raw <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information_raw.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'
segment <- organ_subtype(in_seq = meta[, c('Subtype', 'sub', 'non_nt_percent', 'segment')])
for(i in names(segment)) {
  A <- segment[[i]][, c('sub', 'non_nt_percent', 'segment')] %>% mutate(color='Reserved')
  # A <- table_DF(A$non_nt_percent) %>% arrange(desc(x)) %>% mutate(cum=cumsum(Freq)/nrow(A)*100, sub=i)
  A[A$non_nt_percent > 5, 'color'] <- 'Removed'
  segment[[i]] <- A
}


ambi_letter <- do.call(rbind, segment) %>% as.data.table() %>% 
  mutate(sub=factor(sub, levels=seg_sub_level), x=as.numeric(x))

plot_list <- list()
for(i in names(class)) {
  # plot_list[[i]] <-
  ambi_letter %>% filter(sub %in% class[[i]]) %>% 
    ggplot()+
    geom_hline(yintercept = 5, color='red', linetype='dashed')+
    geom_jitter(mapping=aes(x=sub, y=non_nt_percent, color=color), width=0.05, size=0.75)+
    geom_violin(mapping=aes(x=sub, y=non_nt_percent), width=1, fill='pink')+
    # geom_histogram(mapping=aes(x=non_nt_percent, fill=sub))
    # geom_point(smapping=aes(x=x, y=cum, color=sub), size=2)+
    scale_y_continuous(labels = seq(0, 100, 5), breaks = seq(0, 100, 5))+
    scale_color_manual(values = c('#E26309', '#034FDE'))+
    labs(x='', y='Ambiguous letters proportion', color='')+
    guides(color = guide_legend(override.aes = list(size=5), nrow = 1))+
    theme_bw()+gg_theme+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=24), legend.position="bottom")
}

for(i in names(plot_list)) {
  png(paste('~/Analysis/aiv/merge/0307/result/', i,'_ambiguous_proportion.png', sep=''), width = 9600, height = 5400, res = 600)
  print(plot_list[[i]])
  dev.off()
}



# ambi_letter %>% 
#   ggplot(aes(y = sub, x= non_nt_percent)) +
#   ggridges::geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.00000000000000000000000001)+
#   scale_fill_viridis(name = "Quality score", option = "H") +
#   labs(x='Read quality score') +
#   geom_vline(xintercept = 7, linetype="dotted", 
#              color = "Red", size=1.5)+
#   scale_x_continuous(limits = c(0, 25), breaks = c(seq(0,25,5)))+
#   # guides(fill=guide_legend(title=""))+
#   theme_bw()+gg_theme

lapply(segment, function(x) {boxplot.stats(x$non_nt_percent)$out %>% min()}) %>% unlist()

# alignment length -----------------------------------------------------------
align <- read_as_list('/home/eric/Analysis/aiv/merge/0307/ORF_filter_outlier/', prefix = '_aligned_fill_trimed_remove_outlier.fa$', file_type = 'fasta', as.string = T)
align_df <- data.frame(seq=align %>% unlist() %>% unlist()) %>% tibble::rownames_to_column('header') %>% 
  mutate(nt_l=nchar(str_remove_all(seq, pattern='-')), header=str_remove(header, pattern='.*outlier\\.'))
align_df <- cbind(align_df, header_cleaning(align_df$header, pattern = '\\|'))
colnames(align_df)[4:ncol(align_df)] <- c("Accesion_number","Strain_number","Subtype", 'Clade', "segment", 'Collection_date'
                                          ,"Location","Host", "Header_Host", "Host_type")
align_df <- align_df %>% mutate(gap=str_count(seq, pattern = '-'), 
                                non_atcg=str_count(pattern = '[^atcg-]', seq)) #

segment <- organ_subtype(in_seq = align_df)

align_l <- matrix(lapply(segment, function(x){summary(x[, 'nt_l'])}) %>% unlist(), nrow = length(segment), ncol = 6, byrow = T) %>% 
  as.data.frame() 

rownames(align_l) <- names(segment)
colnames(align_l)[1:6] <- names(segment[[31]][, 'nt_l'] %>% summary)
align_l <- merge(align_l %>% tibble::rownames_to_column('seg'), as.data.frame(ORF_length) %>% tibble::rownames_to_column('seg'), by = 'seg')

plot_align_l <- segment %>% do.call(what=rbind) %>% as.data.table() %>% select(c('nt_l', 'sub', 'segment'))
plot_align_l %>% 
  ggplot()+
  geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_l %>% as.numeric(), fill=segment))+
  # geom_hline(plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color=='ORF'), 
  #            mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=1.5)+
  # geom_jitter(aes(x=1, y=nt_number, color=color), width = 0.25)+
  # geom_violin(aes(x=1, y=nt_number), fill='pink')+
  scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Genome length', title='Genome length distribution', fill='')+
  # facet_wrap(~factor(sub, level=seg_sub_level))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16))


plot_list <- list()
for(i in seg_sub_level) {
  test <- plot_align_l %>% filter(sub %in% i) %>% table_DF() %>% mutate(cum=cumsum(Freq), nt_l=nt_l%>% as.numeric())
  test <- test %>% mutate(cum=cum/(test[nrow(test), 5])*100)
  l1 <- max(test$nt_l)*0.9
  l2 <- max(test$nt_l)*0.95
  l3 <- max(test$nt_l)*0.99
  test[(nrow(test)+1):(nrow(test)+3), 1] <- c(l1, l2, l3)
  test[(nrow(test)-2):(nrow(test)), 2] <- c('l1', 'l2', 'l3')
  print(c(l1, l2, l3))
  plot_list[[i]] <- 
    ggplot()+
    geom_point(test %>% filter(sub == i), mapping=aes(x=nt_l , y=cum, fill='black'), size=3)+
    geom_line(test %>% filter(sub == i), mapping=aes(x=nt_l %>% as.numeric(), y=cum, group=1, color='pink'), linewidth=1.2)+
    geom_vline(test %>% filter(sub %in% c('l1', 'l2', 'l3')), mapping = aes(xintercept = nt_l, color=sub), linetype='dashed', linewidth=1.5)+
    # geom_text(aes(x=max(test$nt_l)*0.9, y=110), label= paste('90% ORF'), hjust=-0.1, size=5, color='pink')+
    # geom_vline(test %>% filter(sub == 'l'), mapping = aes(xintercept = nt_l), color='red', linetype='dashed', linewidth=1.5)+
    # geom_text(aes(x=max(test$nt_l)*0.95, y=110), label= paste('95% ORF'), hjust=-0.1, size=5, color='red')+
    # geom_vline(test %>% filter(sub == 'l'), mapping = aes(xintercept = nt_l), color='maroon', linetype='dashed', linewidth=1.5)+
    # geom_text(aes(x=max(test$nt_l)*0.99, y=110), label= paste('99% ORF'), hjust=-0.1, size=5, color='maroon')+
    # scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))+
    # scale_x_continuous(breaks = c(seq(min(test$nt_l), max(test$nt_l), 50), max(test$nt_l)), 
    #                    limits = c(min(test$nt_l), max(test$nt_l)))+
    labs(x='ORF length', y='Cummulative proportion', title=i, fill='')+
    # facet_wrap(~factor(sub, level=seg_sub_level))+
    theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16), legend.position = 'none')
}



dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_cumm_orf_length.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}



# pd_distb <- list()
for(i in seg_sub_level) {
  plot_list <- list()
  cat(paste('Doing', i, 'calculation\n', sep=' '))
  dis <- readRDS(paste('/home/eric/Analysis/aiv/merge/0307/distance/', i, '_dist_matrix.RData', sep = ''))
  n <- ncol(dis)
  # dis_stat <- data.frame(pd=table(dis) %>% as.data.frame()) %>% mutate(pd.dis=as.numeric(pd.dis))
  # dis_stat[1, 2] <- dis_stat[1, 2]-n
  # dis_stat$pd.Freq <- dis_stat$pd.Freq/2/(n*n-n)
  A <- segment[[i]] %>% select(seq) %>% unlist %>% header_cleaning(pattern = '')
  cons <- c()
  for(j in 1:ncol(A)) {
    cons <- c(cons, (table_DF(A[, j]) %>% arrange(desc(Freq)))[1, 1])
  }
  cons_mat <- matrix(cons, nrow = nrow(A), ncol = length(cons), byrow = T)
  cons_mat <- ifelse(A==cons_mat, yes = 0, no = 1)
  B <- merge(segment[[i]] , data.frame(pd=colSums(dis)/n) %>% tibble::rownames_to_column('header'), by = 'header')
  B$diff_cons <- rowSums(cons_mat)
  
  plot_list[[2]] <- 
    ggplot(B, aes(x=diff_cons, y=pd))+
    geom_point()+
    sm_statCorr(color = '#0f993d', corr_method = 'spearman',linetype = 'dashed')+
    # geom_line(aes(x=pd, y=Freq1, group=1))+
    # geom_histogram(aes(x=pd), bins = 50)+
    labs(title = 'Scatter of Different coordiantes and MPD\nfor each sequence', x='Different coordiantes between concensus', y='Mean pairwise distance\nfor each seqeunce (%)')+
    theme_bw()+gg_theme
  
  plot_list[[3]] <- 
    ggplot(B, aes(x=gap, y=pd))+
    geom_point()+
    sm_statCorr(color = '#0f993d', corr_method = 'spearman',linetype = 'dashed')+
    # geom_line(aes(x=pd, y=Freq1, group=1))+
    # geom_histogram(aes(x=pd), bins = 50)+
    labs(title = 'Scatter of gap numbers and MPD\nfor each sequence', x='Gaps ', y='Mean pairwise distance\nfor each seqeunce (%)')+
    theme_bw()+gg_theme
  
  plot_list[[4]] <- 
    ggplot(B, aes(x=gap, y=pd))+
    geom_point(aes(x=non_atcg, y=pd))+
    sm_statCorr(color = '#0f993d', corr_method = 'spearman',linetype = 'dashed')+
    # geom_line(aes(x=pd, y=Freq1, group=1))+
    # geom_histogram(aes(x=pd), bins = 50)+
    labs(title = 'Scatter of non-ATCG characters and MPD\nfor each sequence', x='Non-ATCG characters', y='Mean pairwise distance\nfor each seqeunce (%)')+
    theme_bw()+gg_theme
  
  pd_plot <- data.frame(pd=colSums(dis)/n) %>% table_DF() %>% arrange(pd) %>% 
    mutate(g=1, Freq1=cumsum(Freq), pd=as.numeric(pd))
  
  plot_list[[1]] <-
    ggplot(pd_plot)+
    geom_point(aes(x=pd, y=Freq1))+
    geom_line(aes(x=pd, y=Freq1, group=1))+
    # geom_histogram(aes(x=pd), bins = 50)+
    labs(title = paste(i, 'Cummulative distribution', sep = ', '), x='Mean pairwise distance\nfor each seqeunce (%)', y='Cummulative Frequence')+
    theme_bw()+gg_theme
  gc()
  # pd_distb[[i]] <- cowplot::plot_grid(plotlist = plot_list)
  # align_l[i, 'k80_pd (%)'] <- sum(dis)/(n*n-n)*100
}

write.csv(align_l[seg_sub_level, c(3, 7)] %>% rename('aligned_ORF_length'='Median'), 
          file = '~/Analysis/aiv/merge/0307/result/aligned_ORF_length.csv')

pd_df <- list()
for(i in paste('H', seq(1, 16, 1), sep = '')) {
  cat(paste('Doing', i, 'calculation\n', sep=' '))
  dis <- readRDS(paste('/home/eric/Analysis/aiv/merge/0307/distance/', i, '_dist_matrix.RData', sep = ''))
  n <- ncol(dis)
  
  pd_df[[i]] <- data.frame(pd=dis[lower.tri(dis, diag = F)], x=i)
  gc()
  # align_l[i, 'k80_pd (%)'] <- sum(dis)/(n*n-n)*100
}

pd_df %>% do.call(what=rbind) %>% 
  ggplot()+
  geom_violin(aes(y=pd, x=x))+
  labs(title = i, x='pairwise distance (%)', y='Proportion')+
  theme_bw()+gg_theme



# stop_ORF_genome ---------------------------------------------------------
# args$in_fasta_dir <- '/home/eric/Analysis/aiv/merge/0307/aligned_fasta/'

# aligned <- read_as_list(path = args$in_fasta_dir, prefix = 'aligned.fa', file_type = 'fasta')
j="NS"
# j="all_clean_H15_aligned" 
# out_list <- list()
# stop_codon_list <- list()
gapthres_stop <- gapthres_orf_length <- gapthres_genome_length <-
  matrix('', nrow = length(seg_sub_level), ncol = length(c((c(seq(40, 80, 5), seq(82, 100, 2)))))) %>% as.data.frame() 
colnames(gapthres_genome_length) <- colnames(gapthres_orf_length) <- colnames(gapthres_stop) <- 
  c(seq(40, 80, 5), seq(82, 100, 2))
z <- 0
for(j in seg_sub_level)  {
  cat(paste('Start processing ', j, '\n', sep=''))
  align <- read.fasta(paste('/home/eric/Analysis/aiv/merge/0307/aligned_fasta/', "gisaid_IRD_merged_", j, '_aligned.fa', sep = ''), as.string = F) %>% 
    do.call(what=rbind) %>% as.data.frame() 
  seg <- j
  z <- z+1
  
  if(grepl(pattern = '(H7|H5)', j))  { #grepl(pattern = '(H7|H5|H9)', j)
    cat(paste('Performing Gap fill procedure', '\n', sep = ''))
    gap <- apply(align[, 1:ncol(align)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% as.data.frame()
    gap$index <- seq(1, nrow(gap))
    
    cons <- c()
    for(i in 1:ncol(align))  {
      test <- table(align[, i])
      cons <- c(cons, which.max(test) %>% names())
    }
    cons <- paste(cons, collapse = '')
    
    atg1 <- gregexpr(cons, pattern = 'atg')[[1]][1]
    
    ## limit in ORF
    fill <- gap#[((atg1)+50):(nrow(gap)-200), ]
    fill <- fill[fill$. >= nrow(align)*0.95, ]
    
    
    fill$in2 <- c(fill$index[-1], 0)
    fill$in2 <- fill$in2-fill$index
    fill[fill$in2 < 15, 'in2'] <- 1
    rownames(fill) <- seq(1, nrow(fill),1)
    
    v <- (fill[!(fill$in2 %in% 1), ] %>% rownames())
    
    if(length(v) == 1)  {
      v <- ifelse(identical(v, character(0)), yes = nrow(fill), no = v)
    }
    row <- c(1, v)
    sub_fill_list <- list()
    x=0
    for(i in 2:(length(row)))  {
      B <- fill[row[i-1]:row[i], ]
      x=x+1
      sub_fill_list[[x]] <- B
      # if(nrow(B)>=5)  {
      #  
      # }
    }
    
    for(i in seq_along(sub_fill_list))  {
      sub_fill <- sub_fill_list[[i]]
      if((sub_fill[nrow(sub_fill), 2]-sub_fill[2, 2]) > 100) {
        cat(paste('Processing region: ', sub_fill[2,2], ' to ', sub_fill[nrow(sub_fill),2], sep = ''))
        test <- align[, (sub_fill[2,2]-5):(sub_fill[nrow(sub_fill),2]+5)]
        l <- ncol(test)
        
        test$fill <- apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '') %>% str_remove_all(pattern = '-')}) #
        test$sub_l <- l-nchar(test$fill)
        
        for(k in 1:nrow(test))  {
          test[k, 'fill'] <- paste(test[k, 'fill'], rep('-', test[k, 'sub_l']) %>% paste(collapse = ''), sep = '')
        }
        
        test <- str_split(test$fill, pattern = '') %>% do.call(what=rbind) %>% as.data.frame()
        
        align[, (sub_fill[2,2]-5):(sub_fill[nrow(sub_fill),2]+5)] <- test
      }else(next)
    }
  }
  
  y <- 0
  for(thres in colnames(gapthres_stop) %>% as.numeric()) {
    y <- y+1
    A <- align
    
    gap <- apply(A[, 1:ncol(A)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% 
      as.data.frame() %>% mutate(index=row_number(), prop=./nrow(A)*100)
    
    gap <- gap[gap$prop < thres, ]
    A <- A[, gap$index]
    # gapthres_genome_length[z, y] <- ncol(A)
    
    cons <- c()
    for(i in 1:ncol(A))  {
      test <- table(A[, i])
      cons <- c(cons, which.max(test) %>% names())
    }
    atg1 <- (paste(cons, collapse = '') %>% gregexpr(pattern = 'atg') %>% unlist())[1]
    
    test <- apply(A[1:nrow(A), atg1:ncol(A)], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
    s <- seq(1, nchar(test[1]), 3)
    e <- seq(3, nchar(test[1]), 3)
    s <- s[1:min(length(s), length(e))]
    e <- e[1:min(length(s), length(e))]
    B <- sapply(test, function(x) {substring(x, s, e)}) %>% t()
    
    if(!seg %in% 'NS') {
      B <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
      stop1 <- which.max(B$.)
      test <- A[, atg1:(stop1*3+atg1-1-3)]
      B <- sapply(apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
                  , function(x) {substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))}) %>% t()
      stop <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
    }
    
    if(seg %in% 'NS') {
      for(i in seq_len(nrow(B))) {
        stop <- B[i, ] %in% c('taa', 'tag', 'tga') %>% which()
        stop <- stop[stop >190] %>% min()
        if(!is.infinite(stop)) {
          B[i, (stop):(ncol(B))] <- '---'
        }
      }
      stop <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
      
      test <- apply(B %>% as.data.frame(), MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
      test <- sapply(test, function(x) {str_sub(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))}) %>% t() %>% as.data.frame()
      
      remove <- ''
      for(i in seq_len(ncol(test))) {
        remove <- c(remove, all(test[, i] == '-'))
      }
      test <- test[, -c(which(remove[-1] %>% as.logical()))]
    }
    
    gapthres_orf_length[z, y] <- ncol(test)
    gapthres_stop[z, y] <- stop %>% sum()
    # gc()
    print(c(thres, ncol(test), stop %>% sum()))
  }
  
  rownames(gapthres_orf_length)[z] <- rownames(gapthres_stop)[z] <- seg
}


j='gisaid_IRD_merged_H1_aligned'
gap_prop <- list()
for(j in names(aligned)) {
  cat(paste('Start processing ', j, '\n', sep=''))
  seg <- str_extract(j, pattern = '[A-Z]+[0-9]+|[A-Z]{2}_a') %>% str_remove(pattern = '_a')
  align <- aligned[[j]] %>% do.call(what=rbind) %>% as.data.frame() 
  
  gap_prop[[seg]] <- data.frame(prop=(ifelse(align=='-', yes = 1, no = 0) %>% colSums())/nrow(align), seg=seg)
}
gc()
gap_prop_plot <- gap_prop %>% do.call(what=rbind) %>% as.data.frame() 

gap_prop_plot %>% mutate(seg=factor(seg, levels=seg_sub_level)) %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  geom_hline(mapping = aes(yintercept = 90), color='red', linetype='dashed', linewidth=1)+
  geom_jitter(aes(x=seg, y=prop*100, color=seg), width = 0.25)+
  # geom_violin(aes(x=seg, y=prop*100), fill='pink')+
  # scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(y='Genome length', title='Genome length distribution before/after filtering', fill='')+
  guides(color = guide_legend(nrow = 2, position = 'bottom', direction = ))+
  theme_bw()+gg_theme+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16))

plot_list <- list()
for(i in names(gap_prop)) {
  A <- gap_prop[[i]]
  A <- A %>% mutate(coord=seq(1, nrow(A))) %>% as.data.table()
  plot_list[[i]] <-
    A %>% 
    ggplot()+
    # geom_line(mapping=aes(x=coord, y=prop*100), linewidth=0.5)+
    geom_point(mapping=aes(x=coord, y=prop*100))+
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))+
    scale_x_continuous(breaks = seq(0, max(A$coord), 250), limits = c(0, max(A$coord)))+
    labs(y='Gap proportion', x='Coordinates',  color='', title=i)+
    # guides(color = guide_legend(nrow = 2, position = 'bottom', direction = ))+
    theme_bw()+gg_theme#+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16))
}

dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_gap_proportion_bycoord.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}



gap_prop_cum <- lapply(gap_prop, function(y) {y %>% table_DF() %>% arrange(x) %>% mutate(cum=cumsum(Freq))}) #/nrow(y)*100
gap_prop_plot <- gap_prop_cum %>% do.call(what=rbind) %>% as.data.frame() 

plot_list <- list()
for(j in names(class)) {
  B <- gap_prop_plot %>% mutate(seg=factor(seg, levels=seg_sub_level), prop=as.numeric(prop)*100) %>% 
    filter(seg %in% class[[j]]) %>% as.data.table()
  plot_list[[j]] <-
    B %>% 
    ggplot()+
    # geom_hline(yintercept = 90, color='red', linetype='dashed')+
    geom_point(mapping=aes(x=cum, y=prop, color=seg), size=2)+
    geom_line(mapping=aes(x=cum, y=prop, group=seg, color=seg), linewidth=1)+
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))+
    scale_x_continuous(breaks = seq(0, max(B$cum), 250), limits = c(0, max(B$cum)))+
    labs(y='Gap proportion', x='',  color='')+
    # guides(color = guide_legend(nrow = 2, position = 'bottom', direction = ))+
    theme_bw()+gg_theme+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16), legend.key.size = unit(2, 'lines'))
}

dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_gap_cum_proportion.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}




out <- cbind(gapthres_stop %>% tibble::rownames_to_column('seg') %>% tidyr::gather(key, value, c(2:(ncol(gapthres_stop)+1))) %>% rename(stop=value), 
             gapthres_orf_length %>% tibble::rownames_to_column('seg') %>% tidyr::gather(key, value, c(2:(ncol(gapthres_orf_length)+1)))%>% select(c('value')) %>% rename(orf_length=value)) 
# write.csv(out, file = '~/Analysis/aiv/merge/0307/stop_orf_length.csv', row.names = F)
out <- fread(file = '~/Analysis/aiv/merge/0307/stop_orf_length.csv')
plot_list <- list()
dir.create('~/Analysis/aiv/merge/0307/plot/')
for(i in unique(out$seg)) {
  B <- out %>% filter(seg %in% i) %>% tidyr::gather(key1, value, c(3:4)) %>% mutate(value=as.numeric(value))
  # B$key <- factor(x = B$key, levels =  c(seq(40, 95, 5), 99, 100))
  B$key <- as.numeric(B$key)
  coeff <- ((B[20:38, 'value'] %>% max())/(B[1:19, 'value'] %>% max() %>% log10())*0.8)
  yint <- B %>% filter(key1 %in% 'orf_length') %>% select(value) %>% unlist() %>% as.numeric() %>% median()
  
  stop <- (B %>% filter(key1 %in% 'stop') %>% select(value) %>% unlist() %>% as.numeric())[14] 
  
  print(c(i, yint, stop))
  # plot_list[[i]]
  png(paste('~/Analysis/aiv/merge/0307/plot/', i, '_stop_orf_genome.png', sep = ''), width = 9600, height = 5400, res = 600)
  print(
    ggplot()+
      geom_point(B[20:38, ], mapping=aes(x=key, y=value, color=key1), size=3)+
      geom_line(B[20:38, ], mapping=aes(x=key, y=value, group=key1, color=key1), linewidth=1.2)+
      geom_point(B[1:19, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), size=3)+
      geom_line(B[1:19, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), linewidth=1.2)+
      geom_hline(yintercept = yint, color='red', linetype='dashed')+
      geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
      geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
      # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
      # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
      labs(x='Gap threshold (%)', title = i, color='')+
      scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
      scale_y_continuous(name = "length", sec.axis = sec_axis(~./coeff, name="log10(Number of stop codons)"))+
      scale_x_continuous( breaks=   c(seq(40, 80, 5), seq(82, 100, 2)), labels =  c(seq(40, 80, 5), seq(82, 100, 2)))+
      theme_bw()+gg_theme
  )
  dev.off()
}


out <- out %>% mutate(stop=stop+1, seg=factor(seg, levels=seg_sub_level)) %>% tidyr::gather(key1, value, 3:4) %>% as.data.table()
class <- list(Internal=seg_sub_level[c(1:3, 20, 30, 31)], HA=seg_sub_level[c(4:19)], N=seg_sub_level[c(21:29)])
plot_list <- list()
for(i in c(out$key1 %>% unique())) {
  for(j in names(class)) {
    B <- out %>% filter(key1==i, seg %in% class[[j]])
    if(i == 'stop') {
      plot_list[[paste(i, j, sep = '_')]] <- 
        ggplot(B)+
        geom_point(mapping=aes(x=key, y=value, color=seg), size=3)+
        geom_line(mapping=aes(x=key, y=value, group=seg, color=seg), linewidth=1.2)+
        # geom_point(B[1:14, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), size=3)+
        # geom_line(B[1:14, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), linewidth=1.2)+
        # geom_hline(yintercept = yint, color='red', linetype='dashed')+
        # geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
        # geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
        # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
        # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
        labs(x='Gap threshold (%)', color='', y='Log10(Numbers of stop codon)')+
        # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
        # scale_y_continuous(name = "length", sec.axis = sec_axis(~./coeff, name="log10(Number of stop codons)"))+
        # scale_y_continuous(name="log10(Number of stop codons)", trans = 'log10')+
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_x_continuous(breaks= c(seq(40, 80, 5), seq(82, 100, 2)), labels =  c(seq(40, 80, 5), seq(82, 100, 2)))+
        theme_bw()+gg_theme+annotation_logticks(sides = "l")+
        theme(legend.text = element_text(size = 16), legend.key.size = unit(2.5, 'lines'))
    }
    if(i == 'orf_length') {
      plot_list[[paste(i, j, sep = '_')]] <-
        ggplot(B)+
        geom_point(mapping=aes(x=key, y=value, color=seg), size=3)+
        geom_line(mapping=aes(x=key, y=value, group=seg, color=seg), linewidth=1.2)+
        # geom_point(B[1:14, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), size=3)+
        # geom_line(B[1:14, ], mapping=aes(x=key, y=log10(value+1) * coeff, group=key1, color=key1), linewidth=1.2)+
        # geom_hline(yintercept = yint, color='red', linetype='dashed')+
        # geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
        # geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
        # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
        # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
        labs(x='Gap threshold (%)', color='', y='identified open reading frame length')+
        # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
        # scale_y_continuous(name = "length", sec.axis = sec_axis(~./coeff, name="log10(Number of stop codons)"))+
        scale_y_continuous(breaks= c(seq(0, max(B$value), 125)), labels =  c(seq(0, max(B$value), 125)), limits = c(0, max(B$value)))+
        scale_x_continuous(breaks= c(seq(40, 80, 5), seq(82, 100, 2)), labels =  c(seq(40, 80, 5), seq(82, 100, 2)))+
        theme_bw()+gg_theme+theme(legend.text = element_text(size = 16), legend.key.size = unit(2.5, 'lines'))
    }
  }
}


dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}


A <- read.fasta(file = '~/Analysis/aiv/merge/0307/ORF_filter_outlier/gisaid_IRD_merged_NS_aligned_fill_trimed_remove_outlier.fa', as.string = T) %>% do.call(what=rbind) %>% 
  as.data.frame() %>% mutate(V1=str_remove_all(V1, pattern = '-+$')) %>% mutate(length=nchar(V1))

B <- sapply(A$V1, function(x) {substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))}) %>% t()

apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()

# processing log ----------------------------------------------------------
# args$in_fasta_dir <- '/home/eric/Analysis/aiv/merge/0307/'
# args$out_fasta_dir <- '/home/eric/Analysis/aiv/merge/0307/'

log <- fread(paste('/home/eric/Analysis/aiv/merge/0307/', 'clean_log.csv', sep = ''), header = T) %>% as.data.frame()
log[is.na(log$Segment), 'Segment'] <- 'NA'
log_sub <- fread(paste('/home/eric/Analysis/aiv/merge/0307/', 'clean_log_sub.csv', sep = ''), header = T) %>% as.data.frame()


log <- log %>% filter(Segment!='all') %>% tidyr::gather(key, value, 2:ncol(log)) %>% 
  mutate(key=factor(key, levels=unique(key)), value=value %>% as.numeric()
         , Segment=factor(Segment, levels=c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS')))

png('~/Analysis/aiv/merge/0307/result/log.png', width = 9600, height = 5400, res = 600)
ggplot()+
  # geom_bar(log %>% filter(Segment=='all'), mapping=aes(x=key, y=value * coeff, fill=Segment), stat='identity', size=3)+
  geom_point(log %>% filter(Segment!='all'), mapping=aes(x=key, y=value, color=Segment), size=3)+
  geom_line(log %>% filter(Segment!='all'), mapping=aes(x=key, y=value, group=Segment, color=Segment), linewidth=1.2)+
  ggrepel::geom_text_repel(data = log %>% filter(key=='raw')
                           , mapping=aes(x=key, y=value, label=Segment, color=Segment)
                           , size=6, show.legend = F, hjust=3, seed = 713,
                           point.padding = 0.4, 
                           nudge_x = -0.1,
                           nudge_y = 4,
                           segment.linetype = 4,
                           segment.curvature = 1e-5,
                           arrow = arrow(length = unit(0.015, "npc")))+
  # geom_hline(yintercept = yint, color='red', linetype='dashed')+
  # geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
  # geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
  # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  labs(x='', color='Segment', y='Number of sequences')+
  # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
  # scale_y_continuous(name = "length", sec.axis = sec_axis(~.*coeff, name="log10(Number of stop codons)"))+
  # scale_x_continuous( breaks=  c(seq(40, 95, 5), 100), labels =  c(seq(40, 95, 5), 100))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16), legend.position="bottom")
dev.off()



log_sub <- log_sub %>% filter(Segment!='all') %>% tidyr::gather(key, value, 2:ncol(log_sub)) %>% 
  mutate(key=factor(key, levels=unique(key)), value=value %>% as.numeric()
         , Segment=factor(Segment, levels=c(seg_sub_level)))

png('~/Analysis/aiv/merge/0307/result/log_HA.png', width = 9600, height = 5400, res = 600)
ggplot()+
  # geom_bar(log %>% filter(Segment=='all'), mapping=aes(x=key, y=value * coeff, fill=Segment), stat='identity', size=3)+
  geom_point(log_sub %>% filter(Segment %in% paste('H', 1:16, sep = '')), mapping=aes(x=key, y=value, color=Segment), size=3)+
  geom_line(log_sub %>% filter(Segment %in% paste('H', 1:16, sep = '')), mapping=aes(x=key, y=value, group=Segment, color=Segment), linewidth=1.2)+
  # ggrepel::geom_text_repel(data = log_sub %>% filter(key=='explicit subtype', Segment %in% paste('H', 1:16, sep = ''))
  #                          , mapping=aes(x=key, y=value, label=Segment, color=Segment)
  #                          , size=6, show.legend = F, hjust=3, seed = 713,
  #                          point.padding = 0.4, 
  #                          nudge_x = -0.1,
  #                          nudge_y = 4,
  #                          segment.linetype = 4,
  #                          segment.curvature = 1e-5,
  #                          arrow = arrow(length = unit(0.015, "npc")))+
  # geom_hline(yintercept = yint, color='red', linetype='dashed')+
  # geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
  # geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
  # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  labs(x='', title = 'Processing log', color='Segment', y='Number of sequences')+
  # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
  # scale_y_continuous(name = "length", sec.axis = sec_axis(~.*coeff, name="log10(Number of stop codons)"))+
  # scale_x_continuous( breaks=  c(seq(40, 95, 5), 100), labels =  c(seq(40, 95, 5), 100))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16))
dev.off()



png('~/Analysis/aiv/merge/0307/result/log_NA.png', width = 9600, height = 5400, res = 600)
ggplot()+
  # geom_bar(log %>% filter(Segment=='all'), mapping=aes(x=key, y=value * coeff, fill=Segment), stat='identity', size=3)+
  geom_point(log_sub %>% filter(Segment %in% paste('N', 1:9, sep = '')), mapping=aes(x=key, y=value, color=Segment), size=3)+
  geom_line(log_sub %>% filter(Segment %in% paste('N', 1:9, sep = '')), mapping=aes(x=key, y=value, group=Segment, color=Segment), linewidth=1.2)+
  # ggrepel::geom_text_repel(data = log_sub %>% filter(key=='explicit subtype', Segment %in% paste('H', 1:16, sep = ''))
  #                          , mapping=aes(x=key, y=value, label=Segment, color=Segment)
  #                          , size=6, show.legend = F, hjust=3, seed = 713,
  #                          point.padding = 0.4, 
  #                          nudge_x = -0.1,
  #                          nudge_y = 4,
  #                          segment.linetype = 4,
  #                          segment.curvature = 1e-5,
  #                          arrow = arrow(length = unit(0.015, "npc")))+
  # geom_hline(yintercept = yint, color='red', linetype='dashed')+
  # geom_text(aes(x=40, y=(yint-100), label= paste('Red dashed line Median ORF length:', yint)), hjust=-0.5, size=5, color='red')+
  # geom_text(aes(x=70, y=(log10(stop+1)*coeff)+50, label= paste('Stop codons:', stop)), size=5, color='red')+
  # geom_point(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  # geom_line(B[1:14, ], mapping=aes(x=key, y=log(value+1), group=key1, color=key1))+
  labs(x='', title = 'Processing log', color='Segment', y='Number of sequences')+
  # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
  # scale_y_continuous(name = "length", sec.axis = sec_axis(~.*coeff, name="log10(Number of stop codons)"))+
  # scale_x_continuous( breaks=  c(seq(40, 95, 5), 100), labels =  c(seq(40, 95, 5), 100))+
  theme_bw()+gg_theme+theme(axis.text.x = element_text(size=16))
dev.off()


# NS length investigation -------------------------------------------------
# same processing procedure gain the high stop codon
# args$in_fasta_dir <- '/home/eric/Analysis/aiv/merge/0307/aligned_fasta/'
# args$out_fasta_dir <- '/home/eric/Analysis/aiv/merge/0307/ORF_fasta/'

aligned <- read_as_list(path ='/home/eric/Analysis/aiv/merge/0307/aligned_fasta/', prefix = 'aligned.fa', file_type = 'fasta')
names(aligned) <- str_extract(names(aligned), pattern = '[A-Z]+[0-9]+|[A-Z]{2}_a') %>% str_remove(pattern = '_a')

# length mutation (under 86 % gap threshold)
stop_coord_list <- list()
stop_seq_list <- list()
for(j in seg_sub_level) {
  cat(paste('Start processing ', j, '\n', sep=''))
  A <- aligned[[j]] %>% do.call(what=rbind) %>% as.data.frame() 
  
  gap <- apply(A[, 1:ncol(A)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% 
    as.data.frame() %>% mutate(index=row_number(), prop=./nrow(A)*100)
  
  gap <- gap[gap$prop <= 86, ]
  A <- A[, gap$index]
  
  cons <- c()
  for(i in 1:ncol(A))  {
    test <- table(A[, i])
    cons <- c(cons, which.max(test) %>% names())
  }
  atg1 <- (paste(cons, collapse = '') %>% gregexpr(pattern = 'atg') %>% unlist())[1]
  
  test <- apply(A[1:nrow(A), atg1:ncol(A)], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
  s <- seq(1, nchar(test[1]), 3)
  e <- seq(3, nchar(test[1]), 3)
  s <- s[1:min(length(s), length(e))]; e <- e[1:min(length(s), length(e))]
  B <- sapply(test, function(x) {substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))}) %>% t()
  stop_coord_list[[j]] <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
  stop_seq_list[[j]] <- apply(B, MARGIN = 1, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
}

for(i in names(stop_coord_list)) {
  stop_coord_list[[i]] <- (stop_coord_list[[i]]/length(aligned[[i]]))*100
}

stop_coord_plot <- do.call(rbind, stop_coord_list) %>% tibble::rownames_to_column('seg') %>% rename('prop'='.') %>% as.data.table()
stop_coord_plot[, c('seg', 'coord')] <- header_cleaning(stop_coord_plot$seg, pattern = '\\.')

plot_list <- list()
for(i in names(class)) {
  if(i == 'N') {title <- 'NA'}
  else{title <- i}
  plot_data <- stop_coord_plot %>% filter(seg %in% class[[i]]) %>% 
    mutate(seg=factor(seg, levels=class[[i]]), 
           coord=as.numeric(coord))
  plot_list[[i]] <-
    ggplot(plot_data) +
    geom_tile(aes(x = coord %>% as.numeric(), y = seg, color = prop, fill = prop), 
              color = "grey", alpha = 0.8, height = 0.75) +
    scale_fill_gradientn(
      colors = c("grey", "white", "blue"),  # Example color scheme
      limits = c(0, 100),  # Ensure fill scale goes from 0 to 100
      breaks = c(0, 25, 50, 75, 100),  # Explicitly define the breaks
      labels = c("0%", "25%", "50%", "75%", "100%"))+
    scale_x_continuous(breaks = c(1, seq(100, max(plot_data$coord), 100), max(plot_data$coord)), 
                       labels = c(1, seq(100, max(plot_data$coord), 100), max(plot_data$coord)))+
    labs(x = "Amino Acid Coordinates", y = "Gene Segment", fill = "Stop Codon (%)", 
         title = paste('Stop Codon Usage Across ', title, ' Gene Segments', sep = '')) +
    theme_bw() +gg_theme + theme(legend.position = "bottom", axis.text.x = element_text(size=16))
  
    # ggplot()+
    # geom_point(mapping=aes(x=coord %>% as.numeric(), y=prop, color=seg), size=2)+
    # geom_line(mapping=aes(x=coord %>% as.numeric(), y=prop, group=seg, color=seg), linewidth=0.5)+
    # geom_hline(yintercept = 100, color='red', linetype='dashed')+
    # labs(x='Amino acid coordinates', y='Stop codon proportion', title = i, color='')+
    # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
    # scale_y_continuous(breaks=  c(seq(0, 100, 5)), labels =  c(seq(0, 100, 5)), limits = c(0, 100))+
    # theme_bw()+gg_theme+theme(legend.position = 'none')
}


dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_stop_codon_coord.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}


stop_seq_plot <- lapply(stop_seq_list, function(x){table_DF(x[,1])}) # %>% do.call(what=rbind) %>% tibble::rownames_to_column('seg') %>% as.data.table()
for(i in names(stop_seq_plot)) {
  A <- stop_seq_plot[[i]] %>% mutate(x=as.numeric(x))
  total <- sum(A$Freq)
  stop_seq_plot[[i]] <- rbind(A[A$x <5, ], data.frame(x='≧5', Freq=A[A$x >=5, 'Freq'] %>% sum())) %>% mutate(Freq=Freq/total*100)
}
stop_seq_plot <- stop_seq_plot %>% do.call(what=rbind) %>% tibble::rownames_to_column('seg') %>% as.data.table() 
stop_seq_plot[, c('seg', 'coord')] <- header_cleaning(stop_seq_plot$seg, pattern = '\\.') 
stop_seq_plot <- stop_seq_plot %>% filter(seg %in% class$Internal) %>% 
  mutate(x=factor(x, levels=rev(c(0:4, '≧5'))), seg=factor(seg, levels=seg_sub_level))


png(paste('~/Analysis/aiv/merge/0307/result/', 'stop_codon_sequence.png', sep = ''), width = 9600, height = 5400, res=600)
stop_seq_plot %>%
  ggplot()+
  geom_bar(mapping=aes(x=seg, y=Freq, fill=x), stat='identity', position='stack')+
  labs(y='Stop codon frequency', title = 'Stop codon frequency in segments', fill='')+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+gg_theme+theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16))
dev.off()

# B <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
# stop1 <- which.max(B$.)
# test <- A[, atg1:(stop1*3+atg1-1-3)]
# test$result <- apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
# test$epi <- str_extract(rownames(test), pattern = 'EPI_ISL_[0-9]+|IRD_[0-9]+')
# out <- matrix(0,0,0) %>% as.data.frame()
# out[seq(1, nrow(test)*2, 2), 1] <- paste('>', rownames(test), sep = '')
# out[seq(2, nrow(test)*2, 2), 1] <- test$result
# write.table(out, file = paste('~/Analysis/aiv/merge/0307/', j, 'NS_test.fa', sep = ''), sep = '\n'
#             , quote = F, col.names = F, row.names = 


# alternative splicing site (under 90 % gap threshold)



# NS2 identification by each sequence ------------------------------------------------------
j="gisaid_IRD_merged_NS_aligned"
cat(paste('Start processing ', j, '\n', sep=''))
A <- aligned[[j]]
header <- names(A)
A <- A %>% do.call(what=rbind) %>% as.data.frame() 
rownames(A) <- header

gap <- apply(A[, 1:ncol(A)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% 
  as.data.frame() %>% mutate(index=row_number(), prop=./nrow(A)*100)
gap <- gap[gap$prop < 90, ]

A <- A[, gap$index]


cons <- c()
for(i in 1:ncol(A))  {
  test <- table(A[, i])
  cons <- c(cons, which.max(test) %>% names())
}
atg1 <- (paste(cons, collapse = '') %>% gregexpr(pattern = 'atg') %>% unlist())[1]

NS2_1end <- (((paste(cons, collapse = '') %>% gregexpr(pattern = 'gtaga') %>% unlist())[1])+1)

NS2_2start <- ((paste(cons, collapse = '') %>% gregexpr(pattern = 'ccaggaca')%>% unlist())[1]+2)

NS2_1 <- apply(A[1:nrow(A), c(atg1:NS2_1end)], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')}) %>% 
  sapply(function(x) {substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))}) %>% t() %>% as.data.frame()
NS2_2 <- apply(A[1:nrow(A), c(NS2_2start:ncol(A))], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')}) %>% 
  sapply(function(x) {substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))}) %>% t() %>% as.data.frame()

paste(NS2_1$V32, NS2_1$V33, sep = '') %>% table_DF() %>% 
  ggplot()+
  geom_bar(mapping=aes(x=x, y=Freq), stat='identity', fill='#93C5F3')+
  geom_text(mapping=aes(x=x, y=Freq, label=Freq), stat='identity', hjust=0.5, vjust=-0.5, size=6)+
  labs(title = "5' terminal alternative splcing site", y='Frequency', x='nt combination')+
  theme_bw()+gg_theme


paste(NS2_2[, 1], NS2_2[, 2], sep = '') %>% table_DF() %>% 
  ggplot()+
  geom_bar(mapping=aes(x=x, y=Freq), stat='identity', fill='#93C5F3')+
  geom_text(mapping=aes(x=x, y=Freq, label=Freq), stat='identity', hjust=0.5, vjust=-0.5, size=6)+
  labs(title="3' terminal alternative splcing site", y='Frequency', x='nt combination')+
  theme_bw()+gg_theme

B <- rbind(paste(NS2_1$V32, NS2_1$V33, sep = '') %>% table_DF(prop = T) %>% mutate(g='NS2_1'), 
           paste(NS2_2[, 1], NS2_2[, 2], sep = '') %>% table_DF(prop = T) %>% mutate(g='NS2_2')) 
B %>% 
  ggplot()+
  geom_bar(mapping=aes(x=g, y=Freq*100, fill=x), stat='identity', position = 'stack')+
  # geom_text(mapping=aes(x=x, y=Freq, label=Freq), stat='identity', hjust=0.5, vjust=-0.5, size=6)+
  labs(x="3' terminal alternative splcing site", y='Frequency')+
  theme_bw()+gg_theme


# tree heatmap ------------------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame() %>% 
  mutate(Strain_number=as.character(Strain_number))

clade_table <- table_DF(meta[grepl(meta$Subtype, pattern='H5'), 'Clade'])
clade_table[, 'Clade'] <- 'GsGD-others'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4'), 'Clade'] <- '2.3.4.4'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4b'), 'Clade'] <- '2.3.4.4b'
clade_table[grep(clade_table$x, pattern = 'nonGsGD'), 'Clade'] <- 'nonGsGD'
clade_table[grep(clade_table$x, pattern = 'IRD_unlabel'), 'Clade'] <- 'IRD_unlabel'

tree_list <- readRDS('/home/eric/Analysis/aiv/merge/0307/iq_tree/internal_NA_tree.RData')

plot_list <- list()
# tree_list <- list()
for(i in c(class$Internal)) {
  i <- 'PB2'
  cat(paste('Loading tree\n', sep = ''))
  tree <- ape::read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i, '_aligned_iqtree.nexus', sep = ''))
  
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
  
  info <- merge(info, segment_groups[[i]], by='Strain_number', all=T)
  rownames(info) <- info$Strain_number
  info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
  info[, 18] <- factor(info[, 18], levels = info[, 18] %>% unique())
  
  # for(j in 18:ncol(info)) {
  #   group <- table_DF(info[, j]) %>% arrange(desc(Freq))
  #   if(nrow(group) > 15) {
  #     group <- group[1:15, ]
  #   }# else(group <- group %>% filter(Freq>100))
  #   rownames(group) <- group$x
  #   info[!info[, j] %in% group$x, j] <- NA
  #   remain_g <- info[, j] %>% unique()
  #   group <- group[remain_g[!is.na(remain_g)], ]
  #   
  #   # info[!info$g %in% group$x, 'clade'] <- 'Minor'
  #   
  # }
  
  
  rownames(info) <- info$Strain_number %>% as.character()
  info <- info[rank, ] # %>% mutate(index=seq(1, nrow(info)))
  info$Strain_number <- factor(info$Strain_number, levels =rev(rank))
  # info$data_source <- str_remove_all(info$Accesion_number, pattern = '_[0-9]+')
  # info[!(info$data_source %in% 'IRD'), 'data_source'] <- 'GISAID'
  
  #change tree header to match the data
  new_header <- info$Strain_number %>% as.character()
  tree$tip.label <- new_header
  #
  
  
  # phylo_raw <- ggtree(tree) 
  # tree_list[[i]] <- phylo_raw
  phylo_raw <- tree_list[[i]]
  clade <- info[, c(1, ncol(info))]
  phylo <- phylo_raw %<+% clade +
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
    labs(x="Serotype", fill="Serotype", color="Serotype") + ylab(NULL)+
    # scale_color_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
    # scale_fill_manual(values = c('MPN1'='#FFD97D','MPN2'='#FF8F54','MPN6'='#FF5356','MPN8'='#B63B3D','H7N9'='#71CE7B','H9N2'='#B9A4ED', 'Other subtype'='#EBF2F7'))+
    scale_color_brewer(palette = 'Set2')+
    scale_fill_brewer(palette = 'Set2')+
    guides(fill = guide_legend(nrow = 1, direction = "vertical"))+
    phylo_theme+theme(legend.position = "bottom") 
  
  
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
  plot_list[[i]] <-
  combined_plot + labs(title = paste(i, "segment phylogenetic tree (iqtree)")) +
    theme(plot.title = element_text(face = "bold", size = 25, hjust=2.5))
  
}


for(i in names(plot_list)) {
  png(paste('/home/eric/Analysis/aiv/merge/0307/result/', i, '_iqtree_lineage_global.png', sep = ''), width=9600, height=5400, res = 600)
  plot_list[[i]] %>% print()
  dev.off()
}
# saveRDS(tree_list, '/home/eric/Analysis/aiv/merge/0307/iq_tree/internal_NA_tree.RData')


# mean partristic distance calculation -------------------------------------
mpd <- data.frame(mpd=1:100)
for(i in class[[3]]) {
  A <- data.frame(mpd=readRDS(paste("~/Analysis/aiv/merge/0307/distance/iq/", i, "_iq_patristic_group.RData", sep = ''))[[1]] %>% unlist,
                  seg=i) %>% mutate(mpd=mpd*100, cut=cut(mpd, c(seq(0, 10, 2.5), 15, 20, 100), right = F))
  A <- hist(x = A[, 1], breaks = round(A$mpd, 0) %>% max(), plot = F)
  mpd <- merge(mpd, data.frame(mpd=A$breaks[-(length(A$breaks))], count=A$counts), by = 'mpd', all = T)
  colnames(mpd)[ncol(mpd)] <- i
}

mpd <- replace(mpd, is.na(mpd), 0)
mpd <- mpd[c(which(rowSums(mpd[,-1]) !=0)), ]
mpd[, -1] <- apply(mpd[, -1], MARGIN = 2, function(x){x/sum(x)*100})


Heatmap(mpd, colorRamp2(c(0, (min(mpd)+0.0000001), 30, 100),  c('grey50',  "#20BCE4", '#E7F8FC',"#E44820"))
        , show_column_dend = F
        , column_names_rot = 45, cluster_columns = F
        , column_title = "sgmRNA subclass"
        , column_title_gp = gpar(fontsize = 32, fontface = "bold")
        , column_names_gp = gpar(fontsize = 18)
        , column_order = (colnames(mpd))
        , cluster_rows = F
        , show_row_dend = T, row_names_side = 'left'
        , row_names_gp = gpar(fontsize = 18)
        , heatmap_legend_param= list(title= 'Log percentage\n'
                                     , labels_gp = gpar(fontsize = 20)
                                     , title_gp = gpar(fontsize = 24)
                                     , legend_width = unit(4, "cm"))
        , width = ncol(mpd)*unit(9, "mm")
        , height = nrow(mpd)*unit(9, "mm"))

mpd %>% do.call(what=rbind) %>% as.data.table() %>% 
  ggplot()+
  # geom_boxplot(aes(x=seg, y=mpd*100%>% log10()))+
  geom_bar(aes(x=seg, fill =cut))+
  labs(x='Gap threshold (%)', color='', y='identified open reading frame length')+
  # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
  # scale_y_continuous(name = "length", sec.axis = sec_axis(~./coeff, name="log10(Number of stop codons)"))+
  # scale_y_continuous(breaks= c(seq(0, 100, 10)), labels =  c(seq(0, 100, 10)), limits = c(0, 100))+
  theme_bw()+gg_theme+theme(legend.text = element_text(size = 16), legend.key.size = unit(2.5, 'lines'))

# Agreement and lineage numbers among fixed threshold ---------------------
plot_list <- list()
for(i in c(class$Internal, class$N)) {
  i <- 'N9'
  A <- readRDS(paste('/home/eric/Analysis/aiv/merge/0307/distance/agree_pd/', i, '_agreement_result.RData', sep = ''))
  
  agree_thres <- do.call(rbind, A$agreement_table) %>% as.data.frame() %>% filter(group != 'bs') %>% 
    mutate(dis=str_remove(dis, pattern = 'd_'), 
           group=str_replace(group, pattern='_groups', replacement='_lineages'),
           group=str_replace(group, pattern='u', replacement='AS'))
  
  coeff <- 1/max(agree_thres$y)
  print(max(agree_thres$y))
  # plot_list[[i]] <-
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_Agreement and lineage numbers among fixed threshold.png', sep = ''), width = 9600, height = 5400, res=600)
  ggplot()+
    geom_bar(agree_thres %>% filter(group %in% c('ft_lineages', 'iq_lineages')), 
             mapping=aes(x=dis, y=y, fill=group), position = 'dodge', stat = 'identity')+
    geom_point(agree_thres %>% filter(group %in% c('AS', 'AS-same', 'AS-distinct')), 
               mapping=aes(x=dis, y=y/coeff, color=group), size=4)+
    geom_line(agree_thres %>% filter(group %in% c('AS', 'AS-same', 'AS-distinct')), 
              mapping=aes(x=dis, y=y/coeff, color=group, group=group), linewidth=1.5)+
    # geom_point(agree_thres %>% filter(group %in% c('ft_groups', 'iq_groups')), mapping=aes(x=dis, y=y/coeff, color=group))+
    scale_y_continuous(name = "Number of lineages",
                       sec.axis = sec_axis(transform=~.*coeff, name="Agreement statistic"))+
    scale_color_manual(values = c('#FF8F54','#5AA977','#8864e0'))+
    labs(x='Fixed MPD threshold', fill='Lineage numbers', color='Agreement', 
         title = paste(i, ', Agreement and lineage numbers among fixed threshold', sep = ''))+
    theme_bw()+gg_theme
  dev.off()
  # plot_list[[i]]
  rm(agree_thres)
}


dir_create('~/Analysis/aiv/merge/0307/result/')
for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_Agreement and lineage numbers among fixed threshold.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}


# MPD percentile removed nodes --------------------------------------------
ft_node_list <- read_as_list(path = '~/Analysis/aiv/merge/0307/distance/ft/', prefix = '_ft_patristic_group_bympd_percentile.RData', file_type = 'RData')
names(ft_node_list) <- str_remove_all(names(ft_node_list), pattern = '_.*')
ft_node_list <- ft_node_list %>% lapply(FUN = function(x){x <- x[[1]]})
ft_node_list <- ft_node_list[which(!names(ft_node_list) %in% 'H5')]


iq_node_list <- read_as_list(path = '~/Analysis/aiv/merge/0307/distance/iq/', prefix = '_iq_patristic_group_bympd_percentile.RData', file_type = 'RData')
names(iq_node_list) <- str_remove_all(names(iq_node_list), pattern = '_.*')
iq_node_list <- iq_node_list %>% lapply(FUN = function(x){x <- x[[1]]})
iq_node_list <- iq_node_list[which(!names(iq_node_list) %in% 'H5')]

lapply(ft_node_list, function(x){length(x)}) %>% unlist
lapply(iq_node_list, function(x){length(x)}) %>% unlist

# if(nrow(mean_all) < 2500) {percentile <- c(seq(0.0005, 0.01, 0.0005))}
# if(2500 < nrow(mean_all) & nrow(mean_all) < 20000) {percentile <- seq(0.00025, 0.005, 0.00025)}
# if(20000 < nrow(mean_all) ) {percentile <- c(seq(0.00005, 0.001, 0.00005))}

small <- which((lapply(iq_node_list, function(x){length(x)}) %>% unlist) < 2500) %>% names
median <- which(((2500 < lapply(iq_node_list, function(x){length(x)}) & 
                    lapply(iq_node_list, function(x){length(x)}) < 20000) %>% unlist)) %>% names
huge <- which((20000 < lapply(iq_node_list, function(x){length(x)}))) %>% names

# small
plot_list <- list()
for(i in small) {
  remove_list <- list()
  for(j in c(seq(0.0005, 0.01, 0.0005))) { #c(seq(0.00005, 0.001, 0.00005)) seq(0.0005, 0.01, 0.0005) seq(0.00025, 0.005, 0.00025)
    remove <- data.frame(seg=i, ft=0, iq=0)
    rownames(remove) <- i
    
    mean_ft <- data.frame(index=1:length(ft_node_list[[i]]), MPD= ft_node_list[[i]] %>% unlist())
    mean_iq <- data.frame(index=1:length(iq_node_list[[i]]), MPD= iq_node_list[[i]] %>% unlist())
    
    # total['PB2', c('ft', 'iq')] <- c(nrow(mean_ft), nrow(mean_iq))
    remove[i, c('ft', 'iq')] <- c(nrow(mean_ft %>% filter(MPD > quantile(mean_ft$MPD, 1-j))), 
                                  nrow(mean_iq %>% filter(MPD > quantile(mean_iq$MPD, 1-j))))
    
    # total_p <- total %>% tidyr::gather(key, value, 2:3) %>% mutate(key=factor(key, levels=c('ft', 'iq')))
    remove_list[[as.character(1-j)]] <- remove %>% tidyr::gather(key, value, 2:3)  %>% mutate(key=factor(key, levels=c('ft', 'iq')))
  }
  
  remove_list[['total']] <- data.frame(seg=i, key=c('ft', 'iq'), value=c(nrow(mean_ft), nrow(mean_iq)))
  remove <- remove_list %>% do.call(what=rbind) %>% as.data.frame() %>% 
    tibble::rownames_to_column('dis') %>% mutate(dis=(dis %>% str_remove(pattern = '\\.[1-2]$')))
  remove[1:40, 1] <- as.numeric(remove[1:40, 1])*100
  
  plot_list[[i]] <-
    ggplot(remove)+
    geom_col(aes(x=dis, y=value, fill=key), position='dodge', color='black')+
    # scale_x_continuous(breaks=unique(remove$dis), labels=unique(remove$dis))+
    scale_y_log10() +
    labs(x='MPD percentile (%)', y='Removed nodes', fill='Method', title=paste(i, ' percentile removed nodes', sep = ''))+
    theme_bw()+gg_theme+theme(axis.text.x = element_text(size=14))
  
}
for(i in median) {
  remove_list <- list()
  for(j in c(seq(0.00025, 0.005, 0.00025))) { #c(seq(0.00005, 0.001, 0.00005)) seq(0.0005, 0.01, 0.0005) seq(0.00025, 0.005, 0.00025)
    remove <- data.frame(seg=i, ft=0, iq=0)
    rownames(remove) <- i
    
    mean_ft <- data.frame(index=1:length(ft_node_list[[i]]), MPD= ft_node_list[[i]] %>% unlist())
    mean_iq <- data.frame(index=1:length(iq_node_list[[i]]), MPD= iq_node_list[[i]] %>% unlist())
    
    # total['PB2', c('ft', 'iq')] <- c(nrow(mean_ft), nrow(mean_iq))
    remove[i, c('ft', 'iq')] <- c(nrow(mean_ft %>% filter(MPD > quantile(mean_ft$MPD, 1-j))), 
                                  nrow(mean_iq %>% filter(MPD > quantile(mean_iq$MPD, 1-j))))
    
    # total_p <- total %>% tidyr::gather(key, value, 2:3) %>% mutate(key=factor(key, levels=c('ft', 'iq')))
    remove_list[[as.character(1-j)]] <- remove %>% tidyr::gather(key, value, 2:3)  %>% mutate(key=factor(key, levels=c('ft', 'iq')))
  }
  
  remove_list[['total']] <- data.frame(seg=i, key=c('ft', 'iq'), value=c(nrow(mean_ft), nrow(mean_iq)))
  remove <- remove_list %>% do.call(what=rbind) %>% as.data.frame() %>% 
    tibble::rownames_to_column('dis') %>% mutate(dis=(dis %>% str_remove(pattern = '\\.[1-2]$')))
  remove[1:40, 1] <- as.numeric(remove[1:40, 1])*100
  
  plot_list[[i]] <-
    ggplot(remove)+
    geom_col(aes(x=dis, y=value, fill=key), position='dodge', color='black')+
    # scale_x_continuous(breaks=unique(remove$dis), labels=unique(remove$dis))+
    scale_y_log10() +
    labs(x='MPD percentile (%)', y='Removed nodes', fill='Method', title=paste(i, ' percentile removed nodes', sep = ''))+
    theme_bw()+gg_theme+theme(axis.text.x = element_text(size=14))
  
}
for(i in huge) {
  remove_list <- list()
  for(j in c(seq(0.00005, 0.001, 0.00005))) { #c(seq(0.00005, 0.001, 0.00005)) seq(0.0005, 0.01, 0.0005) seq(0.00025, 0.005, 0.00025)
    remove <- data.frame(seg=i, ft=0, iq=0)
    rownames(remove) <- i
    
    mean_ft <- data.frame(index=1:length(ft_node_list[[i]]), MPD= ft_node_list[[i]] %>% unlist())
    mean_iq <- data.frame(index=1:length(iq_node_list[[i]]), MPD= iq_node_list[[i]] %>% unlist())
    
    # total['PB2', c('ft', 'iq')] <- c(nrow(mean_ft), nrow(mean_iq))
    remove[i, c('ft', 'iq')] <- c(nrow(mean_ft %>% filter(MPD > quantile(mean_ft$MPD, 1-j))), 
                                  nrow(mean_iq %>% filter(MPD > quantile(mean_iq$MPD, 1-j))))
    
    # total_p <- total %>% tidyr::gather(key, value, 2:3) %>% mutate(key=factor(key, levels=c('ft', 'iq')))
    remove_list[[as.character(1-j)]] <- remove %>% tidyr::gather(key, value, 2:3)  %>% mutate(key=factor(key, levels=c('ft', 'iq')))
  }
  
  remove_list[['total']] <- data.frame(seg=i, key=c('ft', 'iq'), value=c(nrow(mean_ft), nrow(mean_iq)))
  remove <- remove_list %>% do.call(what=rbind) %>% as.data.frame() %>% 
    tibble::rownames_to_column('dis') %>% mutate(dis=(dis %>% str_remove(pattern = '\\.[1-2]$')))
  remove[1:40, 1] <- as.numeric(remove[1:40, 1])*100
  
  plot_list[[i]] <-
    ggplot(remove)+
    geom_col(aes(x=dis, y=value, fill=key), position='dodge', color='black')+
    # scale_x_continuous(breaks=unique(remove$dis), labels=unique(remove$dis))+
    scale_y_log10() +
    labs(x='MPD percentile (%)', y='Removed nodes', fill='Method', title=paste(i, ' percentile removed nodes', sep = ''))+
    theme_bw()+gg_theme+theme(axis.text.x = element_text(size=14))
  
}

for(j in names(plot_list)) {
  png(filename = paste('/home/eric/Analysis/aiv/merge/0307/result/', j, 'removed_node_percentile.png', sep = ''), width = 9600, height = 5400, res = 600)
  plot_list[[j]] %>% print()
  dev.off()
}


# MPD clustering sequence numbers  ----------------------------------------
# iq_grouping <- '~/Analysis/aiv/merge/0307/distance/iq/MP_iq_patristic_group.RData'
# ft_grouping <- '~/Analysis/aiv/merge/0307/distance/ft/MP_ft_patristic_group.RData'
# out <- readRDS(ft_grouping)[[2]]
# 
# rank <- colnames(out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V2=as.numeric(V2), name=paste(V1, V2, sep = '_')) %>%
#   arrange(V2) %>% select(name) %>% unlist()
# 
# out <- out[, c('Strain_number', rank)]
# if(all(out[, ncol(out)] == out[, (ncol(out)-1)])) {out <- out[, -ncol(out)]}
# if(colnames(out)[2] %in% 'd_0.0075') {
#   out <- out %>% select(-c())
# }
# 
# 
# 
# iqtree_out <- readRDS(iq_grouping)[[2]]
# rank <- colnames(iqtree_out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V2=as.numeric(V2), name=paste(V1, V2, sep = '_')) %>%
#   arrange(V2) %>% select(name) %>% unlist()
# 
# iqtree_out <- iqtree_out[, c('Strain_number', rank)]
# if(all(iqtree_out[, ncol(iqtree_out)] == iqtree_out[, (ncol(iqtree_out)-1)])) {iqtree_out <- iqtree_out[, -ncol(iqtree_out)]}
# if(colnames(iqtree_out)[2] %in% 'd_0.0075') {
#   iqtree_out <- iqtree_out %>% select(-c('d_0.0075'))
# }
# 
# 
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
# 
# colnames(FT_IQ_group) <- str_replace_all(string = colnames(FT_IQ_group), pattern = '\n', replacement = '_')

mpd_thres <- list(
  PB2=99.94,
  PB1=99.92,
  PA=99.93,
  NP=99.915,
  MP=99.96,
  NS=99.91,
  N1=99.85,
  N2=99.975,
  N3=99.9,
  N4=99.4,
  N5=99.2,
  N6=99.65,
  N7=99.45,
  N8=99.7,
  N9=99.5
)# %>% do.call(what=rbind)  %>% as.data.frame() %>% tibble::rownames_to_column('seg') %>% rename('thres'='V1')
meta_geno <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta_geno[is.na(meta_geno$segment), 'segment'] <- 'NA'
segment_groups <- list()
for(i in names(mpd_thres)) {
  colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
  iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group_bympd_percentile.RData', sep = ''))[[2]]
  
  colnames(iq_grouping)[-1] <- colnames(iq_grouping)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
    mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
  
  iq_grouping <- iq_grouping[, colN]
  
  group <- table_DF(iq_grouping[, 2])
  # segment_groups[[i]] <- group
  rownames(group) <- group$x
  
  remain_g <- iq_grouping[, 2] %>% unique()
  group <- group[remain_g[!is.na(remain_g)], ]
  group[, 'edit'] <- paste('g', 1:nrow(group), sep = '')
  colnames(iq_grouping)[2] <- 'clades'
  
  for(k in seq_len(nrow(group))) {
    iq_grouping[iq_grouping[, 2] %in% group[k, 'x'], 2] <- group[k, 3]
  }
  segment_groups[[i]] <- iq_grouping
  
  meta_geno <- merge(meta_geno, iq_grouping, by = 'Strain_number', all = T)
  colnames(meta_geno)[ncol(meta_geno)] <- paste(i, '_group', sep = '')
}


# for(i in c(class[c(1,3)] %>% unlist())) {
#   colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
#   png(paste('~/Analysis/aiv/merge/0307/result/', i, 'agreement.png'), width = 9600, height = 5400, res = 600)
#   readRDS(paste('~/Analysis/aiv/merge/0307/distance/', i, '_agreement_result.RData', sep = ''))[['agreement_plot']] %>% print()
#   dev.off()
# }


group_count <- list()
for(i in names(mpd_thres)) {
  colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
  iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group_bympd_percentile.RData', sep = ''))[[2]]
  
  colnames(iq_grouping)[-1] <- colnames(iq_grouping)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
    mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
  
  iq_grouping <- iq_grouping[, colN]
  
  group_count[[i]] <- table_DF(iq_grouping[,2]) %>% mutate(seg=i) %>% select(-c('x'))
}
test <- group_count %>% do.call(what=rbind) %>% as.data.frame()

png('~/Analysis/aiv/merge/0307/result/group_internal_count.png', width = 9600, height = 5400, res = 600)
test %>% filter(seg %in% class$Internal) %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  # geom_hline(plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color=='ORF'), 
  #            mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=3)+
  geom_violin(aes(x=seg, y=Freq ), fill='pink')+
  geom_jitter(aes(x=seg, y=Freq , color=seg), width = 0.25)+
  scale_color_brewer(palette = 'Set2')+
  # scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Group count', title='Group count distribution', fill='')+
  guides(color = guide_legend(nrow = 2, position = 'bottom', direction = ))+
  theme_bw()+gg_theme
dev.off()

png('~/Analysis/aiv/merge/0307/result/group_NA_count.png', width = 9600, height = 5400, res = 600)
test %>% filter(seg %in% class$N) %>% 
  ggplot()+
  # geom_boxplot(aes(x=factor(sub, level=seg_sub_level), y=nt_number, fill=g))+
  # geom_hline(plot_genome_l %>% filter(sub %in% c('PB2', 'PB1', 'PA', 'NP', 'MP', 'NS'), color=='ORF'), 
  #            mapping = aes(yintercept = nt_number), color='red', linetype='dashed', linewidth=3)+
  geom_violin(aes(x=seg, y=Freq ), fill='pink')+
  geom_jitter(aes(x=seg, y=Freq , color=seg), width = 0.25)+
  scale_color_brewer(palette = 'Set3')+
  # scale_y_continuous(breaks = seq(0, 3000, 250), limits = c(0, 3000))+
  labs(x='Internal genes', y='Group count', title='Group count distribution', fill='')+
  guides(color = guide_legend(nrow = 2, position = 'bottom', direction = ))+
  theme_bw()+gg_theme
dev.off()

group_count %>% do.call(what=rbind) %>% as.data.frame() %>% mutate(p=paste('seg', 'Freq'))
A <- data.frame(p=lapply(group_count, function(x){paste(c(x$seg %>% unique(), x$Freq), collapse = '_')}) %>% unlist())
A <- header_cleaning(A$p, pattern = '_')
colnames(A) <- c('seg', paste('g', seq(1, (ncol(A)-1)), sep = '_'))
write.csv(A, file = paste('~/Analysis/aiv/merge/0307/result/group_count.csv'), row.names = F)

within_dis_list <- list()
between_dis_list <- list()
heatmap_list <- list()
for(i in names(mpd_thres)) {
  rank <- read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i, '_aligned_iqtree.tree', sep = ''))$tip.label %>% 
    str_extract(pattern = '[0-9]+\\|H|[0-9]+_H') %>% str_remove(pattern = '\\|H|_H')
  colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
  iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group.RData', sep = ''))[[2]][, colN]
  colnames(iq_grouping)[2] <- 'groups'
  iq_grouping <- iq_grouping %>% group_split(groups) %>% as.list() %>% lapply(as.data.frame)
  dis <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_dist_matrix.RData', sep = ''))[[1]]
  colnames(dis) <- rownames(dis) <- str_extract(colnames(dis), pattern = '[0-9]+\\|H|[0-9]+_H') %>% str_remove(pattern = '\\|H|_H')
  dis <- dis[rank, rank]
  
  # heatmap_list[[i]] <-
  Heatmap(dis, colorRamp2(c(0, max(dis)/2, max(dis)),  c("#20BCE4", '#E7F8FC',"#E44820"))
          , show_column_dend = F, show_row_names = F
          , column_names_rot = 45, cluster_columns = F
          # , column_title = "Distance"
          , column_title_gp = gpar(fontsize = 32, fontface = "bold")
          , column_names_gp = gpar(fontsize = 18)
          , column_order = (colnames(dis))
          , cluster_rows = F, show_column_names = F
          , show_row_dend = T, row_names_side = 'left'
          , row_names_gp = gpar(fontsize = 18)
          , heatmap_legend_param= list(title= 'Pairwised\npatristic\ndistance\n'
                                       , labels_gp = gpar(fontsize = 16)
                                       , title_gp = gpar(fontsize = 20)
                                       , legend_width = unit(4, "cm"))
          , width = 2500*unit(0.1, "mm")
          , height = 2500*unit(0.1, "mm"))
  
  
  # within_dis <- c()
  # x <- 0
  # for(j in seq_along(iq_grouping)) {
  #   x <- x+1
  #   seq <- iq_grouping[[j]][, 1]
  #   sub_dis <- dis[seq, seq]
  #   n <- ncol(sub_dis)
  #   within_dis[x] <- sum(sub_dis)/(n^2/2-n)
  # }
  # within_dis_list[[i]] <- within_dis
  # 
  # between_dis <- c()
  # x <- 0
  # for(j in seq_along(iq_grouping)) {
  #   seq <- iq_grouping[[j]][, 1]
  #   for(k in seq_along(iq_grouping)) {
  #     if(k<=j) {next}
  #     x <- x+1
  #     seq1 <- iq_grouping[[k]][, 1]
  #     sub_dis <- dis[seq, seq1]
  #     n <- ncol(sub_dis)
  #     between_dis[x] <- mean(sub_dis)
  #   }
  # }
  # between_dis_list[[i]] <- between_dis
}


for(i in names(heatmap_list)) {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, 'heatmap.png'), width = 9600, height = 5400, res = 600)
  print(heatmap_list[[i]])
  dev.off()
}


hc_reault <- list()
for(i in names(mpd_thres)) {
  rank <- read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i, '_aligned_iqtree.tree', sep = ''))$tip.label %>% 
    str_extract(pattern = '[0-9]+\\|H|[0-9]+_H') %>% str_remove(pattern = '\\|H|_H')
  colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
  iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group.RData', sep = ''))[[2]][, colN]
  colnames(iq_grouping)[2] <- 'groups'
  # iq_grouping <- iq_grouping %>% group_split(groups) %>% as.list() %>% lapply(as.data.frame)
  dis <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', 'PB2', '_iq_patristic_dist_matrix.RData', sep = ''))[[1]]
  dis <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/ft/', 'PB2', '_ft_patristic_dist_matrix.RData', sep = ''))[[1]]
  colnames(dis) <- rownames(dis) <- str_extract(colnames(dis), pattern = '[0-9]+\\|H|[0-9]+_H') %>% str_remove(pattern = '\\|H|_H')
  dis <- dis[rank, rank]
  max(dis)
  # dis[ceiling(which.max(dis)/nrow(dis)), (which.max(dis)/nrow(dis)-ceiling(which.max(dis)/nrow(dis)-1))*nrow(dis)+1]
  
  hclust_avg <- hclust(dis %>% as.dist(), method = 'complete')
  gs <- iq_grouping[, 2] %>% unique() %>% length()
  cut_avg <- data.frame(hc=cutree(hclust_avg, k = gs)) %>% tibble::rownames_to_column('Strain_number')
  iq_grouping <- merge(iq_grouping, cut_avg, by = 'Strain_number')
  # table(iq_grouping$groups, iq_grouping$hc) %>% print()
  
  result <- list()
  for(j in 1:nrow(iq_grouping)) {
    # sub_ft <-
    iq_gg_mat <- ifelse(iq_grouping[j, 'groups']==iq_grouping[j:nrow(iq_grouping), 'groups'], yes = 1, no = 0)
    hc_gg_mat <- ifelse(iq_grouping[j, 'hc']==iq_grouping[j:nrow(iq_grouping), 'hc'], yes = 1, no = 0)
    
    result[[j]] <- (iq_gg_mat+hc_gg_mat) %>% table_DF() %>% mutate(group=paste(x, Freq, sep = '_')) %>%
      select(group) %>% unlist()
  }
  
  result <- result %>% unlist()
  n <- nrow(iq_grouping)
  grid <- (n*(n-1)/2)
  agree_diff <- (str_remove(result[grepl(pattern='0_', result, )], pattern = '0_') %>% as.numeric() %>% sum())
  disagree <- str_remove(result[grepl(pattern='1_', result, )], pattern = '1_') %>% as.numeric() %>% sum()
  agree_same <- (str_remove(result[grepl(pattern='2_', result, )], pattern = '2_') %>% as.numeric() %>% sum())-n
  
  agreement <- (agree_diff+agree_same)/grid
  agreement %>% print()
  # c(agreement, (agree_same)/grid, (agree_diff)/grid)
  
  hc_reault[[i]] <- iq_grouping
}


result <- data.frame(seg=class[c(1, 3)] %>% unlist) %>% mutate(iq_min=0, ft_min=0)
rownames(result) <- result$seg
for(i in rownames(result)) {
  cat(paste('doing', i, '\n'))
  result[i, 'iq_max'] <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_dist_matrix.RData', sep = ''))[[1]] %>% 
    max()
  result[i, 'ft_max'] <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/ft/', i, '_ft_patristic_dist_matrix.RData', sep = ''))[[1]] %>% max()
}


B <- result %>% tidyr::gather(key, value, 2:3)

result %>% tidyr::gather(key, value, 2:3) %>% 
  ggplot()+
  geom_point(mapping=aes(x=key, y=value, color=seg), size=3)+
  geom_line(mapping=aes(x=key, y=value, group=seg, color=seg), linewidth=1.2)+
  labs(x='', title = paste('Maximum partistic distance between FastTree and IQ-TREE'), y='partistic distance', color='')+
  # scale_color_manual(values = c('#958BE3', '#E3B88B', '#95CE7E'))+
  theme_bw()+gg_theme


# genetic distance --------------------------------------------------------
#intra group
# result_list <- list()
# for(i in c('H5', 'PB2', 'PB1', 'MP', 'PA', 'NP', 'NS', paste('N', seq(1, 9), sep = '')))  {
#   dis <- readRDS(paste('Analysis/aiv/merge/0307/distance/', i, '_dist_matrix.RData', sep = ''))
#   suppressWarnings(group <- read_as_list(paste('/home/eric/Analysis/aiv/merge/0307/group/groups/', i, '/', sep = ''), file_type = 'txt'))
#   for(j in names(group))  {
#     A <- group[[j]][2, 1]
#     A_EPI <- group[[j]][2, 1] %>% str_remove_all(pattern = 'tree tree_1 = [&R]') %>% strsplit(split='\'') %>% do.call(what=rbind) %>% 
#       t() %>% as.data.frame()
#     group[[j]] <- data.frame(epi=A_EPI[grep(pattern = '^EPI|^IRD', A_EPI$V1),], group=j) 
#   }
#   group <- do.call(rbind, group) %>% as.data.frame()
#   group$ui <- str_extract(group$epi, pattern = '.*\\|H[0-9]+') %>% str_remove(pattern = '\\|H[0-9]+')
#   result <- c()
#   for(j in unique(group$group)) {
#     n <- length(A)
#     result[j] <- (dis[group[group$group %in% j, 1], group[group$group %in% j, 1]] %>% sum())/(n*n-n) 
#   }
#   result_list[[i]] <- result
# }
# saveRDS(result_list, 'intra_dist.RData')
# 
# # inter-group
# result_list <- list()
# for(i in c('H5', 'PB2', 'PB1', 'MP', 'PA', 'NP', 'NS', paste('N', seq(1, 9), sep = '')))  {
#   dis <- readRDS(paste('Analysis/aiv/merge/0307/distance/', i, '_dist_matrix.RData', sep = ''))
#   suppressWarnings(group <- read_as_list(paste('/home/eric/Analysis/aiv/merge/0307/group/groups/', i, '/', sep = ''), file_type = 'txt'))
#   for(j in names(group))  {
#     A <- group[[j]][2, 1]
#     A_EPI <- group[[j]][2, 1] %>% str_remove_all(pattern = 'tree tree_1 = [&R]') %>% strsplit(split='\'') %>% do.call(what=rbind) %>% 
#       t() %>% as.data.frame()
#     group[[j]] <- data.frame(epi=A_EPI[grep(pattern = '^EPI|^IRD', A_EPI$V1),], group=j) 
#   }
#   group <- do.call(rbind, group) %>% as.data.frame()
#   group$ui <- str_extract(group$epi, pattern = '.*\\|H[0-9]+') %>% str_remove(pattern = '\\|H[0-9]+')
#   result <- c()
#   for(j in unique(group$group)) {
#     for(k in unique(group$group)) {
#       if(j!=k) {
#         sub_dis <- dis[group[group$group %in% j, 1], group[group$group %in% k, 1]]
#         result[paste(j, k, sep = '_')] <- mean(sub_dis)  
#       }
#     }
#   }
#   result_list[[i]] <- result
# }
# saveRDS(result_list, 'inter_dist.RData')
# 
# 
# dist <- readRDS('inter_dist.RData') %>% lapply(FUN = function(x) {x*100})
# result <- do.call(dist, what=rbind) %>% as.data.frame()
# for(i in seq_along(dist)) {
#   if(length(dist[[i]])==ncol(result)) {next}
#   else{result[i, (length(dist[[i]])+1):ncol(result)] <- ''}
# }
# result$m <- lapply(dist, mean) %>% unlist()
# intra <- t(result) %>% as.data.frame() %>% tibble::rownames_to_column('group') %>% tidyr::gather(key, value, c(2:17))
# intra <- intra[!intra$value %in% '', ] %>% 
#   mutate(int='Inter')
# intra$value <- as.numeric(intra$value)
# 
# 
# dist <- readRDS('intra_dist.RData') %>% lapply(FUN = function(x) {x*100})
# result <- do.call(dist, what=rbind) %>% as.data.frame()
# for(i in seq_along(dist)) {
#   if(length(dist[[i]])==ncol(result)) {next}
#   else{result[i, (length(dist[[i]])+1):ncol(result)] <- ''}
# }
# result$m <- lapply(dist, mean) %>% unlist()
# inter <- t(result) %>% as.data.frame() %>% tibble::rownames_to_column('group') %>% tidyr::gather(key, value, c(2:17))
# inter <- inter[!inter$value %in% '', ] %>% 
#   mutate(int='Intra')
# inter$value <- as.numeric(inter$value)
# 
# A <- rbind(intra, inter) 
# A$int <- factor(A$int, levels = c('Intra', 'Inter')) 
# 
# # png('~/Analysis/aiv/plot/boxplot_genetic_distance.png')
# plot_list <- list()
# for(i in unique(A$key)) {
#   # plot_list[[i]] <- 
#   png('~/Analysis/aiv/merge/0307/plot/boxplot_genetic_distance.png', width = 9600, height = 5400, res = 600)
#   ggplot(A %>% filter(group != 'm'))+
#     # geom_boxplot(aes(x=key, y=value, color=key))+
#     geom_boxplot(aes(x=key, y=value, fill=int),  position = position_dodge(1))+
#     geom_point(aes(x=key, y=value, fill=int),  position = position_dodge(1), size = 2)+
#     geom_text(data = A %>% filter(group %in% 'm', int == 'Intra')
#               , mapping=aes(x=key, y=-1, label=value %>% round(digits = 2)), color='#E68865'
#               , size=6, show.legend = F, vjust=1)+
#     geom_text(data = A %>% filter(group %in% 'm', int == 'Inter')
#               , mapping=aes(x=key, y=max((A %>% filter(group %in% 'm', int == 'Inter'))$value)
#                             , label=value %>% round(digits = 2)), color='#3BAEDA'
#               , size=6, show.legend = F, vjust=-1)+
#     scale_y_continuous(limits = c(-1, (max((A %>% filter(group %in% 'm', int == 'Inter'))$value))+2), 
#                        breaks =  c(0, seq(2.5, (max(A$value))+2, 2.5) ))+
#     scale_color_manual(values=c('#E68865', '#3BAEDA'), labels=c('Inter (within groups)', 'Intra (between groups)'))+
#     labs(x='Segment', y="Genetic distance (%)", color='Inter or\nIntra groups'
#          , title = "Genetic distance (Kimura's 2 parameter)", fill='')+
#     theme_bw()+gg_theme+theme(axis.text.x = element_text(size = 18))
#   dev.off()
# }
# 
# # dev.off()
# # for(i in unique(A$key)) {A$line <- }
# 
# # result_list[[1]]
# dir_create('~/Analysis/aiv/merge/0307/plot/')
# for(i in names(plot_list))  {
#   png(paste('~/Analysis/aiv/merge/0307/plot/', i, '_distance.png', sep = ''), width = 5400, height = 5400, res=600)
#   print(plot_list[[i]])
#   dev.off()
# }
# 
# 
# png('~/Analysis/aiv/merge/0307/plot/parallel_line_genetic_distance.png')
# ggplot(A %>% filter(group %in% 'm'))+
#   geom_point(aes(x=int, y=value, color=key))+
#   geom_line(aes(x=int, y=value, color=key, group=key))+
#   ggrepel::geom_text_repel(data = A %>% filter(group %in% 'm', int == 'Intra')
#                            , mapping=aes(x=int, y=value, label=key, color=key)
#                            , size=6, show.legend = F, hjust=3, 
#                            point.padding = 0.2, 
#                            nudge_x = -.2,
#                            nudge_y = .5,
#                            segment.linetype = 4,
#                            segment.curvature = 1e-20,
#                            arrow = arrow(length = unit(0.015, "npc")))+
#   ggrepel::geom_text_repel(data = A %>% filter(group %in% 'm', int == 'Inter')
#                            , mapping=aes(x=int, y=value, label=key, color=key)
#                            , size=6, show.legend = F, hjust=-3
#                            , point.padding = 0.2, 
#                            nudge_x = .15,
#                            nudge_y = .5,
#                            segment.linetype = 6,
#                            segment.curvature = -1e-20,
#                            arrow = arrow(length = unit(0.015, "npc")))+
#   scale_y_continuous(breaks = c(seq(0, 30, 2.5)))+
#   scale_x_discrete(labels=c('Intra (within groups)', 'Inter (between groups)'))+
#   labs(x='Intra or Inter groups', y="Kimura's 2 parameter genetic distance (%)", color='Segment')+
#   theme_linedraw()+gg_theme+theme(legend.key.size = unit(2.5, 'lines'), )
# dev.off()

intra_dist_list <- list()
inter_dist_list <- list()
wilcoxon_list <- list()
for(j in c(class$Internal)) {
  cat(paste('Doing ', j, ' segment\n', sep = ''))
  A <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/inter_intra_gd/', j, 
                     '_GTR_genetic_distance_aa_ident.RData', sep = ''))
  
  # intra_list <- A$intra
  intra_dist_list[[j]]  <- A$intra %>% unlist() %>% as.numeric()
  # names(intra_list) <- names(intra_list) %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
  #   mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
  # intra_dist_list[[j]] <- intra_list[[paste('d_', mpd_thres[[j]], sep = '')]]
  
  inter_dist_list[[j]]  <- A$inter %>% unlist() %>% as.numeric()
  # inter_list <- A$inter
  # names(inter_list) <- names(inter_list) %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
  #   mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
  # inter_dist_list[[j]] <- inter_list[[paste('d_', mpd_thres[[j]], sep = '')]]
  wilcoxon_list[[j]] <- wilcox.test(intra_dist_list[[j]], inter_dist_list[[j]], exact = FALSE, alternative = 'less')[['p.value']]
  # wilcoxon_list[[j]] <- t.test(intra_dist_list[[j]], inter_dist_list[[j]], exact = FALSE, alternative = 'less')[['p.value']]
}


dis_intra <- do.call(rbind, intra_dist_list) %>% as.data.frame()
for(i in seq_along(intra_dist_list)) {
  if(length(intra_dist_list[[i]])==ncol(dis_intra)) {next}
  else{dis_intra[i, (length(intra_dist_list[[i]])+1):ncol(dis_intra)] <- ''}
}
dis_intra$m <- lapply(intra_dist_list, mean) %>% unlist()
dis_intra <- t(dis_intra) %>% as.data.frame() %>% tibble::rownames_to_column('group')
dis_intra <- dis_intra %>% tidyr::gather(key, value, c(2:ncol(dis_intra)))
dis_intra <- dis_intra[!dis_intra$value %in% '', ] %>% mutate(int='Intra', value=as.numeric(value))


dis_inter <- do.call(rbind, inter_dist_list) %>% as.data.frame()
for(i in seq_along(inter_dist_list)) {
  if(length(inter_dist_list[[i]])==ncol(dis_inter)) {next}
  else{dis_inter[i, (length(inter_dist_list[[i]])+1):ncol(dis_inter)] <- ''}
}
dis_inter$m <- lapply(inter_dist_list, mean) %>% unlist()
dis_inter <- t(dis_inter) %>% as.data.frame() %>% tibble::rownames_to_column('group')
dis_inter <- dis_inter %>% tidyr::gather(key, value, c(2:ncol(dis_inter)))
dis_inter <- dis_inter[!dis_inter$value %in% '', ] %>% mutate(int='Inter', value=as.numeric(value))



A <- rbind(dis_intra, dis_inter) 
A$int <- factor(A$int, levels = c('Intra', 'Inter')) 
A$label <- paste(A$value %>% round(digits = 2), '%', sep = '')
pv <- do.call(rbind, wilcoxon_list) %>% as.data.frame() %>% rename('value'='V1') %>% tibble::rownames_to_column('key') %>% 
  mutate(group='p', int='p', label=value %>% round(digits = 3))
A <- rbind(A, pv[, colnames(A)])
A <- replace(A, is.na(A), 0)

# A$key <- (1-(str_remove(A$key, pattern = 'd_') %>% as.numeric())) %>% as.character()
# A$label <- paste(A$value %>% round(digits = 2), '%', sep = '')


# A$key <- factor(A$key, levels = rev(A$key %>% unique()))
# if(args$segment %in% 'PB2') {
#   A$key <- (1-(str_remove(A$key, pattern = 'd_') %>% as.numeric())) %>% as.character()
# }
# if(args$segment != 'PB2') {
#   A$key <- (str_remove(A$key, pattern = 'd_')) %>% as.character()
# }


plot_data <- A %>% filter(key %in% class$Internal)
plot_data$key <- factor(plot_data$key, levels = class$Internal)
p <- ggplot(plot_data %>% filter(!(group %in% c('m', 'p'))))+
  # geom_boxplot(aes(x=key, y=value, color=key))+
  geom_boxplot(aes(x=key, y=value, fill=int),  position = position_dodge(1))+
  geom_point(aes(x=key, y=value, fill=int),  position = position_dodge(1), size = 2)+
  # annotate(geom="text", x=0.5, y=-2.85, size=7, label="p-value: ", color="black")+
  # geom_text(data = A %>% filter(key %in% class$Internal, group %in% 'p', int == 'Inter')
  #           , mapping=aes(x=key, y=(max((A %>% filter(int == 'Inter'))$value))*1.1
  #                         , label=label), color='#3BAEDA'
  #           , size=4, show.legend = F)+ #0.6*6/length(unique(A$key))
  scale_y_continuous(limits = c(-3, (max((plot_data %>% filter(int == 'Inter'))$value))+10), 
                     breaks =  c(0, seq(5, (max(A$value))+2, 5) ))+
  scale_fill_manual(values=c('#E68865', '#3BAEDA'), labels=c('Intra (within lineage)', 'Inter (between lineage)'))+
  labs(x='', y="Genetic distance (%)", color='Inter or\nIntra lineage', fill='')+
  theme_bw()+gg_theme+theme(axis.title.x = element_blank(), legend.position = 'bottom')


x_positions <- seq_along(class$Internal)
y_positions <- c((plot_data %>% group_split(key) %>% lapply(FUN = function(x){max(x$value)}) %>% unlist())+5) # Adjust these based on your plot scale
p_value <- pv %>% filter(key %in% class$Internal) %>% mutate(key=factor(key, levels = class$Internal)) %>% select(value) %>% unlist

# Apply the function to all entries in your dataset
for (i in x_positions) {
  p <- add_p_annotation(p, p_value[i], x_positions[i], y_positions[i], size = 5)
}

# Print the final plot
print(p)

png(paste('~/Analysis/aiv/merge/0307/result/', 'Internal_thres_GTR_dis.png', sep = ''),
    width = 9600, height = 5400, res = 600)
p
dev.off()



plot_data <- A %>% filter(key %in% class$N)
plot_data$key <- factor(plot_data$key, levels = class$N)
p <- ggplot(plot_data %>% filter(!(group %in% c('m', 'p'))))+
  # geom_boxplot(aes(x=key, y=value, color=key))+
  geom_boxplot(aes(x=key, y=value, fill=int),  position = position_dodge(1))+
  geom_point(aes(x=key, y=value, fill=int),  position = position_dodge(1), size = 2)+
  # annotate(geom="text", x=0.5, y=-2.85, size=7, label="p-value: ", color="black")+
  # geom_text(data = A %>% filter(key %in% class$Internal, group %in% 'p', int == 'Inter')
  #           , mapping=aes(x=key, y=(max((A %>% filter(int == 'Inter'))$value))*1.1
  #                         , label=label), color='#3BAEDA'
  #           , size=4, show.legend = F)+ #0.6*6/length(unique(A$key))
  scale_y_continuous(limits = c(-3, (max((plot_data %>% filter(int == 'Inter'))$value))+10), 
                     breaks =  c(0, seq(5, (max(A$value))+2, 5) ))+
  scale_fill_manual(values=c('#E68865', '#3BAEDA'), labels=c('Intra (within groups)', 'Inter (between groups)'))+
  labs(x='', y="GTR Genetic distance (%)", color='Inter or\nIntra groups'
       , title = paste('NA genes', ", Genetic distance (GTR model)", sep = ''), fill='')+
  theme_bw()+gg_theme+theme(axis.title.x = element_blank(), legend.position = 'bottom')


x_positions <- seq_along(class$N)
y_positions <- c((plot_data %>% group_split(key) %>% lapply(FUN = function(x){max(x$value)}) %>% unlist())+5) # Adjust these based on your plot scale
p_value <- pv %>% filter(key %in% class$N) %>% mutate(key=factor(key, levels = class$N)) %>% select(value) %>% unlist

# Apply the function to all entries in your dataset
for (i in x_positions) {
  p <- add_p_annotation(p, p_value[i], x_positions[i], y_positions[i], size = 5)
}

# Print the final plot
print(p)

png(paste('~/Analysis/aiv/merge/0307/result/', 'NA_thres_GTR_dis.png', sep = ''),
    width = 9600, height = 5400, res = 600)
p
dev.off()

# groups_molecular_rate ---------------------------------------------------
plot_list <- list()
for(seg in c(class$N, class$Internal)) {
  iqtree_out <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', seg, '_iq_patristic_group_bympd_percentile.RData', sep = ''))[[2]]
  
  colnames(iqtree_out)[-1] <- colnames(iqtree_out)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
    mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
  
  iqtree_out <- iqtree_out[, c('Strain_number', paste('d_', mpd_thres[[seg]], sep = ''))]
  
  group <- table_DF(iqtree_out[, 2]) %>% arrange(desc(Freq))
  if(nrow(group) > 12) {
    group <- group[1:12, ]
  }# else(group <- group %>% filter(Freq>100))
  rownames(group) <- group$x
  iqtree_out[!iqtree_out[, 2] %in% group$x, 2] <- NA
  remain_g <- iqtree_out[, 2] %>% unique()
  group <- group[remain_g[!is.na(remain_g)], ]
  group[, 'edit'] <- paste('g', 1:nrow(group), sep = '')
  
  for(k in seq_len(nrow(group))) {
    iqtree_out[iqtree_out[, 2] %in% group[k, 'x'], 3] <- group[k, 3]
  }
  iqtree_out[, 3] <- factor(iqtree_out[, 3], levels = c(paste('g', 1:12, sep = ''), NA))
  
  
  A <- fread(paste('~/Analysis/aiv/merge/0307/timetree/', seg,'_clock_results/rtt.csv', sep = ''))
  root <- A[A$name %in% 'NODE_0000000', 2] %>% unlist()
  
  A$evolve_time <- (A$date-root)
  A$rate <- A$`root-to-tip distance`/A$evolve_time
  A <- A[A$rate >0 & A$rate <0.01, ]
  A <- A[c(grepl(A$name, pattern='NODE_0000000') %>% which, (!grepl(A$name, pattern='NODE')) %>% which), ]
  A$Strain_number <- A$name %>% str_extract(pattern = '_[0-9]+_H|\\|[0-9]+\\|H') %>% str_remove(pattern = '_|\\|') %>% 
    str_remove(pattern = '_H|\\|H')
  A <- merge(A, iqtree_out, by = 'Strain_number')
  
  colnames(A)[ncol(A)] <- 'color'
  
  p1 <- ggplot(A)+
    geom_point(aes(x=date, y=`root-to-tip distance`, color=color), size=1)+
    # geom_density(aes(mean_all$MPD))
    # scale_x_continuous(breaks = seq(min(plot_df$x), max(plot_df$x), 0.05), labels = seq(min(plot_df$x), max(plot_df$x), 0.05))+
    labs(x='Inferred year', y='Root to tip distance', title=paste(seg, ', Root to tip to time regrssion', sep = ''), 
         color='Clades')+
    theme_bw()+gg_theme
  
  
  # p2 <- ggplot(A)+
  #   # geom_histogram(aes(x=rate, y=(after_stat(count)/sum(after_stat(count)))*100), bins=100)+
  #   geom_density(aes(A$rate, y=(after_stat(count)/sum(after_stat(count)))*100, group=color, color=color), )+
  #   # scale_x_continuous(breaks = seq(min(plot_df$x), max(plot_df$x), 0.05), labels = seq(min(plot_df$x), max(plot_df$x), 0.05))+
  #   labs(x='Molecular rate (nt/year)', y='Percentage (%)', title='Molecular rate distribution of N6 sequences', color='groups')+
  #   theme_bw()+gg_theme
  
  
  library(ggridges)
  library(viridis)
  p2 <- 
    ggplot(A, aes(y = factor(color, levels = c(paste('g', 1:12, sep = '')) %>% rev), x= rate,  fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1, rel_min_height = 0.001) +
    scale_fill_viridis(name = "Molecular rate", option = "H") +
    labs(x='Molecular rate (nt/year)', y='Assigned clade', title=paste('Molecular rate distribution of', seg, 'sequences')
         , color='Clades')+
    # geom_vline(xintercept = 7, linetype="dotted", 
    #            color = "Red", size=1.5)+
    # scale_x_continuous(limits = c(0, 25), breaks = c(seq(0,25,5)))+
    # guides(fill=guide_legend(title=""))+
    theme_bw()+gg_theme+theme(title = element_text(size = 30-10),)
  
  
  plot_list[[seg]] <- ggarrange(p1, p2, common.legend = TRUE, legend='bottom')
  # annotate_figure(top = text_grob("", face = "bold", size = 24)) #, color = "red", face = "bold", size = 14)
}


for(i in names(plot_list))  {
  png(paste('~/Analysis/aiv/merge/0307/result/', i, '_molecular_rate.png', sep = ''), width = 9600, height = 5400, res=600)
  print(plot_list[[i]])
  dev.off()
}



# assign genotype 0621 ---------------------------------------------------------
# meta_geno <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
# meta_geno[is.na(meta_geno$segment), 'segment'] <- 'NA'
# segment_groups <- list()
# for(i in names(mpd_thres)) {
#   colN <- c('Strain_number', paste('d_', mpd_thres[i], sep = ''))
#   iq_grouping <- readRDS(paste('~/Analysis/aiv/merge/0307/distance/iq/', i, '_iq_patristic_group_bympd_percentile.RData', sep = ''))[[2]]
#   
#   colnames(iq_grouping)[-1] <- colnames(iq_grouping)[-1] %>% header_cleaning(pattern = '_') %>% mutate(V3=as.numeric(V2)) %>% 
#     mutate(V3=100*(1-V3), name=paste(V1, V3, sep = '_')) %>% select(name) %>% unlist()
#   
#   iq_grouping <- iq_grouping[, colN]
#   
#   group <- table_DF(iq_grouping[, 2])
#   # segment_groups[[i]] <- group
#   rownames(group) <- group$x
#   
#   remain_g <- iq_grouping[, 2] %>% unique()
#   group <- group[remain_g[!is.na(remain_g)], ]
#   group[, 'edit'] <- paste('g', 1:nrow(group), sep = '')
#   
#   
#   for(k in seq_len(nrow(group))) {
#     iq_grouping[iq_grouping[, 2] %in% group[k, 'x'], 2] <- group[k, 3]
#   }
#   colnames(iq_grouping)[2] <- 'group'
#   segment_groups[[i]] <- iq_grouping
#   meta_geno <- merge(meta_geno, iq_grouping, by = 'Strain_number', all = T)
#   colnames(meta_geno)[ncol(meta_geno)] <- paste(i, '_group', sep = '')
# }
# 

# meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
# meta[is.na(meta$segment), 'segment'] <- 'NA'
# meta$p <- paste(meta$Isolate_Id, meta$Location, meta$Collection_Date, sep = '_')
# eight <- list()
# # group_out
# for(i in c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'))  {
#   if((i %in% c('NA')))  {
#     p <- sub(pattern = 'A', replacement = '[0-9]+', x = i)
#     A <- segment_groups[grep(names(segment_groups), pattern=p)] %>% do.call(what=rbind) %>% as.data.frame() %>% 
#       tibble::rownames_to_column('subtype') %>% mutate(subtype=str_remove_all(subtype, pattern = '\\..*'), group=paste(subtype, group, sep='')) %>% 
#       select(-c(subtype))
#     sub_meta <- meta[meta$segment %in% i, ]
#     A <- merge(A, sub_meta %>% mutate(Strain_number=as.character(Strain_number)), by = 'Strain_number', all.x = T)
#     # if(identical(i, 'NA'))  {
#     #   A[, c('group', 's_g')] <- A[, c('s_g')]
#     # }
#   } 
#   if((i %in% c('HA')))  {
#     sub_meta <- meta[meta$segment %in% i, ]
#     A <- data.frame(Strain_number=sub_meta$Strain_number, group=sub_meta$sub)
#     A <- merge(A, sub_meta %>% mutate(Strain_number=as.character(Strain_number)), by = 'Strain_number', all.x = T)
#   } 
#   if(!(i %in% c('HA', 'NA'))) {
#     A <- segment_groups[[i]]
#     sub_meta <- meta[meta$segment %in% i, ]
#     A <- merge(A, sub_meta, by = 'Strain_number', all.x = T)
#   }
#   # A <- A %>% arrange(desc(Collection_Date)) %>% distinct(Isolate_Id, .keep_all = T) # remove duplicate epi
#   eight[[i]] <- A %>% select(-c('seq'))
# }
# 
# 
# genotype <- data.frame(p=meta$p %>% unique())
# for(i in names(eight))  {
#   A <- eight[[i]] %>% distinct(p, .keep_all = T)
#   colnames(A)[2] <- i
#   genotype <- merge(genotype, A[, c('p', i)], by = 'p', all.x=T)
#   genotype[, i] <- replace(genotype[, i], is.na(genotype[, i]), 'X')
# }
# genotype <- genotype[!is.na(genotype$p), ]
# 
# genotype$genotype <- paste(genotype$PB2, genotype$PB1, genotype$PA, genotype$HA,
#                            genotype$NP, genotype$`NA`, genotype$MP, genotype$NS, sep = '_')


A <- table_DF(genotype$genotype) %>% arrange(desc(Freq))
A <- A[!grepl(A$x, pattern='X'), ]
A[, 'HL'] <- 'LPAI'
A[grepl(header_cleaning(A$x, pattern = '_') %>% select(V4) %>% unlist, pattern='H5'), 'HL'] <- 'HPAI'

png('~/Analysis/aiv/merge/0307/result/genotype_top30_bar.png', width = 9600, height = 5400, res = 600)
A %>% top_n(30, Freq) %>% mutate(x=factor(x, levels=c(x))) %>% 
  ggplot()+
  geom_bar(aes(x=Freq, y=x, fill=HL), stat='identity', position='dodge')+
  geom_text(aes(x=Freq, y=x, label=Freq, color=HL), hjust=-0.5, fontface = "bold")+
  labs(fill='', x='', y='Genotype', title=paste('Genotype distribution (Top 30)'), color='')+
  scale_fill_manual(values=c('#66C8E6', '#E68466'))+
  scale_color_manual(values=c('#66C8E6', '#E68466'))+
  # scale_color_brewer(palette = 'Dark2')+
  theme_bw()+gg_theme+theme(axis.text.y = element_text(size = 18-4, family = "mono")) #, color=colors
dev.off()


meta_geno <- merge(genotype, meta %>% distinct(p, .keep_all = T), by = 'p')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()


A <- meta_geno %>% distinct(genotype, .keep_all = T)
A <- table_DF(A$Subtype) %>% arrange(desc(Freq))
A$x <- factor(A$x, levels = A$x)

A$HA <- str_extract(A$x, pattern = 'H[0-9]+')
A$N <- str_extract(A$x, pattern = 'N[0-9]+')
A$HL <- 'LPAI'
A[A$HA %in% c('H5', 'H7'), 'HL'] <- 'HPAI'

A %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=Freq, color=HL))+
  labs(x='Years', y='Subtypes', size='Types of genotype', color='Genotype Frequency', 
       title='Genotype frequency among years and subtypes')+
  scale_size(range = c(3,35), breaks = seq(1, max(A$Freq), length.out =10)%>% round(digits = 0))+
  scale_color_manual(values=c('#66C8E6', '#E68466'))+
  theme_bw()+gg_theme+theme(legend.position="bottom")

A %>% filter(Freq >10) %>% 
  ggplot()+
  geom_bar(aes(x=Freq, y=x, fill=HL), stat='identity', position='dodge')+
  geom_text(aes(x=Freq, y=x, label=Freq, color=HL), hjust=-0.5, fontface = "bold")+
  labs(fill='', x='', y='Genotype', title=paste('Genotype distribution (more than 100)'), color='')+
  scale_fill_manual(values=c('#66C8E6', '#E68466'))+
  scale_color_manual(values=c('#66C8E6', '#E68466'))+
  # scale_color_brewer(palette = 'Dark2')+
  theme_bw()+gg_theme+theme(axis.text.y = element_text(size = 18-4, family = "mono")) #, color=colors



A <- merge(table_DF(meta_geno$Subtype) %>% rename('n_virus'='Freq'), 
           meta_geno %>% distinct(genotype, .keep_all = T) %>% select(Subtype) %>% table_DF() %>% rename('n_genotype'='Freq'), 
           by.x='x', by.y = 'Subtype')

A$genotype_density <- A$n_genotype/A$n_virus*100
A$p_genotype <- A$n_genotype/sum(A$n_genotype)*100
A$HA <- str_extract(A$x, pattern = 'H[0-9]+') %>% factor(levels = class$HA)
A$N <- str_extract(A$x, pattern = 'N[0-9]+') %>% factor(levels = class$N)

png('~/Analysis/aiv/merge/0307/result/genotype_subtype_distb_dot.png', width = 9600, height = 5400, res = 600)
A %>% as.data.table() %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=n_virus, color=n_genotype))+
  labs(x='', y='', color='Types of genotype   ', size='       Numbers of virus')+
  # scale_color_gradient()  
  viridis::scale_color_viridis(limits = c(0, max(A$n_genotype)), end = max(A$n_genotype)/197, 
                               breaks = c(seq(0, max(A$n_genotype), 50), max(A$n_genotype)), 
                               labels = c(seq(0, max(A$n_genotype), 50), max(A$n_genotype)))+
  scale_size(range = c(1,28), breaks = 10^(seq(0, log10(max(A$n_virus)), length.out =5)) %>% round(digits = 0), 
             labels = 10^(seq(0, log10(max(A$n_virus)), length.out =5)) %>% round(digits = 0))+
  # scale_color_manual(values=c('#66C8E6', '#E68466'))+
  theme_bw()+gg_theme+theme(legend.position="bottom")
dev.off()


sub_meta_geno <- meta_geno %>% filter(Collection_Date >= 2014)
A <- merge(table_DF(sub_meta_geno$Subtype) %>% rename('n_virus'='Freq'), 
           sub_meta_geno %>% distinct(genotype, .keep_all = T) %>% select(Subtype) %>% table_DF() %>% rename('n_genotype'='Freq'), 
           by.x='x', by.y = 'Subtype') %>% mutate(y='≥ 2014    ') %>% filter(n_virus > 100)

A$genotype_density <- A$n_virus/A$n_genotype
A$p_genotype <- A$n_genotype/sum(A$n_genotype)*100
A$HA <- str_extract(A$x, pattern = 'H[0-9]+') %>% factor(levels = class$HA)
A$N <- str_extract(A$x, pattern = 'N[0-9]+') %>% factor(levels = class$N)
A$HA_type <- str_extract(A$HA %>% as.character(), pattern = '[0-9]+') %>% as.numeric()
A <- A %>% arrange(HA_type)



sub_meta_geno <- meta_geno %>% filter(Collection_Date < 2014)
B <- merge(table_DF(sub_meta_geno$Subtype) %>% rename('n_virus'='Freq'), 
           sub_meta_geno %>% distinct(genotype, .keep_all = T) %>% select(Subtype) %>% table_DF() %>% rename('n_genotype'='Freq'), 
           by.x='x', by.y = 'Subtype') %>% mutate(y='< 2014    ')  %>% filter(n_virus > 100)

B$genotype_density <- B$n_virus/B$n_genotype
B$p_genotype <- B$n_genotype/sum(B$n_genotype)*100
B$HA <- str_extract(B$x, pattern = 'H[0-9]+') %>% factor(levels = class$HA)
B$N <- str_extract(B$x, pattern = 'N[0-9]+') %>% factor(levels = class$N)
B$HA_type <- str_extract(B$HA %>% as.character(), pattern = '[0-9]+') %>% as.numeric()
B <- B %>% arrange(HA_type)


A %>% filter(x %in% intersect(A$x, B$x)) %>% as.data.table() %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=n_virus, color=n_genotype))+
  labs(x='', y='', color='Types of genotype   ', size='       Numbers of virus')+
  # scale_color_gradient()  
  viridis::scale_color_viridis(limits = c(0, max(A$n_genotype)), end = max(A$n_genotype)/197, 
                               breaks = c(seq(0, max(A$n_genotype), 50), max(A$n_genotype)), 
                               labels = c(seq(0, max(A$n_genotype), 50), max(A$n_genotype)))+
  scale_size(range = c(3,34*max(A$n_virus)/6883), breaks = 10^(seq(0, log10(max(A$n_virus)), length.out =5)) %>% round(digits = 0), 
             labels = 10^(seq(0, log10(max(A$n_virus)), length.out =5)) %>% round(digits = 0))+
  # scale_color_manual(values=c('#66C8E6', '#E68466'))+
  theme_bw()+gg_theme+theme(legend.position="bottom")


B %>% filter(x %in% intersect(A$x, B$x)) %>% as.data.table() %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=n_virus, color=n_genotype))+
  labs(x='', y='', color='Types of genotype   ', size='       Numbers of virus')+
  # scale_color_gradient()  
  viridis::scale_color_viridis(limits = c(0, max(B$n_genotype)), end = max(B$n_genotype)/197, 
                               breaks = c(seq(0, max(B$n_genotype), 50), max(B$n_genotype)), 
                               labels = c(seq(0, max(B$n_genotype), 50), max(B$n_genotype)))+
  scale_size(range = c(3,34*max(B$n_virus)/6883), breaks = 10^(seq(0, log10(max(B$n_virus)), length.out =5)) %>% round(digits = 0), 
             labels = 10^(seq(0, log10(max(B$n_virus)), length.out =5)) %>% round(digits = 0))+
  # scale_color_manual(values=c('#66C8E6', '#E68466'))+
  theme_bw()+gg_theme+theme(legend.position="bottom")



rbind(A %>% filter(x %in% intersect(A$x, B$x)), B %>% filter(x %in% intersect(A$x, B$x))) %>% 
  mutate(x=factor(x, levels=A$x)) %>% as.data.table() %>% 
  ggplot()+
  geom_bar(mapping=aes(x=x, y=n_virus, fill=y), position = 'dodge', stat='identity')+
  labs(x='Subtype', y='Numbers of virus', fill='')+
  # scale_color_gradient()  
  scale_fill_manual(values=c('#66C8E6', '#353692'))+
  theme_bw()+gg_theme


rbind(A %>% filter(x %in% intersect(A$x, B$x)), B %>% filter(x %in% intersect(A$x, B$x))) %>% 
  mutate(x=factor(x, levels=A$x)) %>% as.data.table() %>% 
  ggplot()+
  geom_bar(mapping=aes(x=x, y=genotype_density, fill=y), position = 'dodge', stat='identity')+
  labs(x='Subtype', y='Numbers of virus', fill='')+
  # scale_color_gradient()  
  scale_fill_manual(values=c('#66C8E6', '#353692'))+
  theme_bw()+gg_theme


axis_margin <- 7

p1 <- 
  rbind(A %>% filter(x %in% intersect(A$x, B$x)), B %>% filter(x %in% intersect(A$x, B$x))) %>% 
  mutate(x=factor(x, levels=A$x)) %>% as.data.table() %>% 
  ggplot() +
  geom_col(mapping=aes(y=x, x=n_virus, fill=y), position = 'dodge', stat='identity')+
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values=c('#66C8E6', '#353692'))+
  labs(fill='Years    ', x='Virus strains')+
  guides(fill = FALSE)+
  theme_bw()+gg_theme+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
  )

p2 <- 
  rbind(A %>% filter(x %in% intersect(A$x, B$x)), B %>% filter(x %in% intersect(A$x, B$x))) %>% 
  mutate(x=factor(x, levels=A$x)) %>% as.data.table() %>% 
  ggplot() +
  geom_col(mapping=aes(y=x, x=n_genotype, fill=y), position = 'dodge', stat='identity')+
  scale_fill_manual(values=c('#66C8E6', '#353692'))+
  labs(fill='Years    ', x='Genotypes')+
  theme_bw()+gg_theme+
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
  )


ggpubr::ggarrange(p1, p2, common.legend = TRUE, legend='bottom')
# annotate_figure(top = text_grob(paste("Total nodes and removed nodes under", j, " MPD threshold"), face = "bold", size = 24)) #, color = "red", face = "bold", size = 14)



# treemap internal/NA -----------------------------------------------------
# library
library(treemap)

bygroup_result <- list()
for(i in c(class$Internal, "NA")) {
  bygroup <- eight[[i]] %>% group_split(clades) %>% as.list() %>% lapply(as.data.frame)
  country_list <- list()
  x=1
  for(j in seq_along(bygroup)) {
    x=x+1
    country_list[[paste('L', j, sep = '')]] <- word(bygroup[[j]]$Location, start=1, end=1, sep='/') %>% table_DF(prop = F)
  }
  bygroup_result[[i]] <- 
    do.call(what=rbind, country_list) %>% as.data.frame() %>% rename('cont'='x') %>% tibble::rownames_to_column('x') %>% 
    mutate(x=str_remove(x, pattern='\\..*'))
}

cont_short <- 
  data.frame(full=c("Africa", "Asia", "Europe", "NorthAmerica", "Oceania", "SouthAmerica", "Antarctica"), 
             short=c('AF', 'AS', 'EU', 'NA', 'OC', 'SA', 'AN'))

test <-
  do.call(what=rbind, bygroup_result) %>% as.data.frame() %>% rename('clade'='x') %>% tibble::rownames_to_column('seg') %>%
  mutate(seg=str_remove(seg, pattern='\\..*')) %>% mutate(root='root') %>% filter(seg != 'HA')
test <- test[, c(ncol(test), 1:4)] #%>% filter(seg %in% c('MP', 'PA'))
for(i in 1:nrow(cont_short)) {
  test[test$cont %in% cont_short[i, 1], 'cont'] <- cont_short[i, 2]
}

p_list <- list()
for(i in class$Internal) {
  seg_clade <- test %>% filter(seg == i) %>% mutate(p1=paste(root, seg, sep = '.'), p2=paste(root, seg, clade, sep = '.')) %>%
    select(c('p1', 'p2', 'clade', 'Freq')) %>%
    group_split(p2) %>% as.list() %>% lapply(as.data.frame)
  seg_clade <- data.frame(from=lapply(seg_clade, function(x){x[, 2] %>% unique()})  %>% unlist(),
                          to=lapply(seg_clade, function(x){x[, 3] %>% unique()}) %>% unlist(),
                          Freq=lapply(seg_clade, function(x){x[, 4] %>% sum}) %>% unlist(),
                          last=lapply(seg_clade, function(x){x[, 1] %>% unique()}) %>% unlist())
  seg_clade <- seg_clade %>% top_n(wt = Freq, 12) %>% filter(Freq > (Freq %>% sum())*0.01)
  
  data <- test %>% filter(seg == i, clade %in% seg_clade$to) %>% select(c('clade', 'cont', 'Freq')) %>%  
    group_split(clade) %>% as.list() %>% lapply(as.data.frame)
  for(j in seq_along(data)) {
    data[[j]][, 'percent'] <-  round(data[[j]][, 'Freq'] / sum(data[[j]][, 'Freq']) * 100, digits = 2)
    data[[j]][(nrow(data[[j]])+1), ] <- c(data[[j]][1, 1], 'Total', sum(data[[j]][, 3]), 100)
    data[[j]] <- data[[j]] %>% 
      mutate(Freq=as.numeric(Freq), 
             percent=as.numeric(percent)) %>% 
      arrange(desc(Freq)) %>%  # Arrange within each clade by descending Freq
      mutate(cumulative_Freq = cumsum(Freq),  # Cumulative sum for segment placement
             label_position = cumulative_Freq - (Freq / 2))
  }
  
  # Combine all data back into a single dataframe and adjust labels
  data <- do.call(rbind, data) %>% as.data.frame() %>%
    mutate(label = paste(cont, ', ', percent, '%', sep = ''),
           cont = factor(cont, levels = c("AF", "AS", "EU", "NA", "OC", "AN", "SA", 'Total')))
  data$clade <- factor(data$clade, levels = paste('L', data$clade %>% header_cleaning('L') %>% arrange(as.numeric(V2)) %>% select(V2) %>% unlist(), sep = '') %>% 
                         unique())
  # Remove labels for values below 5% to avoid clutter
  # data[data$percent < 5, 'label'] <- ''
  
  # Plot the data with correct stacking and text positions
  # ggplot(data) +
  #   geom_col(aes(x = clade, y = Freq, fill = forcats::fct_reorder2(cont, Freq, clade, .desc = TRUE))) +
  #   geom_text(aes(x = clade, y = label_position, label = label), size = 3) +
  #   scale_fill_brewer(palette = 'Set2') +
  #   labs(fill = 'Continent', x='Major clades', y='Frequency') +
  #   theme_bw() + gg_theme
  
  
  # Custom labels:
  # p_list[[i]] <- 
  # ggplot(data) +
  #   geom_bar(aes(x = clade, y = Freq, fill = cont), stat = "identity", position = position_dodge(width = 0.9)) +
  #   geom_text(aes(x = clade, y = Freq, label = label, group = cont), 
  #             position = position_dodge(width = 0.9), vjust = 0.5, hjust=-0.1, size = 4, angle=90) +
  #   scale_fill_brewer(palette = 'Set2') +
  #   labs(fill = 'Continent', x = 'Major lineages', y = 'Frequency', title = paste(i, ' major lineages geographycal distribution', sep='')) +
  #   scale_y_continuous(limits = c(0, max(data$Freq)+250))+
  #   scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))+
  #   theme_bw() + gg_theme
  
  n <- (unique(data$clade) %>% length())
  png(filename = paste('~/Analysis/aiv/merge/0307/result/', i, '_clade_geo_distb_hbar.png', sep = ''),width = 9600, height = 5400, res=600)
  ggplot() +
    geom_bar(data %>% filter(cont != 'Total'), 
             mapping=aes(x = percent, y = clade, fill = cont), stat = "identity", position = 'stack') +
    geom_text(data %>% filter(cont == 'Total'), 
              mapping=aes(y = clade, x = 105, label = Freq),
              position = position_dodge(width = 0.9), vjust = 0.5, hjust=0, size = 5, angle=0) +
    geom_vline(xintercept = c(-1, 101), linetype = "solid", color = "black") +
    # Add the text 'Total\nsequence\nnumber' at the top of the numbers
    geom_text(mapping=aes(y = (unique(data$clade) %>% length())+0.5, x = 105, label = "Total\nsequence\nnumber"),
              hjust = 0.1, vjust = 0, size = 5) +
    scale_fill_brewer(palette = 'Set2') +
    labs(fill = 'Continent', y = 'Major lineages', x = 'Continent ratio', title = paste(i, ' major lineages geographycal distribution', sep='')) +
    scale_y_discrete(expand = expansion(c(1/n, 2/n))) +
    scale_x_continuous(limits = c(-5,110), labels = seq(0, 100, 10), breaks =  seq(0, 100, 10))+
    theme_bw() + gg_theme
  dev.off()
  # treemap(data,
  #         index = c("clade", "cont"),
  #         vSize="Freq",
  #         vColor = "clade",
  #         type = "categorical",
  #         fontsize.labels = c(28, 14), # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
  #         fontcolor.labels = c("black", "grey30"), # Color of labels
  #         fontface.labels = c(20, 10), # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
  #         # bg.labels = c("transparent"), # Background color of labels
  #         align.labels = list(
  #           c("left", "top"),
  #           c("right", "bottom")
  #         ), # Where to place labels in the rectangle?
  #         bg.labels = 255,
  #         overlap.labels = 0.5, # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
  #         inflate.labels = F,
  #         palette = "Set3", 
  #         title = paste(i, ' major clades geographycal distribution', sep=''), 
  #         fontsize.title = 30, 
  #         position.legend = "none")
  
}

for(i in names(p_list)) {
  png(filename = paste('~/Analysis/aiv/merge/0307/result/', i, '_clade_geo_distb_hbar.png', sep = ''),width = 9600, height = 5400, res=600)
  print(p_list[[i]])
  dev.off()
}



bygroup <- eight[['NA']] %>% group_split(clades) %>% as.list() %>% lapply(as.data.frame)
# names(bygroup) <- lapply(bygroup, function(x){x[1, 2]})
country_list <- list()
x=1
for(j in seq_along(bygroup)) {
  x=x+1
  country_list[[j]] <- word(bygroup[[j]]$Location, start=1, end=1, sep='/') %>% table_DF(prop = F)
  names(country_list)[j] <- bygroup[[j]]$clades %>% unique()
}

na_result <- 
  do.call(what=rbind, country_list) %>% as.data.frame() %>% rename('cont'='x') %>% tibble::rownames_to_column('x') %>% 
  mutate(x=str_remove(x, pattern='\\..*'), na_sub=str_extract(x, pattern='N[0-9]+'), group=str_extract(x, pattern='L[0-9]+'))

for(i in 1:nrow(cont_short)) {
  na_result[na_result$cont %in% cont_short[i, 1], 'cont'] <- cont_short[i, 2]
}


p_list <- list()
for(i in class$N) {
  sub_clade <- na_result %>% filter(na_sub %in% i) %>% 
    group_split(group) %>% as.list() %>% lapply(as.data.frame)
  A <- lapply(sub_clade, function(x){x[, 3] %>% sum}) %>% unlist()
  
  data <- sub_clade[which((A/sum(A) > 0.01))]
  # lapply(sub_clade, function(x){round((x[, 3]/sum(x[, 3])*100), digit=2)}) %>% unlist()
  for(j in seq_along(data)) {
    data[[j]][, 'percent'] <-  round(data[[j]][, 'Freq'] / sum(data[[j]][, 'Freq']) * 100, digits = 2)
    data[[j]][(nrow(data[[j]])+1), ] <- c(data[[j]][1, 1], 'Total', sum(data[[j]][, 3]), data[[j]][1, 4], data[[j]][1, 5], 100)
    data[[j]] <- data[[j]] %>% 
      mutate(Freq=as.numeric(Freq), 
             percent=as.numeric(percent)) %>% 
      arrange(desc(Freq)) %>%  # Arrange within each clade by descending Freq
      mutate(cumulative_Freq = cumsum(Freq),  # Cumulative sum for segment placement
             label_position = cumulative_Freq - (Freq / 2))
  }
  
  data <- do.call(rbind, data) %>% as.data.frame() %>%
    mutate(label = paste(cont, ', ', percent, '%', sep = ''),
           cont = factor(cont, levels = c("AF", "AS", "EU", "NA", "OC", "AN", "SA", 'Total')))
  
  data$group <- factor(data$group, 
                       levels = c(paste('L', data$group %>% header_cleaning('L') %>% arrange(as.numeric(V2)) %>% select(V2) %>% unlist(), sep = '') %>%  unique()))
  
  # data[data$prop <= 5, 'label'] <- ''
  
  
  # data <- data %>% filter(na_sub %in% 'N1', Freq > (Freq %>% sum())*0.01)
  # Custom labels:
  # p_list[[i]] <- 
  # ggplot(data) +
  #   geom_bar(aes(x = group, y = Freq, fill = cont), stat = "identity", position = position_dodge(width = 0.9)) +
  #   geom_text(aes(x = group, y = Freq, label = label, group = cont), 
  #             position = position_dodge(width = 0.9), vjust = 0.5, hjust=-0.1, size = 4, angle=90) +
  #   scale_fill_brewer(palette = 'Set2') +
  #   labs(fill = 'Continent', x = 'Major lineages', y = 'Frequency', title = paste(i, ' major lineages geographycal distribution', sep='')) +
  #   scale_y_continuous(limits = c(0, max(data$Freq)+250))+
  #   scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))+
  #   theme_bw() + gg_theme
  
  
  num_groups <- unique(data$group) %>% length()
  
  label_positon <- num_groups+0.5
  # Set expansion values
  y_expansion_bottom <- 1 / num_groups
  y_expansion_top <- 1.8 / num_groups
  
  # p_list[[i]] <-
  png(filename = paste('~/Analysis/aiv/merge/0307/result/', i, '_clade_geo_distb_hbar.png', sep = ''),width = 9600, height = 5400, res=600)
  ggplot() +
    geom_bar(data %>% filter(cont != 'Total'), 
             mapping=aes(x = percent, y = group, fill = cont), stat = "identity", position = 'stack') +
    geom_text(data %>% filter(cont == 'Total'), 
              mapping=aes(y = group, x = 105, label = Freq),
              position = position_dodge(width = 0.9), vjust = 0.5, hjust=0, size = 5, angle=0) +
    geom_vline(xintercept = c(-1, 101), linetype = "solid", color = "black") +
    # Add the text 'Total\nsequence\nnumber' at the top of the numbers
    geom_text(mapping=aes(y = label_positon, x = 105, label = "Total\nsequence\nnumber"),
              hjust = 0.1, vjust = 0, size = 5) +
    scale_fill_brewer(palette = 'Set2') +
    labs(fill = 'Continent', y = 'Major lineages', x = 'Continent ratio', title = paste(i, ' major lineages geographycal distribution', sep='')) +
    scale_y_discrete(expand = expansion(c(y_expansion_bottom, y_expansion_top))) +
    scale_x_continuous(limits = c(-5,110), labels = seq(0, 100, 10), breaks =  seq(0, 100, 10))+
    theme_bw() + gg_theme
  dev.off()
  # treemap(data,
  #         index = c('group', "cont"),
  #         vSize="Freq",
  #         vColor = "group",
  #         type = "categorical",
  #         fontsize.labels = c(28, 14), # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
  #         fontcolor.labels = c("black", "grey30"), # Color of labels
  #         fontface.labels = c(20, 10), # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
  #         # bg.labels = c("transparent"), # Background color of labels
  #         align.labels = list(
  #           c("left", "top"),
  #           c("right", "bottom")
  #         ), # Where to place labels in the rectangle?
  #         bg.labels = 255,
  #         overlap.labels = 0.5, # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
  #         inflate.labels = F,
  #         palette = "Set3",
  #         title = paste(i, ' major clades geographycal distribution', sep=''),
  #         fontsize.title = 30,
  #         position.legend = "none"
  # )
}

for(i in names(p_list)) {
  png(filename = paste('~/Analysis/aiv/merge/0307/result/', i, '_clade_geo_distb_hbar.png', sep = ''),width = 9600, height = 5400, res=600)
  print(p_list[[i]])
  dev.off()
}



# assign genotype old ---------------------------------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'

group <- list()
for(i in list.files(path = '/home/eric/Analysis/aiv/merge/0307/group/groups/'))  {
  clade <- read_as_list(paste('/home/eric/Analysis/aiv/merge/0307/group/groups/' ,i , '/', sep = ''), file_type = 'txt')
  for(j in names(clade))  {
    clade[[j]] <- sub(pattern = 'tree tree_1 \\= \\[\\&R\\] ', clade[[j]][2,1], replacement = '') %>%
      strsplit(split='\'') %>% unlist() %>% grep(pattern = '^EPI', value = T) %>% as.data.frame()
  }
  names(clade) <- LETTERS[1:length(clade)]
  clade <- do.call(rbind, clade) %>% as.data.frame() %>% tibble::rownames_to_column('group')
  clade$group <- sub(clade$group, pattern = '\\..*', replacement = '') %>% str_remove(pattern = 'up')
  clade[, c('epi', 'strain_number')] <- (clade$. %>% strsplit(split = '\\|') %>% do.call(what=rbind))[, c(1,2)]
  clade$ui <- paste(clade$epi, clade$strain_number, sep = '|')
  clade$s_g <- paste(i, clade$group, sep = '')
  group[[i]] <- clade
}
genotype <- assign_geno(group_list = group, meta = meta, remove_homo = F)

A <- table_DF(genotype$genotype) %>% arrange(desc(Freq))
A <- A[!grepl(A$x, pattern='X'), ]
A[, 'HL'] <- 'LPAI'
A[grepl(header_cleaning(A$x, pattern = '_') %>% select(V4) %>% unlist, pattern='[A-Z]'), 'HL'] <- 'HPAI'

png('~/Analysis/aiv/merge/0307/plot/genotype.png', width = 9600, height = 5400, res = 600)
A %>% filter(Freq> 100) %>% mutate(x=factor(x, levels=c(x))) %>% 
  ggplot()+
  geom_bar(aes(x=Freq, y=x, fill=HL), stat='identity', position='dodge')+
  geom_text(aes(x=Freq, y=x, label=Freq, color=HL), hjust=-0.5, fontface = "bold")+
  labs(fill='', x='', y='Genotype', title=paste('Genotype distribution (more than 100)'), color='')+
  scale_fill_manual(values=c('#66C8E6', '#E68466'))+
  scale_color_manual(values=c('#66C8E6', '#E68466'))+
  # scale_color_brewer(palette = 'Dark2')+
  theme_bw()+gg_theme+theme(axis.text.y = element_text(size = 18-4, family = "mono")) #, color=colors
dev.off()


meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()


A <- table(meta_geno$genotype, meta_geno$Subtype) %>% as.data.frame() %>% filter(Freq!=0)
# A[, 3:4] <- 

A$HA <- str_extract(A$Var2, pattern = 'H[0-9]+')
A$N <- str_extract(A$Var2, pattern = 'N[0-9]+')

A %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=Freq))+
  labs(x='Years', y='Subtypes', size='Types of genotype', color='Genotype Frequency', 
       title='Genotype frequency among years and subtypes')+
  scale_size(range = c(3,30), breaks = seq(1, max(A$Freq), length.out =4)%>% round(digits = 0))+
  scale_color_gradient(low = '#A8E6FF', high = '#003CB5')+
  theme_bw()+gg_theme+theme(legend.position="bottom")



A <- meta_geno %>% distinct(genotype, .keep_all = T)
A <- table_DF(A$Subtype) %>% as.data.frame() %>% filter(Freq!=0)
# A[, 3:4] <- 

A$HA <- str_extract(A$x, pattern = 'H[0-9]+')
A$N <- str_extract(A$x, pattern = 'N[0-9]+')

A %>% 
  ggplot()+
  geom_point(mapping=aes(x=HA, y=N, size=Freq))+
  labs(x='Years', y='Subtypes', size='Types of genotype', color='Genotype Frequency', 
       title='Genotype frequency among years and subtypes')+
  scale_size(range = c(3,30), breaks = seq(1, 23, length.out =4)%>% round(digits = 0))+
  scale_color_gradient(low = '#A8E6FF', high = '#003CB5')+
  theme_bw()+gg_theme+theme(legend.position="bottom")

A %>% 
  ggplot()+
  geom_bar(aes(x=Freq, y=x, fill=HL), stat='identity', position='dodge')+
  geom_text(aes(x=Freq, y=x, label=Freq, color=HL), hjust=-0.5, fontface = "bold")+
  labs(fill='', x='', y='Genotype', title=paste('Genotype distribution (more than 100)'), color='')+
  scale_fill_manual(values=c('#66C8E6', '#E68466'))+
  scale_color_manual(values=c('#66C8E6', '#E68466'))+
  # scale_color_brewer(palette = 'Dark2')+
  theme_bw()+gg_theme+theme(axis.text.y = element_text(size = 18-4, family = "mono")) #, color=colors




# segment groups on H5 ----------------------------------------------------
tree_list <- read_as_list(path = '~/Analysis/aiv/merge/0307/', file_type = 'tree', prefix = 'aligned.nwk$')
names(tree_list) <- str_extract(names(tree_list), pattern = '[A-Z]+[0-9]+_a|[a-z]+[A-Z]+_a|[A-Z]+_a|NS1true_a') %>% 
  str_remove_all(pattern = '_a')


group <- list()
for(i in list.files(path = '/home/eric/Analysis/aiv/merge/0307/group/groups/'))  {
  clade <- read_as_list(paste('/home/eric/Analysis/aiv/merge/0307/group/groups/' ,i , '/', sep = ''), file_type = 'txt')
  for(j in names(clade))  {
    clade[[j]] <- sub(pattern = 'tree tree_1 \\= \\[\\&R\\] ', clade[[j]][2,1], replacement = '') %>%
      strsplit(split='\'') %>% unlist() %>% grep(pattern = '^EPI|^IRD', value = T) %>% as.data.frame()
  }
  names(clade) <- LETTERS[1:length(clade)]
  clade <- do.call(rbind, clade) %>% as.data.frame() %>% tibble::rownames_to_column('group')
  clade$group <- sub(clade$group, pattern = '\\..*', replacement = '') %>% str_remove(pattern = 'up')
  clade[, c('epi', 'strain_number')] <- (clade$. %>% strsplit(split = '\\|') %>% do.call(what=rbind))[, c(1,2)]
  clade$ui <- paste(clade$epi, clade$strain_number, sep = '|')
  clade$s_g <- paste(i, clade$group, sep = '')
  group[[i]] <- clade
}
genotype <- assign_geno(group_list = group, meta = meta)


clade_table <- table_DF(meta[grepl(meta$Subtype, pattern='H5'), 'Clade'])
clade_table[, 'Clade'] <- 'GsGD-others'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4'), 'Clade'] <- '2.3.4.4'
clade_table[grep(clade_table$x, pattern = '2\\.3\\.4\\.4b'), 'Clade'] <- '2.3.4.4b'
clade_table[grep(clade_table$x, pattern = 'nonGsGD'), 'Clade'] <- 'nonGsGD'
clade_table[grep(clade_table$x, pattern = 'IRD_unlabel'), 'Clade'] <- 'IRD_unlabel'

i='H5'
tree <- tree_list[[i]]

tip_name <- tree$tip.label %>% str_remove_all(pattern = "'")

info <- header_cleaning(tip_name, pattern = '\\|')
colnames(info) <- c("Accesion_number","Strain_number","Subtype", 'Clade', "Segment", 'Collection_date'
                    ,"Location","Host", "Header_Host", "Host_type")

country <- meta[meta$Strain_number %in% info$Strain_number, c('Strain_number', 'Location')]
rownames(country) <- country$Strain_number
country <- country[info$Strain_number, ]
info$Location <- country$Location
info$continent <- info$Location %>% header_cleaning(pattern = '\\/') %>% select(V1) %>% unlist()

info[info$Accesion_number %in% 'EPI_ISL_66116', 'Collection_date'] <- 1970
info$Host_type <- replace(info$Host_type, info$Host_type %in% '', 'Wild')
info$Year <- word(info$Collection_date, sep = '-', 1) %>% as.numeric()

rank <- info$Strain_number
meta_geno <- merge(genotype, meta[meta$Strain_number %in% rank, ], by = 'Isolate_Id')
info <- merge(info, 
              meta_geno[meta_geno$Strain_number %in% rank, c('Strain_number', "PB2", "PB1", "PA",  "HA",  "NP", "NA",  "MP", "NS", 'genotype')], 
              by = 'Strain_number', all = T)
rownames(info) <- info$Strain_number
info <- info[rank, ]

info[, 'NA_new'] <- 'Minor_NA'
info[info$`NA` %in% (table_DF(info$`NA`) %>% filter(Freq>100))[, 1], 'NA_new'] <- 
  info[info$`NA` %in% (table_DF(info$`NA`) %>% filter(Freq>100))[, 1], 'NA']
info[is.na(info$`NA`), 'NA_new'] <- NA
info$NA_new <- factor(info$NA_new, levels = c("N1A", "N1B", "N2B", "N2C", "N3B", "N6A", "N6B", "N6C", "N8A", 'Minor_NA'))

for(j in c("PB2", "PB1", "PA",  "NP", "HA",  "MP", "NS")) {
  info[, j] <- replace(info[, j], info[, j]=='X',NA)
  info[, j] <- factor(info[, j], levels = c(LETTERS[1:14], NA))
}

info$time_group <- 0
OwO <- str_remove(info$Year, pattern = '-.*') %>% as.numeric()

info[OwO < 1996, 'time_group'] <- "before_1995"
info[OwO <= 2013 & OwO >= 1996, 'time_group'] <- "1996~2013"
info[OwO <= 2020 & OwO >= 2014, 'time_group'] <- "2014~2020"
info[OwO >=2021, 'time_group'] <- "after_2021"

time <- c("before_1995", "1996~2013", "2014~2020", "after_2021")
info$time_group <- factor(info$time_group, levels = time)




info$Strain_number <- factor(info$Strain_number, levels = rev(info$Strain_number))
#change tree header to match the data
new_header <- info$Strain_number
tree$tip.label <- new_header
# phylo <- ggtree(tree)

#Heatmap from sequence information
Heatmap_1 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=PB2, color=PB2), linewidth=0.05)+ 
  labs(x="PB2", fill='PB2', color='PB2') + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")

Heatmap_2 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=PB1, color=PB1), linewidth=0.05)+ 
  labs(x="PB1", fill="PB1", color="PB1") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")

Heatmap_3 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=PA, color=PA), linewidth=0.05)+
  labs(x="PA", fill="PA", color="PA") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")

Heatmap_4 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=HA, color=HA), linewidth=0.05)+
  labs(x="HA", fill="HA", color="HA") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")


Heatmap_5 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=NP, color=NP), linewidth=0.05)+
  labs(x="NP", fill="NP", color="NP") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")


Heatmap_6 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=NA_new, color=NA_new), linewidth=0.05)+
  labs(x="NA", fill="NA", color="NA") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme

Heatmap_7 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=MP, color=MP), linewidth=0.05)+
  labs(x="MP", fill="MP", color="MP") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme+theme(legend.position = "none")


Heatmap_8 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=NS, color=NS), linewidth=0.05)+ 
  labs(x="NS", fill="Other_segments", color="Other_segments") + ylab(NULL)+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = 'Set3')+
  phylo_theme	

Heatmap_9 <- ggplot(info)+geom_tile(aes(x='', y=Strain_number,fill=continent, color=continent), linewidth=0.05)+ 
  labs(x="Continent", fill='Continent', color='Continent') + ylab(NULL)+
  scale_color_brewer(palette = 'Accent')+
  scale_fill_brewer(palette = 'Accent')+
  phylo_theme

Heatmap_10 <- ggplot(info)+geom_tile(aes(x="",y=Strain_number, fill=time_group, color=time_group), linewidth=0.05)+
  labs(x="Year", fill="Year", color="Year") + ylab(NULL)+
  scale_color_brewer(palette = 'Blues')+
  scale_fill_brewer(palette = 'Blues')+
  phylo_theme

png(paste('~/Analysis/aiv/merge/0307/plot/H5_genotype_heatmap','.png', sep = ''), width = 9600, height = 5400, res=600)
as.ggplot(Heatmap_1 %>% insert_left(phylo,width = 12) %>% insert_right(Heatmap_2) %>% insert_right(Heatmap_3) %>% 
            insert_right(Heatmap_4) %>% 
            insert_right(Heatmap_5)%>% insert_right(Heatmap_6) %>% 
            insert_right(Heatmap_7)%>% insert_right(Heatmap_8) %>% 
            insert_right(Heatmap_9)%>% insert_right(Heatmap_10))+
  ggtitle(paste(i,"segment phylogenetic tree"))+theme(plot.title = element_text(face = "bold",size = 25 ,hjust = 0.5))
dev.off()










# genotype lasting period -------------------------------------------------
meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]

meta_geno$genotype <- paste(meta_geno$genotype, meta_geno$Clade, sep = '_')

meta_geno_sub <- meta_geno %>% filter(HA %in% c('A', 'B', 'C'), Subtype %in% 'H5N1', Clade != 'IRD_unlabel') %>% arrange(Collection_Date) %>% 
  group_split(genotype) %>% as.list() %>% lapply(as.data.frame)


which((lapply(meta_geno_sub, nrow) %>% unlist()) ==1)

result <- matrix(0,0,0) %>% as.data.frame()
x <- 0
which((lapply(meta_geno_sub, nrow) %>% unlist()) >2)
for(i in which((lapply(meta_geno_sub, nrow) %>% unlist()) >2)) {
  x <- x+1
  A <- meta_geno_sub[[i]]
  # A[A$Clade %in% 'IRD_unlabel', 'Clade'] <- (table_DF(A$Clade) %>% arrange(desc(Freq)))[1, 1]
  
  
  A[, 'Collection_Date'] <- str_remove(A[, 'Collection_Date'], pattern = '-$|--$')
  
  A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}$'), 'Collection_Date'] <- 
    paste(A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}$'), 'Collection_Date'], '-01-01', sep = '')
  
  A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}-[0-9]{2}$'), 'Collection_Date'] <- 
    paste(A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}-[0-9]{2}$'), 'Collection_Date'], '-01', sep = '') 
  
  A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}-[0-9]{2}-[0-9]{2}$'), 'Collection_Date'] <- 
    A[grepl(A[, 'Collection_Date'], pattern='[0-9]{4}-[0-9]{2}-[0-9]{2}$'), 'Collection_Date']
  
  B <- A %>% mutate(p=paste(Collection_Date, Location)) %>% distinct(p, .keep_all = T)
  # mean(as.Date(B[2:(nrow(B)), 'Collection_Date']) -as.Date(B[1:(nrow(B)-1), 'Collection_Date']))
  
  result[x, 'date'] <- (as.Date(A[nrow(A), 'Collection_Date'])-as.Date(A[1, 'Collection_Date']))
  result[x, 'median_time_interval'] <- median(as.Date(B[2:(nrow(B)), 'Collection_Date']) -as.Date(B[1:(nrow(B)-1), 'Collection_Date']))
  result[x, 'n_clade'] <- unique(A$Clade) %>% length()
  result[x, 'n_seq'] <- (A$Clade) %>% length()
  # result[x, 'n_seq'] <- (A$Clade) %>% length()
  result[x, 5:13] <- A[1, 2:10]
}


result <- result %>% filter(date !=0)

ggplot(result)+
  geom_violin(aes(x=HA, y=log(as.numeric(median_time_interval)+1), fill=HA))+
  geom_jitter(aes(x=HA, y=log(as.numeric(median_time_interval)+1)), width = 0.25)+
  theme_bw()+gg_theme


wilcox.test(result[result$HA %in% 'A', 'median_time_interval'] %>% as.numeric(), 
            result[result$HA %in% 'C', 'median_time_interval']%>% as.numeric())

wilcox.test(result[result$HA %in% 'A', 'median_time_interval'] %>% as.numeric(), 
            result[result$HA %in% 'B', 'median_time_interval']%>% as.numeric())


wilcox.test(result[result$HA %in% 'C', 'median_time_interval'] %>% as.numeric(), 
            result[result$HA %in% 'B', 'median_time_interval']%>% as.numeric())

A <- read_as_list('~/Analysis/aiv/merge/0307/test/', sep='', prefix = 'outlier.fa.log$', file_type = 'txt')

lapply(A, function(x){grep(x[, 1], pattern='invariable', value = T)}) %>% unlist()

# genotype clustering -----------------------------------------------------
meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
# meta_geno$Clade_new <- factor(meta_geno$Clade_new , levels = c("2.3.4.4","2.3.4.4b", "GsGD-others","nonGsGD", 'IRD_unlabel'))

# meta_geno <- meta_geno %>% filter(Clade_new %in% c('2.3.4.4'))
meta_geno <- meta_geno %>% arrange(y) %>% filter(y > 1996) %>% distinct(genotype, .keep_all = T)#, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
# meta_geno <- meta_geno[grepl(meta_geno$Subtype, pattern='H5'), ]
meta_geno <- meta_geno[(meta_geno$Collection_Date %>% str_remove_all(pattern = '[0-9]+')=='--'), ]
meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())

A <- meta_geno %>% mutate(p=paste(Collection_Date, genotype)) %>% distinct(p, .keep_all = T) %>% filter(Clade_new %in% c('2.3.4.4', '2.3.4.4b'))

out <- data.frame(g='g', y=A$Collection_Date %>% as.numeric()) %>% filter(y >2010)
# out <- rbind(out, data.frame(g='H5N1', y=meta_geno %>% filter(Subtype=='H5N1') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='H5N8', y=meta_geno %>% filter(Subtype=='H5N8') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4', y=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Collection_Date) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4b', y=meta_geno %>% filter(Clade_new=='2.3.4.4b', y>2016) %>% select(Collection_Date) %>% unlist() %>% as.numeric()))

ggplot(out , aes(x=y, y=g, fill=g))+
  geom_violin()+
  # geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_x_continuous(breaks = seq(2010, 2023, 1))+
  labs(x='Year', y='', fill='', title='Genotype frequency')+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank())


out <- data.frame(g='g', y=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Collection_Date) %>% unlist() %>% as.numeric())
# out <- rbind(out, data.frame(g='H5N1', y=meta_geno %>% filter(Subtype=='H5N1') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='H5N8', y=meta_geno %>% filter(Subtype=='H5N8') %>% select(y) %>% unlist() %>% as.numeric()))
out <- rbind(out, data.frame(g=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Subtype) %>% unlist(), 
                             y=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Collection_Date) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4b', y=meta_geno %>% filter(Clade_new=='2.3.4.4b', y>2016) %>% select(y) %>% unlist() %>% as.numeric()))

ggplot(out, aes(x=y, y=g, fill=g))+
  geom_violin()+
  # geom_jitter(shape=16, position=position_jitter(0.1))+
  theme_bw()+gg_theme





out <- data.frame(g='g', y=meta_geno %>% filter(Clade_new=='2.3.4.4b') %>% select(Collection_Date) %>% unlist() %>% as.numeric())
# out <- rbind(out, data.frame(g='H5N1', y=meta_geno %>% filter(Subtype=='H5N1') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='H5N8', y=meta_geno %>% filter(Subtype=='H5N8') %>% select(y) %>% unlist() %>% as.numeric()))
out <- rbind(out, data.frame(g=meta_geno %>% filter(Clade_new=='2.3.4.4b') %>% select(Subtype) %>% unlist(), 
                             y=meta_geno %>% filter(Clade_new=='2.3.4.4b') %>% select(Collection_Date) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4b', y=meta_geno %>% filter(Clade_new=='2.3.4.4b', y>2016) %>% select(y) %>% unlist() %>% as.numeric()))

ggplot(out, aes(x=y, y=g, fill=g))+
  geom_violin()+
  # geom_jitter(shape=16, position=position_jitter(0.1))+
  theme_bw()+gg_theme




library(ggstream)

meta_geno_sub <- meta_geno %>% filter(Clade_new %in% c('2.3.4.4', '2.3.4.4b', 'GsGD-others')
                                      , Subtype !='H9N2', HA %in% c('A', 'B', 'C'))
A <- (table_DF(meta_geno_sub$y, prop = T) %>% mutate(g='g'))[, c(1, 3, 2)]
B <- table(meta_geno_sub$y, meta_geno_sub$Subtype) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq!=0) %>% 
  mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
B$Var2 <- paste(B$Var2, paste(rep(' ', 29), collapse = ''), sep = '')
B$Freq <- B$Freq/(B$Freq %>% sum())

subtype_stream <- ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title='')+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(breaks = seq(1996, 2023, 1))+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank())



A <- meta_geno %>% filter(Subtype == 'H5N1')
B <- table(A$y, A$Subtype) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq>0) #%>% mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
# B <- table_DF(B$Var2)
# B <- B[B$Var2 %in% (table_DF(A$genotype) %>% filter(Freq > 60) %>% select(x) %>% unlist), ]
B$Var2 <- paste(B$Var2, paste(rep(' ', 29), collapse = ''), sep = '')
B$Freq <- B$Freq/(B$Freq %>% sum())

ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title='')+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(2014, 2023))+
  scale_fill_manual(values = c('#958BE3'))+
  theme_bw()+gg_theme


A <- meta_geno %>% filter(Subtype == 'H5N1')
B <- table(A$y, A$genotype) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq>0) #%>% mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
# B <- table_DF(B$Var2)
B <- B[B$Var2 %in% (table_DF(A$genotype) %>% filter(Freq > 40) %>% select(x) %>% unlist), ]
B$Freq <- B$Freq/(B$Freq %>% sum())

N1_stream <- ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title=paste('H5N1 genotypes', ' (Total genotypes:', table_DF(A$genotype) %>% nrow(), ')', sep = ''), 
       subtitle = 'Frequency more than 40 are included')+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(2014, 2023))+
  scale_fill_brewer(palette = 'Set3')+
  theme_bw()+gg_theme





A <- meta_geno %>% filter(Subtype == 'H5N8')
B <- table(A$y, A$Subtype) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq>0) #%>% mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
# B <- table_DF(B$Var2)
# B <- B[B$Var2 %in% (table_DF(A$genotype) %>% filter(Freq > 60) %>% select(x) %>% unlist), ]
B$Var2 <- paste(B$Var2, paste(rep(' ', 29), collapse = ''), sep = '')
B$Freq <- B$Freq/(B$Freq %>% sum())

N8_sub_stream <- ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title='')+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(2014, 2023))+
  scale_fill_manual(values = c('#95CE7E'))+
  theme_bw()+gg_theme


A <- meta_geno %>% filter(Subtype %in% c('H5N8', 'H5N1'))
B <- table(A$y, A$genotype) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq>0) #%>% mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
B <- B[B$Var2 %in% (table_DF(A$genotype) %>% filter(Freq > 50) %>% select(x) %>% unlist), ]
B$Freq <- B$Freq/(B$Freq %>% sum())
B$g <- str_extract(B$Var2, pattern = 'N[0-9]')

N8_stream <- 
  ggplot(B %>% filter(g=='N8'), aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title=paste('H5N8 genotypes', ' (Total genotypes: 39)', sep = ''), 
       subtitle = 'Frequency more than 50 are included')+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(2014, 2023))+
  scale_fill_brewer(palette = 'Set1')+
  theme_bw()+gg_theme

N1_stream <- 
  ggplot(B %>% filter(g=='N1'), aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25) +
  labs(x='Year', y='', fill='', title=paste('H5N1 genotypes', ' (Total genotypes: 130)', sep = ''), 
       subtitle = 'Frequency more than 50 are included')+
  scale_y_continuous(limits = c(-1, 1))+
  scale_x_continuous(limits = c(2014, 2023))+
  scale_fill_brewer(palette = 'Set1')+
  theme_bw()+gg_theme

# cowplot::plot_grid(plotlist = list(subtype_stream, N1_stream, N8_stream), ncol=1)
cowplot::plot_grid(plotlist = list(N8_sub_stream, N8_stream), ncol=1)
cowplot::plot_grid(plotlist = list(N1_sub_stream, N1_stream), ncol=1)
cowplot::plot_grid(plotlist = list(N8_stream, N1_stream), ncol=1)





meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'
meta_geno <- meta_geno %>% arrange(y) %>% filter(y > 1996) #, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
# meta_geno <- meta_geno[grepl(meta_geno$Collection_Date, pattern='[0-9]+-[0-9]+-[0-9]+'), ]
# meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())

A <- meta_geno %>% mutate(p=paste(Collection_Date, genotype)) %>% filter(HA %in% c('A', 'B', 'C')
                                                                         , Clade_new != c('IRD_unlabel')
                                                                         , Clade_new != c('nonGsGD'))
# A[!A$Clade_new %in% c('2.3.4.4b', '2.3.4.4'), 'Clade_new'] <- 'non-2.3.4.4'
# A[A$Clade_new %in% c('2.3.4.4b', '2.3.4.4'), 'Clade_new'] <- '2.3.4.4'
# A[!A$Clade_new %in% c('2.3.4.4b'), 'Clade_new'] <- 'non-2.3.4.4b'
# A[A$Clade_new %in% c('2.3.4.4b'), 'Clade_new'] <- '2.3.4.4b'
A[, 'time_group'] <- '2013~2016'
A[A$Collection_Date <=2013, 'time_group'] <- '≤2013'
A[A$Collection_Date >2016, 'time_group'] <- '>2016'
A$continent <- header_cleaning(A$Location, pattern = '\\/') %>% select(V1) %>% unlist

C <- table(A$time_group, A$Subtype, A$Clade_new) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
# C[, 3:5] <- header_cleaning(C$Var1, pattern = ' ')
C$Var1 <- factor(C$Var1, levels = c("≤2013",  '2013~2016', ">2016"))
C$Var3 <- factor(C$Var3, levels =  c('GsGD-others', '2.3.4.4', '2.3.4.4b'))

png('~/Analysis/aiv/merge/0307/plot/Genotypes_frequency_in_Clade.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1', 'H5N8'))) +
  geom_col(aes(x = Var3, y = Freq, fill=Var1), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(x='Clade', y='Frequency of genotypes', fill='Time internal', 
       title='Frequency of genotypes in subtype/clade/time interval')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()


C <- table(A$time_group, A$Subtype, A$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
C$Var1 <- factor(C$Var1, levels = c("≤2013",  '2013~2016', ">2016"))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'

png('~/Analysis/aiv/merge/0307/plot/Genotypes_frequency_in_continent.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1', 'H5N8'))) +
  geom_col(aes(x = Var3, y = Freq, fill=Var1), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(x='Continent', y='Frequency of genotypes', fill='Time internal', 
       title='Frequency of genotypes in subtype/continent/time interval')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()


C <- table(A$Clade_new, A$Subtype, A$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
# C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4', '2.3.4.4'))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'
C$Var1 <- factor(C$Var1, levels =  c('GsGD-others', '2.3.4.4', '2.3.4.4b'))

png('~/Analysis/aiv/merge/0307/plot/Genotypes_freq_in_Clade_continent.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1'))) +
  geom_col(aes(x = Var1, y = Freq, fill=Var3), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(fill='Continent', y='Frequency of genotypes', x='Clade', 
       title='H5N1 frequency of genotypes between clade and continent')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A', 'pink', 'lightblue'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  # facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()


%>% filter(Clade_new %in% '2.3.4.4b') #%>% filter(Subtype %in% 'H5N1')

C <- table(B$Clade_new, B$host_type, B$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
# C <- rbind(C, C[73:96, ] %>% mutate(Var3='SouthAmerica', Freq=0))
# C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4', '2.3.4.4'))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
# C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'
C$Var1 <- factor(C$Var1, levels =  c('GsGD-others', '2.3.4.4', '2.3.4.4b'))
ggplot(C %>% filter(Var3 != 'SouthAmerica')) +
  geom_col(aes(x = Var1, y = Freq, fill=Var2), position='stack', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(fill='Host type', y='Genotypes frequency', x='Clade', 
       title='H5N1 genotypes frequency between clade/continent/host type')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A', 'pink', 'lightblue'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var3) +
  theme_bw()+gg_theme







meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'
meta_geno <- meta_geno %>% arrange(y) %>% filter(y > 1996) %>% distinct(genotype, .keep_all = T)#, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
# meta_geno <- meta_geno[grepl(meta_geno$Collection_Date, pattern='[0-9]+-[0-9]+-[0-9]+'), ]
# meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())
A <- meta_geno %>% mutate(p=paste(Collection_Date, genotype)) %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))
# A[!A$Clade_new %in% c('2.3.4.4b'), 'Clade_new'] <- 'non-2.3.4.4b'
A[, 'time_group'] <- '2013~2016'
A[A$Collection_Date <=2013, 'time_group'] <- '≤2013'
A[A$Collection_Date >2016, 'time_group'] <- '>2016'
A$continent <- header_cleaning(A$Location, pattern = '\\/') %>% select(V1) %>% unlist


C <- table(A$time_group, A$Subtype, A$Clade_new) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
C$Var1 <- factor(C$Var1, levels = c("≤2013",  '2013~2016', ">2016"))
C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4b', '2.3.4.4b'))

png('~/Analysis/aiv/merge/0307/plot/Genotypes_type_in_Clade.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1', 'H5N8'))) +
  geom_col(aes(x = Var3, y = Freq, fill=Var1), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(x='Clade', y='Types of new emerging genotypes', fill='Time internal', 
       title='Types of new emerging genotypes in subtype/clade/time interval')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()


C <- table(A$time_group, A$Subtype, A$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
C <- rbind(C, C[73:96, ] %>% mutate(Var3='SouthAmerica', Freq=0))
C$Var1 <- factor(C$Var1, levels = c("≤2013",  '2013~2016', ">2016"))
# C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4', '2.3.4.4'))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'

png('~/Analysis/aiv/merge/0307/plot/Genotypes_type_in_continent.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1', 'H5N8'))) +
  geom_col(aes(x = Var3, y = Freq, fill=Var1), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(x='Continent', y='Types of new emerging genotypes', fill='Time internal', 
       title='Types of new emerging genotypes in subtype/continent/time interval')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()




meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'
meta_geno <- meta_geno %>% arrange(y) %>% filter(y > 1996) %>% distinct(genotype, .keep_all = T)#, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
# meta_geno <- meta_geno[grepl(meta_geno$Collection_Date, pattern='[0-9]+-[0-9]+-[0-9]+'), ]
# meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())
A <- meta_geno %>% mutate(p=paste(Collection_Date, genotype)) %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))
# A[!A$Clade_new %in% c('2.3.4.4b'), 'Clade_new'] <- 'non-2.3.4.4b'
A[, 'time_group'] <- '2013~2016'
A[A$Collection_Date <=2013, 'time_group'] <- '≤2013'
A[A$Collection_Date >2016, 'time_group'] <- '>2016'
A$continent <- header_cleaning(A$Location, pattern = '\\/') %>% select(V1) %>% unlist


C <- table(A$Clade_new, A$Subtype, A$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
C <- rbind(C, C[73:96, ] %>% mutate(Var3='SouthAmerica', Freq=0))
# C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4', '2.3.4.4'))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'
C$Var1 <- factor(C$Var1, levels =  c('GsGD-others', '2.3.4.4', '2.3.4.4b'))

png('~/Analysis/aiv/merge/0307/plot/Genotypes_type_in_Clade_continent.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1'))) +
  geom_col(aes(x = Var1, y = Freq, fill=Var3), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(fill='Continent', y='Types of new emerging genotypes', x='Clade', 
       title='H5N1 new emerging types  of genotypes between clade and continent')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A', 'pink', 'lightblue'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  # facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()


B <- A #%>% filter(Subtype %in% 'H5N1')

C <- table(B$Clade_new, B$host_type, B$continent) %>% as.data.frame(stringsAsFactors=F) #%>% filter(Freq!=0)
# C <- rbind(C, C[73:96, ] %>% mutate(Var3='SouthAmerica', Freq=0))
# C$Var3 <- factor(C$Var3, levels =  c('non-2.3.4.4', '2.3.4.4'))
C[C$Var3 %in% c('NorthAmerica'), 'Var3'] <- 'North\nAmerica'
# C[C$Var3 %in% c('SouthAmerica'), 'Var3'] <- 'South\nAmerica'
C$Var1 <- factor(C$Var1, levels =  c('GsGD-others', '2.3.4.4', '2.3.4.4b'))
ggplot(C) +
  geom_col(aes(x = Var1, y = Freq, fill=Var2), position='stack', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(fill='Host type', y='Types of new emerging genotypes', x='Host_y', 
       title='HPAI new emerging types of genotypes between clade/continent/host type')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A', 'pink', 'lightblue'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  facet_wrap(~ Var3) +
  theme_bw()+gg_theme



A <- meta_geno %>% filter(y >= 1996) %>% filter(Subtype %in% c('H5N1', 'H5N8', 'H7N9'), y >= 1996) %>% distinct(genotype, .keep_all = T) %>% 
  filter(Clade_new %in% c('2.3.4.4', '2.3.4.4b', 'GsGD-others'))
A$time_group <- ifelse(A$Collection_Date > 2013, yes = '>2013', no='≤2013')
A$p <- paste(A$time_group, A$host_type)

B <- table(A$time_group, A$host_type) %>% as.data.frame(stringsAsFactors=as.data.frame(stringsAsFactors=F) %>% filter(Freq!=0) %>% 
                                                          mutate(p=paste(Var1, Var2)) F) %>% filter(Freq!=0) %>% 
  mutate(p=paste(Var1, Var2)) 

C <- table(A$genotype, A$p) %>% as.data.frame(stringsAsFactors=F) %>% filter(Freq!=0)
B <- merge(B, table_DF(C$Var2) %>% rename('Type'='Freq'), by.x = 'p', by.y = 'x')
B <- B[order(B$Var2), ]
B$Var1 <- factor(B$Var1, levels = c("≤2013", ">2013"))
B$div <- (B$Freq/B$Type) %>% round(digits = 2)

ggplot(B) +
  geom_bar(aes(x = Var2, y = Type, fill=Var1), position='dodge', stat='identity', color='black')+
  geom_text(aes(x=Var2, y=0, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7, vjust=1.5)+
  labs(x='Host type', y='Genotypes', fill='Time\ninternal', 
       title='Genotypes in Host type/time interval')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E'))+
  scale_y_continuous(breaks = seq(0, max(B$Freq), 100))+
  theme_bw()+gg_theme


A <- meta_geno %>% filter(Subtype %in% c('H5N1'), y >= 1996)


png('~/Analysis/aiv/merge/0307/plot/Genotypes_type_in_Clade_continent.png', width = 9600, height = 5400, res = 600)
ggplot(C %>% filter(Var2 %in% c('H5N1'))) +
  geom_col(aes(x = Var1, y = Freq, fill=Var3), position='dodge', color='black')+
  # geom_text(aes(x=Var2, y=-2, label=Freq, fill=Var1), position = position_dodge(width = .9), size=7)+
  labs(fill='Continent', y='Types of new emerging genotypes', x='Clade', 
       title='H5N1 frequency of genotypes between clade and continent')+
  scale_fill_manual(values = c('#958BE3', '#95CE7E', '#E9967A', 'pink', 'lightblue'))+
  # scale_y_continuous(breaks = seq(0, 75, 5))+
  # scale_fill_manual(values = c('#877ECE','#CEA77E','#95CE7E'))+
  # facet_wrap(~ Var2) +
  theme_bw()+gg_theme
dev.off()





# To France ---------------------------------------------------------------
meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'

meta_geno <- meta_geno %>% arrange(y) #%>% filter(y > 1996) #%>% distinct(genotype, .keep_all = T)#, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
# meta_geno <- meta_geno[(meta_geno$Collection_Date %>% str_remove_all(pattern = '[0-9]+')=='--'), ]
meta_geno <- meta_geno[grepl(meta_geno$Collection_Date, pattern='[0-9]+-[0-9]+-[0-9]+'), ]
meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())

A <- meta_geno %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))


out <- data.frame(g='Frequency     ', y=A$Collection_Date %>% as.numeric()) #%>% filter(y >2010)

# out <- rbind(out, data.frame(g='H5N1', y=meta_geno %>% filter(Subtype=='H5N1') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='H5N8', y=meta_geno %>% filter(Subtype=='H5N8') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4', y=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Collection_Date) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4b', y=meta_geno %>% filter(Clade_new=='2.3.4.4b', y>2016) %>% select(Collection_Date) %>% unlist() %>% as.numeric()))

frequency_violin <-
  ggplot(out , aes(x=y, y=g, fill=g))+
  geom_violin()+
  # geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_x_continuous(breaks = seq(2003, 2023, 1))+
  scale_fill_manual(values = '#034FDE')+
  labs(x='Year', y='', fill='', title='Genotype frequency')+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank(), axis.title.x = element_blank())


# B <- table(A$Collection_Date %>% round(digits = 1), A$Subtype) %>% as.data.frame(stringsAsFactors=F) %>%
#   filter(Freq!=0) %>% 
#   mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
# # B$Var2 <- paste(B$Var2, paste(rep(' ', 7), collapse = ''), sep = '')
# B$Var2 <- 'Frequency'
# B$Freq1 <- B$Freq/(B$Freq %>% sum())
# 
# frequency_violin <-
#   ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
#   geom_stream(color = 1, lwd = 0.25)+ 
#   labs(x='Year', y='', fill='', title='Genotype frequency')+
#   # scale_fill_manual(values = '#034FDE')+
#   # scale_y_continuous(limits = c(-0.01, 0.01))+
#   scale_x_continuous(breaks = seq(2003, 2023, 1))+
#   theme_bw()+gg_theme+theme(axis.text.y = element_blank())



library(ggstream)

meta_geno_sub <- meta_geno %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))
B <- table(meta_geno_sub$Collection_Date %>% round(digits = 1), meta_geno_sub$Subtype) %>% as.data.frame(stringsAsFactors=F) %>%
  filter(Freq!=0) %>% 
  mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
B$Var2 <- paste(B$Var2, paste('             '), sep = '')
B$Freq1 <- B$Freq/(B$Freq %>% sum())

subtype_stream <-
  ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25)+ 
  labs(x='Year', y='', fill='', title='GsGd Subtypes temporal evoluation')+
  # scale_y_continuous(limits = c(-0.01, 0.01))+
  scale_x_continuous(breaks = seq(2003, 2023, 1))+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank(), axis.title.x = element_blank())


meta_geno_sub <- meta_geno %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))#Clade_new %in% c('2.3.4.4b', '2.3.4.4', 'GsGD-others'))
B <- table(meta_geno_sub$Collection_Date %>% round(digits = 1), meta_geno_sub$Clade_new) %>% as.data.frame(stringsAsFactors=F) %>%
  filter(Freq!=0) %>% 
  mutate(Var1=as.numeric(Var1)) #%>% group_split(Var2) %>% as.list() %>% lapply(as.data.frame)
B$Freq1 <- B$Freq/(B$Freq %>% sum())

clade_stream <-
  ggplot(B, aes(x = Var1 %>% as.numeric(), y = Freq, fill = Var2)) +
  geom_stream(color = 1, lwd = 0.25)+ 
  labs(x='Year', y='', fill='', title='GsGd clade temporal evoluation')+
  # scale_y_continuous(limits = c(-0.01, 0.01))+
  scale_x_continuous(breaks = seq(2003, 2023, 1))+
  scale_fill_brewer(palette = 'Set3')+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank(), axis.title.x = element_blank())



meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id')
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'

meta_geno <- meta_geno %>% arrange(y) %>% filter(y > 1996) %>% distinct(genotype, .keep_all = T)#, Clade_new %in% c('2.3.4.4', '2.3.4.4b')
meta_geno <- meta_geno[grepl(meta_geno$Collection_Date, pattern='[0-9]+-[0-9]+-[0-9]+'), ]
meta_geno$Collection_Date <- lubridate::decimal_date(meta_geno$Collection_Date %>% as.Date())

A <- meta_geno %>% mutate(p=paste(Collection_Date, genotype)) %>% filter(HA %in% c('A', 'B', 'C'), Clade_new != c('IRD_unlabel'))

out <- data.frame(g='Type                   ', y=A$Collection_Date %>% as.numeric())# %>% filter(y >2010)
# out <- rbind(out, data.frame(g='H5N1', y=meta_geno %>% filter(Subtype=='H5N1') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='H5N8', y=meta_geno %>% filter(Subtype=='H5N8') %>% select(y) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4', y=meta_geno %>% filter(Clade_new=='2.3.4.4') %>% select(Collection_Date) %>% unlist() %>% as.numeric()))
# out <- rbind(out, data.frame(g='2.3.4.4b', y=meta_geno %>% filter(Clade_new=='2.3.4.4b', y>2016) %>% select(Collection_Date) %>% unlist() %>% as.numeric()))

type_violin <- 
  ggplot(out , aes(x=y, y=g, fill=g))+
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.1), show.legend=F)+
  scale_x_continuous(breaks = seq(2003, 2023, 1))+
  scale_fill_manual(values = '#A6BEE7')+
  labs(x='Year', y='', fill='', title='Types of new emerging genotype')+
  theme_bw()+gg_theme+theme(axis.text.y = element_blank())#+guides(fill='none')


png('~/Analysis/aiv/merge/0307/plot/Stream_violin.png', width = 9600, height = 5400, res = 600)
cowplot::plot_grid(plotlist = list(clade_stream, subtype_stream, frequency_violin, type_violin), ncol=1)
dev.off()



# result numbers of virus and genotypes -----------------------------------
meta <- fread(file = '~/Analysis/aiv/merge/0307/gisaid_ird_merge_meta_information.csv') %>% as.data.frame()
meta[is.na(meta$segment), 'segment'] <- 'NA'

out <- table_DF(table_DF(meta$Isolate_Id)[, 2]) %>% mutate(x=as.numeric(x)) %>% arrange(x)
out[9, ] <- c('≥9', out[9:nrow(out), 2] %>% sum())
out <- out[1:9, ] %>% mutate(x=factor(x, levels=c(1:8, '≥9')))

png('~/Analysis/aiv/merge/0307/plot/virus_segment_distb.png', width = 9600, height = 5400, res = 600)
ggplot(out)+
  geom_bar(mapping=aes(x=x, y=as.numeric(Freq)), fill='#A6BEE7', stat='identity', color='black')+
  geom_text(mapping=aes(x=x, y=as.numeric(Freq), label=as.numeric(Freq)), vjust=-0.3, size=7)+
  labs(x='Numbers of genome', y='Frequency', title='Virus genome distribution', 
       subtitle = 'Total virus: 57,183')+
  theme_bw()+gg_theme
dev.off()


meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id') %>% 
  filter(Isolate_Id %in% (table_DF(meta$Isolate_Id) %>% filter(Freq>=8) %>% select(x) %>% unlist()))
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno <- meta_geno %>% arrange(y) %>% distinct(genotype, .keep_all = T)
meta_geno$HL <- ifelse(meta_geno$HA %in% c('A', 'B', 'C'), yes = 'HPAI', no = 'LPAI')

out <- table_DF(table_DF(meta_geno$genotype)[, 2]) %>% mutate(x=as.numeric(x)) %>% arrange(x)
out <- table_DF(meta_geno$genotype) %>% arrange(Freq)


out <- table(meta_geno$genotype, meta_geno$y, meta_geno$HL) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0) %>% arrange(Freq)

# genotype_each_year <- table_DF(out$Var2) %>% filter(x >= 1996)
genotype_each_year <- table(out$Var2, out$Var3) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0) %>% arrange(Freq) %>% filter(Var1 >= 1996)

png('~/Analysis/aiv/merge/0307/plot/genotype_year_distb.png', width = 9600, height = 5400, res = 600)
ggplot(genotype_each_year)+
  geom_bar(mapping=aes(y=Var1, x=as.numeric(Freq), fill=Var2), position='dodge', stat='identity', color='black')+
  # geom_text(mapping=aes(y=Var2, x=as.numeric(Freq), label=as.numeric(Freq)), hjust=-0.1, size=7)+
  labs(x='Types of genotype', y='Year', title='Genotypes and their first identitfied year\n(1996~2023, complete virus)', 
       subtitle = 'Total types of genotype: 3219', fill='')+
  scale_fill_brewer(palette = 'Set2')+
  scale_x_continuous(breaks = seq(0, 300, 25))+
  theme_bw()+gg_theme
dev.off()



meta_geno <- merge(genotype, meta %>% distinct(Isolate_Id, .keep_all = T), by = 'Isolate_Id') %>% 
  filter(Isolate_Id %in% (table_DF(meta$Isolate_Id) %>% filter(Freq>=8) %>% select(x) %>% unlist()))
meta_geno[, 'Clade_new'] <- 'nonGsGD'
for(k in seq_len(nrow(clade_table)))  {
  meta_geno[meta_geno$Clade %in% clade_table[k, 'x'], 'Clade_new'] <- clade_table[k, 'Clade']
}
meta_geno[meta_geno$HA %in% 1:16, 'Clade_new'] <- 'nonGsGD'
meta_geno <- meta_geno[!grepl(meta_geno$genotype, pattern='X'), ]
meta_geno$y <- meta_geno$Collection_Date %>% header_cleaning('-') %>% select(V1) %>% unlist()
meta_geno <- meta_geno %>% arrange(y) %>% distinct(genotype, .keep_all = T) %>% 
  filter(HA %in% c('A', 'B', 'C'), Clade_new != 'IRD_unlabel')


out <- table_DF(table_DF(meta_geno$genotype)[, 2]) %>% mutate(x=as.numeric(x)) %>% arrange(x)
out <- table_DF(meta_geno$genotype) %>% arrange(Freq)


out <- table(meta_geno$genotype, meta_geno$y, meta_geno$Clade_new) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0) %>% arrange(Freq)
out <- table(meta_geno$genotype, meta_geno$y, meta_geno$Subtype) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0) %>% arrange(Freq)

genotype_each_year <- table(out$Var2, out$Var3) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0) %>% arrange(Freq)
ggplot(genotype_each_year)+
  geom_bar(mapping=aes(y=Var1, x=as.numeric(Freq), fill=Var2), position='stack', stat='identity', color='black')+
  # geom_text(mapping=aes(y=Var1, x=as.numeric(Freq), label=as.numeric(Freq), fill=Var2), hjust=-0.1, size=7)+
  labs(x='Types of genotype', y='Year', title='Genotypes and their emerging year', 
       subtitle = 'Total types of genotype: 3464')+
  scale_x_continuous(breaks = seq(0, 45, 5))+
  theme_bw()+gg_theme















