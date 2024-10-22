suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidytree))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(TDbook))
suppressPackageStartupMessages(library(aplot))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(smplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(treemap))


seg_sub_level <- c('PB2', 'PB1','PA', paste('H', seq(1:16), sep = ''), 'NP', paste('N', seq(1:9), sep = ''), 'MP', 'NS')


read_as_list <- function(path, sep='auto', prefix='', header='auto', as.string = F, cols=NA, file_type='txt') {
  out_list <- list()
  file_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  # file_type <- sub(str_extract(file_path[1], pattern = '\\..*'), pattern = '.', replacement = '')
  file_path <- as_fs_path(file_path)
  if(file_type == 'nexus')  {
    for(i in file_path) {
      out_list[[i]] <- ape::read.nexus(i)
    }
  }
  if(file_type == 'nwk')  {
    for(i in file_path) {
      out_list[[i]] <- ape::read.tree(i)
    }
  }
  if(file_type == 'fasta')  {
    # out_list <- file_path %>% map(.f = function(path){read.fasta(path, as.string = F, whole.header = T)})
    for(i in file_path) {
      out_list[[i]] <- read.fasta(i, as.string = as.string, whole.header = T)
    }
  }
  if(file_type == 'RData')  {
    # out_list <- file_path %>% map(.f = function(path){readRDS(path)})
    for(i in file_path) {
      out_list[[i]] <- readRDS(i)
    }
  }
  if(file_type == 'csv')  {
    out_list <- file_path %>% map(.f = function(path){read.csv(path, header = T, sep = ',')})
  }
  if(file_type == 'txt')  {
    for(i in file_path) {
      out_list[[i]] <- fread(i, header = header, stringsAsFactors = F, sep = sep) %>% as.data.frame()
    }
  }
  names(out_list) <- file_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)  {
    names(out_list) <- sub(names(out_list), pattern = pattern[i], replacement = '')
  }
  return(out_list)
}




organ_subtype <- function(in_seq)  {
  in_seq$sub <- in_seq$segment
  # colnames(in_seq)[subtype_col] <- 'segment'
  A <- group_split(in_seq, segment) %>% as.list() %>% lapply(as.data.frame)
  names(A) <- lapply(A, function(x) {x[1,'segment']})
  
  # %>% group_split(segment) %>% as.list() %>% lapply(as.data.frame)
  ## NA subtype
  if(any(names(A) %in% 'NA')) {
    B <- A[['NA']] 
    B$sub <- stringr::str_extract(B$Subtype, pattern = 'N.*')
    B <- B[!(is.na(B$sub)), ]
    na_sub <- B %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame)
    names(na_sub) <- lapply(na_sub, function(x) {x[1,'sub']})
    A <- A[-which(names(A) %in% c('NA'))]
  }else(na_sub='')
  
  
  if(any(names(A) %in% 'HA')) {
    B <- A[['HA']]
    B$sub <- stringr::str_extract(B$Subtype, pattern = 'H[0-9]+')
    B <- B[!(is.na(B$sub)), ] 
    
    ha_sub <- B %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame)
    names(ha_sub) <- lapply(ha_sub, function(x) {x[1,'sub']})
    A <- A[-which(names(A) %in% c('HA'))]
  }else(ha_sub='')
  
  test <- c(A, na_sub, ha_sub)
  
  return(test) #[!test == '']
}


table_DF <- function(x, prop=F) {
  A <- table(x) %>% as.data.frame(stringsAsFactors=F)
  if(isTRUE(prop))  {
    A <- merge(table(x) %>% as.data.frame(stringsAsFactors=F), table(x) %>% prop.table() %>% as.data.frame(stringsAsFactors=F), 
               by='x', all=T)
    colnames(A) <- c('x', 'Freq', 'Prop')
    
  }
  return(A)
}

strain_number_extract <- function(header) {
  str_extract(header, '[0-9]+_H|[0-9]+\\|H') %>% str_remove('_H|\\|H')
}


seg_level <- c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS')
seg_sub_level <- c('PB2', 'PB1','PA', paste('H', seq(1:16), sep = ''), 'NP', paste('N', seq(1:9), sep = ''), 'MP', 'NS')


header_cleaning <- function(header_vector, pattern='/') {
  test <- str_split(header_vector, pattern = pattern) %>% do.call(what=rbind) %>% as.data.frame()
  A <- lapply(str_split(header_vector, pattern = pattern), length) %>% unlist() %>% as.data.frame() %>% 
    mutate(index=seq(1, nrow(test)))
  for(i in unique(A$.)) {
    if(identical(i, ncol(test))) {next}
    test[A[A$. %in% i, 'index'], (i+1):ncol(test)] <- ''
  }
  return(test)
}

assign_geno <- function(group_list, meta, remove_homo=T)  {
  eight <- list()
  # group_out
  for(i in c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'))  {
    if((i %in% c('HA', 'NA')))  {
      p <- sub(pattern = 'A', replacement = '[0-9]+', x = i)
      A <- group_list[grep(names(group_list), pattern=p)] %>% do.call(what=rbind) %>% as.data.frame()
      sub_meta <- meta[meta$segment %in% i, ]
      A <- merge(A, sub_meta, by = 'ui', all = T)
      if(identical(i, 'HA'))  {
        A[is.na(A$s_g), c('group', 's_g')] <- A[is.na(A$s_g), c('sub')] %>% str_extract(pattern = '[0-9]+')
      }
      if(identical(i, 'NA'))  {
        A[, c('group', 's_g')] <- A[, c('s_g')]
      }
    } 
    if(!(i %in% c('HA', 'NA'))) {
      A <- group_list[[i]]
      sub_meta <- meta[meta$segment %in% i, ]
      A <- merge(A, sub_meta, by = 'ui', all.x = T)
    }
    A <- A %>% distinct(Isolate_Id, .keep_all = T) # remove duplicate epi
    eight[[i]] <- A[, -c(9)]
  }
  # return(eight)
  
  # A <- Reduce(intersect, lapply(eight, function(x) {x[, 'Isolate_Id']}))
  # eight <- lapply(eight, function(x)  {x[x[, 'Isolate_Id'] %in% A, ]})
  # lapply(eight, function(x)  {table(is.na(x[, 2]))})
  
  genotype <- data.frame(Isolate_Id=meta$Isolate_Id %>% unique())
  for(i in names(eight))  {
    A <- eight[[i]]
    colnames(A)[2] <- i
    genotype <- merge(genotype, A[, c('Isolate_Id', i)], by = 'Isolate_Id', all.x=T)
    genotype[, i] <- replace(genotype[, i], is.na(genotype[, i]), 'X')
  }
  genotype <- genotype[!is.na(genotype$Isolate_Id), ]
  
  genotype$genotype <- paste(genotype$PB2, genotype$PB1, genotype$PA, genotype$HA,
                             genotype$NP, genotype$`NA`, genotype$MP, genotype$NS, sep = '_')
  return(genotype)
}


write_fasta <- 
  function(sequence, header, file_out) {
    out <- matrix(0,0,0) %>% as.data.frame()
    if(identical(stringr::str_extract(header[1], pattern = '>'), '>')) {
      next
    }else{header <- paste('>', header, sep = '')}
    
    out[seq(1, length(sequence)*2, 2), 1] <- header
    out[seq(2, length(sequence)*2, 2), 1] <- sequence
    
    write.table(out, file = file_out, quote = F, col.names = F, row.names = F)
  }

strain_number_extract <- function(header) {
  str_extract(header, '[0-9]+_H|[0-9]+\\|H') %>% str_remove('_H|\\|H')
}


add_p_annotation <- function(p, value, x_pos, y_pos, size=5, linewidth=0.5) {
  p + 
    annotate("text", x = x_pos, y = y_pos + 5, label = format(value, digits = 3, scientific = TRUE), size = size) +
    geom_segment(aes(x = x_pos - 0.2, xend = x_pos + 0.2, y = y_pos, yend = y_pos), linewidth = linewidth) +
    geom_segment(aes(x = x_pos - 0.2, xend = x_pos - 0.2, y = y_pos, yend = y_pos - 4), linewidth = linewidth) +
    geom_segment(aes(x = x_pos + 0.2, xend = x_pos + 0.2, y = y_pos, yend = y_pos - 4), linewidth = linewidth)
}

extract_leaf_nodes <- function(tree) {
  cat("Extracting leaf nodes for each internal node...\n")
  Ntips <- length(tree$tip.label)
  Nnodes <- tree$Nnode
  
  node_leafe <- vector("list", Nnodes)  # Pre-allocate list for performance
  
  for (i in 1:Nnodes) {
    n_node <- Ntips + i
    n_desc <- getDescendants(tree, n_node)
    
    # Filter to include only tips
    leaf_indices <- n_desc[n_desc <= Ntips]
    
    # Extract and clean tip labels
    node_leafe[[i]] <- tree$tip.label[leaf_indices] %>% 
      str_remove_all(pattern = "'")
  }
  
  names(node_leafe) <- as.character(1:Nnodes)
  
  return(node_leafe)
}

gg_theme= theme(
  title = element_text(size = 30-8),
  axis.title.x = element_text(size = 28-4),
  axis.text.x = element_text(size = 18), # -6, angle=90,hjust=0.95,vjust=0.2
  axis.title.y = element_text(size = 28-4),
  axis.text.y = element_text(size = 18),
  plot.subtitle = element_text(size = 24),
  plot.caption = element_text(size = 30),
  legend.title = element_text(size = 20), 
  legend.text = element_text(size = 16+2), 
  legend.key.size = unit(3, 'lines'), 
  # legend.key.height = unit(5, "cm"),
  strip.text = element_text(size = 20), 
  strip.background = element_blank())



phylo_theme <- theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     legend.key.height = unit(0.15, 'cm'),
                     legend.key.width= unit(1.5, 'cm'),
                     axis.title.y = element_text(face = "bold",size = 15),
                     axis.title.x = element_text(face = "bold",size = 15),
                     legend.text = element_text(size=20-8),
                     legend.title = element_text(size=22-8))

# 
custom_color <- c("#FF7751", "#66E4FF", "#9D66FF", "#FF66C1", "#66FFB2", "#FF66F7", "#DFFF66",
                       "#66FF66", "#66B3FF", "#FFC266", "#66FFE4", "#9DFF66",
                       "#FF6685", "#6686FF", "#FFF566", "#66FF86", "#CF66FF", "#FF6666", "#6686FF")

custom_color <- c("#FFB3A1", "#A1E8FF", "#B3A1FF", "#FFA1E4", "#A1FFCF", "#FFA1FA", "#E8FFB3",
                           "#B3FFA1", "#A1CFFF", "#FFD4A1", "#A1FFF3", "#CFFA1A",
                           "#FFA1B3", "#A1A6FF", "#FFF9A1", "#A1FFB3", "#D4A1FF", "#FFA1A1", "#A1B3FF")

custom_color <- c("#DFFF66", "#51D4FF", "#9DFF66", "#FF51B2", "#A1FFCF", "#66B3FF", "#FF7751",
                  "#517FFF", "#66FF86", "#FFB751", "#51FFDA", "#D4A1FF",
                  "#FF5177", "#5171FF", "#FFF551", "#B3FFA1", "#C051FF", "#FF5151")



