for(i in class$Internal) {
  lineage <- eight[[i]]
  rownames(lineage) <- lineage$Strain_number
  
  rank <- ape::read.nexus(paste('~/Analysis/aiv/merge/0307/iq_tree/gisaid_IRD_merged_', i, '_aligned_iqtree.nexus', sep = ''))$tip.label %>% 
    str_remove_all(pattern = "'") %>% str_extract(pattern = '[0-9]+_H') %>% str_remove(pattern = '_H')
  
  lineage <- lineage[rank, ]
  # lineage_split <- lineage[, 2] %>% table_DF() %>% filter(Freq > nrow(lineage)*0.005) %>% select(x) %>% unlist()
  # lineage[!lineage[, 2] %in% lineage_split, 'clades'] <- NA
  
  
  sample_table <- lineage[, 2] %>% table_DF(prop = T) %>% 
    mutate(sample=round(Freq*nrow(lineage)*0.05, digits = 0), all=Freq*nrow(lineage)) %>% 
    tibble::column_to_rownames('x') %>% arrange(desc(all)) %>% mutate(cum=cumsum(Freq), Freq=100*Freq) %>% 
    filter(cum < 0.95)
  print(i)
  print(sample_table)
  
  sample_table <- lineage[, 2] %>% table_DF(prop = T) %>% 
    mutate(sample=round(Freq*nrow(lineage)*0.05, digits = 0), all=Freq*nrow(lineage)) %>% 
    tibble::column_to_rownames('x') %>% arrange(desc(all)) %>% mutate(cum=cumsum(Freq), Freq=100*Freq) %>% 
    filter(Freq >= 1)
  
  print(i)
  print(sample_table)
  # %>% 
  #   filter(cum < 0.95)
  # 
  
  
  align <- read.fasta(paste('~/Analysis/aiv/merge/0307/aligned_fasta/gisaid_IRD_merged_', i, '_aligned.fa', sep = ''), as.string = T) %>%
    do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('header') %>% rename('align'='V1') %>%
    mutate(align=str_replace_all(align, pattern = 'u|p|x|e|i|o|l', replacement = 'n'))
  align$Strain_number <- header_cleaning(align$header, pattern='\\|') %>% select('V2') %>% unlist
  
  align$aa <-
    Biostrings::DNAStringSet(align$align %>% str_remove_all(pattern = '-')) %>%
    Biostrings::translate(if.fuzzy.codon = 'X') %>% as.character()
  
  lineage <- merge(lineage, align, by = 'Strain_number')
  
  lineage$p <- paste(lineage$clades, lineage$aa, sep = '_')
  lineage$new_header <- paste(lineage$new_header, '|', lineage$clades, '_', sep = '')
  rownames(lineage) <- lineage$Strain_number
  lineage <- lineage[rank, ]
  
  subset <- lineage %>% filter(!is.na(clades)) %>%
    distinct(p, .keep_all = T) %>% group_split(clades) %>% as.list() %>% lapply(as.data.frame)
  names(subset) <- lapply(subset, FUN = function(x){x[1, 'clades']}) %>% unlist()
  
  
  sample <- list()
  for(j in rownames(sample_table)) {
    A <- subset[[j]][, c('new_header', 'align')]
    set.seed(713)
    sample[[j]] <- A[sample(x = 1:nrow(A), size = sample_table[j, 2], replace = F), ]
  }
  sample <- do.call(rbind, sample)
  # 
  # 
  # write_fasta(sequence = sample$align, header = sample$new_header
  #             , file_out = paste('/home/eric/Analysis/aiv/merge/0307/sampling_subset/', i, '_5p.fasta', sep = ''))
}


A <- meta[, -c(6, 16, 17, 21, 22, 23)]
colnames(A)
write.table(A, '/home/eric/Analysis/aiv/merge/0307/seqeunce_meta_information_preprocessed.csv', sep = ',', quote = F, row.names = F)

B <- meta_raw[, -c(6, 16, 17, 21)]
colnames(B)
write.table(B, '/home/eric/Analysis/aiv/merge/0307/seqeunce_meta_information_raw.csv', sep = ',', quote = F, row.names = F)

colN <- c('Isolate_Id', 'header', 'Subtype', 'Clade', 'Location', 'Collection_Date', 'Host', "scentific_name", "H_order" , "H_family",
          'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'genotype', 'complete', 'multiples_geno')
C <- meta_geno[, colN] %>% arrange(Collection_Date)
colnames(C)[16] <- 'NA_g'
write.table(C, '/home/eric/Analysis/aiv/merge/0307/virus_genotype_meta_information.csv', sep = '\t', quote = F, row.names = F)

A <- fread('/home/eric/Analysis/aiv/merge/0307/virus_genotype_meta_information.csv')


A <- organ_subtype(meta)
out_fasta_dir <- "~/Analysis/aiv/merge/0307/raw_fasta/"
dir.create(out_fasta_dir)
for(i in names(A))  {
  B <- A[[i]]
  write_fasta(sequence = B$seq, header = B$new_header, file_out = paste(out_fasta_dir, i, '_raw.fa', sep = ''))
}

meta$new_header <- 
  paste(meta$Isolate_Id, meta$Strain_number, 
        meta$Subtype %>% str_remove(pattern = 'A/'), 
        meta$Clade, 
        meta$segment, meta$Collection_Date,
        paste(word(meta$Location, 1, sep = '/'), meta$country, sep = '/'), 
        meta$Host, meta$strain_host, meta$host_type, sep = '|') 



align <- read_as_list('/home/eric/Analysis/aiv/merge/0307/ORF_filter_outlier/', 
                      prefix = '_aligned_fill_trimed_remove_outlier.fa$', file_type = 'fasta', as.string = T)


# PB1, PA, NP, MP sampling  -----------------------------------------------
outdir <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0922/'
dir.create(outdir)
for(i in c(class$Internal[-c(1, 6)])) {
  a <- align[grep(names(align), pattern=i)] %>% do.call(what=cbind) %>% as.data.frame()
  colnames(a) <- 'sequence'
  a <- cbind(header_cleaning(rownames(a), pattern = '\\|') %>% select(V2), a) %>% rename('Strain_number'='V2') %>%
    mutate(sequence=as.character(sequence))
  rownames(a) <- 1:nrow(a)
  a <- a %>% mutate(sequence=str_replace_all(sequence, pattern = 'u|p|x|e|i|o|l', replacement = 'n'))
  
  a$aa <-
    Biostrings::DNAStringSet(a$sequence %>% str_remove_all(pattern = '-')) %>%
    Biostrings::translate(if.fuzzy.codon = 'X') %>% as.character()
  
  
  info <- eight[[i]] %>% mutate(new_header=paste(new_header, clades, sep = '|'))
  seg <- merge(info[, c('Strain_number', 'clades', 'new_header')], a) %>% 
    mutate(p=paste(clades, aa)) %>% distinct(p, .keep_all = T)
  sampling <- table_DF(info$clades) %>% mutate(sample=round(Freq*0.05, digits = 0)) %>% tibble::column_to_rownames('x')
  seg <- seg %>% group_split(clades) %>% as.list() %>% lapply(as.data.frame)
  names(seg) <- lapply(seg, function(x){x[1, 2]}) %>% unlist()
  out <- list()
  for(s in names(seg)) {
    subseg <- seg[[s]]
    set.seed(713)
    out[[s]] <- subseg[sample(x=1:nrow(subseg), size=sampling[s, 'sample'], replace = F), ]
  }
  out <- out %>% do.call(what=rbind) %>% as.data.frame()
  write_fasta(sequence = out$sequence, header = out$new_header, file_out = paste(outdir, i, '_5p_sampling.fa', sep = ''))
}



# internal sampling  -----------------------------------------------
outdir <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/0922/'
dir.create(outdir)
for(i in c(class$Internal)) {
  a <- align[grep(names(align), pattern=i)] %>% do.call(what=cbind) %>% as.data.frame()
  colnames(a) <- 'sequence'
  a <- cbind(header_cleaning(rownames(a), pattern = '\\|') %>% select(V2), a) %>% rename('Strain_number'='V2') %>%
    mutate(sequence=as.character(sequence))
  rownames(a) <- 1:nrow(a)
  a <- a %>% mutate(sequence=str_replace_all(sequence, pattern = 'u|p|x|e|i|o|l', replacement = 'n'))
  
  a$aa <-
    Biostrings::DNAStringSet(a$sequence %>% str_remove_all(pattern = '-')) %>%
    Biostrings::translate(if.fuzzy.codon = 'X') %>% as.character()
  
  seg <- merge(eight[[i]][, c('Strain_number', 'new_header')], a) %>%
    distinct(aa, .keep_all = T)
  # set.seed(713)
  # out <- seg[sample(x=1:nrow(seg), size=nrow(info)*0.05, replace = F), ]
  
  write_fasta(sequence = out$sequence, header = out$new_header, file_out = paste(outdir, i, '_5p_sampling.fa', sep = ''))
}



# identicle aa ------------------------------------------------------------
outdir <- '/home/eric/Analysis/aiv/merge/0307/sampling_subset/1015/'
dir.create(outdir)
for(i in c(class$Internal)) {
  a <- align[grep(names(align), pattern=i)] %>% do.call(what=cbind) %>% as.data.frame()
  colnames(a) <- 'sequence'
  a <- cbind(header_cleaning(rownames(a), pattern = '\\|') %>% select(V2), a) %>% rename('Strain_number'='V2') %>%
    mutate(sequence=as.character(sequence))
  rownames(a) <- 1:nrow(a)
  a <- a %>% mutate(sequence=str_replace_all(sequence, pattern = 'u|p|x|e|i|o|l', replacement = 'n'))
  
  a$aa <-
    Biostrings::DNAStringSet(a$sequence %>% str_remove_all(pattern = '-')) %>%
    Biostrings::translate(if.fuzzy.codon = 'X') %>% as.character()
  
  
  seg <- merge(eight[[i]][, c('Strain_number', 'new_header')], a) %>%
    distinct(aa, .keep_all = T)
  # remove <- c()
  # for(i in 1:nrow(seg)) {
  #   same <- grep(seg$aa, pattern=str_replace_all(seg[i, 'aa'], pattern = '\\*', replacement = '\\\\*'))
  #   # same <- grep(seg$aa, pattern=str_replace_all(seg[2162, 'aa'], pattern = '\\*', replacement = '\\\\*'))
  #   if(length(same) != 1){
  #     remove <- c(remove, (seg[same, 'aa'] %>% nchar() %>% which.min()))
  #   }else{next}
  # }
  # set.seed(713)
  # out <- seg[sample(x=1:nrow(seg), size=nrow(info)*0.05, replace = F), ]
  
  write_fasta(sequence = seg$sequence, header = seg$new_header, file_out = paste(outdir, i, '_aa_sampling.fa', sep = ''))
}
