.libPaths(c("/home/eric/R/x86_64-pc-linux-gnu-library/4.3/", "/usr/lib/R/library/"))
source('~/R/aiv/function.R')
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
args$meta <- '~/Analysis/aiv/gisaid_all/meta/'
args$gisaid_seq_fasta <- '~/Analysis/aiv/gisaid_all/all.fasta'
args$ird_seq_fasta <- '/home/eric/Analysis/aiv/ird/IRD_Sequence.fa'
args$length_down_limit <- 0.8
args$known_outlier <- 0.8
args$out_fasta_prefix <- 'gisaid_IRD_merged_'
args$out_fasta_dir <- '~/Analysis/aiv/merge/0307/'
cat(paste('End of sequence preprocessing', sep = ''))

meta <- read_as_list(path = args$meta , prefix = '.csv', header = T, file_type = 'txt') %>% 
  # lapply(function(x){select(x, Isolate_Id, Subtype, Clade, Location, Host, Collection_Date, Submission_Date)}) %>% 
  do.call(what=rbind) %>% distinct(Isolate_Id, .keep_all = T) %>% mutate()

ncbi <- meta[, grep(colnames(meta), pattern = 'INSDC')] %>% unlist %>% str_split(pattern = ', ') %>% unlist %>% unique()

meta <- meta %>% select(Isolate_Id, Subtype, Clade, Location, Host, Collection_Date, Submission_Date)

seq <- read.fasta(args$gisaid_seq_fasta, as.string = T) %>% 
  do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('header')

seq$epi <- seq_name <- str_extract(seq$header, pattern = 'EPI_ISL_[0-9]+')
seq$segment <- str_extract(seq$header, pattern = '\\|EPI_ISL_.*') %>% sub(pattern='\\|EPI', replacement='EPI') %>% 
  str_remove(pattern = 'EPI_ISL_[0-9]*\\|') %>% 
  str_remove(pattern = '\\|[0-9]')



# organize header
inform <- str_split(seq$header, pattern = '\\|') %>% do.call(what=rbind) %>% as.data.frame()
# inform$V1 <- paste(inform$V1, '/', sep = '')
# inform$V1 <- sub(inform$V1, pattern = '.*?\\/', replacement='A/', )

colnames(inform)[1:5] <- c("Strain_name", "Date", "Strain_number", "Isolate_Id", 'segment')
seq <- cbind(seq[, -c(4)], inform)

seq <- merge(seq, meta, by = 'Isolate_Id', all.x = T) %>% as.data.frame() %>% 
  mutate(strain_host=word(Strain_name, sep = '/', 2) %>% tolower(), 
         new_header = '', 
         Location=str_remove_all(Location, pattern = ' '),
         ui= paste(Isolate_Id, Strain_number, sep = '|')) %>% rename(seq=V1) %>% select(-c(V6))

seq <- replace(seq, seq=='', 'na')

seq[seq$Location %in% '/Kurdistan', 'Location'] <- 'Asia/Kurdistan'

# seq$country <- word(seq$Location %>% str_remove_all(pattern = ' '), start = 2, sep = '\\/')

gis_country <- seq$Location %>% header_cleaning(pattern = '\\/') %>% select(V1, V2)
gis_country[gis_country$V2 %in% 'Kurdistan', 'V1'] <- 'Asia'

continent <- table(gis_country$V1, gis_country$V2) %>% as.data.frame(stringsAsFactors=F) %>% 
  filter(Freq!=0, Var2!='', Var1!='')
continent[(nrow(continent)+1), 1:2] <- c('Africa', 'Lesotho')




country_df <- data.frame(ird=c('USA', 'Russia', 'VietNam', 'HongKong', 'SouthKorea', 'Iran', 'Laos', 'Bolivia', 'GazaStrip', 'Czechoslovakia', 'CotedIvoire', 'DemocraticRepublicoftheCongo', 'NewZealand', 'SouthAfrica', 'UnitedArabEmirates'),
                         gis=c('UnitedStates', 'RussianFederation', 'Vietnam', 'HongKong(SAR)', 'Korea,Republicof', 'Iran,IslamicRepublicof', "Lao,People'sDemocraticRepublic", 'Bolivia,PlurinationialStateof', 'PalestinianTerritory', 'CzechRepublic', "Coted'Ivoire", 'Congo,theDemocaticRepublicof', 'NewZealand', 'SouthAfrica', 'UnitedArabEmirates'))

# processing IRD sequences: import
ird <- read.fasta(args$ird_seq_fasta, as.string = T, whole.header = T) %>% 
  do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('header') %>% 
  mutate(header=str_remove_all(header, pattern=' ')
         , index = word(header, 1, sep = '/')) %>% rename(seq=V1)

# header organization, remove one non-rules virus (mixed)
ird_header <- header_cleaning(header_vector = ird$header, pattern = '\\/A\\/')[, -3] %>% 
  mutate(V1=str_replace_all(V1, replacement = '|', pattern='/')) %>% 
  mutate(header=paste(V1, V2, sep = '|') %>% str_replace(pattern = '[0-9]\\|N[0-9]', replacement = ',mix'))  # Fix the mixed subtype sequence (label in header)

ird_header <- header_cleaning(header_vector = ird_header$header, pattern = '\\|') %>% filter(V9 == '') %>% select(-c(V9)) %>% 
  mutate(host= word(V8, 1, sep = '/'), county= word(V8, 2, sep = '/')) # remove the mixed subtype sequence (label in header)


colnames(ird_header) <- c('segment_code1', 'segment_code2', 'Avian', 'segment'
                          , 'Subtype', 'country', 'Date', 'virus_code', 'host', 'county')

# ird_header <- ird_header[!ird_header$Date %in% '--' ,] # remove the unclear date


# merge header with sequence
ird <- ird[ird$index %in% ird_header$segment_code1, -3]
ird <- cbind(ird, ird_header[, c('segment_code1', 'segment_code2' , 'virus_code', 'segment', 'Subtype'
                                 , 'country', 'county', 'Date', 'host')])
if(all(ird$segment_code1==word(ird$header, 1, sep = '/')) != TRUE) stop('Mismatching')


# extract country and identified continent by GISAID data
for(i in seq_len(nrow(country_df)))  {
  ird[ird$country %in% country_df[i, 1], 'country'] <- country_df[i, 2]
}


ird$continent <- ''
for(i in seq_len(nrow(continent)))  {
  ird[ird$country %in% continent[i, 2], 'continent'] <- continent[i, 1]
}


# segment organization
ird$segment <- str_extract(ird$segment, pattern = '([A-Z]{2})[0-9]|([A-Z]{2})')
ird$sub <- ird$segment


# identifed IRD unique virus
ird_only <- ird[ird$virus_code %in% ird[ird$segment_code2 %in% setdiff(ird$segment_code2, ncbi), 'virus_code'], ] %>% 
  group_split(virus_code) %>% as.list() %>% lapply(as.data.frame)

ird_name <- c()
x <- 0
for(i in seq_along(ird_only))  { #
  x <- x+1
  set.seed(x)
  test <- paste('IRD', paste(sample(c(0:9), 8, replace = T), collapse = ''), sep='_')
  while(nrow(ird_name %in% test %>% table())==2) {
    x <- x+1
    set.seed(x)
    test <- paste('IRD', paste(sample(c(0:9), 8, replace = T), collapse = ''), sep='_')
    # stop(paste('No.', i, ', seed is ', x, ', Repeat IRD names', sep = ''))
  }
  ird_name <- c(ird_name, test)
  names(ird_name)[i] <- x
}
ird_name <- data.frame(id=ird_name) %>% tibble::rownames_to_column('seed')

names(ird_only) <- ird_name$id


ird_seq <- ird_only %>% do.call(what=rbind) %>% as.data.frame() %>% tibble::rownames_to_column('Isolate_Id') %>% 
  mutate(Isolate_Id=str_remove(Isolate_Id, pattern='\\..*'), 
         Location=paste(continent, country, county, sep = '/'), 
         Clade='IRD_unlabel', Submission_Date='IRD_unlabel', 
         new_header='', 
         ui=paste(Isolate_Id, segment_code1, sep = '|')) 


# A <- organ_subtype(in_seq = ird_seq) %>% do.call(what=rbind) %>% as.data.frame()
ird_seq <- ird_seq[, c('Isolate_Id', "header", 'seq', 'Isolate_Id', 'segment_code1', 'Date', 
                       'segment_code1', 'segment', 'Subtype', 'Clade', 'Location', 'host', 
                       'Date', 'Date', 'host', 'new_header', 'ui')]

colnames(ird_seq) <- colnames(seq)

# str_extract(result$Location, pattern='^/.*/') %>% table




### processing the merge sequences
result <- rbind(seq, ird_seq)
result <- result %>% mutate(nt_number=nchar(seq)) %>% mutate(non_nt_number=str_count(pattern = '[^atcg]', seq)) %>% 
  mutate(non_nt_percent= non_nt_number/nt_number*100) 


not_process <- colnames(result) %in% c('seq') %>% which()
result[, -c(not_process)] <- apply(result[, -c(not_process)], MARGIN = 2
                                   , function(x) {
                                     str_remove_all(x, pattern = "GPS.*") %>%
                                       str_remove_all(pattern = "; GPS.*") %>%
                                       str_remove_all(pattern = "\\(.*E\\)") %>%
                                       str_remove_all(pattern = "<.*") %>% 
                                       str_remove_all(pattern = " ") %>% 
                                       str_replace_all(pattern = "'", replacement = '_') %>% 
                                       str_replace_all(pattern = "\\(|\\)|\\;|\\+|\\:|\\,", replacement = '_') # ? # 
                                     # str_replace_all(pattern = "'", replacement = '_') %>% 
                                     # str_replace_all(pattern = "GPS.*\\|", replacement = '|') %>% 
                                     # str_replace_all(pattern = "\\;|\\+|\\:|\\,", replacement = '_') %>%
                                     # str_replace_all(pattern = "\\(|\\)", replacement = '_')
                                   })
result[, 'Subtype'] <- str_replace_all(result[, 'Subtype'], pattern = 'H01', replacement = 'H1')

write.csv(organ_subtype(result) %>% do.call(what=rbind) %>% as.data.frame(), 
          file = paste(args$out_fasta_dir, 'gisaid_ird_merge_meta_information_raw', '.csv', sep = ''), 
          quote = F, row.names = F)

log <- table(result$segment) %>% as.data.frame(stringsAsFactors=F)
colnames(log) <- c('Segment', 'raw')
log[9, ] <- c('all', sum(log[, 2]))

log$Segment <- factor(log$Segment, levels = c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS' ,'all'))

####### QC

## step1: ambiguous subtype
result$Subtype <- str_remove(string = result$Subtype, pattern = 'A / ')
result <- result[grepl(result$Subtype, pattern= 'H[0-9]+N[0-9]+'), ]
result <- result[!grepl(result$Subtype, pattern= '^mix.*|^M.*'), ]
result[, 'Subtype'] <- str_replace_all(result[, 'Subtype'], pattern = 'H01', replacement = 'H1')


log <- merge(log, table_DF(result$segment), by.x = 'Segment', by.y = 'x', all = T)
colnames(log)[ncol(log)] <- c('explicit subtype')
log[9, ncol(log)] <- c(nrow(result))


log_sub <- organ_subtype(in_seq = result) %>% 
  lapply(function(x) {nrow(x)}) %>% do.call(what=rbind) %>% as.data.frame() %>%
  tibble::rownames_to_column('Segment')

colnames(log_sub) <- c('Segment', 'explicit subtype')


## step2: complete geography and temporal information 
# collection data contain at least year
result <- result[grepl(result$Collection_Date, pattern= '[0-9]{4}.*'), ]
# Location without country
result$country <- word(result$Location %>% str_remove_all(pattern = ' '), start = 2, sep = '\\/')
result <- result[!is.na(result$country), ]
# A <- result[is.na(result$country), ]

log <- merge(log, table_DF(result$segment), by.x = 'Segment', by.y = 'x', all = T)
colnames(log)[ncol(log)] <- c('explicit geography and temporal')
log[log$Segment %in% 'all', ncol(log)] <- nrow(result)

log_sub <- merge(log_sub, organ_subtype(in_seq = result) %>% 
                   lapply(function(x) {nrow(x)}) %>% do.call(what=rbind) %>% as.data.frame() %>%
                   tibble::rownames_to_column('Segment')
                 , by = 'Segment', all = T)

colnames(log_sub)[ncol(log_sub)] <- c('explicit geography and temporal')




## step3: length filtering, current use the 80% length median
cat(paste('Box threshold is ', args$length_down_limit, '\n', sep = ''))

A <- organ_subtype(in_seq = result)
out_list <- list()
out_list_remove <- list()
for(i in names(A))  {
  B <- A[[i]]
  # up_lim <- boxplot.stats(B$nt_number, coef = args$box_cut %>% as.numeric())$stats %>% max()
  # down_lim <- boxplot.stats(B$nt_number, coef = args$box_cut %>% as.numeric())$stats %>% min()
  out_list[[i]] <- B[B$nt_number > median(B$nt_number)*0.80, ]
  out_list_remove[[i]] <- B[B$nt_number <= median(B$nt_number)*0.80, ] # | B$nt_number >= median(B$nt_number)*1.05
}

result <- do.call(rbind, out_list)


log <- merge(log, table_DF(result$segment), by.x = 'Segment', by.y = 'x', all = T)
colnames(log)[ncol(log)] <- c('80% median')
log[log$Segment %in% 'all', ncol(log)] <- nrow(result)

log_sub <- merge(log_sub, organ_subtype(in_seq = result) %>% 
                   lapply(function(x) {nrow(x)}) %>% do.call(what=rbind) %>% as.data.frame() %>%
                   tibble::rownames_to_column('Segment')
                 , by = 'Segment', all = T)

colnames(log_sub)[ncol(log_sub)] <- c('80% median')


### step4: ambiguous characters
result <- result[result$non_nt_percent < 5, ]


log <- merge(log, table_DF(result$segment), by.x = 'Segment', by.y = 'x', all = T)
colnames(log)[ncol(log)] <- c('95% explicit letters')
log[log$Segment %in% 'all', ncol(log)] <- nrow(result)

log_sub <- merge(log_sub, organ_subtype(in_seq = result) %>% 
                   lapply(function(x) {nrow(x)}) %>% do.call(what=rbind) %>% as.data.frame() %>%
                   tibble::rownames_to_column('Segment')
                 , by = 'Segment', all = T)

colnames(log_sub)[ncol(log_sub)] <- c('95% explicit letters')




not_process <- colnames(result) %in% c('seq') %>% which()
result[, -c(not_process)] <- apply(result[, -c(not_process)], MARGIN = 2
                   , function(x) {
                     str_remove_all(x, pattern = "GPS.*") %>%
                       str_remove_all(pattern = "; GPS.*") %>%
                       str_remove_all(pattern = "\\(.*E\\)") %>%
                       str_remove_all(pattern = "<.*") %>% 
                       str_remove_all(pattern = " ") %>% 
                       str_replace_all(pattern = "'", replacement = '_') %>% 
                       str_replace_all(pattern = "\\(|\\)|\\;|\\+|\\:|\\,", replacement = '_') # ? # 
                     # str_replace_all(pattern = "'", replacement = '_') %>% 
                     # str_replace_all(pattern = "GPS.*\\|", replacement = '|') %>% 
                     # str_replace_all(pattern = "\\;|\\+|\\:|\\,", replacement = '_') %>%
                     # str_replace_all(pattern = "\\(|\\)", replacement = '_')
                   })
# seq[, c('Host', 'strain_host')]
result[, c('Host', 'strain_host')] <- apply(result[, c('Host', 'strain_host')]
                                         , MARGIN = 2
                                         , function(x) {tolower(x) %>% 
                                             str_remove_all(pattern = '_$|^_|\\?') %>% 
                                             str_replace_all(pattern = '-', replacement = '_')}) %>% as.data.frame()
result[, c('Host', 'strain_host')] <- replace(result[, c('Host', 'strain_host')], is.na(result[, c('Host', 'strain_host')])
                                           , values = 'wild_bird')

B <- table_DF(c(result$Host, result$strain_host)) %>% mutate(raw_type='Wild') 
B[grep(B$x, pattern = '(domes|^chicken|^duck|^goose|^turkey$)'), 'raw_type'] <- 'Domestic'

result$host_type <- 'Wild'
result[result$Host %in% B[B$raw_type %in% 'Domestic', 'x'], 'host_type'] <- 'Domestic'
result[result$strain_host %in% B[B$raw_type %in% 'Domestic', 'x'], 'host_type'] <- 'Domestic'



result$new_header <- paste(result$Isolate_Id, result$Strain_number, 
                           result$Subtype %>% str_remove(pattern = 'A/'), 
                           result$Clade, 
                           result$segment, result$Collection_Date,
                        paste(word(result$Location, 1, sep = '/'), result$country, sep = '/'), 
                        result$Host, result$strain_host, result$host_type, sep = '|') 


A <- organ_subtype(result)
args$out_fasta_dir <- "~/Analysis/aiv/merge/0307/"
dir.create(args$out_fasta_dir)
for(i in names(A))  {
  B <- A[[i]]
  out <- matrix(0,0,0) %>% as.data.frame()
  out[seq(1, nrow(B)*2, 2), 1] <- paste('>', B$new_header, sep = '')
  out[seq(2, nrow(B)*2, 2), 1] <- B$seq
  
  write.table(out, file = paste(args$out_fasta_dir, args$out_fasta_prefix, i, '.fa', sep = '')
              , quote = F, col.names = F, row.names = F)
}



write.csv(log, paste(args$out_fasta_dir, 'clean_log.csv', sep = ''), quote = F, row.names = F)
write.csv(log_sub, paste(args$out_fasta_dir, 'clean_log_sub.csv', sep = ''), quote = F, row.names = F)



write.csv(result, file = paste(args$out_fasta_dir, 'gisaid_ird_merge_meta_information', '.csv', sep = '')
          , quote = F, row.names = F)
