library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(elsa)
library(rgdal)
library(spdep)
library(sf)

# Mon May 19 16:20:42 2025 ------------------------------


# setwd("D://梁雋承/研究所/集大成HH_as_case_2014to2023/LISA/")


# ---------------------------------------------------------------------
# 開始分析
gb.map.nb <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/LISA/data/neighbor.csv")
# NOnb <- gb.map.nb[gb.map.nb$neigh1 == 0,]
# OwO <- NOnb$id
alltype <- read.csv("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/data/H5HPAI_genotypenum_2014_to_2021_in_grid.csv")


# fix id 20276
WTF <- gb.map.nb[20276,2]
WTF1 <- strsplit(WTF,split=":")
WTF2 <- unlist(WTF1)
gb.map.nb[20276,2] <- WTF2[[1]]
gb.map.nb[20276,3] <- WTF2[[2]] %>% as.integer()


gb.map.nb$neigh1 <- gb.map.nb$neigh1 %>% as.integer()

str(gb.map.nb)

################## EU ###############
EU <- alltype[alltype$Europe==1, ]
EUnb <- gb.map.nb[gb.map.nb$id %in% EU$Id, ]

# i=636
# EUnb[i, 1]
# EUnb[i, 2:9]
for (i in 1:nrow(EUnb)) {
  nb_list <- EUnb[i, 2:9] %>% as.list()
  nb_list <- nb_list %>% unlist()
  nb_list_new <- nb_list[nb_list %in% EUnb$id]
  #identical(nb_list_new, structure(integer(0), names = character(0))) 判斷 nb_list_new 是否為 named logical(0)
  if (identical(nb_list_new, structure(integer(0), names = character(0)))) {
    EUnb[i, 2:9] <- NA
  } else {
    EUnb[i, 2:(1+length(nb_list_new))] <- nb_list_new
    EUnb[i, (1+length(nb_list_new)+1):9] <- NA
  }
}

is.na(EUnb$neigh1) %>% table()


#把孤立網格(沒有臨格)先丟掉，不然 localmoran() 不能跑
# 已確認過無臨格皆無case (不過後續分析不影響，因為會先丟掉)
EUNOnb <- EUnb[is.na(EUnb$neigh1),]
EU_OwO <- EUNOnb$id


EUnb1 <- EUnb[!(EUnb$id %in% EU_OwO),]
EU1 <- EU[!(EU$Id %in% EU_OwO),]

#給予 fack_ID 使結果對應回原 ID
new_col_num <- ncol(EU1)+1
EU1[, new_col_num] <- 1:nrow(EU1)
colnames(EU1)[new_col_num] <- "fack_ID"

AwA <- list()

# 生成 fack_id 的臨格表
# i = 1
for (i in 1:nrow(EUnb1)) {
  A <- EUnb1[i,2:9]
  B <- as.list(A)
  C <- B[!(B %in% "NA")]
  D <- as.integer(C)
  for (j in seq_along(D)) {
    E <- D[[j]]
    D[[j]] <- EU1[EU1$Id %in% E, 'fack_ID']
  }
  AwA[[i]] <- D
}

class(AwA) <- "nb"
OAO <- nb2listw(AwA, zero.policy = T)
QAQ <- EU1$type_count


LISA <- localmoran(QAQ, OAO,
                   zero.policy = T,
                   alternative = "two.sided")

LISA1 <- as.data.frame(LISA)
LISA1[,6] <- 0
colnames(LISA1)[6] <- 'class'
for (i in 1:nrow(EU1)) {
  value1 <- EU1[i,'type_count']
  diff1 <- value1 - mean(EU1$type_count)
  Z <- LISA1[i,4]
  P <- LISA1[i,5]
  # p 顯著性/ diff 本身是高是低 / Z 與自身同質異質
  if (P <= 0.05 & diff1 > 0 & Z > 0) {
    LISA1[i,6] <- 'HH'
  }
  if (P <= 0.05 & diff1 < 0 & Z > 0) {
    LISA1[i,6] <- 'LL'
  }
  if (P <= 0.05 & diff1 > 0 & Z < 0) {
    LISA1[i,6] <- 'HL'
  }
  if (P <= 0.05 & diff1 < 0 & Z < 0) {
    LISA1[i,6] <- 'LH'
  }
  if (P > 0.05) {
    LISA1[i,6] <- 'nonsig'
  }
}

EU1[,6:11] <- 0
colnames(EU1)[6:11] <- colnames(LISA1)

for (i in 1:nrow(LISA1)) {
  EU1[EU1$fack_ID %in% i,6:11] <- LISA1[i,1:6]
}

LISA_EU <- merge(x = EU, y = EU1, 
                     by.x = "Id", by.y = "Id", all = T)

na <- is.na(LISA_EU$class)
# table(na)

LISA_EU[na, 'class'] <- 'independent_grid'

write.csv(LISA_EU, file = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/LISA/result/LISA_EU.csv", row.names = FALSE)

# ---------------------------------------------------------------------------------------------------
################## AS ###############
AS <- alltype[alltype$Asia==1, ]
ASnb <- gb.map.nb[gb.map.nb$id %in% AS$Id, ]

# i=636
# ASnb[i, 1]
# ASnb[i, 2:9]
for (i in 1:nrow(ASnb)) {
  nb_list <- ASnb[i, 2:9] %>% as.list()
  nb_list <- nb_list %>% unlist()
  nb_list_new <- nb_list[nb_list %in% ASnb$id]
  #identical(nb_list_new, structure(integer(0), names = character(0))) 判斷 nb_list_new 是否為 named logical(0)
  if (identical(nb_list_new, structure(integer(0), names = character(0)))) {
    ASnb[i, 2:9] <- NA
  } else {
    ASnb[i, 2:(1+length(nb_list_new))] <- nb_list_new
    ASnb[i, (1+length(nb_list_new)+1):9] <- NA
  }
}

is.na(ASnb$neigh1) %>% table()


#把孤立網格(沒有臨格)先丟掉，不然 localmoran() 不能跑
# 已確認過無臨格皆無case (不過後續分析不影響，因為會先丟掉)
ASNOnb <- ASnb[is.na(ASnb$neigh1),]
AS_OwO <- ASNOnb$id


ASnb1 <- ASnb[!(ASnb$id %in% AS_OwO),]
AS1 <- AS[!(AS$Id %in% AS_OwO),]

#給予 fack_ID 使結果對應回原 ID
new_col_num <- ncol(AS1)+1
AS1[, new_col_num] <- 1:nrow(AS1)
colnames(AS1)[new_col_num] <- "fack_ID"

AwA <- list()

# 生成 fack_id 的臨格表
# i = 1
for (i in 1:nrow(ASnb1)) {
  A <- ASnb1[i,2:9]
  B <- as.list(A)
  C <- B[!(B %in% "NA")]
  D <- as.integer(C)
  for (j in seq_along(D)) {
    E <- D[[j]]
    D[[j]] <- AS1[AS1$Id %in% E, 'fack_ID']
  }
  AwA[[i]] <- D
}

class(AwA) <- "nb"
OAO <- nb2listw(AwA, zero.policy = T)
QAQ <- AS1$type_count


LISA <- localmoran(QAQ, OAO,
                   zero.policy = T,
                   alternative = "two.sided")

LISA1 <- as.data.frame(LISA)
LISA1[,6] <- 0
colnames(LISA1)[6] <- 'class'
for (i in 1:nrow(AS1)) {
  value1 <- AS1[i,'type_count']
  diff1 <- value1 - mean(AS1$type_count)
  Z <- LISA1[i,4]
  P <- LISA1[i,5]
  # p 顯著性/ diff 本身是高是低 / Z 與自身同質異質
  if (P <= 0.05 & diff1 > 0 & Z > 0) {
    LISA1[i,6] <- 'HH'
  }
  if (P <= 0.05 & diff1 < 0 & Z > 0) {
    LISA1[i,6] <- 'LL'
  }
  if (P <= 0.05 & diff1 > 0 & Z < 0) {
    LISA1[i,6] <- 'HL'
  }
  if (P <= 0.05 & diff1 < 0 & Z < 0) {
    LISA1[i,6] <- 'LH'
  }
  if (P > 0.05) {
    LISA1[i,6] <- 'nonsig'
  }
}

AS1[,6:11] <- 0
colnames(AS1)[6:11] <- colnames(LISA1)

for (i in 1:nrow(LISA1)) {
  AS1[AS1$fack_ID %in% i,6:11] <- LISA1[i,1:6]
}

LISA_AS <- merge(x = AS, y = AS1, 
                 by.x = "Id", by.y = "Id", all = T)

na <- is.na(LISA_AS$class)
# table(na)

LISA_AS[na, 'class'] <- 'independent_grid'

write.csv(LISA_AS, file = "/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/LISA/result/LISA_AS.csv", row.names = FALSE)











