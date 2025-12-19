library(readr)
library(readxl)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(viridis)
library(aplot)
library(cowplot)
library(ggplotify)
library(vegan)


# Wed May 21 11:23:32 2025 ------------------------------


# input shp
library(sf)

Global_Fish_100km <- st_read('/media/dyclab/新增磁碟區/Jeffery/data/shp/Global_Fish_100km/Global_Fish_100km.shp')



# shpfile <- Global_Fish_100km
# csv_data <- all_clean
# lon_col <- "Longitude"
# lat_col <- "Latitude"
# geom_col <- 'Id'

library('rlang') # for sym()
library(dplyr)
# count point in polygon function
count_point_in_polygon <- function(shpfile, csv_data, lon_col, lat_col, geom_col){
  # transform csv to geometry (in new column)
  points_sf <- st_as_sf(csv_data, coords = c(lon_col, lat_col), crs = 4326)
  
  #change project system in case of wrong position
  shpfile <- st_transform(shpfile, crs = 4326) 
  points_sf <- st_transform(points_sf, crs = 4326) 
  
  sf_use_s2(FALSE)  # 禁用 s2 引擎, 從球面計算轉換為平面計算
  
  #join
  joined_data <- st_join(points_sf, shpfile)
  # is.na(joined_data$Id) %>% table()
  
  sf_use_s2(TRUE)
  
  joined_data_in_shp <- joined_data[!(joined_data$Id %in% NA), ]
  
  # 自訂函數，使用 group_split() 並以變數方式傳入欄位名稱
  split_by_column <- function(data, group_column_var) {
    group_column_sym <- sym(group_column_var)   # 將欄位名稱轉換為符號
    
    data %>%
      group_by(!!group_column_sym) %>%   # 使用 !! 解包符號
      group_split()
  }
  
  A <- split_by_column(data = joined_data_in_shp, group_column = geom_col)
  shpfile$point_count <- 0
  
  for (i in seq_along(A)) {
    AAA <- A[[i]]
    AAA_Id <- AAA[1, geom_col] %>% as.data.frame()
    AAAA_Id <- AAA_Id[1,1]
    shpfile[shpfile$Id %in% AAAA_Id, 'point_count'] <- nrow(AAA)
  }
  return(shpfile)
}


# FAO data
FAO_AIV <- read_csv('/media/dyclab/新增磁碟區/Jeffery/data/FAO/epidemiology-raw-data_202410081208.csv')

#some data unclean (move on cell right)
FAO_AIV_QAQ <- FAO_AIV[is.na(FAO_AIV$Latitude), ]
FAO_AIV_OwO <- FAO_AIV[!(is.na(FAO_AIV$Latitude)), ]

FAO_AIV_QAQ1 <- FAO_AIV_QAQ
FAO_AIV_QAQ1[, 9:17] <- FAO_AIV_QAQ[, 10:18]
FAO_AIV_QAQ1$Human.deaths <- NA

FAO_AIV <- rbind(FAO_AIV_OwO, FAO_AIV_QAQ1)

# Pick up NICD data and fix it
non_NICD <- FAO_AIV[!(FAO_AIV$Diagnosis.status == 'NICD'), ]
unclean <- FAO_AIV[FAO_AIV$Diagnosis.status == 'NICD', ]

clean_of_NICD <- unclean
clean_of_NICD[, 12:15] <- unclean[, 13:16]

NICD_correct <- rbind(non_NICD, clean_of_NICD)


# is.na(NICD_correct$Latitude) %>% table()
# is.na(NICD_correct$Longitude) %>% table()
# class(NICD_correct$Latitude)
# class(NICD_correct$Longitude)


NICD_correct$Longitude <- NICD_correct$Longitude %>% as.numeric()

NICD_correct_position <- NICD_correct[!(is.na(NICD_correct$Longitude)), ]

#create ob m.y. colume
OwO <- str_split(NICD_correct_position$Observation.date..dd.mm.yyyy., pattern = '/')

NICD_correct_position$OB_year <- NA
NICD_correct_position$OB_month <- NA


unclean$OB_year <- NA
unclean$OB_month <- NA

# debugQAQ <- NICD_correct[i,]

# is.na(NICD_correct_position$Latitude) %>% table()

for (i in seq_along(OwO)) {
  error_test_OwO <- try(OwO[[i]][[3]])
  if (inherits(error_test_OwO, 'try-error')) {
    unclean <- rbind(unclean, NICD_correct_position[i,])
  } else{
    NICD_correct_position[i, 'OB_year'] <- OwO[[i]][[3]]
    NICD_correct_position[i, 'OB_month'] <- OwO[[i]][[2]]
  }
}

# is.na(NICD_correct_position$Event.ID) %>% table()




str(NICD_correct_position)
class(NICD_correct_position)
NICD_correct_position <- NICD_correct_position %>% as.data.frame()


NICD_correct_position$OB_year <- NICD_correct_position$OB_year %>% as.integer()
NICD_correct_position$OB_month <- NICD_correct_position$OB_month %>% as.integer()

# table(NICD_correct_position$OB_year)
# class(NICD_correct_position)
# is.na(NICD_correct_position$Event.ID) %>% table()
is.na(NICD_correct_position$OB_year) %>% table()


# test <- NICD_correct_position[(is.na(NICD_correct_position$OB_year)), ]


NICD_correct_position_with_y <- NICD_correct_position[!(is.na(NICD_correct_position$OB_year)), ]




NICD_correct_position_to_2023 <- NICD_correct_position_with_y[NICD_correct_position_with_y$OB_year < 2024, ] %>% as.data.frame()

class(NICD_correct_position_to_2023)
# nrow(NICD_correct_position) - nrow(NICD_correct_position_to_2023)
# class(NICD_correct_position)
# table(NICD_correct_position_to_2023$OB_year == 2023)
# class(NICD_correct_position_to_2023)
is.na(NICD_correct_position_to_2023$Event.ID) %>% table()


# fao_2023 <- NICD_correct_position_to_2023[NICD_correct_position_to_2023$OB_year %in% 2023, ] %>% as.data.frame()
# fao_before_2023 <- NICD_correct_position_to_2023[NICD_correct_position_to_2023$OB_year < 2023, ] %>% as.data.frame()
# 
# # nrow(fao_2023) + nrow(fao_before_2023)
# # nrow(NICD_correct_position_to_2023)
# # 
# 
# 
# fao_2023_06 <- fao_2023[fao_2023$OB_month < 7, ] %>% as.data.frame()
# 
# fao_to_2023_06 <- rbind(fao_before_2023, fao_2023_06)

fao_2014_to_2021 <- NICD_correct_position_to_2023[(NICD_correct_position_to_2023$OB_year <= 2021)&(NICD_correct_position_to_2023$OB_year >= 2014), ] %>% as.data.frame()



# H5 HPAI / Domestic / Confirm cases / outbreak / 2014-2021

fao_2014_to_2021_H5_HPAI <- fao_2014_to_2021[grepl(x = fao_2014_to_2021$Serotype, pattern = "H5")&grepl(x = fao_2014_to_2021$Serotype, pattern = "HPAI"), ]

# unique(fao_2014_to_2021_H5_HPAI$Serotype)

fao_2014_to_2021_H5_HPAI_domestic <- fao_2014_to_2021_H5_HPAI[fao_2014_to_2021_H5_HPAI$Animal.type %in% "Domestic",]

# unique(fao_2014_to_2021_H5_HPAI_domestic$Animal.type)

fao_2014_to_2021_H5_HPAI_domestic_confirm <- fao_2014_to_2021_H5_HPAI_domestic[fao_2014_to_2021_H5_HPAI_domestic$Diagnosis.status %in% 'Confirmed', ]

# unique(fao_2014_to_2021_H5_HPAI_domestic_confirm$Diagnosis.status)

table(fao_2014_to_2021_H5_HPAI_domestic_confirm$Species)

remove_Species <- c("Canine", "Cats", "Dogs", "Ferret", "Swine", 'Unspecified Env. Sample')


fao_2014_to_2021_H5_HPAI_domestic_confirm_farm <- fao_2014_to_2021_H5_HPAI_domestic_confirm[!(fao_2014_to_2021_H5_HPAI_domestic_confirm$Species %in% remove_Species), ]




FAO_H5_HPAI_outbreak_poultry_farm <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                   csv_data = fao_2014_to_2021_H5_HPAI_domestic_confirm_farm, 
                                   lon_col = "Longitude", 
                                   lat_col = "Latitude", 
                                   geom_col = 'Id')



FAO_H5_HPAI_outbreak_poultry_farm1 <- FAO_H5_HPAI_outbreak_poultry_farm %>% as.data.frame()

colnames(FAO_H5_HPAI_outbreak_poultry_farm1)
FAO_H5_HPAI_outbreak_poultry_farm2 <- FAO_H5_HPAI_outbreak_poultry_farm1[, c('Id', 'point_count')]


FAO_H5_HPAI_outbreak_poultry_farm2$point_count <- ifelse(test = FAO_H5_HPAI_outbreak_poultry_farm2$point_count > 0,
                                                        yes = 1,
                                                        no = 0)


colnames(FAO_H5_HPAI_outbreak_poultry_farm2) <- c('Id', 'H5_HPAI_farm')


write.csv(FAO_H5_HPAI_outbreak_poultry_farm2, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/variable_arrange_for_second_stage/整理好的/H5_HPAI_farm.csv', row.names = F)







