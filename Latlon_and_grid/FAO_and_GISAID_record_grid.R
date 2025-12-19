library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)


# Sun May  4 22:58:18 2025 ------------------------------

Location_more_than_3_fixed <- read_csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/FAO_and_GISAID_record_grid_data/whole_genotype_with_latlon.csv')

# ----------------------------------------------------------------------------------------

Location_more_than_3_fixed$C_year <- NA
OwO <- str_split(string = Location_more_than_3_fixed$Collection_Date, pattern = '-')
for (i in seq_along(OwO)) {
  Location_more_than_3_fixed[i, 'C_year'] <- OwO[[i]][[1]]
}




Location_more_than_3_fixed_before_2013 <- Location_more_than_3_fixed[Location_more_than_3_fixed$C_year <= 2013, ]


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


GISAID_with_uncomplete <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                csv_data = Location_more_than_3_fixed_before_2013, 
                                lon_col = "Longitude", 
                                lat_col = "Latitude", 
                                geom_col = 'Id')

# table(Location_more_than_3_fixed$Complete)
# class(Location_more_than_3_fixed$Complete)

Location_more_than_3_fixed_before_2013_completed <- Location_more_than_3_fixed_before_2013[Location_more_than_3_fixed_before_2013$Complete > 0, ]


GISAID_with_complete <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                                 csv_data = Location_more_than_3_fixed_before_2013_completed, 
                                                 lon_col = "Longitude", 
                                                 lat_col = "Latitude", 
                                                 geom_col = 'Id')

# ------------------------------------------------------------------------------------
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
# is.na(NICD_correct_position$OB_year) %>% table()


NICD_correct_position_with_ym <- NICD_correct_position[!(is.na(NICD_correct_position$OB_year)), ]




NICD_correct_position_to_2023 <- NICD_correct_position_with_ym[NICD_correct_position_with_ym$OB_year < 2024, ] %>% as.data.frame()

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

fao_before_2013 <- NICD_correct_position_to_2023[NICD_correct_position_to_2023$OB_year <= 2013, ] %>% as.data.frame()


FAO_grid <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                               csv_data = fao_before_2013, 
                                               lon_col = "Longitude", 
                                               lat_col = "Latitude", 
                                               geom_col = 'Id')

# ----------------------------------------------------------------



FAO_grid <- FAO_grid %>% as.data.frame()
GISAID_with_complete <- GISAID_with_complete %>% as.data.frame()
GISAID_with_uncomplete <- GISAID_with_uncomplete %>% as.data.frame()

sum(FAO_grid$point_count)
sum(GISAID_with_complete$point_count)
sum(GISAID_with_uncomplete$point_count)

# FAO_grid
# GISAID_with_complete
# GISAID_with_uncomplete : all of gisaid data with lat/lon (location more than 3) ()


# make new complete_grid_data
complete_grid_data <- read_csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/data/complete_grid_data.csv')
complete_grid_data1 <- complete_grid_data[, 1:3]



FAO_grid1 <- FAO_grid
FAO_grid1[FAO_grid1$point_count > 0, 'point_count'] <- 1
FAO_grid2 <- FAO_grid1[, c(1,3)]
colnames(FAO_grid2)[2] <- 'FAO_AIV_before_2013'


GISAID_with_complete1 <- GISAID_with_complete
GISAID_with_complete1[GISAID_with_complete1$point_count > 0, 'point_count'] <- 1
GISAID_with_complete2 <- GISAID_with_complete1[, c(1,3)]
colnames(GISAID_with_complete2)[2] <- 'GISAID_complete_before_2013'


GISAID_with_uncomplete1 <- GISAID_with_uncomplete
GISAID_with_uncomplete1[GISAID_with_uncomplete1$point_count > 0, 'point_count'] <- 1
GISAID_with_uncomplete2 <- GISAID_with_uncomplete1[, c(1,3)]
colnames(GISAID_with_uncomplete2)[2] <- 'GISAID_all_before_2013'



complete_grid_data2 <- merge(x = complete_grid_data1, y = FAO_grid2, by = 'Id')
complete_grid_data3 <- merge(x = complete_grid_data2, y = GISAID_with_complete2, by = 'Id')
complete_grid_data4 <- merge(x = complete_grid_data3, y = GISAID_with_uncomplete2, by = 'Id')


complete_grid_data4$Group_2 <- 0

complete_grid_data4[complete_grid_data4$FAO_AIV_before_2013 %in% 1 | complete_grid_data4$GISAID_all_before_2013 %in% 1, 'Group_2'] <- 1

complete_grid_data4$Group_3 <- complete_grid_data4$GISAID_complete_before_2013



# combine LISA result


LISA_AS <- read_csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/LISA/result/LISA_AS.csv')
LISA_EU <- read_csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/LISA/result/LISA_EU.csv')

LISA_AS_HH <- LISA_AS[LISA_AS$class %in% 'HH', ]
LISA_AS_HH_Id <- LISA_AS_HH$Id

complete_grid_data4$AS_HH <- 0
complete_grid_data4[complete_grid_data4$Id %in% LISA_AS_HH_Id, 'AS_HH'] <- 1


LISA_EU_HH <- LISA_EU[LISA_EU$class %in% 'HH', ]
LISA_EU_HH_ID <- LISA_EU_HH$Id

complete_grid_data4[, 'EU_HH'] <- 0
complete_grid_data4[complete_grid_data4$Id %in% LISA_EU_HH_ID, 'EU_HH'] <- 1



write.csv(complete_grid_data4, '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/Latlon_and_grid/FAO_and_GISAID_record_grid_data/complete_grid_data_20250528.csv', row.names = F)






# Mon May  5 09:45:25 2025 ------------------------------
# H5 HPAI outbreak farm 2004-2013 -----------------------------------------
# before this, please run the code from 101-205 (FAO data ~ fao_before_2013)

# fao_before_2013

fao_2000_to_2013 <- fao_before_2013[fao_before_2013$OB_year >= 2000, ]

fao_2000_to_2013$Serotype %>% table()

fao_2000_to_2013_H5HPAI <- fao_2000_to_2013[grepl(x = fao_2000_to_2013$Serotype, pattern = 'H5')&grepl(x = fao_2000_to_2013$Serotype, pattern = 'HPAI'), ]

fao_2000_to_2013_H5HPAI$Serotype %>% table()




fao_2000_to_2013_H5HPAI$Animal.type %>% table()

fao_2000_to_2013_H5HPAI_Domestic <- fao_2000_to_2013_H5HPAI[fao_2000_to_2013_H5HPAI$Animal.type %in% 'Domestic',]

fao_2000_to_2013_H5HPAI_Domestic$Species %>% table()

species_need_to_remove <- c('Cats', 'Dogs', 'Donkey', 'Goats', 'Swine')

fao_2000_to_2013_H5HPAI_Domestic_bird <- fao_2000_to_2013_H5HPAI_Domestic[!(fao_2000_to_2013_H5HPAI_Domestic$Species %in% species_need_to_remove), ]

fao_2000_to_2013_H5HPAI_Domestic_bird$Diagnosis.status %>% table()

fao_2000_to_2013_H5HPAI_Domestic_bird_confirmedCase <- fao_2000_to_2013_H5HPAI_Domestic_bird[fao_2000_to_2013_H5HPAI_Domestic_bird$Diagnosis.status %in% 'Confirmed', ]


FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                   csv_data = fao_2000_to_2013_H5HPAI_Domestic_bird_confirmedCase, 
                                   lon_col = "Longitude", 
                                   lat_col = "Latitude", 
                                   geom_col = 'Id')


FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm1 <- FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm %>% as.data.frame()

colnames(FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm1)[3] <- 'FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm'

FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm2 <- FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm1[, -2]
FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm2$FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm <- 
  ifelse(test = FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm2$FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm >=1, yes = 1, no = 0)



write.csv(FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm2, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/Latlon_and_grid/H5_HPAI_outbreak_farm/FAO_2000_to_2013_H5HPAI_outbreak_poultryfarm.csv', row.names = F)







