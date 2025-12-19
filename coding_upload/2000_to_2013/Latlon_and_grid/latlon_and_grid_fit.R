library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)


# Sun May  4 22:11:05 2025 ------------------------------

# fit 2004-2013 genotype data into 100*100 km2 grid
Location_more_than_3_fixed <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/data/Location_more_than_3_fixed.csv')

# class(Location_more_than_3_fixed$C_year)
Location_more_than_3_fixed_2000_to_2013 <- Location_more_than_3_fixed[(Location_more_than_3_fixed$C_year >= 2000)&(Location_more_than_3_fixed$C_year <= 2013), ]

# input shp
library(sf)

Global_Fish_100km <- st_read('/media/dyclab/新增磁碟區/Jeffery/data/shp/Global_Fish_100km/Global_Fish_100km.shp')

shpfile <- Global_Fish_100km 
csv_data <- Location_more_than_3_fixed_2000_to_2013
lon_col = "Longitude"
lat_col = "Latitude"
geom_col = 'Id'

library('rlang') # for sym()

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


test1 <- count_point_in_polygon(shpfile = Global_Fish_100km, 
                                csv_data = Location_more_than_3_fixed_2000_to_2013, 
                                lon_col = "Longitude", 
                                lat_col = "Latitude", 
                                geom_col = 'Id')
sum(test1$point_count)

# count_point_in_polygon_result <- test1[, -2]

complete_grid_data <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/data/complete_grid_data.csv')
ID_continent <- complete_grid_data[, c(1,2,3)]

H5HPAI_2000_to_2013_in_grid <- merge(x = ID_continent, y = test1, by = 'Id')
H5HPAI_2000_to_2013_in_grid <- H5HPAI_2000_to_2013_in_grid %>% as.data.frame()
H5HPAI_2000_to_2013_in_grid <- H5HPAI_2000_to_2013_in_grid[, -4]

write.csv(H5HPAI_2000_to_2013_in_grid, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/Latlon_and_grid/data/H5HPAI_2000_to_2013_in_grid.csv', row.names = F)





#---------------------------------------------------------------------
# fit 2014-2023 genotype data into 100*100 km2 grid, count type num
Location_more_than_3_fixed <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/data/Location_more_than_3_fixed.csv')

Location_more_than_3_fixed_2000_to_2013 <- Location_more_than_3_fixed[(Location_more_than_3_fixed$C_year >= 2000)&(Location_more_than_3_fixed$C_year <= 2013), ]


# input shp
library(sf)

Global_Fish_100km <- st_read('/media/dyclab/新增磁碟區/Jeffery/data/shp/Global_Fish_100km/Global_Fish_100km.shp')

# shpfile <- Global_Fish_100km 
# csv_data <- Location_more_than_3_fixed_after_2014
# lon_col = "Longitude"
# lat_col = "Latitude"
# geom_col = 'Id'
# feature_type = 'Genotype'
# 
# i=7874

library('rlang') # for sym()
library(dplyr)

# count point in polygon function
count_point_in_polygon_with_feature_type <- function(shpfile, csv_data, lon_col, lat_col, geom_col, feature_type){
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
  
  shpfile1 <- shpfile %>% as.data.frame()
  shpfile2 <- shpfile1[, geom_col] %>% as.data.frame()
  colnames(shpfile2) <- geom_col
  shpfile2$type_count <- 0
  for (i in 1:nrow(shpfile2)) {
    target_OwO <- shpfile2[i, geom_col]
    joined_data_in_shp1 <- joined_data_in_shp %>% as.data.frame()
    judge_string <- joined_data_in_shp1[, geom_col]
    sub_data <- joined_data_in_shp1[judge_string %in% target_OwO, ]
    sub_data_table <- table(sub_data[, feature_type]) %>% as.data.frame()
    shpfile2[i, 'type_count'] <- nrow(sub_data_table)
  }
  return(shpfile2)
}


ID_typenum <- count_point_in_polygon_with_feature_type(shpfile = Global_Fish_100km, 
                                         csv_data = Location_more_than_3_fixed_2000_to_2013, 
                                         lon_col = "Longitude", 
                                         lat_col = "Latitude", 
                                         geom_col = 'Id',
                                         feature_type = 'Genotype')


complete_grid_data <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/Latlon_and_grid/data/complete_grid_data.csv')
ID_continent <- complete_grid_data[, c(1,2,3)]

H5HPAI_genotypenum_2000_to_2013_in_grid <- merge(x = ID_continent, y = ID_typenum, by = 'Id')


write.csv(H5HPAI_genotypenum_2000_to_2013_in_grid, file = '/media/dyclab/新增磁碟區/Jeffery/HH_analysis/2000_to_2013/Latlon_and_grid/data/H5HPAI_genotypenum_2000_to_2013_in_grid.csv', row.names = F)
























