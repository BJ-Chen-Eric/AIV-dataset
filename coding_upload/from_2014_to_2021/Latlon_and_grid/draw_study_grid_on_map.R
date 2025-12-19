library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(sf)


# Sun Apr 20 19:26:17 2025 ------------------------------

# combine four plot together

Global_Fish_100km <- st_read('/media/dyclab/新增磁碟區/Jeffery/data/shp/Global_Fish_100km/Global_Fish_100km.shp')


complete_grid_data <- read.csv('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/FAO_and_GISAID_record_grid_data/complete_grid_data_20250519.csv')


global_map <- st_read('/media/dyclab/新增磁碟區/Jeffery/HH_analysis/genotype_data_visualization1/QGIS/圓餅圖shp/1996_to_2013.shp')


#EU
EU <- complete_grid_data[complete_grid_data$Europe == 1,]

EU_map <- merge(x = EU, y = Global_Fish_100km, by = 'Id', all.x = T)

str(EU_map)

EU_map$Group_2 <- EU_map$Group_2 %>% as.character()
EU_map$Group_3 <- EU_map$Group_3 %>% as.character()




EU_map1 <- st_as_sf(EU_map)




EU_outer_border <- st_union(EU_map1)
# ggplot(data = EU_outer_border) +
#   geom_sf()

EU_outer_border1 <- EU_outer_border %>% st_make_valid()

EU_outer_border_clip <- st_intersection(global_map, EU_outer_border1)





EU_map1$Group_2 <- factor(EU_map1$Group_2, levels = c(1, 0))



ggplot() +
  geom_sf(data = EU_map1, aes(fill = Group_2), color = 'white') +
  geom_sf(data = EU_outer_border_clip, fill = NA, color = 'black', linewidth = 0.1) +
  scale_fill_manual(name = "Category", 
                    values = c('1' = '#138b00', 
                               '0' = 'white'), 
                    labels = c(
                      'others',
                      'Group 2 grids'
                    )  # 自定義圖例名稱
  )+
  theme_bw()+
  labs(title = "Europe group 2 study grid", x = "Lon", y = "Lat")+
  theme(axis.text = element_text(size = 50),
        axis.title = element_text(size = 60),
        plot.title = element_text(size = 70),
        legend.title = element_text(size = 70),
        legend.text = element_text(size = 45),
        legend.key.width  = unit(50, "pt")) +
  guides(fill = guide_legend(
    reverse = T,
    keywidth = unit(2, "cm"),  # 調整圖例鍵的寬度
    keyheight = unit(3, "cm")  # 調整圖例鍵的高度
  ))






the_plot_EU2 <- ggplot() +
  geom_sf(data = EU_map1, aes(fill = Group_2), color = 'white') +
  geom_sf(data = EU_outer_border_clip, fill = NA, color = 'black', linewidth = 0.2) +
  scale_fill_manual(name = "Category", 
                    values = c('1' = '#138b00', 
                               '0' = 'white'), 
                    labels = c(
                      'Group 1 grids',
                      'others'
                    )  # 自定義圖例名稱
  )+
  theme_bw()+
  labs(title = "", x = "Lon", y = "Lat")
  # theme(axis.text = element_text(size = 50),
  #       axis.title = element_text(size = 60),
  #       plot.title = element_text(size = 70),
  #       legend.title = element_text(size = 70),
  #       legend.text = element_text(size = 45),
  #       legend.key.width  = unit(50, "pt")) +
  # guides(fill = guide_legend(
  #   reverse = T,
  #   keywidth = unit(2, "cm"),  # 調整圖例鍵的寬度
  #   keyheight = unit(3, "cm")  # 調整圖例鍵的高度
  # ))

EU_map1$Group_3 <- factor(EU_map1$Group_3, levels = c(1, 0))

the_plot_EU3 <- ggplot() +
  geom_sf(data = EU_map1, aes(fill = Group_3), color = 'white') +
  geom_sf(data = EU_outer_border_clip, fill = NA, color = 'black', linewidth = 0.2) +
  scale_fill_manual(name = "Category", 
                    values = c('1' = '#138b00', 
                               '0' = 'white'), 
                    labels = c(
                      'Group 2 grids',
                      'others'
                    )  # 自定義圖例名稱
  )+
  theme_bw()+
  labs(title = "", x = "Lon", y = "Lat")
  # theme(axis.text = element_text(size = 50),
  #       axis.title = element_text(size = 60),
  #       plot.title = element_text(size = 70),
  #       legend.title = element_text(size = 70),
  #       legend.text = element_text(size = 45),
  #       legend.key.width  = unit(50, "pt")) +
  # guides(fill = guide_legend(
  #   reverse = T,
  #   keywidth = unit(2, "cm"),  # 調整圖例鍵的寬度
  #   keyheight = unit(3, "cm")  # 調整圖例鍵的高度
  # ))






#AS

AS <- complete_grid_data[complete_grid_data$Asia == 1,]

AS_map <- merge(x = AS, y = Global_Fish_100km, by = 'Id', all.x = T)

str(AS_map)

AS_map$Group_2 <- AS_map$Group_2 %>% as.character()
AS_map$Group_3 <- AS_map$Group_3 %>% as.character()



library(spData)
library(ggplot2)
library(sf)


AS_map1 <- st_as_sf(AS_map)



AS_outer_border <- st_union(AS_map1)
# ggplot(data = AS_outer_border1) +
#   geom_sf()

AS_outer_border1 <- AS_outer_border %>% st_make_valid()

AS_outer_border_clip <- st_intersection(global_map, AS_outer_border1)


lon_center = 100
lon_edge = lon_center + 180
matrix_center = matrix(c(lon_edge, lon_edge, 90, -90), ncol = 2)
line_center = sf::st_linestring(matrix_center)
line_center_sfc = st_sfc(line_center, crs = st_crs(AS_outer_border_clip))
polygon_cent = st_buffer(line_center_sfc, dist = 1)
# plot(AS_type_num_map1$geom)
# plot(polygon_cent, add = TRUE)

AS_outer_border_clip1 = sf::st_difference(AS_outer_border_clip, polygon_cent)
# plot(AS_type_num_map2$geom)


crs = "+proj=eqc +lon_0=100 +lat_ts=0 +datum=WGS84"
# different projection could be used
# "+proj=moll +lon_0=+100"
# '+proj=merc +lon_0=100'
# '+proj=laea +lon_0=100 +lat_0=30'
# "+proj=eqc +lon_0=100 +lat_ts=0 +datum=WGS84"

AS_outer_border_clip2 = st_transform(AS_outer_border_clip1, crs = crs)










lon_center = 100
lon_edge = lon_center + 180
matrix_center = matrix(c(lon_edge, lon_edge, 90, -90), ncol = 2)
line_center = sf::st_linestring(matrix_center)
line_center_sfc = st_sfc(line_center, crs = st_crs(AS_map1))
polygon_cent = st_buffer(line_center_sfc, dist = 1)
plot(AS_map1$geom)
plot(polygon_cent, add = TRUE)

AS_map2 = sf::st_difference(AS_map1, polygon_cent)
plot(AS_map2$geom)


crs = "+proj=eqc +lon_0=100 +lat_ts=0 +datum=WGS84"
# different projection could be used
# "+proj=moll +lon_0=+100"
# '+proj=merc +lon_0=100'
# '+proj=laea +lon_0=100 +lat_0=30'
# "+proj=eqc +lon_0=100 +lat_ts=0 +datum=WGS84"

AS_map3 = st_transform(AS_map2, crs = crs)



AS_map3$Group_2 <- factor(AS_map3$Group_2, levels = c(1, 0))






the_plot_AS2 <- ggplot() +
  geom_sf(data = AS_map3, aes(fill = Group_2), color = 'white') +
  geom_sf(data = AS_outer_border_clip2, fill = NA, color = 'black', linewidth = 0.2) +
  scale_fill_manual(name = "Category", 
                    values = c('1' = '#138b00', 
                               '0' = 'white'), 
                    labels = c(
                      'Group 1 grids',
                      'others'  
                    )  # 自定義圖例名稱
  )+
  theme_bw()+
  labs(title = "", x = "Lon", y = "Lat")
  # theme(axis.text = element_text(size = 50),
  #       axis.title = element_text(size = 60),
  #       plot.title = element_text(size = 70),
  #       legend.title = element_text(size = 70),
  #       legend.text = element_text(size = 45),
  #       legend.key.width  = unit(50, "pt")) +
  # guides(fill = guide_legend(
  #   reverse = T,
  #   keywidth = unit(2, "cm"),  # 調整圖例鍵的寬度
  #   keyheight = unit(3, "cm")  # 調整圖例鍵的高度
  # ))

AS_map3$Group_3 <- factor(AS_map3$Group_3, levels = c(1, 0))

the_plot_AS3 <- ggplot() +
  geom_sf(data = AS_map3, aes(fill = Group_3), color = 'white') +
  geom_sf(data = AS_outer_border_clip2, fill = NA, color = 'black', linewidth = 0.2) +
  scale_fill_manual(name = "Category", 
                    values = c('1' = '#138b00', 
                               '0' = 'white'), 
                    labels = c(
                      'Group 2 grids',
                      'others'  
                    )  # 自定義圖例名稱
  )+
  theme_bw()+
  labs(title = "", x = "Lon", y = "Lat")





library(cowplot)
sp_plot <- plot_grid(the_plot_EU2, the_plot_EU3, the_plot_AS2, the_plot_AS3, nrow = 2, ncol = 2)

sp_plot1 <- ggdraw(sp_plot) +
  draw_label("a", x = 0.025, y = 0.98, size = 10, fontface = 'bold') +  # 標記 'a' 在 plot1
  draw_label("b", x = 0.5, y = 0.98, size = 10, fontface = 'bold')+ # 標記 'b' 在 plot2
  draw_label("c", x = 0.025, y = 0.5, size = 10, fontface = 'bold')+
  draw_label("d", x = 0.5, y = 0.5, size = 10, fontface = 'bold')


ggsave("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/study_grid_visulization/2014_to_2021.png", plot = sp_plot1, width = 15, height = 9, dpi = 300, bg = "white")

library(svglite)
ggsave("/media/dyclab/新增磁碟區/Jeffery/HH_analysis/from_2014_to_2021/Latlon_and_grid/study_grid_visulization/2014_to_2021.svg", plot = sp_plot1, width = 15, height = 9, units = "in")


