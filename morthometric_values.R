library(dggridR)
library(dplyr)
library(terra)
library(raster)
library(sf)
library(stars)
library(tibble)

#Функция вычисления морфометрических величин по гексагональным сеткам 
morphometric_values = function(path, res, cellsize, radius) {
  raster_obj = raster(path) #Чтение растра
  ext_coord_dem = extent(raster_obj) #Определение координат ограничивающего прямоугольника
  
  #Построение дискретной глобальной сеточной системы
  dggs = dgconstruct(projection = "ISEA", res = res)
  
  #Получение координат отдельных ячеек в заданном охвате
  grid = dgrectgrid(dggs,   
                     minlat = ext_coord_dem[3], minlon = ext_coord_dem[1],
                     maxlat = ext_coord_dem[4], maxlon = ext_coord_dem[2], cellsize = cellsize)
  
  centroid = st_centroid(grid) #Определение координат центроидов ячеек
  
  raster_obj2 = read_stars(path) #Чтение растра
  
  #Извлечение значений высот из ЦМР и запись их в ячейки гексагональной сетки
  hex_values = st_extract(raster_obj2,
                          centroid,
                          bilinear = TRUE) 
  colnames(hex_values) = c("vus", "geometry")
  
  #Построение буфура для поиска соседних ячеек
  zone = st_buffer(centroid$geometry, dist = radius)
  zone = as.data.frame(zone)
  
  #Нахождение пересечения буфера ячейки с соседними ячейками
  ints = list()
  for (i in zone$geometry) {
    ints = c(ints, st_intersects(i, centroid$geometry))
  }
  ints1 = cbind(ID = 1:length(ints), ints)
  ints1 = as.data.frame(ints1)
  hex_values1 = cbind(ID = 1:length(hex_values$vus), hex_values) #Фрейм с id ячеек и номерами ячеек
  
  #Присвоение высот ячейкам в соответствии с их идентификатором
  d = list()
  for (a in ints1$ints) {
    temp_list = c()
    for (b in a) {
      temp_list = c(temp_list, filter(hex_values1, ID == b)$vus)
    }
    d = c(d, list(temp_list))
  }
  new_table = mutate(ints1, high = d)
  new_table = cbind(new_table, grid$geometry) #Фрейм со значениями id, высот и координат ячеек
  
  #Вычисление крутизны и эскпозиции
  slope = c()
  aspect = c()
  for (high in new_table$high) {
    di = (high[2] - high[4] + high[4] - high[6])/2
    dj = (high[1] - high[4] + high[4] - high[7])/2
    dk = (high[3] - high[4] + high[4] - high[5])/2
    if (res %% 2 == 0) {
      dy = di + dj*sin(pi/6) - dk*sin(pi/6)
      dx = dj*cos(pi/6) + dk*cos(pi/6)
    }
    if (res %% 2 == 1) {
      dy = di*cos(pi/6) + dj*cos(pi/6)  
      dx = dk + dj*sin(pi/6) - di*sin(pi/6)  
    }
    slope = c(slope, (atan(sqrt(dx^2 + dy^2))*180/pi))
    aspect = c(aspect, (atan2(dy, -dx))*180/pi)
  }
  
  #Вычисление кривизны
  curvature = c()
  for (high in new_table$high) {
    di2 = 2*high[4] - high[6] - high[2]
    dj2 = 2*high[4] - high[3] - high[5]
    dk2 = 2*high[4] - high[1] - high[7]
    if (res %% 2 == 0) {
      dy2 = di2 + dj2*sin(pi/6) - dk2*sin(pi/6)
      dx2 = dj2*cos(pi/6) + dk2*cos(pi/6)
    }
    if (res %% 2 == 1) {
      dy2 = di2*cos(pi/6) + dj2*cos(pi/6)  
      dx2 = dk2 + dj2*sin(pi/6) - di2*sin(pi/6)  
    }
    curvature = c(curvature, (sqrt(dx2^2 + dy2^2))*180/pi/(cellsize*111000))
  }
  
  #Добавление полученных значений в полученный ранее фрейм
  data_end = add_column(new_table, slope = slope, aspect = aspect, curvature = curvature, .before = 'geometry')
  #Удаление строк с нулевыми значениями
  data_end = data_end[rowSums(is.na(data_end)) == 0, ]
  
  #Перевод значений экспозиции в шкалу, имеющую диапазон от 0° до 360°
  aspect2 = c()
  for (as in data_end$aspect) {
    if (as < 0) {
      aspect2 = c(aspect2, as + 360)
    }
    else {
      aspect2 = c(aspect2, as)
    }
  }
  
  #Создание итогового геодатафрейма
  df2 = data_end[,-5]
  df2 =  add_column(df2, aspect = aspect2, .before = 'geometry')
  new_df = subset(df2, select = -c(ints, high))
  new_df$ID = as.integer(new_df$ID)
  new_df$slope = as.numeric(new_df$slope)
  new_df$aspect = as.numeric(new_df$aspect)
  new_df$curvature = as.numeric(new_df$curvature)
  return(new_df)
}

#Построение по ЦМР GMTED2010
path = "D:\\semestr\\Kursovay\\RESULTS\\Yablonevy\\Yablonevy_GMTED_last.tif"
res = 21
cellsize = 0.000622
radius = 84
gdf = morphometric_values(path, res, cellsize, radius)
st_write(gdf, 'D:\\semestr\\Kursovay\\RESULTS\\Y_GMTED2010.shp')

#Построение по исходной ЦМР NASADEM
path = 'D:\\semestr\\Kursovay\\Obrezka\\nasadem.tif'
res = 24
cellsize = 0.000115
radius = 15.8
gdf = morphometric_values(path, res, cellsize, radius)
st_write(gdf, 'D:\\semestr\\Kursovay\\RESULTS\\NASADEM.shp')


