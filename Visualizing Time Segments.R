#Cluster Time Segments from Initial Partitioning Model

set.seed(1)

library(dplyr)
library(lubridate)
library(ggplot2)
library(coda)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("~/Documents/Snail Kite Project/Cat-Dirichlet Gibbs Sampler")
source('gibbs functions2.R')


#####################################################
#### Import Original Model Input and Breakpoints ####
#####################################################

setwd("~/Documents/Snail Kite Project/Data")
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
obs<- read.csv("Occupancy Matrix for all Obs and Locs.csv", header = T, sep = ",")
obs1.breakpts<- read.csv("ID1 Breakpoints (5 km).csv", header = T, sep = ",")
obs1.breakpts=obs1.breakpts[,1]
obs12.breakpts<- read.csv("ID12 Breakpoints (5 km).csv", header = T, sep = ",")
obs12.breakpts=obs12.breakpts[,1]
obs19.breakpts<- read.csv("ID19 Breakpoints (5 km).csv", header = T, sep = ",")
obs19.breakpts=obs19.breakpts[,1]
obs27.breakpts<- read.csv("ID27 Breakpoints (5 km).csv", header = T, sep = ",")
obs27.breakpts=obs27.breakpts[,1]


obs1=data.frame(obs) %>% filter(id == 1) %>% dplyr::select(-id)
obs12=data.frame(obs) %>% filter(id == 12) %>% dplyr::select(-id)
obs19=data.frame(obs) %>% filter(id == 19) %>% dplyr::select(-id)
obs27=data.frame(obs) %>% filter(id == 27) %>% dplyr::select(-id)

obs1$time1=1:nrow(obs1)
obs12$time1=1:nrow(obs12)
obs19$time1=1:nrow(obs19)
obs27$time1=1:nrow(obs27)




############################################
#### Summarize Results From First Model ####
############################################

nloc=ncol(obs)-1 #remove time1


obs1.seg=get.summary.stats(obs1.breakpts,obs1,nloc)
obs12.seg=get.summary.stats(obs12.breakpts,obs12,nloc)
obs19.seg=get.summary.stats(obs19.breakpts,obs19,nloc)
obs27.seg=get.summary.stats(obs27.breakpts,obs27,nloc)




####################################################
#### Map and Visualize Time Segments from Model ####
####################################################

assign.time.seg<- function(obs, breakpts, dat) {
 seg.n=apply(obs,1,sum)
 time.seg=rep(1:(length(breakpts)+1), seg.n)
 dat$time.seg<- time.seg
 
 dat
}

dat1=dat %>% filter(id==1)
dat12=dat %>% filter(id==12)
dat19=dat %>% filter(id==19)
dat27=dat %>% filter(id==27)

dat1<- assign.time.seg(obs1.seg, obs1.breakpts, dat1)
dat12<- assign.time.seg(obs12.seg, obs12.breakpts, dat12)
dat19<- assign.time.seg(obs19.seg, obs19.breakpts, dat19)
dat27<- assign.time.seg(obs27.seg, obs27.breakpts, dat27)


### Map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N


#Paths of individual segments
ggplot(data = fl) +
  geom_sf() +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)), ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat1, aes(x = utmlong, y = utmlat, color = time.seg)) +
  scale_color_viridis_c("Time Segment") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()


#Paths for segments 1-10
ggplot(data = fl) +
  geom_sf() +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)), ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat1[as.numeric(dat1$time.seg) <= 10,], aes(x = utmlong, y = utmlat,
                                                               color  = time.seg)) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()


#plot centroids for each cell of grid
# 5 km w 1 cell buffer
grid_5<- raster(extent(min(dat$utmlong), max(dat$utmlong),
                       min(dat$utmlat), max(dat$utmlat)) + 10000)
res(grid_5)<- 5000
proj4string(grid_5)<- CRS("+init=epsg:32617")
grid_5[]<- 0

#calc centroids of cells
grid.cell.locs<- coordinates(grid_5) %>% data.frame()
#write.csv(grid.cell.locs, "Centroids of 5km grid cells.csv", row.names = F)

#Create grid cell borders
borders_5<- rasterToPolygons(grid_5, dissolve = F)
borders_5f<- fortify(borders_5)

#Make grid as DF
grid_5f<- as.data.frame(grid_5, xy = TRUE)

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_tile(data=grid_5f, aes(x=x, y=y), fill="transparent") +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = grid.cell.locs, aes(x, y), color = "red", size=0.25, alpha=0.5) +
  labs(x = "Easting", y = "Northing") +
  theme_bw()



#plotting connections among grid cells for ID1
names(grid.cell.locs)<- c("grid.x", "grid.y")
grid.cell.locs$grid.cell<- 1:length(grid_5)

dat1<- left_join(dat1,grid.cell.locs, by="grid.cell")
dat12<- left_join(dat12,grid.cell.locs, by="grid.cell")
dat19<- left_join(dat19,grid.cell.locs, by="grid.cell")
dat27<- left_join(dat27,grid.cell.locs, by="grid.cell")


ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_tile(data=grid_5f, aes(x=x, y=y), fill="transparent") +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat1, aes(grid.x, grid.y, color=time.seg), size=1) +
  scale_color_viridis_c("Time Segment") +
  labs(x = "Easting", y = "Northing") +
  theme_bw()



### Determine mode within each time segment and connect those with heavy path by time seg
## Use lighter weighted path to connect locations within the same time segment

mode.fun<- function(x) {
  names(table(x))[table(x)==max(table(x))]
  }

attr.cntrs1<- matrix(NA, length(unique(dat1$time.seg)), 4)
attr.cntrs12<- matrix(NA, length(unique(dat12$time.seg)), 4)
attr.cntrs19<- matrix(NA, length(unique(dat19$time.seg)), 4)
attr.cntrs27<- matrix(NA, length(unique(dat27$time.seg)), 4)


id.attr.cntr=function(dat,attr.cntrs) {
  for (i in 1:length(unique(dat$time.seg))) {
    tmp=mode.fun(dat[dat$time.seg == i, "grid.cell"])
    tmp=ifelse(length(tmp) > 1, tmp[1], tmp) #arbitrary selection of 1st cell in set if multiple vals
  
    attr.cntrs[i,1]=as.numeric(tmp)
    attr.cntrs[i,2]=grid.cell.locs[grid.cell.locs$grid.cell==attr.cntrs[i,1], "grid.x"]
    attr.cntrs[i,3]=grid.cell.locs[grid.cell.locs$grid.cell==attr.cntrs[i,1], "grid.y"]
    attr.cntrs[i,4]=max(table(dat[dat$time.seg == i, "grid.cell"]))
  }
  colnames(attr.cntrs)<- c("grid.cell","grid.x","grid.y","size")
  attr.cntrs<- data.frame(attr.cntrs)
  attr.cntrs$time.seg<- 1:length(unique(dat$time.seg))
  
  attr.cntrs
}


attr.cntrs1=id.attr.cntr(dat1, attr.cntrs1)
attr.cntrs12=id.attr.cntr(dat12, attr.cntrs12)
attr.cntrs19=id.attr.cntr(dat19, attr.cntrs19)
attr.cntrs27=id.attr.cntr(dat27, attr.cntrs27)


ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = attr.cntrs1, aes(x=grid.x, y=grid.y)) +
  geom_point(data = attr.cntrs1, aes(x=grid.x, y=grid.y, color = time.seg), pch=21, stroke=1.25,
             size=attr.cntrs1$size/nrow(dat1)*100) +
  scale_color_viridis_c("Time Segment") +
  labs(x = "Easting", y = "Northing") +
  theme_bw()


dat.size1=matrix(NA,nrow(dat1),1)
dat.size12=matrix(NA,nrow(dat12),1)
dat.size19=matrix(NA,nrow(dat19),1)
dat.size27=matrix(NA,nrow(dat27),1)

n.per.tseg=function(dat,dat.size) {
  for (i in 1:length(unique(dat$time.seg))) {
    tmp=table(dat[dat$time.seg == i, "grid.cell"]) %>% data.frame()
    ind=which(dat$time.seg == i)
    for (j in ind){
      dat.size[j,]=tmp[which(tmp$Var1 == dat$grid.cell[j]),2]
    }
  }
  dat.size
}

dat1$size<- n.per.tseg(dat1, dat.size1) %>% as.vector()
dat12$size<- n.per.tseg(dat12, dat.size12) %>% as.vector()
dat19$size<- n.per.tseg(dat19, dat.size19) %>% as.vector()
dat27$size<- n.per.tseg(dat27, dat.size27) %>% as.vector()





ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat1, aes(x=grid.x, y=grid.y), color="gray60", size=0.25) +
  geom_point(data = dat1, aes(x=grid.x, y=grid.y, color = as.numeric(time.seg)), pch=21,stroke=1.25,
             size=dat1$size/nrow(dat1)*100, alpha=0.7) +
  geom_path(data = attr.cntrs1, aes(x=grid.x, y=grid.y)) +
  geom_point(data = attr.cntrs1, aes(x=grid.x, y=grid.y, color = time.seg), pch=21, stroke=1.25,
             size=attr.cntrs1$size/nrow(dat1)*100) +
  scale_color_viridis_c("Time Segment (n=49)") +
  labs(x = "Easting", y = "Northing", title = "ID 1") +
  theme_bw()





### 3D viz

library(plotly)

Sys.setenv("plotly_username"="joshcullen")
Sys.setenv("plotly_api_key"="qjYzTwjSylvlZnYmQzMk")

#ID 1
time.scatter1 = plot_ly() %>%
  add_sf(data=fl, x=~x, y=~y, z=0, type="scatter3d", color=I("black")) %>%
  add_trace(data=attr.cntrs1, x=~grid.x, y=~grid.y, z=0, type="scatter3d", mode="markers",
            color=I("gray40")) %>%
  add_trace(data=attr.cntrs1, x=~grid.x, y=~grid.y, z=~time.seg, type="scatter3d", mode="lines+markers", color=~size, line = list(color = 'black', width = 1)) %>%
  layout(title="Attraction Centers of ID 1",
         scene=list(xaxis = list(title ="Easting", range=c(min(dat$utmlong-50000),
                                                           max(dat$utmlong+50000))),
                    yaxis = list(title = "Northing", range=c(min(dat$utmlat-50000),
                                                             max(dat$utmlat+50000))),
                    zaxis = list(title = "Time Segment", range=c(0,nrow(attr.cntrs1)+1)),
                    camera = list(eye = list(x = -1, y = -1.5, z = 1.5),
                                  up = list(x=1, y=0, z=1), 
                                  center=list(x=0, y=0, z=0))
                    ))


#ID 12
time.scatter12 = plot_ly() %>%
  add_sf(data=fl, x=~x, y=~y, z=0, type="scatter3d", color=I("black")) %>%
  add_trace(data=attr.cntrs12, x=~grid.x, y=~grid.y, z=0, type="scatter3d", mode="markers",
            color=I("gray40")) %>%
  add_trace(data=attr.cntrs12, x=~grid.x, y=~grid.y, z=~time.seg, type="scatter3d", mode="lines+markers", color=~size, line = list(color = 'black', width = 1)) %>%
  layout(title="Attraction Centers of ID 12",
         scene=list(xaxis = list(title ="Easting", range=c(min(dat$utmlong-50000),
                                                           max(dat$utmlong+50000))),
                    yaxis = list(title = "Northing", range=c(min(dat$utmlat-50000),
                                                             max(dat$utmlat+50000))),
                    zaxis = list(title = "Time Segment", range=c(0,nrow(attr.cntrs12)+1)),
                    camera = list(eye = list(x = -1, y = -1.5, z = 1.5),
                                  up = list(x=1, y=0, z=1), 
                                  center=list(x=0, y=0, z=0))
         ))


#ID 19
time.scatter19 = plot_ly() %>%
  add_sf(data=fl, x=~x, y=~y, z=0, type="scatter3d", color=I("black")) %>%
  add_trace(data=attr.cntrs19, x=~grid.x, y=~grid.y, z=0, type="scatter3d", mode="markers",
            color=I("gray40")) %>%
  add_trace(data=attr.cntrs19, x=~grid.x, y=~grid.y, z=~time.seg, type="scatter3d", mode="lines+markers", color=~size, line = list(color = 'black', width = 1)) %>%
  layout(title="Attraction Centers of ID 19",
         scene=list(xaxis = list(title ="Easting", range=c(min(dat$utmlong-50000),
                                                           max(dat$utmlong+50000))),
                    yaxis = list(title = "Northing", range=c(min(dat$utmlat-50000),
                                                             max(dat$utmlat+50000))),
                    zaxis = list(title = "Time Segment", range=c(0,nrow(attr.cntrs19)+1)),
                    camera = list(eye = list(x = -1, y = -1.5, z = 1.5),
                                  up = list(x=1, y=0, z=1), 
                                  center=list(x=0, y=0, z=0))
         ))


#ID 27
time.scatter27 = plot_ly() %>%
  add_sf(data=fl, x=~x, y=~y, z=0, type="scatter3d", color=I("black")) %>%
  add_trace(data=attr.cntrs27, x=~grid.x, y=~grid.y, z=0, type="scatter3d", mode="markers",
            color=I("gray40")) %>%
  add_trace(data=attr.cntrs27, x=~grid.x, y=~grid.y, z=~time.seg, type="scatter3d", mode="lines+markers", color=~size, line = list(color = 'black', width = 1)) %>%
  layout(title="Attraction Centers of ID 27",
         scene=list(xaxis = list(title ="Easting", range=c(min(dat$utmlong-50000),
                                                           max(dat$utmlong+50000))),
                    yaxis = list(title = "Northing", range=c(min(dat$utmlat-50000),
                                                             max(dat$utmlat+50000))),
                    zaxis = list(title = "Time Segment", range=c(0,nrow(attr.cntrs27)+1)),
                    camera = list(eye = list(x = -1, y = -1.5, z = 1.5),
                                  up = list(x=1, y=0, z=1), 
                                  center=list(x=0, y=0, z=0))
         ))






id1_3d<- api_create(time.scatter1)
id12_3d<- api_create(time.scatter12)
id19_3d<- api_create(time.scatter19)
id27_3d<- api_create(time.scatter27)

id1_3d













