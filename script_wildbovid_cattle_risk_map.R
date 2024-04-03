# R code for creating the hotspot maps of the risk of diseases transmission 
# between wild bovid and domestic livestock
# Mar 2024
# include: Figure 4.1, 4.2, 4.4, Table 4.1-4.3 from the manuscript:
# 

#rm(list = ls(all.names = T))

library(sf)
library(sp)
library(rgdal)
library(raster)
library(tidyverse)
library(scales)
library(lattice)
library(rasterVis)
library(xlsx)
library(CoordinateCleaner) # for cleaning  coordinate

set.seed(111)

# Gridded Livestock of the World – 2015 (GLW4) ---------------------------------------------------
# map: https://dataverse.harvard.edu/file.xhtml?fileId=6769711&version=1.0
# article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10687097/

#import rasters:
# 1) cattle density, already cropped to Thailand
c<-raster("./cattle_cr_Da_2015.tif")
plot(c)

# 2) bovid's habitat suitability [3 species: gaur, banteng and wild buffalo]
habit<- list.files(pattern = "_wm_")
habit

bg<-raster(habit[[4]]) #gaur
bj<-raster(habit[[6]]) #banteng
ba<-raster(habit[[2]]) #wild buffalo

habit2<-raster::stack(bg,bj,ba)
habit2
plot(habit2)

#check maximum,minimum values of habitat suitability
mini<-list()
maxi<-list()

for ( i in 1:3) {
  mini[[i]] <- min(na.omit(values(habit2[[i]])))
  
  maxi[[i]] <- max(na.omit(values(habit2[[i]])))
}
mini
maxi

#rescale to 0-1
for ( i in 1:3) {
  values(habit2[[i]]) <-( values(habit2[[i]]) - mini[[i]])/ (maxi[[i]]-mini[[i]])
}

habit2

#disaggregate from 1 cell/10 km^2 to 1 cell/km^2
c_dis<- raster::disaggregate(c, fact = 10)
c_dis

#plot map
plot(c_dis, main = 'Cattle density',cex.main = 1.5)
hist(c_dis, main = 'Cattle density',cex.main = 1.5)

# rescale
# X2 = (X1 - Xminimum)/(Xmaximum - Xminimum)
c_rmax = cellStats(c_dis, "max")
c_rmin = cellStats(c_dis, "min")

c_re<-(c_dis-c_rmin)/(c_rmax-c_rmin)
# or just
c_re<-c_dis/c_rmax

plot(c_re, main = 'Cattle density', cex.main = 1.5)

# resample with gaur's map
c_re2<-raster::resample(x=c_re, y=bg,method="bilinear")
extent(c_re2)==extent(bg)

#plot resample map
my.colors = colorRampPalette(c("grey95","lightblue", "yellow","orangered", "red"))
plot(c_re2,main = 'Cattle density resample', cex.main = 1.5,
     col=my.colors(255)) 

#write rater (.tif) file:
#writeRaster(c_re2, 
#            filename = "cattle_rescale01_Da_2015", 
#            bylayer=TRUE, 
#            format='GTiff',
#            overwrite = TRUE)

# Bivariate map  --------------------------------------------
#import cattle density
#c_re2<- raster("cattle_rescale01_Da_2015.tif")

png("risk_matrices.png",width = 7, height = 8,units = "cm",res=300)
col.matrix <- bivariatemaps::colmat( nquantiles=4,
                                     upperleft  = rgb(255,230,15, maxColorValue=255), 
                                     upperright = rgb(130,0,80, maxColorValue=255), 
                                     bottomleft= "grey", 
                                     bottomright= rgb(0,150,235, maxColorValue=255), 
                                     xlab="Habitat suitability", 
                                     ylab="Livestock density")

dev.off()

#Adjust bivariate.map fx 
# add:   brks2 <- unique(brks)
#change: brks to brks2

bivariate.map_u<- 
  function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
    
    quanmean<-getValues(rasterx)
    temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
    brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
    brks2 <- unique(brks)
    r1<-within(temp, quantile <- cut(quanmean, breaks = brks2, labels = 2:length(brks2),include.lowest = TRUE))
    quantr<-data.frame(r1[,2]) 
    
    quanvar<-getValues(rastery)
    temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
    brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
    brks2 <- unique(brks)
    r2<-within(temp, quantile <- cut(quanvar, breaks = brks2, labels = 2:length(brks2),include.lowest = TRUE))
    quantr2<-data.frame(r2[,2])
    
    as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
    col.matrix2<-colormatrix
    cn<-unique(colormatrix)
    
    for(i in 1:length(col.matrix2)){
      ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
    cols<-numeric(length(quantr[,1]))
    
    for(i in 1:length(quantr[,1])){
      a<-as.numeric.factor(quantr[i,1])
      b<-as.numeric.factor(quantr2[i,1])
      cols[i]<-as.numeric(col.matrix2[b,a])}
    r<-rasterx
    r[1:length(r)]<-cols
    return(r)
  }

bivmap<-list()
habit3<-unstack(habit2)

#x= habitat, y= livestock ไม่งั้นจะผิดสีตอน plot
for ( i in 1:length(habit3)) { 
  bivmap[[i]]<- bivariate.map_u(habit3[[i]], c_re2,
                                colormatrix=col.matrix, 
                                nquantiles=4) 
}
bivmap2<-raster::stack(bivmap)
# Figure 1: Bivariate maps::cattle density and habitat suitability -------------------------

nam_plot<-c("Gaur","Banteng","Wild water buffalo")

png("risk_map_3sp_test.png",width = 20, height = 10,units = "cm",res=600)

plot(bivmap2,frame.plot=F,axes=F,box=F,add=F,legend=F,
       col=as.vector(col.matrix), main=nam_plot, cex.main = 1.2,nc=3) 

dev.off()

# Overlap area ---------------------

# import binary map 
b<-intersect(list.files(pattern = "tha"),
             list.files(pattern = "bin"))
b
#select only gaur, baneteng and buffalo
b2<-raster::stack(b[[2]],b[[3]],b[[1]])
plot(b2)

cellStats(c_re2, "mean")
summary(c_re2,maxsamp=1000000000)
summary(b2,maxsamp=1000000000)
summary(habit2,maxsamp=1000000000)

#setting up a mean of cattle density
threshold <- cellStats(c_re2, "mean")

ch <- c_re2 > threshold # High cattle density : > mean
cl <- c_re2 <= threshold # Low cattle density : <= mean

plot(ch)
plot(cl) 

# binary map = 1 (suitable) assign to habit1
habit1 <- b2 == 1
habit1

#binary map = 0 (unsuitable) assign to habit0
habit0 <- b2 == 0
habit0
# overlap and recode
# high habitat suitability & high cattle density = 4
overlap <- habit1 & ch 
hh <- reclassify(overlap, cbind(1, 4))
hh
# high habitat suitability & low cattle density = 3
overlap2 <- habit1 & cl
hl <- reclassify(overlap2, cbind(1, 3))
hl
# low habitat suitability & high cattle density = 2
overlap3 <- habit0 & ch  
lh <- reclassify(overlap3, cbind(1, 2))
lh
# low habitat suitability & low cattle density = 1
overlap4 <- habit0 & cl  
ll <- reclassify(overlap4, cbind(1, 1))
ll

m<-hh+hl+lh+ll

names(m)<-c( "Gaur","Banteng","Buffalo")
m

# Figure 2: Overlap area of High cattle - High habitat suitability ------------
par(mfcol = c(1, 3))

nam_plot<-c( "Gaur","Banteng","Wild water buffalo")

png("high_risk_3sp.png",width = 20, height = 10,units = "cm",res=600)
raster::plot(hh, frame.plot = FALSE, axes = FALSE, box = FALSE, 
             add = FALSE, legend = F, 
             col = adjustcolor(as.vector(col.matrix), 
                               alpha = 1),
             main = nam_plot,
             cex.main = 1.2, nc=3)
dev.off()

# Save Raster:: run this if we want to save High-High overlap raster 
#for(i in 1:nlayers(m)) {
#  writeRaster(m[[i]], 
#              paste0(nam[[i]], "_overlap_allrisk"), 
#              format="GTiff",overwrite=TRUE)
#}

# plot biviate maps with four risk colors (levels)
# hh hl lh ll map, 4 colours
plot(m, frame.plot = FALSE, axes = FALSE, box = FALSE, add = FALSE, legend = F, 
       col = adjustcolor(as.vector(col.matrix), alpha = 1),
       main = nam_plot, cex.main = 1.2,nc=3) 

# end bivariate plots #

# > import protected areas raster
#unstack high-high
high<-unstack(overlap)
high

pa_th <- raster("PA_thai_re.tif")            
pa_th <- resample(pa_th, high[[1]], method="ngb")

extent(pa_th)==extent(high[[1]])

plot(pa_th, main ='PA Thai')

# Calculate high-high areas (overlap areas) within inside and outside protected areas

zonlist<-list()
zonalboth<-list()
check<-list()

nam<-c("Gaur","Banteng","Buffalo")

for (i in 1:length(high)){
  
  values(high[[i]])[values(high[[i]]) > 0] = 1
  
  # Now croping by presence or absence 
  pres <- high[[i]]
  plot(pres) #high cattle - high habitat
  unique(values(pres))
  
  values(pres)[values(pres) < 1] = NA
  
  absences <- high[[i]]
  plot(absences) # the other area than high-high
  
  values(absences)[values(absences) > 0] = NA
  
  # Now crossing with protected areas and presences
  zonalpres<- data.frame(zonal(area(pres, na.rm=TRUE), pa_th, fun='sum'))
  
  colnames(zonalpres) <- c('Protected area code', 'Area sqkm')
  
  # Area of presence per IUCN type
  
  zonalpres$codef <- as.factor(zonalpres$`Protected area code`)
  
  zonalpres$sp_condition <- paste0(nam[[i]],'_high')
  
  # rename IUCN PA category
  zonaltidy <-zonalpres %>%  mutate(label = fct_recode(codef,
                                                       "IUCN PA category Iab" = "1" ,
                                                       "IUCN PA category II" ="2" ,
                                                       "IUCN PA category V" = "5",
                                                       "Not Applicable"="7",
                                                       "Not protected" = "8")) %>% 
    select(-c(codef)) %>% arrange(desc(-`Protected area code`))%>%
    as.data.frame()%>%
    mutate(percentage_sep = round((`Area sqkm`/ sum(`Area sqkm`) * 100),digits=2)) 
  
  # Crossing with absences
  # Now crossing with protected areas and ABSENCES
  zonalabsent <- data.frame(zonal(area(absences, na.rm=TRUE), pa_th, fun='sum') )
  
  colnames(zonalabsent) <- c('Protected area code', 'Area sqkm')
  
  # Area of absence per IUCN type
  # 8 == unprotected
  
  zonalabsent$codef <- as.factor(zonalabsent$`Protected area code`)
  
  zonalabsent$sp_condition <- paste0(nam[[i]],'_low')
  
  zonaltidya<-zonalabsent %>%  mutate(label = fct_recode(codef,
                                                         "IUCN PA category Iab" = "1" ,
                                                         "IUCN PA category II" ="2" ,
                                                         "IUCN PA category V" = "5",
                                                         "Not Applicable"="7",
                                                         "Not protected" = "8")) %>% 
    select(-c(codef)) %>% arrange(desc(-`Protected area code`))%>%
    as.data.frame()%>%
    mutate(percentage_sep = round((`Area sqkm`/ sum(`Area sqkm`) * 100),digits=2)) 
  
  zonalboth[[i]]<- rbind(zonaltidy, zonaltidya) 
  str(zonalboth[[i]])
  
  # calculate total amount of area inside protected areas (km^2)
  zonlist[[i]] <-
    zonalboth[[i]]%>%
    as.data.frame()%>%
    mutate(percentage_tot = (`Area sqkm`/ sum(`Area sqkm`) * 100)) 
  
  df<-bind_rows(zonlist)
  str(df)
  write.xlsx(df,'Table_Overlap_Thai_test.xlsx', row.names = FALSE)
}

# Table 4.1 : Overlap area ------
df$sp_condition<-as.factor(df$sp_condition)
str(df)

print(n=50,
      df|> 
        group_by(sp_condition,label) |>
        dplyr::summarise(Area = sum(`Area sqkm`),
                         Percent = `percentage_sep`,
                         Percent_tot = `percentage_tot`))
df|> group_by(sp_condition) |>
  dplyr::summarise(Area = sum(`Area sqkm`))

#inside PA
df2<-df |>  
  filter(!label %in% c("Not protected"),
         sp_condition %in% c("Gaur_high","Banteng_high","Buffalo_high")) |>
  group_by(sp_condition) |>
  dplyr::summarise(inside = sum(`Area sqkm`)) 
df2
#outside PA
df3 <-df |>  filter(label %in% c("Not protected"),
                    sp_condition %in% c("Gaur_high","Banteng_high","Buffalo_high")) |>
  select(c(`Area sqkm`))|>
  rename("outside"=`Area sqkm`)

df3
#thailand area square kilometer
tha<-514410
df4 <- df|>    
  group_by(sp_condition) |>
  filter(sp_condition %in% c("Gaur_high","Banteng_high","Buffalo_high")) |>
  dplyr::summarise(Area = sum(`Area sqkm`)) |>
  mutate(percent = round((Area / tha * 100),digits=2)) |>
  select(-c("sp_condition"))

df5<-cbind(df2,df3,df4)
df5

# import the disease occurrence 
d <- read.csv("disease_cattleonly.csv")
names(d)
str(d)

d$lab_diag<-as.factor(d$lab_diag)
d$animal<-as.factor(d$animal)
d$lab_code<-d$lab_diag

str(d)

table(d$lab_diag)

d<-d %>%  mutate(lab_code = fct_recode(lab_code,
                                       "1" =  "Tuberculosis",
                                       "2" =  "Hemorrhagic Septicemia",
                                       "3" =  "LSD",
                                       "4" =  "FMD" ,
                                       "5" =  "Brucellosis"))
str(d)

d2<-clean<-cc_dupl(
  d,
  lon = "x",
  lat = "y",
  species = "lab_diag",
  additions = NULL,
  value = "clean",
  verbose = TRUE)

str(d2)
View(d2)
table(d2$lab_diag)

coordinates(d2) <- ~ x + y
str(d2)

# Create an sf object
d3 <- st_as_sf(d2, coords = c("x","y"), 
               crs = crs("+init=epsg:4326"))
str(d3)

# > Count point in the overlap area
#merge variables and ENM 

risk <- 
  cbind(raster::extract(x = m, y=d3,
                              field = "lab_code", 
                              bind=T,df=T),d2)|>
  drop_na()

head(risk)
print (risk|> group_by(lab_diag) |>
         dplyr::summarise(
           N = n()))
table(risk$lab_diag)

write.xlsx(risk,'risk_overlap_cattle.xlsx')

risk1<- st_as_sf(risk, coords = c("x", "y"), 
               crs = crs("+init=epsg:4326"))

n <- c("Tuberculosis","Hemorrhagic Septicemia","LSD", "FMD", "Brucellosis" )

test<-risk1%>% filter(lab_diag  %in%  n[1]) 

test

# disease occurrences within high and low cattle density
# 0 = low cattle density; < mean 
# 1 = high cattle density; >= mean

count<- list()
count2 <- list()
count3<-list()
#extract raster and bind with disease occurrence data
for (i in 1:length(n)) {
  count[[i]]<- raster::extract(x = ch, 
                               y= risk1 %>% filter(lab_diag  %in%  n[i]), 
                               field = "lab_diag", bind=T,df=T)
  count2[[i]]<- raster::extract(x = cl, 
                                y= risk1 %>% filter(lab_diag  %in%  n[i]), 
                                field = "lab_diag", bind=T,df=T)
}

for (i in 1:length(n)) {
  count3[[i]]<- raster::extract(x = overlap, 
                                y= risk1 %>% filter(lab_diag  %in%  n[i]), 
                                field = "lab_diag", bind=T,df=T)
  
}                   

count3

# disease occurrences within high and low cattle density

count
tb<-table(count[[1]]$layer)
hs<-table(count[[2]]$layer)
lsd<-table(count[[3]]$layer)
fmd <-table(count[[4]]$layer)
bru <-table(count[[5]]$layer)

a  <-rbind(tb,hs,lsd,fmd,bru)
a

# 1 = high
tb2<-table(count2[[1]]$layer)
hs2<-table(count2[[2]]$layer)
lsd2<-table(count2[[3]]$layer)
fmd2 <-table(count2[[4]]$layer)
bru2 <-table(count2[[5]]$layer)

a2  <-rbind(tb2,hs2,lsd2,fmd2,bru2)
a2

#Calculate high and low cattle density

unique(values(ch) )

values(ch)[values(ch) < 1] = NA
values(cl)[values(cl) < 1] = NA

plot(ch,col='red')
plot(cl, col='blue')

# whole raster layer
ar_ch<- area(ch, na.rm=TRUE)
ar_cl<- area(cl, na.rm=TRUE)

# size of the cell in km2
ar_ch2 <- sum(values(ar_ch), na.rm = TRUE) ### applies a correction for latitude, to km2
ar_cl2 <- sum(values(ar_cl), na.rm = TRUE)
# Total area
ar_ch2+ar_cl2

area_high <- sum(values(area(ch, na.rm=TRUE)), na.rm=TRUE)
area_low <- sum(values(area(cl, na.rm=TRUE)), na.rm=TRUE)

# Almost equal, so estimate is good enough
sum(area_high,area_low)== ar_ch2+ar_cl2
round((area_high + area_low), digits = 2) -  round(ar_ch2+ar_cl2, digits=2)

#table 4.2 outbreak per event ------
a
a<-as.data.frame(a)
# low cattle density area
ar_cl2
sum(a$`0`)/ar_cl2
# high cattle density area
ar_ch2 
sum(a$`1`)/ar_ch2 

#total area
sum(a)/ (ar_ch2+ar_cl2)

# 
#Calculate high - high area 
m2 <- unstack(m)
m2

m3<-list()
for (i in 1:length(m2)){
  names(m2)[i]
  m3[[i]]<-tapply(area(m2[[i]]), m2[[i]][], sum)
}

m3

# Table 4.3 :: outbreak events of bovine infectious disease occurrences in the high risk areas-------

risk|> 
  filter(Gaur %in% 4)|>
  group_by(Gaur, lab_diag)|>
  dplyr::summarise(N = n())

risk|> 
  filter(Banteng %in% 4)|>
  group_by(Banteng, lab_diag)|>
  dplyr::summarise(N = n())

risk|> 
  filter(Buffalo %in% 4)|>
  group_by(Buffalo, lab_diag)|>
  dplyr::summarise(N = n())

# Figure 4.4 ::Interface areas of three wild bovid vs land use types. ------------

#high-high (overlap) areas
plot(overlap)

# import land use types, Thailand
l<-raster("./land_cgls_thai.tif")
plot(l)

# Calculate overlap area by land use type
# Overlap == high cattle density - high habitat suitability 

zonlist<-list()
zonalboth<-list()
check<-list()

# Calculate overlapped area in the land use

for (i in 1:length(high)){
  
  names(high)[i]
  values(high[[i]])[values(high[[i]]) > 0] = 1
  
  # Now croping by presence or absence 
  pres <- high[[i]]
  
  unique(values(pres) )
  
  values(pres)[values(pres) < 1] = NA
  
  absences <- high[[i]]
  values(absences)[values(absences) > 0] = NA
  
  plot(absences,col='red')
  
  plot(pres, col='green')# add=TRUE,
  
  # whole raster layer
  ar<- area(high[[i]], na.rm=TRUE)
  
  # size of the cell in km2
  ar2 <- sum(values(ar), na.rm = TRUE) ### applies a correction for latitude, to km2
  
  # Total area
  ar2 
  
  area_present <- sum(values(area(pres, na.rm=TRUE)), na.rm=TRUE)
  
  area_absences<- sum(values(area(absences, na.rm=TRUE)), na.rm=TRUE)
  
  # Almost equal, so estimate is good enough
  sum(area_present,area_absences)==ar2
  round((area_present + area_absences), digits = 2) -  round(ar2, digits=2)
  
  # Now crossing with protected areas and presences
  
  zonalpres<- data.frame(zonal(area(pres, na.rm=TRUE), l, fun='sum'))
  
  colnames(zonalpres) <- c('land type', 'Area sqkm')
  
  # Checking 
  sum(zonalpres$`Area sqkm`) - area_present
  
  zonalpres$codef <- as.factor(zonalpres$`land type`)
  
  zonalpres$sp_condition <- paste0(nam[[i]])
  
  zonaltidy <-zonalpres %>%  mutate(label = fct_recode(codef,
                                                       "Closed evergreen, needle leaf" = "111",
                                                       "Closed evergreen, broad leaf" = "112", 
                                                       "Closed deciduous broad leaf" = "114",
                                                       "Closed unknown" = "116", 
                                                       "Open deciduous  broad leaf" = "124",
                                                       "Open unknown" = "126",
                                                       "Shrubs" = "20",   # Shrubs
                                                       "Herbaceous vegetation" = "30",
                                                       "Cropland" = "40", # Cropland
                                                       "Urban" = "50", # Urban
                                                       "Bare" = "60", # Bare / sparse vegetation
                                                       "Permanent water bodies" = "80", # Permanent water bodies
                                                       "Herbaceous wetland" = "90",   
                                                       "Open Sea" = "200")) %>% 
    arrange(desc(-`land type`))%>%
    as.data.frame()%>%
    mutate(percentage_sep = round((`Area sqkm`/ sum(`Area sqkm`) * 100),digits=2)) 
  
  
  zonalboth[[i]]<- zonaltidy
  
  str(zonalboth[[i]])
  
  # Checking small divergence due to distortion
  check[[i]]<-sum(zonalboth[[i]]$`Area sqkm`) - ar2
  
  # calculate total amount of area inside protected areas (km^2)
  zonlist[[i]] <-
    zonalboth[[i]]%>%
    as.data.frame()%>%
    mutate(percentage_tot = (`Area sqkm`/ sum(`Area sqkm`) * 100)) 
  
  ov<-bind_rows(zonlist)
  str(ov)
  xlsx::write.xlsx(ov,'Table_Overlap_landuse_Thai_test.xlsx', row.names = FALSE)
}

str(ov)

# bar plot
ggplot(ov, aes(x=reorder(label,`Area sqkm`), y=`Area sqkm`,fill=sp_condition))+
  geom_bar(stat='identity', position=position_dodge())+
  # facet_wrap(.~sp_condition,scales = "free_y") +
  
  theme_bw()+
  labs(title = "Interfce areas by land use types",
       x = "Land use",
       y=expression(Area (km^2)))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(text = element_text(size = 10))  +
  theme(axis.text.x = element_text(angle = 90, size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_brewer(palette= "Blues",
                    name  ="Species")+
  coord_flip()

# --- END :) ---#