#### 
#### Crawford Drury*, Roberta E. Martin, David E. Knapp, Joseph Heckler, Joshua Levy, Ruth Gates, Gregory P. Asner
#### *Analysis: Crawford Drury - crawford.drury@gmail.com

library(ggnewscale);library(readxl);library(janitor);library(cowplot);library(egg)
library(caret);library(gbm)
library(plotrix)
library(factoextra);library(vegan)
library(rsq)
library(purrr)
library(Hmisc)
library(betareg)
library(doMC);registerDoMC(cores=1)
library(tidyverse);library(ggdendro)
setwd("~/Projects/CAO/data/")

######################################### READ IN AND FORMAT DATA FOR SUBSTRATE MODEL ########################################################################################## #####
library(sf);library(spdep);library(raster);library(rgdal)
rm(list = ls())
reef_numbers<-as.numeric(c("13","25","42","44"))
for (reef in reef_numbers){
     substrate<-read_csv("./points/substrate2.csv")
     rb_data<-stack(paste0("patch",reef,"_20170930_test3_rb_img"))
     depth<-stack(paste0("patch",reef,"_depth"))
     coords<-filter(substrate,Site==paste0(reef))%>%dplyr::select(x,y)
     meta<-filter(substrate,Site==paste0(reef))%>%dplyr::select(Site,Substrate)
     depth_values<-as.data.frame(raster::extract(depth,coords,method="simple"));names(depth_values)[1]<-"depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords,method="simple"))
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     assign(paste("substrate",reef,sep="_"),cbind(meta,coords,depth_values,bn_values))
}
substrate_data<-do.call(rbind, lapply(ls(pattern = "substrate_"), get))
rm(list=ls(pattern="reef"));rm(list=ls(pattern="rb_"));rm(list=ls(pattern="depth"));rm(list=ls(pattern="meta"));rm(list=ls(pattern="sums"));rm(list=ls(pattern="bn"));rm(list=ls(pattern="coords"))
saveRDS(substrate_data,"./output/processed_spectra/substrate_data_bn")
detach(package:spdep);detach(package:sf);detach(package:raster);detach(package:rgdal)

######################################### FIT SUBSTRATE ######################################################################################################################## ##### 
rm(list = ls())
data<-readRDS("./output/processed_spectra/substrate_data_bn")
substrate_trim<-data%>%select(-x,-y)%>%mutate(Site=as.factor(Site))
set.seed(3839)
train_control <- trainControl(method="cv",number=10,search="grid")
#grid<-expand.grid(interaction.depth=5,n.trees=1000,shrinkage=0.01,n.minobsinnode=5)
#grid<-expand.grid(interaction.depth=seq(1,5),n.trees=seq(100,1000,100),shrinkage=c(0.005,0.01,0.02,0.05),n.minobsinnode=seq(1,5))
grid<-expand.grid(interaction.depth=4,n.trees=300,shrinkage=0.02,n.minobsinnode=1)
seeds <- vector(mode = "list", length = 11) #nrepeats * nfolds +1
for(i in 1:11) seeds[[i]]<- sample.int(n=nrow(grid), nrow(grid)) #last value is at least grid length
seeds[[11]]<-sample.int(1000, 1)

sub_model<-caret::train(Substrate ~., data=substrate_trim, method="gbm",
                         trControl=train_control,
                         preProcess=c("center","scale"),
                         tuneGrid=grid,
                         metric="Kappa",
                         verbose=FALSE)

sub_model
quartz()
plot(sub_model)
confusionMatrix(sub_model)
sub_model$bestTune
sub_model$resample
mean(sub_model$resample$Accuracy);sd(sub_model$resample$Accuracy)
prediction<-predict(sub_model,newdata=substrate_trim,type=c("raw")); cm<-confusionMatrix(prediction,as.factor(substrate_trim$Substrate));cm
saveRDS(sub_model,"./output/model_fits/0_substrate_model_final")
list<-varImp(sub_model)
table<-as.data.frame(list$importance)%>%rownames_to_column(var="Wavelength")%>%filter(grepl("Nanometers",Wavelength))%>%
     mutate(max=max(Overall))%>%mutate(Overall=Overall/max)%>%select(-max)%>%arrange(desc(Overall))%>%
     separate(Wavelength,into=c("trash","Wavelength"),sep=1)%>%select(-trash)%>%
     separate(Wavelength,into=c("Wavelength","trash"),sep=3)%>%select(-trash)%>%rename('Relative VI'=2)%>%arrange(Wavelength)
write.table(table,"../revision/WebTable1.txt",quote=FALSE,row.names = FALSE)
#sub_model<-readRDS("./output/model_fits/0_substrate_model_final")



######################################### CREATE SUBSTRATE LAYERS ############################################################################################################## ##### 
library(sf);library(spdep);library(raster);library(rgdal)
rm(list = ls())
sub_model<-readRDS("./output/model_fits/0_substrate_model_final")
reef_numbers<-as.numeric(c("13","25","42","44"))

for (reef in reef_numbers){
     depth<-stack(paste0("patch",reef,"_depth"))
     rb_data<-stack(paste0("patch",reef,"_20170930_test3_rb_img"))
     extent<-extent(rb_data)
     x<-seq(extent@xmin+.2,extent@xmax+.2,.4)
     y<-seq(extent@ymin+.2,extent@ymax+.2,.4)
     total<-expand.grid(x,y);colnames(total)[1]<-"x";colnames(total)[2]<-"y"
     coords<-total%>%dplyr::select(x,y)%>%rownames_to_column()%>%rename(index=rowname)
     depth_values<-as.data.frame(raster::extract(depth,coords[2:3],method="simple"));names(depth_values)[1]<-"depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords[2:3],method="simple"))
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     substrate<-cbind(coords,depth_values,bn_values)%>%mutate(Site=as.factor(reef))
     substrate<-substrate[complete.cases(substrate), ]
     prediction<-as.data.frame(predict(sub_model,newdata=substrate,type=c("raw")));colnames(prediction)[1]<-"prediction"
     output<-cbind(substrate[,1:3],as.data.frame(prediction))
     saveRDS(output,paste("./output/substrate_prediction",reef,sep="_"))
     coordinates(output)=~x+y
     proj4string(output)=CRS("+proj=utm +zone=4N +datum=WGS84")
     raster::shapefile(output,paste("./output/substrate_prediction",reef,sep="_"),overwrite=TRUE)
}
detach(package:spdep);detach(package:sf);detach(package:raster);detach(package:rgdal)


######################################### F1 - SPECTRA ######################################################################################################################### #####
raw<-readRDS("./output/processed_spectra/substrate_data_bn")
wavelength_labels<-readRDS("./output/misc/wavelength_labels");colnames(raw)[6:57]<-wavelength_labels
data<-raw%>%select(-x,-y,-depth,-Site)

sub_mean<-data%>%group_by(Substrate)%>%summarise_all(funs(mean))
sub_SD<-data%>%group_by(Substrate)%>%summarise_all(funs(sd))
{
     a<-as.data.frame(t(((sub_SD[1,2:53]))))
     b<-as.data.frame(t(((sub_SD[2,2:53]))))
     c<-as.data.frame(t(((sub_SD[3,2:53]))))
     d<-as.data.frame(t(((sub_SD[4,2:53]))))
     e<-as.data.frame(t(((sub_SD[5,2:53]))))
     
     
     a1<-cbind(as.data.frame(t(sub_mean[1,2:53])),a);colnames(a1)[1]<-"Mean";colnames(a1)[2]<-"SD"
     b1<-cbind(as.data.frame(t(sub_mean[2,2:53])),b);colnames(b1)[1]<-"Mean";colnames(b1)[2]<-"SD"
     c1<-cbind(as.data.frame(t(sub_mean[3,2:53])),c);colnames(c1)[1]<-"Mean";colnames(c1)[2]<-"SD"
     d1<-cbind(as.data.frame(t(sub_mean[4,2:53])),d);colnames(d1)[1]<-"Mean";colnames(d1)[2]<-"SD"
     e1<-cbind(as.data.frame(t(sub_mean[5,2:53])),e);colnames(e1)[1]<-"Mean";colnames(e1)[2]<-"SD"
     
     a2<-a1%>%rownames_to_column()%>%mutate(upper=(Mean+SD))%>%mutate(lower=(Mean-SD))%>%mutate(Substrate="Dead");colnames(a2)[1]<-"Wavelength"
     b2<-b1%>%rownames_to_column()%>%mutate(upper=(Mean+SD))%>%mutate(lower=(Mean-SD))%>%mutate(Substrate="Montipora");colnames(b2)[1]<-"Wavelength"
     c2<-c1%>%rownames_to_column()%>%mutate(upper=(Mean+SD))%>%mutate(lower=(Mean-SD))%>%mutate(Substrate="Porites");colnames(c2)[1]<-"Wavelength"
     d2<-d1%>%rownames_to_column()%>%mutate(upper=(Mean+SD))%>%mutate(lower=(Mean-SD))%>%mutate(Substrate="Sand");colnames(d2)[1]<-"Wavelength"
     e2<-e1%>%rownames_to_column()%>%mutate(upper=(Mean+SD))%>%mutate(lower=(Mean-SD))%>%mutate(Substrate="Water");colnames(e2)[1]<-"Wavelength"
}

data<-rbind(a2,b2,c2,d2,e2)
data<-data%>%mutate(Substrate= replace(Substrate,Substrate=="Dead","Hardbottom"))%>%filter(Substrate!="Water")
rm(a,b,c,d,e,a1,b1,c1,d1,e1,a2,b2,c2,d2,e2)

plot1b<-ggplot(data)+
     geom_line(aes(as.numeric(Wavelength),Mean,color=Substrate))+
     geom_ribbon(aes(as.numeric(Wavelength),ymin=lower,ymax=upper,fill=Substrate),alpha=.5)+
     theme_classic(base_size=6)+
     scale_y_continuous(lim=c(0,.23),breaks=seq(0,.2,.05))+
     scale_x_continuous(breaks=seq(425,675,75))+
     xlab("Wavelength")+
     ylab("bn Reflectance")+
     scale_fill_manual(values=c("dimgray","blue","green","tan","yellow"),name="")+
     scale_color_manual(values=c("dimgray","blue","green","tan","yellow"),name="")+
     theme(legend.position ="none",
           legend.justificatio=c("left","top"),
           legend.key.size = unit(0.15,"cm"),  
           legend.box.margin=margin(0,0,-10,15),
           plot.background = element_rect(colour = "gray35", fill="white", size=1),
           axis.title.x=element_blank())

######################################### F1 - SUBSTRATE MAP ################################################################################################################### #####
library(sf);library(ggrepel);library(ggsn);library(cowplot);library(tidyverse);library(patchwork);library(jpeg);library(ggnewscale)
raw_sub<-st_intersection(st_read("./output/substrate_prediction_13.shp"),st_read("./points/clip.shp"))
shp_sub<-fortify(as.data.frame(cbind(raw_sub,st_coordinates(raw_sub))))%>%filter(prediction!="Water")

r13_substrate<-ggplot()+
     geom_raster(data = shp_sub, aes(x = X, y = Y,fill=prediction))+
     coord_fixed(ratio = 1,xlim=c(624650,624980),ylim=c(2372480,2372760))+
     theme_classic(base_size=8)+
     scale_fill_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),
           #legend.position=c(0.8,0.7),
           legend.position="none",
           axis.line=element_blank(),
           axis.ticks=element_blank(),
           legend.key.size=unit(0.25,"cm"),
           panel.border = element_rect(colour = "gray35", fill=NA, size=1))+
     annotate("text",x=624645,y=2372486,label="Water not shown",size=2,hjust=0)+
     annotate("text",x=624990,y=2372486,label=" Best-Fit Accuracy: 99.7%\nCV-Accuracy: 89.1%",size=2,hjust=1)+
     annotate("text",label="Reef 13",x=624645,y=2372760,size=2.5,color="black",hjust=0);r13_substrate

######################################### F1 - LOADING ######################################################################################################################### #####

oahu<-st_transform(st_read("~/Projects/Oahu Map/Layers/coast_n83.shp"),4326)
cropped<-st_crop(oahu,xmin=-158.8,xmax=-157.4,ymin=20.8,ymax=21.8)
fringe<-st_read("~/Projects/Oahu Map/Layers/Fringing Reef.shp")
patch<-st_read("~/Projects/Oahu Map/Layers/Patches.shp")%>%mutate(type="Reef")
raw_points<-st_transform(st_read("~/Projects/Oahu Map/Layers/4reef_numbers.shp"),crs=4326)%>%mutate(type="Study Sites")
kbay_points<-cbind(as.data.frame(st_coordinates(raw_points)),raw_points$Reef)%>%rename(Reef=3)
x1<-624688
x2<-x1+15
y1<-2372570
y2<-y1+15

library(raster);library(RStoolbox)
reef13_dimac<-raster::stack("~/Projects/CAO/data/patch13_20170905_dimac.tif")
extent <- raster::extent(x1-5,x2+5,y1-5,y2+5)
dimac_trim <- raster::crop(reef13_dimac, extent)
substrate_points<-read_tsv("./points/substrate.txt")
known_points<-read_csv("./points/phenotype_A.csv")%>%mutate(Species=str_to_title(Species))%>%mutate(Phenotype=str_to_title(Phenotype))
dimac_color<-RStoolbox::ggRGB(dimac_trim,r=1,b=3,g=2,ggObj=FALSE,stretch="none")
reef13_overall<-stack("~/Projects/CAO/data/patch13_20170930_dimac_match")
dimac_big<-RStoolbox::ggRGB(reef13_dimac,r=1,b=3,g=2,ggObj=FALSE,stretch="none")
depth_data<-fortify(raster::crop(stack("./patch13_depth"),extent))%>%rename(Depth=3)
detach(package:raster);detach(package:RStoolbox)

inset<-ggplotGrob(ggplot()+
                       geom_sf(data=cropped)+
                       theme_minimal(base_size=8)+
                       theme(axis.text.x=element_blank(),
                             axis.title.y=element_blank(),
                             axis.title.x=element_blank(),
                             axis.text.y=element_blank(),
                             panel.grid.minor=element_blank(),
                             panel.grid.major=element_blank(),
                             panel.border = element_rect(colour = "gray35", fill=NA, size=1),
                             panel.background = element_rect(color="white"))+
                       annotate("rect",xmin=-157.86,xmax=-157.78,ymin=21.425,ymax=21.51,color="blue",alpha=0.25,size=0.25)+
                       annotate("text",x=-158.05,y=21.52,label="Oʻahu",size=2.5))

main<-ggplot()+
     geom_sf(data=patch,aes(fill=type),colour="NA",show.legend="polygon")+
     geom_sf(data=fringe,fill="turquoise4",colour="NA")+
     geom_sf(data=raw_points,aes(color=type),size=1,show.legend="point")+
     geom_sf(data=cropped)+
     geom_text_repel(data=kbay_points,aes(X,Y,label=Reef),size=3,seed=3,nudge_x=0.01)+
     scale_fill_manual(values=c("turquoise4"),
                       guide = guide_legend(override.aes = list(linetype = c("blank"), 
                                                                shape = NA)))+
     scale_color_manual(values=c("orange"))+
     theme_classic(base_size=8)+
     coord_sf(ylim=c(21.425,21.51),xlim=c(-157.86,-157.78))+
     theme(#legend.position="none",
           legend.position=c(0.01,0.07),
           legend.justification=c("left","bottom"),
           legend.background = element_rect(fill=alpha('blue', 0)),
           legend.key.size=unit(0.3,"cm"),
           legend.title=element_blank(),
           axis.text.x=element_blank(),
           axis.title.y=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.line=element_blank(),
           panel.border = element_rect(colour = "gray35", fill=NA, size=1),
           legend.spacing.y=unit(-0.2,"line"),
           legend.box.margin=margin(1,1,1,1))+
     scalebar(x.min=-157.78,x.max=-157.86,y.min=21.426,y.max=21.511,transform=TRUE,dist=1,dist_unit="km",
              height = .01, st.dist = 0.03,
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2,
              location="bottomleft")+
     annotation_custom(inset,xmin=-157.815,xmax=-157.777,ymin=21.483,ymax=21.516);main

######################################### F1 - REEF 13 ######################################################################################################################### #####
#points<-bind_rows(known_points%>%select(Species,Phenotype,LAT,LONG)%>%clean_names()%>%rename(substrate=1,x=3,y=4),substrate_points%>%clean_names()%>%select(substrate,x,y)%>%mutate(phenotype=NA)%>%select(substrate,phenotype,x,y))%>%filter(x>x1&x<x2&y>y1&y<y2)
#write.table(points,"~/Desktop/example.txt",sep="\t",row.names=FALSE,quote=FALSE)
points<-read_tsv("./points/example.txt")

reso<-ggplotGrob(ggplot()+
                      geom_raster(aes(x,y,fill=fill),data=dimac_color,color=NA)+
                      scale_color_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
                      scale_fill_identity()+
                      #coord_fixed(ratio = 1,xlim=c(x1+10,x1+15),ylim=c(y1+1,y1+6))+
                      coord_fixed(ratio = 1,xlim=c(x1,x2),ylim=c(y1,y2))+
                      theme_minimal(base_size=6)+
                      theme(panel.border=element_rect(colour="white",fill=NA,size=1),panel.ontop = TRUE,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank(),legend.position="none")+
                      scale_x_continuous(breaks=seq(x1+5,x1+10,0.4))+scale_y_continuous(breaks=seq(y1+1,y1+6,0.4))+
                      #scale_x_continuous(breaks=seq(x1,x2,0.4))+scale_y_continuous(breaks=seq(y1,y2,0.4))+
                      geom_point(aes(x,y,color=substrate),data=points,size=1.5)+
                      annotate("text",label="Known/\nTraining Points",x=x1,y=y1+13,size=2.5,color="white",hjust=0))

dimac<-ggplot() +
     geom_point(aes(x,y,color=substrate),data=points,size=2)+
     geom_raster(data = dimac_big,aes(x = x, y = y,fill=fill),color=NA)+
     scale_fill_identity()+
     scale_color_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     coord_fixed(xlim=c(624650,624980),ylim=c(2372480,2372760))+
     theme_classic(base_size=8)+
     theme(axis.text.x=element_blank(),
           axis.title.y=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.line=element_blank(),
           axis.ticks=element_blank(),
           panel.border=element_rect(colour="gray35",fill=NA,size=1),
           legend.position=c(0.82,0.16),
           legend.key.size=unit(0.2,"cm"),
           legend.spacing.y = unit(0.1, "cm"),
           legend.background = element_rect(fill="white"))+
     scalebar(x.min=624650,x.max=624980,y.min=2372480,y.max=2372760,transform=FALSE,dist=25,dist_unit="m",
              height = .01, st.dist = 0.025,
              box.fill = c("black", "white"),
              box.color = "white", border.size = .1,st.size=2,st.color="white",
              location="bottomleft")+
     annotation_custom(reso,xmin=624843,xmax=Inf,ymin=2372618,ymax=Inf)+
     annotate("text",label="Reef 13",x=624645,y=2372760,size=2.5,color="white",hjust=0)+
     annotate("rect",xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=NA,color="white");dimac

######################################### F1 - TIMELINE ######################################################################################################################## #####
detach(package:sf);detach(package:ggsn)
library(tidyverse);library(lubridate);library(scales);library(readxl);library(janitor)

#downloaded year temp data from NDBC, concatenated, rest of processing happens below
#https://www.ndbc.noaa.gov/station_page.php?station=mokh1 
#'historical data'

rawtemp<-read_table("../../K_Bay_NDBC/data/master.txt")%>%rename(YY=1)%>%filter(YY!="#yr")%>%dplyr::select(YY,MM,DD,hh,mm,WTMP)%>%filter(WTMP<50)%>%unite(date,YY,MM,DD,sep="-")%>%unite(time,hh,mm,sep=":")%>%
     unite(date,date,time,sep=" ")%>%mutate(date=ymd_hm(date))%>%rename(temp=2)%>%mutate(temp=as.numeric(temp))%>%
     group_by(date(date))%>%summarise(mean=mean(temp))%>%rename(date=1)%>%
     mutate(a=month(date))%>%
     mutate(b=day(date))%>%
     unite(day,a,b,sep="-")%>%mutate(dummy_year=1999)%>%unite(dummy_date,dummy_year,day,sep="-")%>%mutate(year=as.factor(year(date)))%>%
     mutate(dummy_date=ymd(dummy_date)) #set up dummy date for plotting all years on same axes using different colors, is simply the day of the year with all observations set as year 1999

plotdata<-rawtemp%>%
     filter(year==2015|year==2016|year==2017|year==2018|year==2019)%>%
     mutate(color=case_when(mean>28.5~"red",mean<=28.5~"black"))
example<-read_xlsx("./points/dummydata.xlsx")%>%clean_names()%>%
     mutate(visual_status=as.factor(visual_status))%>%rownames_to_column()%>%mutate(rowname=as.numeric(rowname))%>%mutate(date=ymd(date))%>%filter(rowname!=9)%>%filter(rowname!=4)
example$visual_status<-factor(example$visual_status,levels=c("Healthy","Bleached"))

timeline<-ggplot()+
     geom_point(aes(date,y,fill=visual_status),data=example,pch=21,size=5)+
     geom_text(aes(date,y,label=historical_phenotype),data=example,color="black",size=2)+
     geom_line(aes(date,mean,color=mean),data=plotdata)+
     scale_color_gradient2(low="black",mid="blue",high="red",name="Temperature",midpoint=25,guide="none")+
     scale_fill_manual(values=c("gray","white"),name="Visual Status")+
     geom_hline(yintercept=28.5,linetype="dotted")+
     scale_x_date(date_breaks="6 months",minor_breaks=waiver(),labels=date_format("%b-%y"))+
     ylab("Temp (°C)")+
     theme_classic(base_size=8)+
     theme(legend.position="none",axis.title.x=element_blank(),
           legend.spacing.y=unit(0.1,"cm"))+
     coord_cartesian(ylim=c(27,31.5))+
     annotate("text",x=as_date("2015-2-5"),label="MMM +1°C",y=28.25,size=2,alpha=0.8,hjust=1,vjust=1,fontface="italic")+
     annotate("text",x=as_date("2015-10-10"),label="Historical Phenotype",y=31.5,size=2,hjust=1,vjust=0.5,fontface="italic")+
     annotate("text",x=as_date("2016-9-15"),label="Recovery",y=31.5,size=2,hjust=1,vjust=0.5,fontface="italic")+
     annotate("text",x=as_date("2017-9-15"),label="GAO Flight",y=31.5,size=2,hjust=1,vjust=0.5,fontface="italic")+
     #annotate("text",x=as_date("2019-4-15"),label="Metabolite Sampling",y=31.5,size=2,hjust=1,vjust=0.5,fontface="italic")+
     annotate("text",x=as_date("2019-10-15"),label="Validation",y=31.5,size=2,hjust=1,vjust=0.5,fontface="italic")+
     
     annotate("text",x=as_date("2015-11-15"),label="Visually Healthy",y=30.5,size=2,hjust=0,vjust=0.5,fontface="italic")+
     annotate("text",x=as_date("2015-11-15"),label="Visually Bleached",y=29.5,size=2,hjust=0,vjust=0.5,fontface="italic")+

     annotate("text",x=as_date("2016-10-25"),y=30,label="}")+
     annotate("text",x=as_date("2016-11-5"),y=30,label="Visually Indistinguishable",size=2,fontface="italic",hjust=0,vjust=0.7)

detach(package:lubridate);detach(package.scales)

w=7.2
h=3.31
quartz(w=w,h=h) #7.2 x 3
test<-ggplotGrob(plot1b)
out<-r13_substrate+annotation_custom(test,xmin=624828,xmax=624988,ymin=2372618,ymax=2372768)
timeline / (main|dimac|out) + plot_layout(heights=c(1.01,3))+
     plot_annotation(tag_levels = c('a'), tag_prefix = '(', tag_suffix = ')') +
     theme(plot.tag.position = c(0.15, 0.95),plot.tag = element_text(size = 8, hjust = 0, vjust = 0))

######################################### PHENOTYPE SPECTRA #################################################################################################################### #####
rm(list = ls())
setwd("~/Projects/CAO/data")
reef_numbers<-as.numeric(c("13","25","42","44"))
for (reef in reef_numbers){
     rawedits<-rbind(read.csv("./points/phenotype_A.csv"),read.csv("./points/phenotype_B.csv"))
     edits<-rawedits%>%select(ID,Site,Species,Phenotype,Tag,Date,LAT,LONG)%>%rename(x=LAT,y=LONG)%>%filter(Date=="A")
     rb_data<-stack(paste0("patch",reef,"_20170905_test3_rb_img"))
     depth<-stack(paste0("patch",reef,"_depth"))
     coords<-filter(edits,Site==paste0(reef))%>%select(x,y)
     meta<-filter(edits,Site==paste0(reef))%>%select(ID,Site,Species,Phenotype,Tag,Date,x,y)
     depth_values<-as.data.frame(raster::extract(depth,coords,method="simple"));names(depth_values)[1]<-"Depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords,method="simple"))
     wavelength_labels<-readRDS("wavelength_labels");colnames(rb_values)[1:52]<-wavelength_labels
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     assign(paste("data",reef,"A",sep="_"),cbind(meta,depth_values,bn_values))
     rm(list=ls(pattern="depth"))
     rm(list=ls(pattern="coords"))
     rm(list=ls(pattern="rb"))
     rm(list=ls(pattern="meta"))
     
     edits<-rawedits%>%select(ID,Site,Species,Phenotype,Tag,Date,LAT,LONG)%>%rename(x=LAT,y=LONG)%>%filter(Date=="B")
     rb_data<-stack(paste0("patch",reef,"_20170930_test3_rb_img"))
     depth<-stack(paste0("patch",reef,"_depth"))
     coords<-filter(edits,Site==paste0(reef))%>%select(x,y)
     meta<-filter(edits,Site==paste0(reef))%>%select(ID,Site,Species,Phenotype,Tag,Date,x,y)
     depth_values<-as.data.frame(raster::extract(depth,coords,method="simple"));names(depth_values)[1]<-"Depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords,method="simple"))
     wavelength_labels<-readRDS("wavelength_labels");colnames(rb_values)[1:52]<-wavelength_labels
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     assign(paste("data",reef,"B",sep="_"),cbind(meta,depth_values,bn_values))
     rm(list=ls(pattern="depth"))
     rm(list=ls(pattern="coords"))
     rm(list=ls(pattern="rb"))
     rm(list=ls(pattern="edits"))
     rm(list=ls(pattern="meta"))
}

newpoints_data<-do.call(rbind, lapply(ls(pattern = "data"), get))
B13<-newpoints_data%>%filter(Site=="13"&Date=="B")
B25<-newpoints_data%>%filter(Site=="25"&Date=="B")
A44<-newpoints_data%>%filter(Site=="44"&Date=="A")
A42<-newpoints_data%>%filter(Site=="42"&Date=="A")

allpoints<-rbind(B13,B25,A44,A42)
svm_newpoints<-allpoints%>%select(-Date,-Site,-x,-y,-Tag,-ID)
svm_metadata<-allpoints%>%select(Date,Site,x,y,Tag,ID)

saveRDS(svm_newpoints,"./output/processed_spectra/phenotype_data_bn")
saveRDS(svm_metadata,"./output/processed_spectra/phenotype_metadata")
######################################### REFLECTANCE WILCOX ################################################################################################################### #####
rm(list = ls())
reef_numbers<-as.numeric(c("13","25","42","44"))
data<-readRDS("./output/processed_spectra/phenotype_data_bn")
data<-data[-117,];set.seed(3839)

sort<-data%>%arrange(desc(Species))%>%
     arrange(desc(Phenotype))%>%
     select(-Depth)%>%
     gather(wavelength,R,-Species,-Phenotype)%>%
     mutate(wavelength=as.factor(wavelength))

prep<-data%>%select(-Depth)%>%
     unite(Sp_Ph,Species,Phenotype,sep="_")%>%
     mutate(Sp_Ph=as.factor(Sp_Ph))%>%
     rename_all(function(x) paste0("A", x))%>%
     rename(Sp_Ph=ASp_Ph)

#write.table(prep,"~/Desktop/data.txt",sep="\t",quote=FALSE,row.names=FALSE)
combins<-combn(levels(prep$Sp_Ph),2);combins<-cbind(combins[,1],combins[,6])
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

list<-as.data.frame(colnames(prep[,2:53]))

model_A421.73<-map(.x=params_list,.f=~wilcox.test(formula=A421.73~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A426.73<-map(.x=params_list,.f=~wilcox.test(formula=A426.73~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A431.74<-map(.x=params_list,.f=~wilcox.test(formula=A431.74~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A436.75<-map(.x=params_list,.f=~wilcox.test(formula=A436.75~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A441.76<-map(.x=params_list,.f=~wilcox.test(formula=A441.76~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A446.76<-map(.x=params_list,.f=~wilcox.test(formula=A446.76~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A451.77<-map(.x=params_list,.f=~wilcox.test(formula=A451.77~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A456.78<-map(.x=params_list,.f=~wilcox.test(formula=A456.78~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A461.79<-map(.x=params_list,.f=~wilcox.test(formula=A461.79~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A466.79<-map(.x=params_list,.f=~wilcox.test(formula=A466.79~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A471.8<-map(.x=params_list,.f=~wilcox.test(formula=A471.8~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A476.81<-map(.x=params_list,.f=~wilcox.test(formula=A476.81~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A481.82<-map(.x=params_list,.f=~wilcox.test(formula=A481.82~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A486.82<-map(.x=params_list,.f=~wilcox.test(formula=A486.82~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A491.83<-map(.x=params_list,.f=~wilcox.test(formula=A491.83~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A496.84<-map(.x=params_list,.f=~wilcox.test(formula=A496.84~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A501.85<-map(.x=params_list,.f=~wilcox.test(formula=A501.85~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A506.85<-map(.x=params_list,.f=~wilcox.test(formula=A506.85~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A511.86<-map(.x=params_list,.f=~wilcox.test(formula=A511.86~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A516.87<-map(.x=params_list,.f=~wilcox.test(formula=A516.87~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A521.88<-map(.x=params_list,.f=~wilcox.test(formula=A521.88~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A526.88<-map(.x=params_list,.f=~wilcox.test(formula=A526.88~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A531.89<-map(.x=params_list,.f=~wilcox.test(formula=A531.89~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A536.9<-map(.x=params_list,.f=~wilcox.test(formula=A536.9~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A541.91<-map(.x=params_list,.f=~wilcox.test(formula=A541.91~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A546.91<-map(.x=params_list,.f=~wilcox.test(formula=A546.91~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A551.92<-map(.x=params_list,.f=~wilcox.test(formula=A551.92~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A556.93<-map(.x=params_list,.f=~wilcox.test(formula=A556.93~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A561.94<-map(.x=params_list,.f=~wilcox.test(formula=A561.94~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A566.94<-map(.x=params_list,.f=~wilcox.test(formula=A566.94~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A571.95<-map(.x=params_list,.f=~wilcox.test(formula=A571.95~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A576.96<-map(.x=params_list,.f=~wilcox.test(formula=A576.96~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A581.97<-map(.x=params_list,.f=~wilcox.test(formula=A581.97~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A586.97<-map(.x=params_list,.f=~wilcox.test(formula=A586.97~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A591.98<-map(.x=params_list,.f=~wilcox.test(formula=A591.98~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A596.99<-map(.x=params_list,.f=~wilcox.test(formula=A596.99~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A602.<-map(.x=params_list,.f=~wilcox.test(formula=A602.~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A607.<-map(.x=params_list,.f=~wilcox.test(formula=A607.~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A612.01<-map(.x=params_list,.f=~wilcox.test(formula=A612.01~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A617.02<-map(.x=params_list,.f=~wilcox.test(formula=A617.02~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A622.03<-map(.x=params_list,.f=~wilcox.test(formula=A622.03~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A627.03<-map(.x=params_list,.f=~wilcox.test(formula=A627.03~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A632.04<-map(.x=params_list,.f=~wilcox.test(formula=A632.04~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A637.05<-map(.x=params_list,.f=~wilcox.test(formula=A637.05~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A642.06<-map(.x=params_list,.f=~wilcox.test(formula=A642.06~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A647.06<-map(.x=params_list,.f=~wilcox.test(formula=A647.06~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A652.07<-map(.x=params_list,.f=~wilcox.test(formula=A652.07~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A657.08<-map(.x=params_list,.f=~wilcox.test(formula=A657.08~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A662.09<-map(.x=params_list,.f=~wilcox.test(formula=A662.09~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A667.09<-map(.x=params_list,.f=~wilcox.test(formula=A667.09~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A672.1<-map(.x=params_list,.f=~wilcox.test(formula=A672.1~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))
model_A677.11<-map(.x=params_list,.f=~wilcox.test(formula=A677.11~Sp_Ph,data=subset(prep,Sp_Ph%in%.x)))


wilcox_pvals <- do.call(cbind,list(t(data.frame(map(.x=model_A421.73,.f="p.value"))),
                                   t(data.frame(map(.x=model_A426.73,.f="p.value"))),
                                   t(data.frame(map(.x=model_A431.74,.f="p.value"))),
                                   t(data.frame(map(.x=model_A436.75,.f="p.value"))),
                                   t(data.frame(map(.x=model_A441.76,.f="p.value"))),
                                   t(data.frame(map(.x=model_A446.76,.f="p.value"))),
                                   t(data.frame(map(.x=model_A451.77,.f="p.value"))),
                                   t(data.frame(map(.x=model_A456.78,.f="p.value"))),
                                   t(data.frame(map(.x=model_A461.79,.f="p.value"))),
                                   t(data.frame(map(.x=model_A466.79,.f="p.value"))),
                                   t(data.frame(map(.x=model_A471.8,.f="p.value"))),
                                   t(data.frame(map(.x=model_A476.81,.f="p.value"))),
                                   t(data.frame(map(.x=model_A481.82,.f="p.value"))),
                                   t(data.frame(map(.x=model_A486.82,.f="p.value"))),
                                   t(data.frame(map(.x=model_A491.83,.f="p.value"))),
                                   t(data.frame(map(.x=model_A496.84,.f="p.value"))),
                                   t(data.frame(map(.x=model_A501.85,.f="p.value"))),
                                   t(data.frame(map(.x=model_A506.85,.f="p.value"))),
                                   t(data.frame(map(.x=model_A511.86,.f="p.value"))),
                                   t(data.frame(map(.x=model_A516.87,.f="p.value"))),
                                   t(data.frame(map(.x=model_A521.88,.f="p.value"))),
                                   t(data.frame(map(.x=model_A526.88,.f="p.value"))),
                                   t(data.frame(map(.x=model_A531.89,.f="p.value"))),
                                   t(data.frame(map(.x=model_A536.9,.f="p.value"))),
                                   t(data.frame(map(.x=model_A541.91,.f="p.value"))),
                                   t(data.frame(map(.x=model_A546.91,.f="p.value"))),
                                   t(data.frame(map(.x=model_A551.92,.f="p.value"))),
                                   t(data.frame(map(.x=model_A556.93,.f="p.value"))),
                                   t(data.frame(map(.x=model_A561.94,.f="p.value"))),
                                   t(data.frame(map(.x=model_A566.94,.f="p.value"))),
                                   t(data.frame(map(.x=model_A571.95,.f="p.value"))),
                                   t(data.frame(map(.x=model_A576.96,.f="p.value"))),
                                   t(data.frame(map(.x=model_A581.97,.f="p.value"))),
                                   t(data.frame(map(.x=model_A586.97,.f="p.value"))),
                                   t(data.frame(map(.x=model_A591.98,.f="p.value"))),
                                   t(data.frame(map(.x=model_A596.99,.f="p.value"))),
                                   t(data.frame(map(.x=model_A602.,.f="p.value"))),
                                   t(data.frame(map(.x=model_A607.,.f="p.value"))),
                                   t(data.frame(map(.x=model_A612.01,.f="p.value"))),
                                   t(data.frame(map(.x=model_A617.02,.f="p.value"))),
                                   t(data.frame(map(.x=model_A622.03,.f="p.value"))),
                                   t(data.frame(map(.x=model_A627.03,.f="p.value"))),
                                   t(data.frame(map(.x=model_A632.04,.f="p.value"))),
                                   t(data.frame(map(.x=model_A637.05,.f="p.value"))),
                                   t(data.frame(map(.x=model_A642.06,.f="p.value"))),
                                   t(data.frame(map(.x=model_A647.06,.f="p.value"))),
                                   t(data.frame(map(.x=model_A652.07,.f="p.value"))),
                                   t(data.frame(map(.x=model_A657.08,.f="p.value"))),
                                   t(data.frame(map(.x=model_A662.09,.f="p.value"))),
                                   t(data.frame(map(.x=model_A667.09,.f="p.value"))),
                                   t(data.frame(map(.x=model_A672.1,.f="p.value"))),
                                   t(data.frame(map(.x=model_A677.11,.f="p.value")))))


row.names(wilcox_pvals) <- unlist(map(.x = params_list, .f = ~ paste0(.x, collapse = "")))

colnames(wilcox_pvals) <- names(prep)[2:53]
output<-as.data.frame(t(wilcox_pvals))%>%
     rownames_to_column()
names(output) <- c("Wavelength","Montipora","Porites")
output$Wavelength <- gsub('A', '', output$Wavelength)

#m_adj<-p.adjust(output$Montipora, method = "bonferroni")
#p_adj<-p.adjust(output$Porites, method = "bonferroni")
#final<-as.data.frame(cbind(output$Wavelength,m_adj,p_adj));names(final) <- c("Wavelength","Montipora","Porites")

write.table(output,"./output/misc/wv_significance.txt",quote=FALSE,row.names = FALSE,sep="\t")
######################################### FIT MONTIPORA PHENOTYPE ############################################################################################################## #####
#setwd("~/nfs_fs02/CAO")
#library(gbm)
#library(caret)
#library(tidyverse)
#library(doMC)
#registerDoMC(cores=1)

#data<-readRDS("./phenotype_data_bn")
#data<-data[-117,]%>%filter(Species=="Montipora")%>%select(-Species)
#grid <-  expand.grid(interaction.depth = seq(1,10,1),
#                     n.trees = seq(5,5000,5),
#                     shrinkage = c(0.001,0.0025,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1),
#                     n.minobsinnode = c(1,15,1))
#nrow(grid)
#set.seed(3839)
#seeds <- vector(mode = "list", length = 11) #nrepeats * nfolds +1
#for(i in 1:11) seeds[[i]]<- sample.int(n=nrow(grid), nrow(grid)) #last value is at least grid length
#seeds[[11]]<-sample.int(1000, 1)

#montipora_gbm<- caret::train(Phenotype ~., data = data, method="gbm",
#                             trControl=trainControl(method="cv",number=5,seeds=seeds),
#                             preProcess=c("center","scale"),
#                             metric="Kappa",
#                             tuneGrid=grid,
#                             bag.fraction=0.75,
#                             verbose=FALSE)
#saveRDS(montipora_gbm,"0_montipora_gbm_model")
######################################### MONTIPORA PHENOTYPE DETAILS ########################################################################################################## #####
rm(list = ls())
model<-readRDS("./output/model_fits/0_montipora_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
data<-raw[-117,]%>%filter(Species=="Montipora")%>%select(-Species)
wavelength_labels<-readRDS("./output/misc/wavelength_labels")

#model
#confusionMatrix(model)
model$bestTune
model$resample
prediction<-predict(model,newdata=data,type=c("raw")); cm<-confusionMatrix(prediction,as.factor(data$Phenotype)); cm
resample<-mean(model$resample$Accuracy);resample
sdresample<-sd(model$resample$Accuracy);sdresample
kresample<-mean(model$resample$Kappa);kresample
list<-varImp(model);list$importance


raw<-readRDS("./output/processed_spectra/phenotype_data_bn")%>%bind_cols(readRDS("./output/processed_spectra/phenotype_metadata"))
data<-raw[-117,]%>%filter(Species=="Montipora")%>%select(Species,Phenotype,Site,Depth,x,y)%>%bind_cols(.,prediction)%>%rename(Prediction=7)
write.table(data,"~/Desktop/m_pheno.txt",sep="\t",quote=FALSE,row.names=FALSE)

######################################### FIT PORITES PHENOTYPE ################################################################################################################ #####
#setwd("~/nfs_fs02/CAO")
#library(gbm)
#library(caret)
#library(tidyverse)
#library(doMC)
#registerDoMC(cores=5)

#data<-readRDS("./phenotype_data_bn")
#data<-data%>%filter(Species=="Porites")%>%select(-Species)
#grid <-  expand.grid(interaction.depth = seq(1,10,1),
#                     n.trees = seq(5,5000,5),
#                     shrinkage = c(0.05,0.06,0.07,0.08,0.085,0.0875,0.09,0.0925,0.095,0.1),
#                     n.minobsinnode = c(5,20,1))
#nrow(grid)
#set.seed(3839)
#seeds <- vector(mode = "list", length = 11) #nrepeats * nfolds +1
#for(i in 1:11) seeds[[i]]<- sample.int(n=nrow(grid), nrow(grid)) #last value is at least grid length
#seeds[[11]]<-sample.int(1000, 1)

#porites_gbm<- caret::train(Phenotype ~., data = data, method="gbm",
#                           trControl=trainControl(method="cv",number=5,seeds=seeds),
#                           preProcess=c("center","scale"),
#                           metric="Kappa",
#                           tuneGrid=grid,
#                           bag.fraction=0.75,
#                           verbose=FALSE)
#saveRDS(porites_gbm,"0_porites_gbm_model")
######################################### PORITES PHENOTYPE DETAILS ############################################################################################################ #####
model<-readRDS("./output/model_fits/0_porites_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
data<-raw[-117,]%>%filter(Species=="Porites")%>%select(-Species)

#model
#confusionMatrix(model)
model$bestTune
model$resample
prediction<-predict(model,newdata=data,type=c("raw")); cm<-confusionMatrix(prediction,as.factor(data$Phenotype)); cm
resample<-mean(model$resample$Accuracy);resample
sdresample<-sd(model$resample$Accuracy);sdresample
kresample<-mean(model$resample$Kappa);kresample
list<-varImp(model);list$importance

raw<-readRDS("./output/processed_spectra/phenotype_data_bn")%>%bind_cols(readRDS("./output/processed_spectra/phenotype_metadata"))
data<-raw[-117,]%>%filter(Species=="Porites")%>%select(Species,Phenotype,Site,Depth,x,y)%>%bind_cols(.,prediction)%>%rename(Prediction=7)
write.table(data,"~/Desktop/p_pheno.txt",sep="\t",quote=FALSE,row.names=FALSE)

#### TABLE 2
cm

######################################### SEPARABILITY INDICES ################################################################################################################# #####
rm(list = ls())
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
metadata<-readRDS("./output/processed_spectra/phenotype_metadata")
data<-bind_cols(raw,metadata)[-117,]

fulldata<-data%>%select(-Date,-Site,-x,-y,-Tag,-ID,-Depth)
sep_mean_montipora<-fulldata%>%group_by(Species,Phenotype)%>%summarise_all(funs(mean))%>%filter(Species=="Montipora")
sep_SD_montipora<-fulldata%>%group_by(Species,Phenotype)%>%summarise_all(funs(sd))%>%filter(Species=="Montipora")
sep_mean_porites<-fulldata%>%group_by(Species,Phenotype)%>%summarise_all(funs(mean))%>%filter(Species=="Porites")
sep_SD_porites<-fulldata%>%group_by(Species,Phenotype)%>%summarise_all(funs(sd))%>%filter(Species=="Porites")

sep_index_montipora<-as.data.frame(t((1.96*(sep_SD_montipora[1,3:52]+sep_SD_montipora[2,3:52]))/abs(sep_mean_montipora[1,3:52]-sep_mean_montipora[2,3:52])))
sep_index_porites<-as.data.frame(t((1.96*(sep_SD_porites[1,3:52]+sep_SD_porites[2,3:52]))/abs(sep_mean_porites[1,3:52]-sep_mean_porites[2,3:52])))

figdata_montipora<-rownames_to_column(sep_index_montipora)%>%mutate(Species="Montipora");colnames(figdata_montipora)[1]<-"Wavelength";colnames(figdata_montipora)[2]<-"Separability"
figdata_porites<-rownames_to_column(sep_index_porites)%>%mutate(Species="Porites");colnames(figdata_porites)[1]<-"Wavelength";colnames(figdata_porites)[2]<-"Separability"
figdata<-rbind(figdata_montipora,figdata_porites)

sep<-ggplot(figdata)+
     geom_point(aes(as.numeric(Wavelength),Separability,color=Species,group=Species))+
     geom_line(aes(as.numeric(Wavelength),Separability,color=Species,group=Species))+
     theme_classic(base_size = 8)+
     scale_color_manual(values=c("blue","green2"))+
     xlab("Wavelength (nm)")+
     ylab("Phenotype Separability")+
     scale_x_continuous(breaks=seq(425,675,50))+
     theme(legend.position = c(0.5, .95),
           legend.text=element_text(size=8),
           legend.title=element_blank(),
           legend.key.size = unit(0.25,"cm"))+
     guides(color = guide_legend(nrow = 1));sep

######################################### PHENOTYPE LAYERS ##################################################################################################################### ##### 
library(sf);library(spdep);library(raster);library(rgdal)
rm(list = ls())
model<-readRDS("./output/model_fits/0_montipora_gbm_model")

reef_numbers<-as.numeric(c("13","25","42","44"))
for (reef in reef_numbers){ #for montipora
     depth<-stack(paste0("patch",reef,"_depth"))
     rb_data<-stack(paste0("patch",reef,"_20170930_test3_rb_img"))
     extent<-extent(rb_data)
     x<-seq(extent@xmin+.2,extent@xmax+.2,.4)
     y<-seq(extent@ymin+.2,extent@ymax+.2,.4)
     total<-expand.grid(x,y);colnames(total)[1]<-"x";colnames(total)[2]<-"y"
     coords<-total%>%dplyr::select(x,y)%>%rownames_to_column()%>%rename(index=rowname)
     depth_values<-as.data.frame(raster::extract(depth,coords[,2:3],method="simple"));names(depth_values)[1]<-"Depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords[,2:3],method="simple"))
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     wavelength_labels<-readRDS("./output/misc/wavelength_labels");colnames(bn_values)[1:52]<-wavelength_labels
     merge<-cbind(coords[1],depth_values,bn_values)                                           
     merge<-merge[complete.cases(merge), ]
     species<-readRDS(paste("./output/substrate_prediction",reef,sep="_"))
     phenotype_full<-inner_join(species,merge,by="index"); colnames(phenotype_full)[4]<-"Species"
     
     phenotype_trim<-phenotype_full%>%filter(Species=="Montipora")%>%filter(Depth<=5)
     prediction<-as.data.frame(predict(model,newdata=phenotype_trim,type=c("prob")))
     pre_output<-cbind(phenotype_trim[,1:5],as.data.frame(prediction))
     output<-pre_output%>%dplyr::filter(nonbleached>=.5)                                       #filter to only nonbleached
     saveRDS(pre_output,paste("./output/phenotype_prediction_montipora_all",reef,sep="_"))
     coordinates(pre_output)=~x+y
     proj4string(pre_output)=CRS("+proj=utm +zone=4N +datum=WGS84")
     raster::shapefile(pre_output,paste("./output/moran/montipora_moran",reef,sep="_"),overwrite=TRUE) #output of everything for moran test
     saveRDS(output,paste("./output/phenotype_prediction_montipora",reef,sep="_"))             #output of only nonbleached predictions for maps
     coordinates(output)=~x+y
     proj4string(output)=CRS("+proj=utm +zone=4N +datum=WGS84")
     raster::shapefile(output,paste("./output/phenotype_prediction_montipora",reef,sep="_"),overwrite=TRUE)
}

model<-readRDS("./output/model_fits/0_porites_gbm_model")

reef_numbers<-as.numeric(c("13","25","42","44"))
for (reef in reef_numbers){ #for porites
     depth<-stack(paste0("patch",reef,"_depth"))
     rb_data<-stack(paste0("patch",reef,"_20170930_test3_rb_img"))
     extent<-extent(rb_data)
     x<-seq(extent@xmin+.2,extent@xmax+.2,.4)
     y<-seq(extent@ymin+.2,extent@ymax+.2,.4)
     total<-expand.grid(x,y);colnames(total)[1]<-"x";colnames(total)[2]<-"y"
     coords<-total%>%dplyr::select(x,y)%>%rownames_to_column()%>%rename(index=rowname)
     depth_values<-as.data.frame(raster::extract(depth,coords[,2:3],method="simple"));names(depth_values)[1]<-"Depth"
     rb_values<-as.data.frame(raster::extract(rb_data,coords[,2:3],method="simple"))
     sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
     bn_values<-rb_values/sums$sumsq
     wavelength_labels<-readRDS("./output/misc/wavelength_labels");colnames(bn_values)[1:52]<-wavelength_labels
     merge<-cbind(coords[1],depth_values,bn_values)                                           
     merge<-merge[complete.cases(merge), ]
     species<-readRDS(paste("./output/substrate_prediction",reef,sep="_"))
     phenotype_full<-inner_join(species,merge,by="index"); colnames(phenotype_full)[4]<-"Species"
     
     #porites output
     phenotype_trim<-phenotype_full%>%filter(Species=="Porites")%>%filter(Depth<=5)
     prediction<-as.data.frame(predict(model,newdata=phenotype_trim,type=c("prob")))
     pre_output<-cbind(phenotype_trim[,1:5],as.data.frame(prediction))
     output<-pre_output%>%filter(nonbleached>=.5)                                              #filter to only nonbleached
     saveRDS(pre_output,paste("./output/phenotype_prediction_montipora_all",reef,sep="_"))
     coordinates(pre_output)=~x+y
     proj4string(pre_output)=CRS("+proj=utm +zone=4N +datum=WGS84")
     raster::shapefile(pre_output,paste("./output/moran/porites_moran",reef,sep="_"),overwrite=TRUE) #output of everything for moran test
     saveRDS(output,paste("./output/phenotype_prediction_montipora",reef,sep="_"))             #output of only nonbleached predictions for maps
     coordinates(output)=~x+y
     proj4string(output)=CRS("+proj=utm +zone=4N +datum=WGS84")
     raster::shapefile(output,paste("./output/phenotype_prediction_porites",reef,sep="_"),overwrite=TRUE)
}
detach(package:spdep);detach(package:sf);detach(package:raster);detach(package:rgdal)


######################################### F2 - MONTIPORA SPECTRA *FIGURE*  ##################################################################################################### #####
#rm(list = ls())
mmodel<-readRDS("./output/model_fits/0_montipora_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
mdata<-raw[-117,]%>%filter(Species=="Montipora")%>%select(-Species)
wavelength_labels<-readRDS("./output/misc/wavelength_labels")
mlist<-varImp(mmodel,scale=FALSE);#mlist$importance
mbar<-mlist$importance%>%rownames_to_column()%>%filter(rowname!="Depth")%>%select(-rowname)%>%bind_cols(as.data.frame(wavelength_labels))%>%mutate(wavelength_labels=as.numeric(as.character(wavelength_labels)))%>%
     mutate(y=0.01)

mplotdata<-mdata%>%rownames_to_column()%>%rename(sample=rowname)%>%select(-Depth)%>%gather(wavelength,reflectance,-Phenotype,-sample)%>%group_by(sample)%>%mutate(order=sample(1:100000,1,replace=F))%>%arrange(desc(order))
mdata%>%group_by(Phenotype)%>%tally()
mplotdata$Phenotype <- factor(mplotdata$Phenotype,levels = c("bleached","nonbleached"),labels=c("Bleached (n=46)","Nonbleached (n=74)"))

mspec<-ggplot(mplotdata)+
     #annotate("rect",xmin=461,xmax=516,ymin=0,ymax=0.23,fill="lightgray")+
     #annotate("rect",xmin=561,xmax=617,ymin=0,ymax=0.23,fill="lightgray")+
     #annotate("rect",xmin=652,xmax=677,ymin=0,ymax=0.23,fill="lightgray")+
     geom_path(aes(as.numeric(wavelength),reflectance,color=Phenotype,group=sample),size=0.25)+
     theme_classic(base_size = 8)+
     scale_x_continuous(breaks=seq(425,675,50))+
     ylab("Visually Healthy\nM.capitata R")+
     theme(legend.position = c(0.2,0.8),
           legend.text=element_text(size=6),
           legend.key.size = unit(0.15,"cm"),
           legend.title=element_text(size=6),
           axis.title=element_text(size=6),,
           axis.text.y=element_text(angle=90,hjust=0.5),
           axis.text.x=element_blank(),
           axis.title.x=element_blank())+
     scale_colour_manual(values=c("orange","blue"),name=c("Historical Phenotype"))+
     #scale_fill_manual(values=c("orange","blue"))+
     scale_y_continuous(lim=c(0,.23),breaks=seq(0,.2,0.1))+
     #new_scale_color()+
     #geom_path(aes(wavelength_labels,y,color=Overall),size=2,data=bar,lineend="square")+
     #scale_color_gradient(low="gray", high="purple")+
     #scale_color_viridis_c()+
     annotate("line", x = 461:516, y=0,size=0.5)+
     annotate("line", x = 561:617, y=0,size=0.5)+
     annotate("line", x = 652:677, y=0,size=0.5)+
     #guides(color=FALSE)+
     annotate("text",x=510,y=0.02,label="GBM Best-Fit Accuracy: 100%, CV-Accuracy: 69.9%",size=2);mspec

######################################### F2 - PORITES SPECTRA *FIGURE* ######################################################################################################## #####
pmodel<-readRDS("./output/model_fits/0_porites_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
pdata<-raw[-117,]%>%filter(Species=="Porites")%>%select(-Species)

plist<-varImp(pmodel,scale=FALSE);plist$importance
pbar<-plist$importance%>%rownames_to_column()%>%filter(rowname!="Depth")%>%select(-rowname)%>%bind_cols(as.data.frame(wavelength_labels))%>%mutate(wavelength_labels=as.numeric(as.character(wavelength_labels)))%>%
     mutate(y=0.01)

pplotdata<-pdata%>%rownames_to_column()%>%rename(sample=rowname)%>%select(-Depth)%>%gather(wavelength,reflectance,-Phenotype,-sample)%>%group_by(sample)%>%mutate(order=sample(1:1000,1,replace=T))%>%arrange(desc(order))
data%>%group_by(Phenotype)%>%tally()
pplotdata$Phenotype <- factor(pplotdata$Phenotype,levels = c("bleached","nonbleached"),labels=c("Bleached (n=93)","Nonbleached (n=61)"))


pspec<-ggplot(pplotdata)+
     geom_path(aes(as.numeric(wavelength),reflectance,color=Phenotype,group=sample),size=0.25)+
     theme_classic(base_size = 8)+
     scale_x_continuous(breaks=seq(425,675,50))+
     xlab("Wavelength")+
     #ylab(expression(atop(NA,atop("Visually Healthy", paste(italic("P. compressa "),"R")))))+
     ylab("Visually Healthy\nP.compressa R")+
     theme(legend.position = c(0.2,0.8),
           legend.text=element_text(size=6),
           legend.key.size = unit(0.15,"cm"),
           legend.title=element_text(size=6),
           axis.title=element_text(size=6),
           axis.text.y=element_text(angle=90,hjust=0.5))+
     scale_colour_manual(values=c("orange","blue"),name=c("Historical Phenotype"))+
     scale_y_continuous(lim=c(0,.23),breaks=seq(0,.2,0.1))+
     #new_scale_color()+
     #geom_path(aes(wavelength_labels,y,color=Overall),size=2,data=bar,lineend="square")+
     #scale_color_gradient(low="gray", high="purple")+
     #scale_color_viridis_c()+
     #guides(color=FALSE)+
     annotate("text",x=510,y=0.02,label="GBM Best-Fit Accuracy: 94.8%, CV-Accuracy: 64.3%",size=2);pspec
######################################### F2 - ACCURACY DETAILS ################################################################################################################ #####
#rm(list = ls())
model<-readRDS("./output/model_fits/0_montipora_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
metadata<-readRDS("./output/processed_spectra/phenotype_metadata")
data<-bind_cols(raw,metadata)[-117,]%>%filter(Species=="Montipora")%>%select(-Species)
prediction<-as.data.frame(predict(model,newdata=data[,2:54],type=c("raw")))%>%rename(prediction=1); cm<-confusionMatrix(prediction$prediction,as.factor(data$Phenotype)); cm

m<-bind_cols(data%>%select(Site,Phenotype),prediction)%>%mutate(Correct=case_when(Phenotype==prediction~1,Phenotype!=prediction~0))%>%
     group_by(Site)%>%
     summarise(Mean = mean(Correct, na.rm=TRUE))%>%mutate(Species="Montipora")

model<-readRDS("./output/model_fits/0_porites_gbm_model")
raw<-readRDS("./output/processed_spectra/phenotype_data_bn")
metadata<-readRDS("./output/processed_spectra/phenotype_metadata")
data<-bind_cols(raw,metadata)[-117,]%>%filter(Species=="Porites")%>%select(-Species)
prediction<-as.data.frame(predict(model,newdata=data[,2:54],type=c("raw")))%>%rename(prediction=1); cm<-confusionMatrix(prediction$prediction,as.factor(data$Phenotype)); cm

p<-bind_cols(data%>%select(Site,Phenotype),prediction)%>%mutate(Correct=case_when(Phenotype==prediction~1,Phenotype!=prediction~0))%>%
     group_by(Site)%>%
     summarise(Mean = mean(Correct, na.rm=TRUE))%>%mutate(Species="Porites")

plot<-bind_rows(m,p)%>%mutate(Site=as.factor(Site))
accuracy<-ggplot(plot)+geom_bar(aes(Site,Mean,fill=Species),stat="Identity",position="dodge")+
     ylab("Best-Fit Accuracy")+
     scale_fill_manual(values=c("blue","green2"),
                       guide=guide_legend(nrow=2,
                                        label.theme = element_text(angle = 90,size=6,face="italic"),
                                        label.position="top",
                                        label.hjust = 0.5,
                                        label.vjust = 0.5,),labels=c("M.capitata","P.compressa"))+
     scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
     theme_classic(base_size=8)+
     theme( legend.key.size = unit(0.2,"cm"),
            legend.position="right",
            legend.spacing.x=unit(0.1,"cm"),
            legend.spacing.y=unit(0.05,"cm"),
            legend.title=element_blank(),
            legend.box.margin=margin(0,0,0,-8),
            axis.text.y=element_text(angle=90,hjust=0.5),
            axis.title=element_text(size=6))+
     xlab("Reef");accuracy

######################################### F2 - COVER STATISTICS ################################################################################################################ #####
library(sf)
#r13<-st_read("./output/substrate_prediction_13.shp")%>% st_set_geometry(NULL)%>%filter(prediction!="Water")%>%group_by(prediction)%>%tally()%>%mutate(Reef='13')
#r25<-st_read("./output/substrate_prediction_25.shp")%>% st_set_geometry(NULL)%>%filter(prediction!="Water")%>%group_by(prediction)%>%tally()%>%mutate(Reef='25')
#r42<-st_read("./output/substrate_prediction_42.shp")%>% st_set_geometry(NULL)%>%filter(prediction!="Water")%>%group_by(prediction)%>%tally()%>%mutate(Reef='42')
#r44<-st_read("./output/substrate_prediction_44.shp")%>% st_set_geometry(NULL)%>%filter(prediction!="Water")%>%group_by(prediction)%>%tally()%>%mutate(Reef='44')

#m13<-st_read("./output/phenotype_prediction_montipora_13.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='13',Species="Montipora")
#m25<-st_read("./output/phenotype_prediction_montipora_25.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='25',Species="Montipora")
#m42<-st_read("./output/phenotype_prediction_montipora_42.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='42',Species="Montipora")
#m44<-st_read("./output/phenotype_prediction_montipora_44.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='44',Species="Montipora")
#p13<-st_read("./output/phenotype_prediction_porites_13.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='13',Species="Porites")
#p25<-st_read("./output/phenotype_prediction_porites_25.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='25',Species="Porites")
#p42<-st_read("./output/phenotype_prediction_porites_42.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='42',Species="Porites")
#p44<-st_read("./output/phenotype_prediction_porites_44.shp")%>% st_set_geometry(NULL)%>%tally()%>%mutate(Reef='44',Species="Porites")
#x<-rbind(rbind(r13,r25,r42,r44)%>%rename(substrate=1)%>%mutate(type="all"), rbind(m13,m25,m42,m44,p13,p25,p42,p44)%>%rename(substrate=3)%>%mutate(type="non"))%>%rename(reef=Reef)
saveRDS(x,"./output/cover_summary")
x<-readRDS("./output/cover_summary")
#proportion by site
out<-x%>%spread(type,n)%>%mutate(bl=all-non)%>%gather(type,n,-reef,-substrate)%>%filter(n!="NA")%>%filter(type!="all"|(substrate!="Montipora"&substrate!="Porites"))%>%
     group_by(reef)%>%
     mutate(substrate=case_when(substrate=="Dead"~"Hardbottom",TRUE ~as.character(substrate)))%>%
     mutate(type=case_when(type=="all"~"",
                           type=="non"~"NB",
                           type=="bl"~"B"))%>%
     mutate(substrate=case_when(substrate=="Montipora"~"M.capitata",
                                substrate=="Porites"~"P.compressa",
                                TRUE~as.character(substrate)))%>%
     unite(group,substrate,type,sep=" ")%>%
     group_by(reef)%>%mutate(reef_sum=sum(n))%>%
     mutate(prop=n/reef_sum)

out$group<-factor(out$group,levels=c("Sand ","Hardbottom ","M.capitata B","M.capitata NB","P.compressa B","P.compressa NB"),
                  labels=c("Sand","Hardbottom","M.capitata (B)","M.capitata (NB)","P.compressa (B)","P.compressa (NB)"))

cover<-ggplot(out)+geom_bar(aes(reef,prop,fill=group),stat="identity")+
     theme_classic(base_size=8)+
     ylab("Predicted Proportion of Reef Area")+xlab("Reef")+
     scale_fill_manual(values=c("tan","dimgray","#ccccff","blue","#ccffcc","green"),name="",
                       guide=guide_legend(nrow=3,byrow=TRUE,
                                          label.theme = element_text(angle = 90,size=6),
                                          label.position="top",
                                          label.hjust = 0,
                                          label.vjust = 0.5))+
     scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
     theme(legend.key.size=unit(0.2,"cm"),
           legend.box.margin=margin(0,0,0,-8),
           legend.spacing.x=unit(0.1,"cm"),
           legend.spacing.y=unit(0.05,"cm"),
           axis.text.y=element_text(angle=90,hjust=0.5),
           axis.title=element_text(size=6));cover

library(cowplot)     
quartz(w=5.5,h=2.2)
plots<-align_plots(pspec,accuracy,cover,align="h",axis="b")
plot_grid(plot_grid(mspec,plots[[1]],ncol=1,rel_heights=c(0.86,1),labels=c("(a)","(b)"),label_size=8),plots[[2]],plots[[3]],nrow=1,rel_widths=c(3,1,1.05),labels=c("","(c)","(d)"),label_size=8,label_x=c(0,-0.08,-0.08))

#coral cover by site
out%>%filter(group=="Hardbottom"|group=="Sand")%>%group_by(reef)%>%mutate(prop_sum_by_reef=1-sum(prop))%>%select(reef,prop_sum_by_reef)%>%distinct()

#proportion reef total area that is thermally tolerant biomass
x%>%spread(type,n)%>%mutate(bl=all-non)%>%gather(type,n,-reef,-substrate)%>%filter(n!="NA")%>%filter(type!="all"|(substrate!="Montipora"&substrate!="Porites"))%>%
     group_by(reef)%>%
     mutate(substrate=case_when(substrate=="Dead"~"Hardbottom",TRUE ~as.character(substrate)))%>%
     mutate(type=case_when(type=="all"~"",
                           type=="non"~"Nonbleached",
                           type=="bl"~"Bleached"))%>%
     group_by(reef)%>%mutate(reef_sum=sum(n))%>%
     filter(type=="Nonbleached")%>%group_by(reef)%>%mutate(pheno_sum=sum(n))%>%
     mutate(prop=pheno_sum/reef_sum)%>%
     select(reef,prop)%>%distinct()

#proportion of corals (only) that are thermally tolerant by reef
x%>%spread(type,n)%>%mutate(bl=all-non)%>%gather(type,n,-reef,-substrate)%>%filter(n!="NA")%>%filter(type!="all"|(substrate!="Montipora"&substrate!="Porites"))%>%
     group_by(reef)%>%
     mutate(substrate=case_when(substrate=="Dead"~"Hardbottom",TRUE ~as.character(substrate)))%>%
     mutate(type=case_when(type=="all"~"",
                           type=="non"~"Nonbleached",
                           type=="bl"~"Bleached"))%>%
     filter(substrate!="Sand",substrate!="Hardbottom")%>%
     group_by(reef)%>%mutate(reef_sum=sum(n))%>%
     filter(type=="Nonbleached")%>%
     group_by(reef)%>%
     mutate(pheno_sum=sum(n))%>%
     mutate(prop=pheno_sum/reef_sum)%>%select(reef,prop)%>%distinct()
     
#proportion of each species that is thermally tolerant at each reef
x%>%spread(type,n)%>%mutate(bl=all-non)%>%gather(type,n,-reef,-substrate)%>%filter(n!="NA")%>%filter(type!="all"|(substrate!="Montipora"&substrate!="Porites"))%>%
     group_by(reef)%>%
     mutate(substrate=case_when(substrate=="Dead"~"Hardbottom",TRUE ~as.character(substrate)))%>%
     mutate(type=case_when(type=="all"~"",
                           type=="non"~"Nonbleached",
                           type=="bl"~"Bleached"))%>%
     filter(substrate!="Sand",substrate!="Hardbottom")%>%
     group_by(reef,substrate)%>%mutate(reef_sp_sum=sum(n))%>%
     filter(type=="Nonbleached")%>%
     mutate(prop=n/reef_sp_sum)

write.table(out,"~/Desktop/area.txt",sep="\t",row.names=FALSE,quote=FALSE)
######################################### SPATIAL STATISTICS ################################################################################################################### ##### 
species<-c("porites","montipora")
reef_numbers<-as.numeric(c("25"))
#run on HPC

#result <- data.frame(matrix(nrow = 8, ncol = 8))
#colnames(result) <- c("reef", "species","moran_1","p_value_1","moran_5","p_value_5","moran_10","p_value_10")
#result$reef<-c(13,13,25,25,42,42,44,44)
#result$species<-rep(c("porites","montipora"),4)
#for (i in reef_numbers){
#     for (j in species){
#          data<-as(st_read(paste("./output/",paste(j,"moran",i,sep="_"),".shp",sep="")),"Spatial")
#          coord<-coordinates(data)
#          dist<-dnearneigh(coord,0,1) #within 1m  
#          weight<-nb2listw(dist, style="W",zero.policy=T) 
#          out<-moran.mc(data$nnblchd,weight,nsim=2,zero.policy=T)
#          index <- which(result[,1] == i & result[,2] == j)
#          result[index,3]<-out$statistic
#          result[index,4]<-out$p.value
#          dist<-dnearneigh(coord,0,5) #within 5m  
#          weight<-nb2listw(dist, style="W",zero.policy=T) 
#          out<-moran.mc(data$nnblchd,weight,nsim=2,zero.policy=T)
#          index <- which(result[,1] == i & result[,2] == j)
#          result[index,5]<-out$statistic
#          result[index,6]<-out$p.value
#          dist<-dnearneigh(coord,0,10) #within 10m  
#          weight<-nb2listw(dist, style="W",zero.policy=T) 
#          out<-moran.mc(data$nnblchd,weight,nsim=2,zero.policy=T)
#          index <- which(result[,1] == i & result[,2] == j)
#          result[index,7]<-out$statistic
#          result[index,8]<-out$p.value
#     }
#}

moran<-bind_rows(readRDS("./output/montipora_moran_results")%>%filter(species=="montipora"),readRDS("./output/porites_moran_results")%>%filter(species=="porites"))
######################################### FORMAT SPECTRA FOR VALIDATION ######################################################################################################## #####
library(sf);library(spdep);library(raster);library(rgdal)
rm(list = ls())
reef<-13

rawedits<-read_excel("Field_Scores.xlsx",sheet="final")
edits<-rawedits%>%dplyr::select(id,lat,long,field_assessment,score)%>%rename(x=lat,y=long)
rb_data<-stack(paste0("patch",reef,"_20170905_test3_rb_img"))
depth<-stack(paste0("patch",reef,"_depth"))
coords<-edits%>%dplyr::select(x,y)
meta<-edits%>%dplyr::select(id,field_assessment,score)
depth_values<-as.data.frame(raster::extract(depth,coords,method="simple"));names(depth_values)[1]<-"Depth"
rb_values<-as.data.frame(raster::extract(rb_data,coords,method="simple"))
wavelength_labels<-readRDS("./output/misc/wavelength_labels");colnames(rb_values)[1:52]<-wavelength_labels
sums<-as.data.frame(sqrt(rowSums(rb_values[,1:52]^2)));colnames(sums)=c('sumsq')                                                 
bn_values<-rb_values/sums$sumsq
assign(paste("validation",reef,sep="_"),cbind(meta,depth_values,bn_values))
spectra<-validation_13%>%dplyr::select(-id,-field_assessment,-score)
meta<-validation_13%>%dplyr::select(id,field_assessment,score)
saveRDS(spectra,"./output/processed_spectra/validation_spectra")
saveRDS(meta,"./output/processed_spectra/validation_metadata")
detach(package:spdep);detach(package:sf);detach(package:raster);detach(package:rgdal)

######################################### F3 - PHENOTYPE MAPS ################################################################################################################## #####
library(sf);library(ggsn)

raw_mcap<-st_intersection(st_read("./output/phenotype_prediction_montipora_13.shp"),st_read("./points/clip.shp"))
shp_mcap<-fortify(as.data.frame(cbind(raw_mcap,st_coordinates(raw_mcap))))%>%#mutate(breaks=cut(nnblchd,breaks=c(0.5,0.9,0.99,0.999,0.9999,1),labels=c("0.5","0.9","0.99","0.999","0.999")))
     mutate(breaks=ntile(nnblchd,6))
table(shp_mcap$breaks)
as.vector(shp_mcap%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))

validation_points<-read_csv("./points/validation.csv")
metabolite_points<-read_xlsx("./metabolites/metabolites_coordinates.xlsx")
mcap<-ggplot()+
     geom_raster(data = shp_mcap, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(624650,624980),ylim=c(2372480,2372760))+
     theme_classic(base_size=8)+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           axis.line=element_blank(),
           axis.ticks=element_blank(),
           legend.position=c(0.8,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Prediction Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5))+
     geom_point(aes(x=lat,y=long),data=validation_points,shape=21,color="orange",size=1)+
     #geom_point(aes(x=LAT,y=LONG),data=metabolite_points,shape=21,color="blue",size=1)+
     scale_fill_viridis_c(direction=-1,name="Prediction\nNonbleached\n(2017)",labels=c("Moderate","","","","","High (~1)"))+
     annotate("text",x=624690,y=2372473,label="Validation Points (n=54)",size=2,color="orange")+
     annotate("segment",x=624700,y=2372483,xend=624712.4,yend=2372531,color="orange")+
     #annotate("text",x=624820,y=2372473,label="Metabolite Points (n=17)",size=2,color="blue")+
     #annotate("segment",x=624795,y=2372483,xend=624743,yend=2372505,color="blue")+
     annotate("segment",x=625000,y=2372710,xend=624712.4,yend=2372531,color="orange",linetype="dashed")+
     annotate("segment",x=625000,y=2372500,xend=624712.4,yend=2372531,color="orange",linetype="dashed")+
     annotate("text",label="Reef 13 - M.capitata",x=624655,y=2372770,size=2,color="black",hjust=0);mcap


raw_pcomp<-st_intersection(st_read("./output/phenotype_prediction_porites_13.shp"),st_read("./points/clip.shp"))
shp_pcomp<-fortify(as.data.frame(cbind(raw_pcomp,st_coordinates(raw_pcomp))))%>%
     mutate(breaks=ntile(nnblchd,6))
as.vector(shp_pcomp%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))

pcomp<-ggplot()+
     geom_raster(data = shp_pcomp, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(624650,624980),ylim=c(2372480,2372760))+
     theme_classic(base_size=8)+
     scale_fill_viridis_c(direction=-1,labels=c("Moderate","","","","","High (0.82+)"))+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           axis.line=element_blank(),
           axis.ticks=element_blank(),
           legend.position=c(0.8,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Predicted Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5))+
     scalebar(x.min=624650,x.max=624980,y.min=2372480,y.max=2372760,transform=FALSE,dist=25,dist_unit="m",
              height = .02, st.dist = 0.03,
              location="bottomleft",
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2)+
     annotate("text",label="Reef 13 - P.compressa",x=624655,y=2372770,size=2,color="black",hjust=0);pcomp

######################################### F3 - FIELD VALIDATION ################################################################################################################ ##### 
#rm(list = ls())
newdata<-readRDS("./output/processed_spectra/validation_spectra")
model<-readRDS("./output/model_fits/0_montipora_gbm_model")
meta<-readRDS("./output/processed_spectra/validation_metadata")
pheno_prediction<-predict(model,newdata=newdata,type="prob")

validation_data<-bind_cols(meta,pheno_prediction)%>%mutate(prediction=case_when(bleached>0.5 ~"Bleached",
                                                                                bleached<=0.5 ~"Nonbleached"))%>%
     filter(field_assessment!="skip")%>%mutate(score=as.numeric(score))%>%
     mutate(match=case_when(prediction==field_assessment~1,prediction!=field_assessment~0))%>%
     mutate(field_assessment=as.factor(field_assessment),prediction=as.factor(prediction))

cm<-confusionMatrix(validation_data$prediction,validation_data$field_assessment);cm
cm$byClass[3] #bleached predicted accuracy
cm$byClass[4] #nonbleacehd predicted accuracy
heats<-as.data.frame(cm$table)%>%mutate(sum=sum(Freq),prop=Freq/sum)
heats$Reference<-factor(heats$Reference,levels=c("Bleached","Nonbleached"),labels=c("B","NB"))
heats$Prediction<-factor(heats$Prediction,levels=c("Bleached","Nonbleached"),labels=c("B","NB"))

a<-ggplot(heats)+geom_tile(aes(Prediction,Reference,fill=Freq))+
     geom_text(aes(Prediction,Reference,label=Freq),size=2)+
     scale_fill_gradient(low="white",high="turquoise4")+
     #coord_fixed(ratio=1)+
     theme_classic(base_size=6)+
     theme(legend.position="none",
           plot.title=element_text(size=6))+
     ggtitle("Validation (2019)");a

wilcox.test(score~prediction,data=validation_data)
validation_data$prediction<-factor(validation_data$prediction,levels=c("Bleached","Nonbleached"),labels=c("B","NB"))

b<-ggplot(validation_data,aes(y=score,x=prediction,color=prediction,fill=prediction))+
     geom_boxplot(fill="white",width=0.25)+
     geom_dotplot(binaxis="y",stackdir="center",dotsize=.75)+
     theme_classic(base_size=6)+
     scale_y_continuous(limits=c(0,5),breaks=seq(1,5))+
     scale_fill_manual(values=c("gray","turquoise4"),name="")+
     scale_color_manual(values=c("gray","turquoise4"),name="")+
     theme(legend.position="none",
           axis.text.y=element_text(angle=90,hjust=0.5))+
     ylab("Visual\nBleaching Score")+xlab("Prediction")+
     annotate("text", x = 0.5, y=0,size=2,label="Wilcoxon p<0.001",fontface='italic',hjust=0);b


quartz(w=7.2,h=2.6)
plot_grid(pcomp,mcap,plot_grid(a,b,ncol=1,align="v",axis="lr",rel_heights=c(2,3),labels=c("(c)","(d)"),label_size=8),nrow=1,rel_widths=c(3,3,1.1),labels=c("(a)","(b)","",""),label_size=8)

######################################### IN SITU ANALYSIS ##################################################################################################################### #####
rm(list = ls())
rawdata<-read_excel("1_FieldSpectra.xlsx",sheet="Working Phenotype")%>%clean_names()
data<-rawdata%>%filter(phenotype!="NA")%>%select(species,site,tag,phenotype,contains('x'))#%>%group_by(tag,site)%>%summarise_at(vars(contains('x')),mean)%>%left_join(.,(rawdata%>%select(tag,species,site,phenotype)%>%distinct()),by=c("tag","site"))%>%select(tag,species,site,phenotype,everything())
data%>%ungroup()%>%filter(species=="Montipora")%>%select(tag)%>%distinct()%>%tally() #n colonies
data%>%ungroup()%>%filter(species=="Porites")%>%select(tag)%>%distinct()%>%tally() #n colonies

m_input<-data%>%ungroup()%>%filter(species=="Montipora")%>%filter()%>%select(phenotype,contains('x'))
ggplot(m_input%>%rownames_to_column((var="row"))%>%gather(wv,ref,-phenotype,-row)%>%separate(wv,into=c("trash","wv"),sep=1)%>%select(-trash)%>%mutate(wv=as.numeric(wv)))+
     geom_line(aes(wv,ref,group=row,color=phenotype))

set.seed(3839)
grid<-expand.grid(interaction.depth=c(1,2,3,4),n.trees=seq(1,1000,50),shrinkage=c(0.002,0.001,0.01,0.005),n.minobsinnode=c(2,3))
seeds <- vector(mode = "list", length = 11) #nrepeats * nfolds +1
for(i in 1:11) seeds[[i]]<- sample.int(n=nrow(grid), nrow(grid)) #last value is at least grid length
seeds[[11]]<-sample.int(1000, 1)
model<- caret::train(phenotype ~., data=m_input, method="gbm",
                     trControl=trainControl(method="cv",number=5,search="grid"),
                     preProcess=c("center","scale"),
                     tuneGrid=grid,
                     metric="Accuracy",
                     bag.fraction=1,
                     verbose=FALSE);model
confusionMatrix(model)
plot(model)
model$bestTune
model$resample
mean(model$resample$Accuracy);sd(model$resample$Accuracy)
prediction<-predict(model,newdata=m_input,type=c("raw"));cm<-confusionMatrix(prediction,as.factor(m_input$phenotype));cm
insitu_prediction<-cbind(as.data.frame(prediction),data%>%filter(species=="Montipora")%>%select(species,site,tag,phenotype))%>%rename(insitu_prediction=5)%>%mutate(tag=as.factor(tag))
remote_model<-readRDS("./output/model_fits/0_montipora_gbm_model")
remote_data<-cbind(readRDS("./output/processed_spectra/phenotype_metadata"),readRDS("./output/processed_spectra/phenotype_data_bn"))[-117,]

remote_prediction<-predict(remote_model,newdata=remote_data,type=c("raw"))
m_eval<-cbind(remote_data%>%select(Tag,Site,Species,Phenotype)%>%mutate(Tag=as.factor(Tag)),remote_prediction)%>%clean_names()%>%
     inner_join(.,insitu_prediction,by=c("tag","site"))%>%select(-species.x,-species.y,-prediction)%>%
     mutate(insitu_prediction=str_to_lower(insitu_prediction))%>%
     mutate(correct=case_when(remote_prediction==phenotype&remote_prediction==insitu_prediction~1,TRUE~0))%>%group_by(tag)%>%sample_n(1)%>%
     ungroup()%>%mutate(species="M")%>%group_by(species)%>%count(correct)


p_input<-data%>%ungroup()%>%filter(species=="Porites")%>%select(phenotype,contains('x'))
ggplot(p_input%>%rownames_to_column((var="row"))%>%gather(wv,ref,-phenotype,-row)%>%separate(wv,into=c("trash","wv"),sep=1)%>%select(-trash)%>%mutate(wv=as.numeric(wv)))+
     geom_line(aes(wv,ref,group=row,color=phenotype))

set.seed(3839)
grid<-expand.grid(interaction.depth=c(1,2,3,4),n.trees=seq(1,500,5),shrinkage=c(0.002,0.001,0.01,0.005),n.minobsinnode=c(1,2,3))
seeds <- vector(mode = "list", length = 11) #nrepeats * nfolds +1
for(i in 1:11) seeds[[i]]<- sample.int(n=nrow(grid), nrow(grid)) #last value is at least grid length
seeds[[11]]<-sample.int(1000, 1)
model<- caret::train(phenotype ~., data=p_input, method="gbm",
                     trControl=trainControl(method="cv",number=5,,seeds=seeds),
                     preProcess=c("center","scale"),
                     tuneGrid=grid,
                     metric="Accuracy",
                     bag.fraction=0.75,
                     verbose=FALSE);model
confusionMatrix(model)
plot(model)
model$bestTune
mean(model$resample$Accuracy);sd(model$resample$Accuracy)
prediction<-predict(model,newdata=p_input,type=c("raw"));cm<-confusionMatrix(prediction,as.factor(p_input$phenotype)); cm
insitu_prediction<-cbind(as.data.frame(prediction),data%>%filter(species=="Porites")%>%select(species,site,tag,phenotype))%>%rename(insitu_prediction=5)%>%mutate(tag=as.factor(tag))
remote_model<-readRDS("./output/model_fits/0_porites_gbm_model")
remote_data<-cbind(readRDS("./output/processed_spectra/phenotype_metadata"),readRDS("./output/processed_spectra/phenotype_data_bn"))[-117,]

remote_prediction<-predict(remote_model,newdata=remote_data,type=c("raw"))
p_eval<-cbind(remote_data%>%select(Tag,Site,Species,Phenotype)%>%mutate(Tag=as.factor(Tag)),remote_prediction)%>%clean_names()%>%
     inner_join(.,insitu_prediction,by=c("tag","site"))%>%select(-species.x,-species.y,-prediction)%>%
     mutate(insitu_prediction=str_to_lower(insitu_prediction))%>%
     mutate(correct=case_when(remote_prediction==phenotype&remote_prediction==insitu_prediction~1,TRUE~0))%>%group_by(tag)%>%sample_n(1)%>%
     ungroup()%>%mutate(species="P")%>%group_by(species)%>%count(correct)


rbind(m_eval,p_eval)%>%group_by(correct)%>%summarise(sum=sum(n))%>%mutate(total=sum(sum),prop=sum/total)
#total 4 of 39 wrong = 89.7% agreement and correct


######################################### GAOxIN SITU REGRESSION ############################################################################################################### #####
rm(list = ls())
raw_insitu<-read_excel("1_FieldSpectra.xlsx",sheet="Working Phenotype")%>%clean_names()%>%filter(phenotype!="NA")%>%select(species,site,tag,phenotype,contains('x'))
raw_remote<-readRDS("./output/processed_spectra/phenotype_data_bn")%>%clean_names()
metadata<-readRDS("./output/processed_spectra/phenotype_metadata")%>%clean_names()
wavelength_labels<-readRDS("./output/misc/wavelength_labels")

remote<-cbind(metadata,raw_remote)%>%select(site,species,phenotype,tag,contains('x'),-x)%>%mutate(tag=as.factor(tag));colnames(remote)[5:56]<-wavelength_labels
insitu<-raw_insitu%>%group_by(site,species,phenotype,tag)%>%summarise_all(funs(mean))%>%mutate(species=as.factor(species),phenotype=as.factor(phenotype),tag=as.factor(tag));colnames(insitu)[5:56]<-wavelength_labels

data<-remote%>%gather(wv,remote,-site,-tag,-species,-phenotype)%>%inner_join(.,insitu%>%gather(wv,insitu,-site,-tag,-species,-phenotype),by=c("site","species","tag","wv"))%>%select(-phenotype.x)

list<-(data%>%select(tag)%>%distinct)$tag
spp<-data%>%select(tag,species)%>%distinct()
magic_for(print,silent=TRUE)

for (i in list){
     summary<-summary(lm(insitu~remote,data=data%>%filter(tag==i)))
     print(summary$adj.r.squared)}
out<-magic_result_as_dataframe()%>%rename(tag=i)%>%inner_join(.,spp,by="tag")%>%rename(r2=2);out

sqrt(mean(out$r2)) #correlatoin coefficient
out%>%group_by(species)%>%summarise(r2=mean(r2))

x = expression(paste("M.capitata ",R^2,"=0.919"))
y = expression(paste("P.compressa ",R^2,"=0.948"))

quartz()
a<-ggplot(data%>%filter(species=="Montipora"))+geom_abline(slope=1,color="gray")+
     geom_line(aes(remote,insitu,group=tag))+
     geom_point(aes(remote,insitu,color=species),alpha=0.3,size=1)+
     scale_y_continuous(lim=c(0,.25),breaks=seq(0,.25,.05))+
     scale_x_continuous(lim=c(0,.25),breaks=seq(0,.25,.05))+
     ylab("in situ Brightness Normalized Reflectance")+
     xlab("Remotely Sensed Brightness Normalized Reflectance (GAO)")+
     annotate("text", x = 0, y=0,hjust=-0.1,size=3,label=x)+
     scale_color_manual(values=c("blue","green"))+
     scale_fill_manual(values=c("blue","green"))+
     theme_classic(base_size=8)+
     theme(legend.position="none");a
b<-ggplot(data%>%filter(species=="Porites"))+geom_abline(slope=1,color="gray")+
     geom_line(aes(remote,insitu,group=tag))+
     geom_point(aes(remote,insitu,color=species),alpha=0.3,size=1)+
     scale_y_continuous(lim=c(0,.25),breaks=seq(0,.25,.05))+
     scale_x_continuous(lim=c(0,.25),breaks=seq(0,.25,.05))+
     ylab("in situ Brightness Normalized Reflectance")+
     xlab("Remotely Sensed Brightness Normalized Reflectance (GAO)")+
     annotate("text",  x = 0, y=0,hjust=-0.1,size=3,label=y)+
     scale_color_manual(values=c("green"))+
     scale_fill_manual(values=c("green"))+
     theme_classic(base_size=8)+
     theme(legend.position="none");b

plot_grid(a,b,nrow=1,align="h",axis="b",labels=c("A","B"),label_size=8,label_x=c(0.12,0.12))

######################################### CHLOROPHYLL   ######################################################################################################################## #####
rawdata<-read_excel("ChlxPhenotypeData.xlsx")%>%clean_names()%>%rename(chl_a=2,chl_c=3)
data<-rawdata%>%filter(chl_a>0)%>%filter(chl_c>0)%>%filter(phenotype!="NA")%>%mutate(phenotype=str_to_lower(phenotype))%>%
     mutate(site=as.factor(site),phenotype=as.factor(phenotype),species=as.factor(species),chl_a=as.numeric(chl_a),chl_c=as.numeric(chl_c))

summary(glm(chl_a~phenotype,data=data%>%filter(species=="M")))
summary(glm(chl_c~site*phenotype,data=data%>%filter(species=="M")))
summary(glm(chl_a~site*phenotype,data=data%>%filter(species=="P")))
summary(glm(chl_c~site*phenotype,data=data%>%filter(species=="P")))

library(lme4)
model_fullA <- lmer(chl_a ~species*phenotype+(1|site),data=data)
model_limA <- lmer(chl_a ~ species + (1|site),data=data)
anova(model_fullA,model_limA)
model_fullC <- lmer(chl_c ~species*phenotype+(1|site),data=data)
model_limC <- lmer(chl_c ~ species + (1|site),data=data)
anova(model_fullC,model_limC)










######################################### REEF 25 SUPPLEMENTAL ################################################################################################################# #####
library(sf);library(ggrepel);library(ggsn);library(cowplot);library(tidyverse);library(patchwork);library(jpeg);library(ggnewscale)
raw_sub<-st_read("./output/substrate_prediction_25.shp")
shp_sub<-fortify(as.data.frame(cbind(raw_sub,st_coordinates(raw_sub))))%>%filter(prediction!="Water")
points<-read_tsv("./points/example.txt")

w=7.2
h=3.17
quartz(w=w,h=h)
r25_substrate<-ggplot()+
     geom_point(aes(x,y,color=substrate),data=points)+
     geom_raster(data = shp_sub, aes(x = X, y = Y,fill=prediction),show.legend = FALSE)+
     coord_fixed(ratio = 1,xlim=c(621910,622105),ylim=c(2373480,2373665))+
     theme_minimal(base_size=8)+
     scale_fill_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     scale_color_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),
           #legend.position=c(0.8,0.7),
           legend.position=c(0.15,0.15),
           legend.key.size=unit(0.25,"cm"))+
     annotate("text",x=624645,y=2372486,label="Water not shown",size=2,hjust=0)+
     #annotate("text",x=624990,y=2372486,label=" Best-Fit Accuracy: 99.7%\nCV-Accuracy: 89.1%",size=2,hjust=1)+
     annotate("text",label="Reef 25",x=621920,y=2373665,size=2.5,color="black",hjust=0);r25_substrate

library(sf);library(ggsn)

raw_mcap<-st_read("./output/phenotype_prediction_montipora_25.shp")
shp_mcap<-fortify(as.data.frame(cbind(raw_mcap,st_coordinates(raw_mcap))))%>%#mutate(breaks=cut(nnblchd,breaks=c(0.5,0.9,0.99,0.999,0.9999,1),labels=c("0.5","0.9","0.99","0.999","0.999")))
     mutate(breaks=ntile(nnblchd,6))
table(shp_mcap$breaks)
as.vector(shp_mcap%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))


mcap<-ggplot()+
     geom_raster(data = shp_mcap, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(621910,622105),ylim=c(2373480,2373665))+
     theme_minimal(base_size=8)+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.8,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Prediction Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5))+
     scalebar(x.min=621910,x.max=621940,y.min=2373485,y.max=2373765,transform=FALSE,dist=25,dist_unit="m",
              height = .01, st.dist = 0.03,
              location="bottomleft",
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2)+
     scale_fill_viridis_c(direction=-1,name="Prediction\nNonbleached\n(2017)",labels=c("Moderate","","","","","High (~1)"));mcap


raw_pcomp<-st_read("./output/phenotype_prediction_porites_25.shp")
shp_pcomp<-fortify(as.data.frame(cbind(raw_pcomp,st_coordinates(raw_pcomp))))%>%
     mutate(breaks=ntile(nnblchd,6))
as.vector(shp_pcomp%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))

pcomp<-ggplot()+
     geom_raster(data = shp_pcomp, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(621910,622105),ylim=c(2373480,2373665))+
     theme_minimal(base_size=8)+
     scale_fill_viridis_c(direction=-1,labels=c("Moderate","","","","","High (0.82+)"))+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.75,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Predicted Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5));pcomp

quartz(w=7.2,h=2.45)
r25_substrate + mcap + pcomp

######################################### REEF 42 SUPPLEMENTAL ################################################################################################################# #####
library(sf);library(ggrepel);library(ggsn);library(cowplot);library(tidyverse);library(patchwork);library(jpeg);library(ggnewscale)
raw_sub<-st_read("./output/substrate_prediction_42.shp")
shp_sub<-fortify(as.data.frame(cbind(raw_sub,st_coordinates(raw_sub))))%>%filter(prediction!="Water")
points<-read_tsv("./points/example.txt")

r42_substrate<-ggplot()+
     geom_point(aes(x,y,color=substrate),data=points)+
     geom_raster(data = shp_sub, aes(x = X, y = Y,fill=prediction),show.legend = FALSE)+
     coord_fixed(ratio = 1,xlim=c(621619,621775),ylim=c(2375340,2375506))+
     theme_minimal(base_size=8)+
     scale_fill_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     scale_color_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),
           #legend.position=c(0.8,0.7),
           legend.position=c(0.85,0.1),
           legend.key.size=unit(0.25,"cm"))+
     annotate("text",x=624645,y=2372486,label="Water not shown",size=2,hjust=0)+
     #annotate("text",x=624990,y=2372486,label=" Best-Fit Accuracy: 99.7%\nCV-Accuracy: 89.1%",size=2,hjust=1)+
     annotate("text",label="Reef 42",x=621619,y=2375506,size=2.5,color="black",hjust=0);r42_substrate

library(sf);library(ggsn)

raw_mcap<-st_read("./output/phenotype_prediction_montipora_42.shp")
shp_mcap<-fortify(as.data.frame(cbind(raw_mcap,st_coordinates(raw_mcap))))%>%#mutate(breaks=cut(nnblchd,breaks=c(0.5,0.9,0.99,0.999,0.9999,1),labels=c("0.5","0.9","0.99","0.999","0.999")))
     mutate(breaks=ntile(nnblchd,6))
table(shp_mcap$breaks)
as.vector(shp_mcap%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))


mcap<-ggplot()+
     geom_raster(data = shp_mcap, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(621619,621775),ylim=c(2375340,2375506))+
     theme_minimal(base_size=8)+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.8,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Prediction Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5))+
     scalebar(x.min=621619,x.max=621649,y.min=2375332,y.max=2375500,transform=FALSE,dist=25,dist_unit="m",
              height = .01, st.dist = 0.03,
              location="bottomleft",
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2)+
     scale_fill_viridis_c(direction=-1,name="Prediction\nNonbleached\n(2017)",labels=c("Moderate","","","","","High (~1)"));mcap

raw_pcomp<-st_read("./output/phenotype_prediction_porites_42.shp")
shp_pcomp<-fortify(as.data.frame(cbind(raw_pcomp,st_coordinates(raw_pcomp))))%>%
     mutate(breaks=ntile(nnblchd,6))
as.vector(shp_pcomp%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))

pcomp<-ggplot()+
     geom_raster(data = shp_pcomp, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(621619,621775),ylim=c(2375340,2375506))+
     theme_minimal(base_size=8)+
     scale_fill_viridis_c(direction=-1,labels=c("Moderate","","","","","High (0.82+)"))+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.75,0.92),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Predicted Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5));pcomp

dev.size()
quartz(w=7.2,h=2.72)
r42_substrate + mcap + pcomp





######################################### REEF 44 SUPPLEMENTAL ################################################################################################################# #####
library(sf);library(ggrepel);library(ggsn);library(cowplot);library(tidyverse);library(patchwork);library(jpeg);library(ggnewscale)
raw_sub<-st_read("./output/substrate_prediction_44.shp")
shp_sub<-fortify(as.data.frame(cbind(raw_sub,st_coordinates(raw_sub))))%>%filter(prediction!="Water")
points<-read_tsv("./points/example.txt")

r44_substrate<-ggplot()+
     geom_point(aes(x,y,color=substrate),data=points)+
     geom_raster(data = shp_sub, aes(x = X, y = Y,fill=prediction),show.legend = FALSE)+
     coord_fixed(ratio = 1,xlim=c(620834,621150),ylim=c(2375300,2375527))+
     theme_minimal(base_size=8)+
     scale_fill_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     scale_color_manual(values=c("dimgray","blue","green","tan"),name="Substrate")+
     theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),
           #legend.position=c(0.8,0.7),
           legend.position=c(0.85,0.13),
           legend.key.size=unit(0.2,"cm"))+
     annotate("text",x=624645,y=2372486,label="Water not shown",size=2,hjust=0)+
     #annotate("text",x=624990,y=2372486,label=" Best-Fit Accuracy: 99.7%\nCV-Accuracy: 89.1%",size=2,hjust=1)+
     annotate("text",label="Reef 44",x=620834,y=2375527,size=2.5,color="black",hjust=0);r44_substrate

library(sf);library(ggsn)

raw_mcap<-st_read("./output/phenotype_prediction_montipora_44.shp")
shp_mcap<-fortify(as.data.frame(cbind(raw_mcap,st_coordinates(raw_mcap))))%>%#mutate(breaks=cut(nnblchd,breaks=c(0.5,0.9,0.99,0.999,0.9999,1),labels=c("0.5","0.9","0.99","0.999","0.999")))
     mutate(breaks=ntile(nnblchd,6))
table(shp_mcap$breaks)
as.vector(shp_mcap%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))


mcap<-ggplot()+
     geom_raster(data = shp_mcap, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(620834,621150),ylim=c(2375300,2375527))+
     theme_minimal(base_size=8)+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.15,0.85),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Prediction Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5))+
     scalebar(x.min=620834,x.max=620834,y.min=2375305,y.max=2375505,transform=FALSE,dist=25,dist_unit="m",
              height = .02, st.dist = 0.05,
              location="bottomleft",
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2)+
     scale_fill_viridis_c(direction=-1,name="Prediction\nNonbleached\n(2017)",labels=c("Moderate","","","","","High (~1)"));mcap

raw_pcomp<-st_read("./output/phenotype_prediction_porites_44.shp")
shp_pcomp<-fortify(as.data.frame(cbind(raw_pcomp,st_coordinates(raw_pcomp))))%>%
     mutate(breaks=ntile(nnblchd,6))
as.vector(shp_pcomp%>%group_by(breaks)%>%summarise(max=max(nnblchd))%>%select(max))

pcomp<-ggplot()+
     geom_raster(data = shp_pcomp, aes(x = X, y = Y,fill=breaks))+
     coord_fixed(ratio = 1,xlim=c(620834,621150),ylim=c(2375300,2375527))+
     theme_minimal(base_size=8)+
     scale_fill_viridis_c(direction=-1,labels=c("Moderate","","","","","High (0.82+)"))+
     theme(axis.text.x=element_blank(),
           axis.title.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position=c(0.15,0.85),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6))+
     guides(fill = guide_colorbar(title = "Predicted Probability\nNonbleached (2017)",
                                  direction="horizontal",
                                  label.position = "bottom",
                                  title.position = "top", 
                                  title.vjust = 1,
                                  label.vjust=0.5,
                                  barwidth = 4,
                                  barheight = 0.5));pcomp

dev.size()
quartz(w=7.2,h=2.1)
r44_substrate + mcap + pcomp




