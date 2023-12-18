### Fly RF test

### Load Packages
{
  library("SciViews")
  library("markovchain")
  library("pracma")
  library("expm")
  library("diagram")
  library("devtools")
  #library("flowcatchR")
  library("gganimate")
  library("ggpubr")
  library("rstatix")
  library("multcomp")
  library(magick)
  library("viridis")    
  library(transformr)
  library(openxlsx)
  library(tidyverse)
  library(janitor)
  library(readxl)
  library(reshape2)
  library(ggplot2)
  library("gridExtra")
  library(openxlsx)
  library(depmixS4)
  options(java.parameters = "-Xmx8g")
  library("viridis")    
  library(tidyverse)
  library(janitor)
  library(readxl)
  library(reshape2)
  library(ggplot2)
  library(ggpointdensity)
  library(mclust)
  #library(rhdf5)
  library(openxlsx)
  library(data.table)
  library(gdata)
  library(ggpubr)
  library(skimr)
  library(ggstatsplot)
  library(mclust)
  library("xlsx")
  library(depmixS4)
  library(qdapRegex)
  library(dplyr)
  library(cycleRtools)
  library("sur")
  library(ranger)
  library(boot)
  library(zoo)
  library(scales)
  library(spm)
  library(caret)
}

###Set Seed for reproduciblity
set.seed(123)

####Groups for analysis
merged.wrmt.final= NULL
numbers <-function(x)(as.numeric(x)) ##Use to convert characters to numbers for xy coords
n= 10
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),] ## Use to randomly remove rows for RF training

####Import the data for both controls and predicted data points 
# Choose folder with files, change pattern to relfect file type
### Choose Folder for Mac
choose.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}

filename <- choose.dir()

if(is.na(filename)) stop("No folder", call. = F)

#Import Data
files.wrmt<- list.files(filename , pattern = ".xlsx", full.names = TRUE) %>%
  set_names()

#Ground Truth Group
merged.wrmt.GT= map_dfr(files.wrmt, read_excel,col_types = "text", .id = "file.ID") %>%
  mutate_at(c("x","y", "timebin","animalNumber"),numbers)

# Ground Truth 
{
  coolio.3.GT= as.data.frame(merged.wrmt.GT)
  Worm.test.GT = coolio.3.GT %>%
    mutate(file= qdapRegex::ex_between(file.ID,"Fixed/", "x")) %>% # change to reflect the path before the file name to the file type
    mutate(worm_ID = paste0(animalNumber,".",file)) %>%
    mutate(file= qdapRegex::ex_between(file.ID,"Fixed/", "_")) %>%
    mutate(plateNumber = case_when(grepl(10, file) ~ 10,
                                   grepl(11, file) ~ 11,
                                   grepl(12, file) ~ 12,
                                   grepl(13, file) ~ 13,
                                   grepl(14, file) ~ 14,
                                   grepl(15, file) ~ 15,
                                   grepl(16, file) ~ 16,
                                   grepl(17, file) ~ 17,
                                   grepl(18, file) ~ 18,
                                   grepl(19, file) ~ 19,
                                   grepl(20, file) ~ 20,
                                   grepl(21, file) ~ 21,
                                   grepl(22, file) ~ 22,
                                   grepl(23, file) ~ 23,
                                   grepl(24, file) ~ 24,
                                   grepl(25, file) ~ 25,
                                   grepl(26, file) ~ 26,
                                   grepl(1, file) ~ 1,
                                   grepl(2, file) ~ 2,
                                   grepl(3, file) ~ 3,
                                   grepl(4, file) ~ 4,
                                   grepl(5, file) ~ 5,
                                   grepl(6, file) ~ 6,
                                   grepl(7, file) ~ 7,
                                   grepl(8, file) ~ 8,
                                   grepl(9, file) ~ 9)) %>%
    mutate(Sex = case_when(grepl('Male', file) ~ 'Male',
                           grepl('Herm', file) ~ 'Herm')) %>%
    mutate(food = case_when(grepl('BF', file) ~ 'Big Food',
                            grepl('SF', file) ~ 'Small Food')) %>%
    mutate(geno = case_when(grepl("OsmTra",file) ~ 'Osm-5::Tra-2',
                            grepl("FFPP",file) ~ 'PPdfr',
                            grepl("Tph", file) ~ 'Tph',
                            grepl("TPH", file) ~ 'Tph',
                            grepl("Pdfr", file) ~ 'Pdfr',
                            grepl("PDFR", file) ~ 'Pdfr',
                            grepl('Cat2',file) ~ 'Cat2',
                            grepl('Cat',file) ~ 'Cat2',
                            grepl('NTs',file) ~ 'Him5',
                            grepl('P',file) ~ 'Him5',
                            TRUE~ 'Him5')) %>%
    mutate(pheromone = case_when(grepl('Ascr', file) ~ 'Ascr#3',
                                 grepl('Ctrl', file) ~ 'Control',
                                 TRUE~ 'None')) %>% 
    mutate(spread = case_when(grepl('Dot', file) ~ 'Dot',
                              grepl('Even', file) ~ 'Even'
    )) %>% 
    mutate(Sex_Reverse = case_when(grepl('Fem', file) ~ 'Feminine',
                                   grepl('Masc', file) ~ 'Masculine',
                                   TRUE ~ 'WT'
    )) %>%
    mutate(WPP = case_when(grepl('Single', file) ~ 'Single',
                           TRUE ~ 'Multi'
    )) %>%
    mutate(Date= qdapRegex::ex_between(worm_ID,"_", ".")) %>%
    
    mutate(group = case_when(Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & pheromone=='Ascr#3' & spread =='Dot' ~ 'HermAscrDot',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & pheromone=='Ascr#3' & spread =='Dot' ~ 'MaleAscrDot',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & pheromone=='Ascr#3' & spread =='Even' ~ 'HermAscrEven',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & pheromone=='Ascr#3' & spread =='Even' ~ 'MaleAscrEven',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & pheromone=='Control' & spread =='Dot' ~ 'HermCtrlDot',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & pheromone=='Control' & spread =='Dot' ~ 'MaleCtrlDot',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & pheromone=='Control' & spread =='Even' ~ 'HermCtrlEven',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & pheromone=='Control' & spread =='Even' ~ 'MaleCtrlEven',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & Sex_Reverse== 'Feminine' ~ 'HermFem',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & Sex_Reverse== 'Masculine' ~ 'HermMasc',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & Sex_Reverse== 'Feminine' ~ 'MaleFem',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & Sex_Reverse== 'Masculine' ~ 'MaleMasc',
                             Sex == 'Herm' & geno == 'Him5' & Sex_Reverse== 'Feminine' ~ 'HermFem',
                             Sex == 'Herm' & geno == 'Him5' & Sex_Reverse== 'Masculine' ~ 'HermMasc',
                             Sex == 'Male' & geno == 'Him5' & Sex_Reverse== 'Feminine' ~ 'MaleFem',
                             Sex == 'Male' & geno == 'Him5' & Sex_Reverse== 'Masculine' ~ 'MaleMasc',
                             Sex == 'Herm' & food == 'Big Food' & geno == 'Him5' ~ 'HermFoodNBF',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' ~ 'HermFoodNSF',
                             #Sex == 'Herm' & geno == 'Him5' ~ 'HermFoodNSF',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & WPP=='Multi' ~ 'HermFoodNSF',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Him5' & WPP=='Single' ~ 'HermSingle',
                             Sex == 'Male' & food == 'Big Food' & geno == 'Him5' ~ 'MaleFoodNBF',
                             #Sex == 'Male' & food == 'Small Food' & geno == 'Him5' ~ 'MaleFoodNSF',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & WPP=='Multi' ~ 'MaleFoodNSF',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Him5' & WPP=='Single' ~ 'MaleSingle',
                             Sex == 'Male' & geno == 'Him5' ~ 'MaleFoodNSF',
                             Sex == 'Herm' & food == 'Big Food' & geno == 'Tph' ~ 'HermTph',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Tph' ~ 'HermTph',
                             Sex == 'Male' & food == 'Big Food' & geno == 'Tph' ~ 'MaleTph',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Tph' ~ 'MaleTph',
                             Sex == 'Herm' & food == 'Big Food' & geno == 'Pdfr' ~ 'HermPdfr',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Pdfr' ~ 'HermPdfr',
                             Sex == 'Male' & food == 'Big Food' & geno == 'Pdfr' ~ 'MalePdfr',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Pdfr' ~ 'MalePdfr',
                             Sex == 'Herm' & food == 'Big Food' & geno == 'Cat2' ~ 'HermCat2',
                             Sex == 'Herm' & food == 'Small Food' & geno == 'Cat2' ~ 'HermCat2',
                             Sex == 'Male' & food == 'Big Food' & geno == 'Cat2' ~ 'MaleCat2',
                             Sex == 'Male' & food == 'Small Food' & geno == 'Cat2' ~ 'MaleCat2',
                             Sex == 'Male' & food == 'Small Food' & geno == 'PPdfr' ~ 'MalePPdfr',
                             Sex == 'Herm' & geno == 'Pdf-1' ~ 'HermPdf1',
                             Sex == 'Herm' & geno == 'Pdf-2' ~ 'HermPdf2',
                             Sex == 'Male' & geno == 'Pdf-1' ~ 'MalePdf1',
                             Sex == 'Male' & geno == 'Pdf-2' ~ 'MalePdf2',
                             Sex == 'Male' & geno == 'Osm-5::Tra-2' ~ 'MaleSensFem'
    )) %>% mutate_at(c("x","y","radius"),numbers) 
  
  }

 class(Worm.test.GT$radius)
##### Calculate Various Features of Tracks Ground_Truth

#####sample rate of video
tim<-24  #seconds/frame
#########
subsamp<-1  #subsample rate frames

tracks.df<-as.data.frame(Worm.test.GT)

 tracks.df <- tracks.df %>%
   mutate(y=y/21276.6, x=x/21276.6) #%>% #1040pix/cm
#   mutate(y = 2.748-y) %>% dplyr::select(-`spread`,-`food`) #%>% na.omit() 

#tracks.df<- tracks.df %>% select(-...1,-`Unnamed: 0`) 

tracks2.GT<-tracks.df %>% select(-spread) %>% 
  # group_by(Date,group,plateNumber,animalNumber) %>% 
  # mutate(x= ((lag(x,2)+lag(x)+x+lead(x)+lead(x,4))/9),
  #        y= ((lag(y,2)+lag(y)+y+lead(y)+lead(y,4))/9)
  # )%>% ungroup() %>% 
  mutate(Date= qdapRegex::ex_between(worm_ID,"_", ".")) %>%
  group_by(Date,group,plateNumber,animalNumber) %>% 
  mutate(linspeed=sqrt((lead(x)-x)^2+(lead(y)-y)^2)) %>% 
  #####converts to mm/s #####
mutate(linspeed=(linspeed/tim)*10) %>% 
  mutate(first_d = (lead(linspeed,1)-linspeed)/tim) %>% 
  # mutate(first_d = abs(first_d)) %>% 
  mutate(dy= y - lag(y))%>%
  mutate(dx= x - lag(x)) %>% 
  mutate(q1= abs(dy/dx))%>%
  mutate(q2= abs(dx/dy))%>%
  mutate(angle = case_when(
    dy==0 & dx>0 ~180,
    dy == 0 & dx < 0 ~ 0,
    dx == 0 & dy > 0 ~ 90,
    dx == 0 & dy < 0 ~ 270,
    dx < 0 & dy < 0 ~ (270 + (180 / pi) * atan(q2)),
    dx > 0 & dy < 0 ~ (180 + (180 / pi) * atan(q1)),
    dx > 0 & dy > 0 ~ (90 + (180 / pi) * atan(q2)),
    dx < 0 & dy > 0 ~ ((180 / pi) * atan(q1))
  )) %>% 
  mutate(dthet= (angle - lead(angle))) %>% 
  mutate(theta= case_when
         (abs(dthet) >180 ~ abs(dthet)-360,
           abs(dthet) <180 ~ dthet)) %>% 
  mutate(ang= (abs(theta))) %>% na.omit() %>% 
  mutate(dt_ang=((lead(ang,1)-ang))/tim,               #average angular acceleration
         dt_ang=abs(dt_ang)) %>% na.omit() %>%         #get rid of random NAs                 
  mutate(d=linspeed*tim,
         rad_vel=(linspeed*sin(ang*(pi/180))), #tangential velocity - when absent straight line (vsin(o))
         tan_vel=(linspeed*cos(ang*(pi/180))),          #radial velocity - total linear velocity of curvilinear motion v= vparallel+vperp or linspeed=tan_vel+rad_vel
         rad_accel=lead(rad_vel)-rad_vel,
         tan_accel=lead(tan_vel)-tan_vel) %>% na.omit() %>% 
  ungroup() %>% group_by(Date,group,animalNumber,plateNumber) %>% 
  mutate(checksum=(rad_vel)+abs(tan_vel),
         mnra=mean(rad_accel),
         mnta=mean(tan_accel),
         mnrv=mean(rad_vel),
         mntv=mean(tan_vel),
         var_rad_accel=(rad_accel-mnra)^2,
         var_tan_accel=(tan_accel-mnta)^2,
         var_rad_vel=(rad_vel-mnrv)^2,
         var_tan_vel=(tan_vel-mntv)^2) %>%          #instantaneous rad_accel linear velocity^2/r
  mutate(turn= case_when(abs(dt_ang)>90 ~1,
                         abs(dt_ang)<=90 ~0)) %>% 
  mutate(r_calc = (lead(linspeed,1)-linspeed)) %>%
  mutate(reversal= case_when((r_calc)<0 ~1,
                             (dt_ang)>=0 ~0)) %>% 
  na.omit() %>% 
  filter(linspeed<500) %>% 
  ungroup() 

tracks2.GT.test= tracks2.GT %>% filter(animalNumber!=5)


### Rando forest 

# Remove points for double check 

graph.Pred= tracks2.GT %>% filter(Date=='2023-12-15') %>% select(-ground_truth)

tracks2.test= tracks2.GT %>% filter(Date!='2023-12-15')



{
  rparam_test<-tracks2.test %>% 
    dplyr::select(worm_ID,Date,timebin,group,animalNumber,plateNumber,`ground_truth`,radius,linspeed,first_d,ang,rad_vel,rad_accel,tan_vel,tan_accel,rad_accel,tan_accel,var_rad_vel,var_tan_vel,var_rad_accel,var_tan_accel) %>% 
    na.omit()
  
  
  
  #encode temporal dependence through centered median smoothing
  system.time(rparam<-rparam_test %>% dplyr::ungroup() %>% group_by(worm_ID) %>% 
                mutate(ang5=rollapply(ang,5,median,align='center',fill=NA),
                       ang10=rollapply(ang,10,median,align='center',fill=NA),
                       ang15=rollapply(ang,15,median,align='center',fill=NA),
                       ang30=rollapply(ang,30,median,align='center',fill=NA),
                       ang60=rollapply(ang,60,median,align='center',fill=NA),
                       ang120=rollapply(ang,120,median,align='center',fill=NA),
                       
                       
                       accel5=rollapply(first_d,5,median,align='center',fill=NA),
                       accel10=rollapply(first_d,10,median,align='center',fill=NA),
                       accel15=rollapply(first_d,15,median,align='center',fill=NA),
                       accel30=rollapply(first_d,30,median,align='center',fill=NA),
                       accel60=rollapply(first_d,60,median,align='center',fill=NA),
                       accel120=rollapply(first_d,120,median,align='center',fill=NA),
                       
                       ls5=rollapply(linspeed,5,median,align='center',fill=NA),
                       ls10=rollapply(linspeed,10,median,align='center',fill=NA),
                       ls15=rollapply(linspeed,15,median,align='center',fill=NA),
                       ls30=rollapply(linspeed,30,median,align='center',fill=NA),
                       ls60=rollapply(linspeed,60,median,align='center',fill=NA),
                       ls120=rollapply(linspeed,120,median,align='center',fill=NA),
                       
                       tvel5=rollapply(tan_vel,5,median,align='center',fill=NA),
                       tvel10=rollapply(tan_vel,10,median,align='center',fill=NA),
                       tvel15=rollapply(tan_vel,15,median,align='center',fill=NA),
                       tvel30=rollapply(tan_vel,30,median,align='center',fill=NA),
                       tvel60=rollapply(tan_vel,60,median,align='center',fill=NA),
                       tvel120=rollapply(tan_vel,120,median,align='center',fill=NA),
                       
                       rvel5=rollapply(rad_vel,5,median,align='center',fill=NA),
                       rvel10=rollapply(rad_vel,10,median,align='center',fill=NA),
                       rvel15=rollapply(rad_vel,15,median,align='center',fill=NA),
                       rvel30=rollapply(rad_vel,30,median,align='center',fill=NA),
                       rvel60=rollapply(rad_vel,60,median,align='center',fill=NA),
                       rvel120=rollapply(rad_vel,120,median,align='center',fill=NA),
                       
                       taccel5=rollapply(tan_accel,5,median,align='center',fill=NA),
                       taccel10=rollapply(tan_accel,10,median,align='center',fill=NA),
                       taccel15=rollapply(tan_accel,15,median,align='center',fill=NA),
                       taccel30=rollapply(tan_accel,30,median,align='center',fill=NA),
                       taccel60=rollapply(tan_accel,60,median,align='center',fill=NA),
                       taccel120=rollapply(tan_accel,120,median,align='center',fill=NA),
                       
                       raccel5=rollapply(rad_accel,5,median,align='center',fill=NA),
                       raccel10=rollapply(rad_accel,10,median,align='center',fill=NA),
                       raccel15=rollapply(rad_accel,15,median,align='center',fill=NA),
                       raccel30=rollapply(rad_accel,30,median,align='center',fill=NA),
                       raccel60=rollapply(rad_accel,60,median,align='center',fill=NA),
                       raccel120=rollapply(rad_accel,120,median,align='center',fill=NA),
                       
                       vtvel5=rollapply(var_tan_vel,5,median,align='center',fill=NA),
                       vtvel10=rollapply(var_tan_vel,10,median,align='center',fill=NA),
                       vtvel15=rollapply(var_tan_vel,15,median,align='center',fill=NA),
                       vtvel30=rollapply(var_tan_vel,30,median,align='center',fill=NA),
                       vtvel60=rollapply(var_tan_vel,60,median,align='center',fill=NA),
                       vtvel120=rollapply(var_tan_vel,120,median,align='center',fill=NA),
                       
                       vrvel5=rollapply(var_rad_vel,5,median,align='center',fill=NA),
                       vrvel10=rollapply(var_rad_vel,10,median,align='center',fill=NA),
                       vrvel15=rollapply(var_rad_vel,15,median,align='center',fill=NA),
                       vrvel30=rollapply(var_rad_vel,30,median,align='center',fill=NA),
                       vrvel60=rollapply(var_rad_vel,60,median,align='center',fill=NA),
                       vrvel120=rollapply(var_rad_vel,120,median,align='center',fill=NA),
                       
                       vtaccel5=rollapply(var_tan_accel,5,median,align='center',fill=NA),
                       vtaccel10=rollapply(var_tan_accel,10,median,align='center',fill=NA),
                       vtaccel15=rollapply(var_tan_accel,15,median,align='center',fill=NA),
                       vtaccel30=rollapply(var_tan_accel,30,median,align='center',fill=NA),
                       vtaccel60=rollapply(var_tan_accel,60,median,align='center',fill=NA),
                       vtaccel120=rollapply(var_tan_accel,120,median,align='center',fill=NA),
                       
                       vraccel5=rollapply(var_rad_accel,5,median,align='center',fill=NA),
                       vraccel10=rollapply(var_rad_accel,10,median,align='center',fill=NA),
                       vraccel15=rollapply(var_rad_accel,15,median,align='center',fill=NA),
                       vraccel30=rollapply(var_rad_accel,30,median,align='center',fill=NA),
                       vraccel60=rollapply(var_rad_accel,60,median,align='center',fill=NA),
                       vraccel120=rollapply(var_rad_accel,120,median,align='center',fill=NA),
                       
                       radius5=rollapply(radius,5,median,align='center',fill=NA),
                       radius10=rollapply(radius,10,median,align='center',fill=NA),
                       radius15=rollapply(radius,15,median,align='center',fill=NA),
                       radius30=rollapply(radius,30,median,align='center',fill=NA),
                       radius60=rollapply(radius,60,median,align='center',fill=NA),
                       radius120=rollapply(radius,120,median,align='center',fill=NA)
                       
                ) %>%
                na.omit() )
  
  
  #############################################
  
  rparam_mod<-rparam[,c(7:91)] #these are column ranges to select the ones that would be useful for the model you might have to change these for your set up
  
  
  ### Check control ratio of LS and plot######## Comment out if not using###############
  # tabpropherm= as.matrix(table(rparam$`ground_truth`))
  # 
  # prop.plotH= as.data.frame(prop.table(tabpropherm)) %>% tibble::rownames_to_column() %>%
  #   rename_at('rowname',~'state') %>% tibble::rownames_to_column() %>% rename_at('rowname',~'sex')
  # 
  # prop.plotH[prop.plotH==1 |prop.plotH==2]<-'Herm'
  # 
  # tabpropmale= as.matrix(table(rparam$`ground_truth`))
  # 
  # prop.plotM.Greg= as.data.frame(prop.table(tabpropmale)) %>% tibble::rownames_to_column() %>%
  #   rename_at('rowname',~'state') %>% tibble::rownames_to_column() %>% rename_at('rowname',~'sex')
  # 
  # prop.plotM[prop.plotM==1 |prop.plotM==2]<-'Male'
  # 
  # prop.plot=rbind(prop.plotM,prop.plotH)
  # 
  # ggplot(prop.plot, aes(x=sex,y=V1,fill=state))+
  #   geom_bar(position = 'stack',stat = 'identity')+
  #   labs(x='Group', y='Proportion of Time Spent in Locomotor State', title='Proportion of time Roaming and Dwelling') +
  #   theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))
  
  ####SUBSET for TRAINING
  set.seed(123)
  train.idx <- sample(nrow(rparam_mod),  4/5*nrow(rparam_mod))
  feat.train <- rparam_mod[train.idx,  ]
  feat.test <- rparam_mod[-train.idx, ]
  
  rg.feat <- ranger(as.factor(`ground_truth`) ~ ., data = feat.train, importance = "impurity_corrected")
  rg.feat$variable.importance
  
  imp.a=rg.feat$variable.importance %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
    rename_at('.',~'importance') %>% rename_at('rowname',~'observation') #%>% top_n(5)
  
  imp.plot.a=ggplot(imp.a,aes(x=reorder(observation,importance),y=importance))+
    geom_bar(stat = 'identity')+
    coord_flip()+
    labs(y='importance',x='observations', title='Random Forest Importance') +
    theme(plot.title = element_text(hjust=0.5, size=15, face='bold'))
  #imp.plot.a
  
  ####Check trained model against remaining untrained subset
  set.seed(123)
  rg.feat <- ranger(as.factor(`ground_truth`) ~ ., data = feat.train)
  pred.feat <- predict(rg.feat, data = feat.test)
  tab<-as.matrix(table(feat.test$`ground_truth`, pred.feat$predictions))
  prop.table(tab, 1)
  plt<-as.data.frame(prop.table(tab, 1))
  colnames(plt)<- c("Ground_truth","Prediction","Freq")
  plt$Freq<-round(plt$Freq,digits=3)
  cm.a<-ggplot(plt, aes(Ground_truth,Prediction, fill= Freq)) +
    geom_tile() + geom_text(aes(label=Freq)) +
    scale_fill_gradient2(low=muted("blue"),high=("orange"),midpoint = 0.5 )+
    labs(x = "Ground Truth",y = "Prediction",title = "Confusion Matrix")
  #cm.a
  
  
  #Test model against full training data set
  rparam2<-rparam[,c(7:91)]
  ranger(as.factor(`ground_truth`) ~ ., data = rparam2)
  
  set.seed(123)
  rg.feat <- ranger(as.factor(`ground_truth`) ~ ., data = feat.train,mtry = 10)
  
  
  ###########################################################################
  ###APPEND APPROPRIATE SEX DESIGNATION (e.g. rg.feat_m/h)   #################
  rg.feat_all<-rg.feat
  ###########################################################################
  pred.feat <- predict(rg.feat, data = rparam2)
  tab<-as.matrix(table(rparam2$`ground_truth`, pred.feat$predictions))
  table(rparam2$`ground_truth`, pred.feat$predictions)
  prop.table(tab, 1)
  rpred<-pred.feat$predictions
  tab
  
  
  # percent.table(rparam2$GT_sex, pred.feat$predictions)
  rparam2<-cbind(rparam,rpred=rpred)
  
  
  
  
  ###########################################################################
  ###APPEND APPROPRIATE SEX DESIGNATION (e.g. rg.feat_m/h)   #################
  rparam_all<-rparam2
  ###########################################################################
  }

### Check Results
imp.plot.a
cm.a

count= tracks2.test %>% count(ground_truth)

### Predict
{
  rparam_test<-graph.Pred %>% 
    dplyr::select(worm_ID,Date,timebin,group,animalNumber,plateNumber,radius,linspeed,first_d,ang,rad_vel,rad_accel,tan_vel,tan_accel,rad_accel,tan_accel,var_rad_vel,var_tan_vel,var_rad_accel,var_tan_accel) %>% 
    na.omit()
  
  ​
  ​
  #encode temporal dependence through centered median smoothing
  system.time(rparam<-rparam_test %>% dplyr::ungroup() %>% group_by(worm_ID) %>% 
                mutate(ang5=rollapply(ang,5,median,align='center',fill=NA),
                       ang10=rollapply(ang,10,median,align='center',fill=NA),
                       ang15=rollapply(ang,15,median,align='center',fill=NA),
                       ang30=rollapply(ang,30,median,align='center',fill=NA),
                       ang60=rollapply(ang,60,median,align='center',fill=NA),
                       ang120=rollapply(ang,120,median,align='center',fill=NA),
                       
                       
                       accel5=rollapply(first_d,5,median,align='center',fill=NA),
                       accel10=rollapply(first_d,10,median,align='center',fill=NA),
                       accel15=rollapply(first_d,15,median,align='center',fill=NA),
                       accel30=rollapply(first_d,30,median,align='center',fill=NA),
                       accel60=rollapply(first_d,60,median,align='center',fill=NA),
                       accel120=rollapply(first_d,120,median,align='center',fill=NA),
                       
                       ls5=rollapply(linspeed,5,median,align='center',fill=NA),
                       ls10=rollapply(linspeed,10,median,align='center',fill=NA),
                       ls15=rollapply(linspeed,15,median,align='center',fill=NA),
                       ls30=rollapply(linspeed,30,median,align='center',fill=NA),
                       ls60=rollapply(linspeed,60,median,align='center',fill=NA),
                       ls120=rollapply(linspeed,120,median,align='center',fill=NA),
                       
                       tvel5=rollapply(tan_vel,5,median,align='center',fill=NA),
                       tvel10=rollapply(tan_vel,10,median,align='center',fill=NA),
                       tvel15=rollapply(tan_vel,15,median,align='center',fill=NA),
                       tvel30=rollapply(tan_vel,30,median,align='center',fill=NA),
                       tvel60=rollapply(tan_vel,60,median,align='center',fill=NA),
                       tvel120=rollapply(tan_vel,120,median,align='center',fill=NA),
                       
                       rvel5=rollapply(rad_vel,5,median,align='center',fill=NA),
                       rvel10=rollapply(rad_vel,10,median,align='center',fill=NA),
                       rvel15=rollapply(rad_vel,15,median,align='center',fill=NA),
                       rvel30=rollapply(rad_vel,30,median,align='center',fill=NA),
                       rvel60=rollapply(rad_vel,60,median,align='center',fill=NA),
                       rvel120=rollapply(rad_vel,120,median,align='center',fill=NA),
                       
                       taccel5=rollapply(tan_accel,5,median,align='center',fill=NA),
                       taccel10=rollapply(tan_accel,10,median,align='center',fill=NA),
                       taccel15=rollapply(tan_accel,15,median,align='center',fill=NA),
                       taccel30=rollapply(tan_accel,30,median,align='center',fill=NA),
                       taccel60=rollapply(tan_accel,60,median,align='center',fill=NA),
                       taccel120=rollapply(tan_accel,120,median,align='center',fill=NA),
                       
                       raccel5=rollapply(rad_accel,5,median,align='center',fill=NA),
                       raccel10=rollapply(rad_accel,10,median,align='center',fill=NA),
                       raccel15=rollapply(rad_accel,15,median,align='center',fill=NA),
                       raccel30=rollapply(rad_accel,30,median,align='center',fill=NA),
                       raccel60=rollapply(rad_accel,60,median,align='center',fill=NA),
                       raccel120=rollapply(rad_accel,120,median,align='center',fill=NA),
                       
                       vtvel5=rollapply(var_tan_vel,5,median,align='center',fill=NA),
                       vtvel10=rollapply(var_tan_vel,10,median,align='center',fill=NA),
                       vtvel15=rollapply(var_tan_vel,15,median,align='center',fill=NA),
                       vtvel30=rollapply(var_tan_vel,30,median,align='center',fill=NA),
                       vtvel60=rollapply(var_tan_vel,60,median,align='center',fill=NA),
                       vtvel120=rollapply(var_tan_vel,120,median,align='center',fill=NA),
                       
                       vrvel5=rollapply(var_rad_vel,5,median,align='center',fill=NA),
                       vrvel10=rollapply(var_rad_vel,10,median,align='center',fill=NA),
                       vrvel15=rollapply(var_rad_vel,15,median,align='center',fill=NA),
                       vrvel30=rollapply(var_rad_vel,30,median,align='center',fill=NA),
                       vrvel60=rollapply(var_rad_vel,60,median,align='center',fill=NA),
                       vrvel120=rollapply(var_rad_vel,120,median,align='center',fill=NA),
                       
                       vtaccel5=rollapply(var_tan_accel,5,median,align='center',fill=NA),
                       vtaccel10=rollapply(var_tan_accel,10,median,align='center',fill=NA),
                       vtaccel15=rollapply(var_tan_accel,15,median,align='center',fill=NA),
                       vtaccel30=rollapply(var_tan_accel,30,median,align='center',fill=NA),
                       vtaccel60=rollapply(var_tan_accel,60,median,align='center',fill=NA),
                       vtaccel120=rollapply(var_tan_accel,120,median,align='center',fill=NA),
                       
                       vraccel5=rollapply(var_rad_accel,5,median,align='center',fill=NA),
                       vraccel10=rollapply(var_rad_accel,10,median,align='center',fill=NA),
                       vraccel15=rollapply(var_rad_accel,15,median,align='center',fill=NA),
                       vraccel30=rollapply(var_rad_accel,30,median,align='center',fill=NA),
                       vraccel60=rollapply(var_rad_accel,60,median,align='center',fill=NA),
                       vraccel120=rollapply(var_rad_accel,120,median,align='center',fill=NA),
                       
                       radius5=rollapply(radius,5,median,align='center',fill=NA),
                       radius10=rollapply(radius,10,median,align='center',fill=NA),
                       radius15=rollapply(radius,15,median,align='center',fill=NA),
                       radius30=rollapply(radius,30,median,align='center',fill=NA),
                       radius60=rollapply(radius,60,median,align='center',fill=NA),
                       radius120=rollapply(radius,120,median,align='center',fill=NA)
                       
                ) %>%
                na.omit() )
  
  
  
  #############################################
  
  #### Prediction
  ​
  rparam2<-rparam[,c(7:90)]
  set.seed(123)
  pred.feat <- predict(rg.feat_all, data = rparam2) #Change rgfeat to right sex
  ​
  rpred<-pred.feat$predictions
  # percent.table(rparam2$GT_sex, pred.feat$predictions)
  rparam2<-cbind(rparam,rpred=rpred)
  ​
  ​
  ​
  ###########################################################################
  ###APPEND APPROPRIATE SEX DESIGNATION (e.g. rg.feat_m/h)   #################
  rparam_allpred<-rparam2
  
  ####### Add State to new plot before rbind
  # rparam_allpred= rparam_allpred %>% mutate(State = case_when(grepl('Herm.Dwell', rpred) ~ 1,
  #                                                             grepl('Herm.Roam', rpred) ~ 2,
  #                                                             grepl('Male.Dwell',rpred)~3,
  #                                                             grepl('Male.Roam',rpred)~4,
  #                                                             grepl('Male.Tail Chase',rpred)~5
  # )) %>% as.data.frame()
  
}

check=unique(rparam_allpred$worm_ID) %>% as.data.frame()

Test.EP= rparam_allpred %>% filter(worm_ID=='1./Male2SF_2023-12-15.') 

### Autocorrelation
data <- tracks2.test$ground_truth #load your categorical series here. e.g. scan(file.choose(), what="character")
Tlen <- length(data)
states <- unique(data)
datanum <- match(data, states)
#for numerical variables
acf(datanum, lag.max = 22000)
################################

#for categorical variables
states <- unique(data)
nostates <- length(states)
datanum <- match(data, states)

#Binarization
bincodes <- diag(1,nostates)
databin <- bincodes[datanum,]
bd<-as.data.frame(bincodes)
db<-as.data.frame(databin)
#relative frequencies
hatpi <- colMeans(databin)
hp<-as.data.frame(hatpi)

maxlag <- 22000
hatbivprob <- array(0,c(nostates,nostates,maxlag))
hbp<-as.data.frame(hatbivprob)

for(k in c(1:maxlag)){#for each lag
  for(i in c(1:nostates)){#for each lagged vector representing a category
    for(j in c(1:nostates)){#for each vector representing a category
      hatbivprob[i,j,k] <- mean(databin[(k+1):Tlen,i]*databin[1:(Tlen-k),j])
    }
  }
}
#compare
indprob <- hatpi %*% t(hatpi)

#Cramer's v
cramer <- rep(0,maxlag)
for(k in c(1:maxlag)){
  cramer[k] <- sqrt(sum((hatbivprob[,,k]-indprob)^2/indprob)/(nostates-1))
}
dat_cramer<-as.data.frame(cramer)
tau<-min(which(dat_cramer <0.9))
#tau<-18

plot(cramer, type="h", xlab = "k", ylab = "Cramer's   v(k)", lwd=4, ylim=c(-1,1))

  

## Function to smooth the data
Smooth.blip<-function(x){ 
  require(zoo)
  require(dplyr)
  y<- rollapply(x$rpred, tau, function(x) names(which.max(table(x))),align ='center', fill= TRUE) #smooth the data by "looking" at a rolling window of 4 points (or based off of autocorrelation value) and find the state that occurs the most and make that the new state
  rpred.smooth<-y
  z<- cbind(x,rpred.smooth)
  return(z)
}


All.smooth.Pred= Test.EP %>% group_split(worm_ID) %>% lapply(Smooth.blip)

All.smooth=rbindlist(All.smooth.Pred) %>% ungroup()



### Plots for Ethogram


#### No Sex LS for simplicity #####################################
#### Ground Truth Plots 

graph.ETH= tracks2.GT %>% filter(worm_ID=='1./Male2SF_2023-12-15.') %>% as.data.frame() %>% 
  mutate(State.GT = case_when(grepl('0', ground_truth) ~ 1,
                           grepl('1', ground_truth) ~ 2,
                           grepl('2',ground_truth)~3,
                           grepl('3',ground_truth)~4,
                           grepl('4',ground_truth)~5))

Test.Eth= All.smooth %>% mutate(State.Pred = case_when(grepl('0', rpred.smooth) ~ 1,
                                                 grepl('1', rpred.smooth) ~ 2,
                                                 grepl('2',rpred.smooth)~3,
                                                 grepl('3',rpred.smooth)~4,
                                                 grepl('4',rpred.smooth)~5))



{
pH.LV<-ggplot(graph.ETH ,aes(x=`timebin`,y=`linspeed`))+
  geom_line()+
  xlab('Time (s)')+
  ylab('Linear Speed (mm/s)')+
  scale_x_continuous(limits = c(0, 18714))
#labs(title = 'Herm State Ground Truth')+
#theme(plot.title = element_text(hjust = 0.5))
pH.LV

pH.Rad<-ggplot(graph.ETH ,aes(x=`timebin`,y=`radius`))+
  geom_line()+
  xlab('Time (s)')+
  ylab('Radius (mm)')+
  scale_x_continuous(limits = c(0,18714))
#labs(title = 'Herm State Ground Truth')+
#theme(plot.title = element_text(hjust = 0.5))
pH.Rad
pH.GT<-ggplot(graph.ETH,aes(x=`timebin`,y=as.numeric(State.GT)))+
  geom_line()+
  xlab('Time (s)')+
  ylab('Ground Truth')+
  scale_y_continuous(limits = c(1, 3))+
  scale_x_continuous(limits = c(0,18714))
#labs(title = 'Herm State Ground Truth')+
#theme(plot.title = element_text(hjust = 0.5))
pH.GT

pH.RF<-ggplot(Test.Eth,aes(x=`timebin`,y=State.Pred))+
  geom_line()+
  xlab('Time (s)')+
  ylab('RF Pred')+
  scale_y_continuous(limits = c(1, 3))+
  scale_x_continuous(limits = c(0,18714))

#labs(title = 'Herm State Ground Truth')+
#theme(plot.title = element_text(hjust = 0.5))
pH.RF

}

# Herm GT vs RF

# Herm GT vs RF
pH.GTRF2<-ggarrange(pH.LV,pH.Rad,pH.GT,pH.RF, nrow = 4, common.legend = F, legend = 'right', heights = c(10, 10))

pH.GTRF2

check.GT= graph.ETH %>% select(timebin,ground_truth,State.GT)
check.Pred= Test.Eth %>% select(timebin,rpred,State.Pred)

check= merge.data.table(graph.ETH,Test.Eth, by= 'timebin') %>% select(`State.GT`,`State.Pred`) %>% 
  na.omit()

tab<-as.matrix(table(check))

prop.table(tab, 1)
plt<-as.data.frame(prop.table(tab, 1))
colnames(plt)<- c("Ground_truth","Prediction","Freq")
plt$Freq<-round(plt$Freq,digits=3)
cm.a<-ggplot(plt, aes(Ground_truth,Prediction, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient2(low=muted("blue"),high=("orange"),midpoint = 0.5 )+
  labs(x = "Ground Truth",y = "Prediction",title = "Confusion Matrix for New Data")
cm.a

## Export Data
file= paste(choose.dir(),"RF_Fly_results_2023-12-15.xlsx",sep = "")
openxlsx::write.xlsx(All.smooth, file=file)
  

