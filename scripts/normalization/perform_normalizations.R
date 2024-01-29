#clear environment
rm(list = ls())


library(knitr) #used to print the results
library(chemCal) #for error propigation
library(tidyverse) #used for everything
library(progress) #for progress bar

# import data -------------------------------------------------------------

#load raw combined data
iso.df <- readRDS("Rdata/raw_data.rds")%>%
  select(IRMS.ID,Name,Weight.mg,Species,Facility,Type,Matrix.1,Matrix.2,isotope.raw,isotope.expected, isotope.precision)

# figure out all of the possible normalization combinations ---------------

#create a list of all possible standards
name.list <- iso.df %>%
  filter(Type!="Linearity") %>% #we don't want the working standards
  pull(Name) %>%
  unique() %>%
  as.character()

#make an empty list to store the standard/QC combinations
combos.list <- list()

onepoint.list <- list(standards=combn(name.list,1))
twopoint.list <- list(standards=combn(name.list,2))
threepoint.list <- list(standards=combn(name.list,3))
fourpoint.list <- list(standards=combn(name.list,4))

#make a list of lists
combos.list <- list(one.point=onepoint.list,
                 two.point=twopoint.list,
                 three.point=threepoint.list,
                 four.point=fourpoint.list)

combos.list[1]

# perform a normalization for each combination ----------------------------

method.list <- unique(names(combos.list))

#make a list of facilities
facility.list <- iso.df %>%
  pull(Facility) %>%
  unique()

#start a name counter at 0
name_counter <- 0

#count the number of combinations
for (combo in combos.list){
  for (col in 1:ncol(combo$standards)){
    name_counter <- name_counter+1
  }
}

n_iter <- name_counter

#make an empty list to store the normalization results
results.list <- list()

#set teo tracking variables to zero
method_counter <- 0
norm_counter <- 0

#initialize progress bar
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [:elapsed || :eta]",
                       total = n_iter,complete = "=",incomplete = "-",current = ">",
                       clear = FALSE, width = 100, show_after=0)
pb$tick(0)

#compute the normalizations
for(combo in combos.list){# for every QC/standard combination
  method_counter <- method_counter+1 #advance the method counter
  for (col in 1:ncol(combo$standards)){ #for every group of standards
    pb$tick() #advance the progress bar
    norm_counter <- norm_counter+1 #advance the normalization counter
    for (facility in facility.list){
      for (myspecies in c("N","C")){
        
        temp.df <- iso.df %>%
          filter(Species==myspecies)%>%
          filter(Facility==facility)
        
        method <- names(combos.list)[method_counter] #return the method
        standards <- combo$standards[,col] #list of standards
        
        #get expected values from the standards
        x <- temp.df %>%
          filter(Name %in% standards) %>%
          pull(isotope.expected)
        
        #get observed values of the standards
        y <- temp.df %>%
          filter(Name %in% standards) %>%
          pull(isotope.raw)
        
        for (name in unique(temp.df$Name)){
          
          z <- temp.df %>%
            filter(Name==name) %>%
            pull(isotope.raw)
        
        #perform different method for one-point calibration
          if (length(standards)==1){
            
            lm <- lm(x-y~1)
            slope <- 1
            slope.confi95 <- NA
            intercept <- coefficients(lm)
            intercept.confi95 <- intercept-confint(lm, level=0.95)[1,1]
            
            isotope.norm <- mean(z+intercept)
            isotope.SE <- summary(lm)$sigma
            isotope.confi95 <- intercept.confi95
            one.point.mean <- mean(x)
            
          } else {
            lm <- lm(y~x)
            
            slope <- coefficients(lm)[2]
            slope.confi95 <- slope-confint(lm, level=0.95)[2,1]
            intercept <- coefficients(lm)[1]
            intercept.confi95 <- intercept-confint(lm, level=0.95)[1,1]
            
            prediction <- inverse.predict(lm, z)
            
            isotope.norm <- prediction$Prediction
            isotope.SE <- prediction$`Standard Error`
            isotope.confi95 <- prediction$Confidence
            
            one.point.mean <- NA
             
          }
          
          new.data <- data.frame(Name=name,
                                 Species=myspecies,
                                 Method=method,
                                 Facility=facility,
                                 Normalization=norm_counter,
                                 Slope=slope,
                                 slope.confi95=slope.confi95,
                                 Intercept=intercept,
                                 intercept.confi95=intercept.confi95,
                                 isotope.max=max(y),
                                 isotope.min=min(y),
                                 isotope.range=max(y)-min(y),
                                 standards=paste(standards,collapse = '_'),
                                 isotope.norm=isotope.norm,
                                 isotope.SE=isotope.SE,
                                 isotope.confi95=isotope.confi95,
                                 one.point.mean)
          
          results.list <- append(results.list,list(new.data))
          
        }
      }
    } 
  }
}

pb$terminate()

results.df <- bind_rows(results.list)

sample_data.df <- iso.df %>%
  select(Name, Type, Matrix.1, Matrix.2, Species, isotope.expected, isotope.precision) %>%
  unique()

results.df <- full_join(sample_data.df, results.df) %>%
  mutate(isotope.deviation=isotope.norm-isotope.expected,
         deviation_significant=ifelse(abs(isotope.deviation)>(isotope.confi95+isotope.precision), TRUE, FALSE)) %>%
  arrange(Name, Species, Facility, Method)

results.df$SRM_logical <- mapply(grepl, results.df$Name, results.df$standards)

results.df <- results.df %>%
  mutate(Type=ifelse(SRM_logical,"SRM",
                     ifelse(Type=="Linearity",NA,"QC"))) %>%
  select(-c("SRM_logical"))

# assess matrix -----------------------------------------------------------

matrix.df <- results.df %>%
  drop_na(Type) %>%
  filter(Type=="SRM") %>%
  group_by(Normalization) %>%
  summarise(SRM.matrix.1 = ifelse(length(unique(Matrix.1))==1,
                                  paste(unique(Matrix.1), collapse = "/"),
                                  "Mixed"),
            SRM.matrix.2 = ifelse(length(unique(Matrix.2))==1,
                                  paste(unique(Matrix.2), collapse = "/"),
                                  "Mixed"),) %>%
  ungroup()

matrix.df %>%
  count(SRM.matrix.1)

matrix.df %>%
  count(SRM.matrix.2)

normalizations.df <- left_join(results.df,matrix.df) %>%
  rename("QC.matrix.1"="Matrix.1",
         "QC.matrix.2"="Matrix.2")

matrix.df <- normalizations.df %>%
  filter(Type=="QC") %>%
  mutate(Matrix.1=ifelse(QC.matrix.1==SRM.matrix.1,
                  "Matched",paste(SRM.matrix.1,"/",QC.matrix.1)),
         Matrix.2=ifelse(QC.matrix.2==SRM.matrix.2,
                         "Matched",paste(SRM.matrix.2,"/",QC.matrix.2)))

normalizations.df <- left_join(normalizations.df, matrix.df)

normalizations.df %>%
  count(Matrix.1)

normalizations.df %>%
  count(Matrix.2)

# determine extrapolation -------------------------------------------------

normalizations.df <- normalizations.df %>%
  filter(!is.na(Type)) %>%
  mutate(Extrapolation=ifelse(isotope.expected<isotope.min | isotope.expected>isotope.max,
                              "Extrapolated","Interpolated"))
normalizations.df %>%
  count(Extrapolation)

# save the normalizations -------------------------------------------------

saveRDS(normalizations.df,file="Rdata/normalized_data.rds")

write.csv(normalizations.df,"export/normalized_data.csv",row.names=FALSE)
