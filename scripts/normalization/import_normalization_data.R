#clear environment
rm(list = ls())

library(tidyverse) #used for everything
library(readxl) #used to import data


# import sample data ------------------------------------------------------

samples.df <- read.csv('raw/input_data.csv') %>% #import sample data
  pivot_longer(c(d15N.expected, d13C.expected, d15N.precision, d13C.precision),
               names_to = c("Species",".value"), names_sep = "\\.") %>% #convert to long data format
  drop_na(expected) %>% #remove any samples that don't have an expected value (aka blanks)
  rename("isotope.expected"="expected", "isotope.precision"="precision")

samples.df$Species <- factor(samples.df$Species,levels=c("d15N","d13C"),labels=c("N","C")) #convert species column to factor

samples.df <- samples.df[c(1:3,ncol(samples.df)-1,ncol(samples.df),6:ncol(samples.df)-2)] #reorganize the dataframe

str(samples.df) #check data structure

# import UNM CSI data -----------------------------------------------------

csi.df <- read.csv("raw/UNMCSI_rawdata.csv") #import IRMS data

csi.df <- csi.df %>%
  select(c("Row", "Identifier.1","Method", "Amount","Peak.Nr","Ampl..28","d.15N.14N","Ampl..44","d.13C.12C")) %>% #select useful columns
  filter(Peak.Nr==4 | Peak.Nr==5) %>% #reutrn the sample peaks for N and C
  mutate(dilution=as.numeric(str_extract(Method, "\\d+"))) %>%
  select(-Method) %>%
  mutate_if(is.character,as.factor) %>% #convert character data to factors
  mutate(Row=as.character(Row))

csi.df$Peak.Nr <- factor(csi.df$Peak.Nr,levels=c(4,5),labels=c("N","C")) #convert peak number to species column

mycolnames <- c("IRMS.ID", "Name","Weight.mg","Species","Amplitude", "isotope.raw", "dilution") #these will be our new column names

nitrogen.df <- csi.df %>% #make a dataframe of the nitrogen data
  filter(Peak.Nr=="N") %>% #return the nitrogen data
  select(-c("Ampl..44","d.13C.12C")) %>% #remove the extranious carbon columns
  drop_na()

colnames(nitrogen.df) <- mycolnames #use our column names

#repeat for carbon
carbon.df <- csi.df %>%
  filter(Peak.Nr=="C") %>%
  select(-c("Ampl..28","d.15N.14N")) %>%
  drop_na()

colnames(carbon.df) <- mycolnames

csi.df <- rbind(nitrogen.df,carbon.df) #create a new dataframe of the nitrogen and carbon data in a long format

csi.df$Facility <- factor("UNM CSI")

csi.df$Amplitude <- csi.df$Amplitude/1000

# import csi linearity data -----------------------------------------------

csi_refgas.df <- read.csv("raw/UNMCSI_refgas.csv") %>%
  filter(Type=="Linearity")

csi_refgas.df$Name <- "Reference Gas"
csi_refgas.df$Species <- "N"

csi_refgas.df <- csi_refgas.df %>%
mutate(IRMS.ID=paste0(row_number(),"_wg"),
       dilution=100) %>%
  select(c("IRMS.ID","Name","Amount","Species","Ampl..28","d.15N.14N", "dilution"))

colnames(csi_refgas.df) <- mycolnames

csi_refgas.df$Facility <- factor("UNM CSI")

csi_refgas.df$Amplitude <- csi_refgas.df$Amplitude/1000

# import elementar data ---------------------------------------------------

epa.df <- read.csv("raw/USEPA_rawdata.csv") #import IRMS data

mg.df <- epa.df %>%
  select(sample.id,N.pct,C.pct) %>%
  pivot_longer(c("N.pct","C.pct"),names_to="Species",values_to="percent") %>%
  mutate(Species=factor(Species,levels=c("N.pct","C.pct"),labels=c("N","C"))) %>%
  group_by(sample.id, Species) %>%
  summarize(percent.mean = mean(percent, na.rm=TRUE),
            percent.95CI = diff(t.test(percent,conf.level = 0.95)$conf.int)/2) %>%
  ungroup() %>%
  rename("Name"="sample.id")

epa.df <- epa.df %>%
  select(c("IRMS.id","sample.id","weight.mg","N.height.nA","d15N.permil","C.height.nA","d13C.permil", "dilution")) %>% #select useful columns
  mutate(across(sample.id,as.factor),
         across(c("weight.mg","N.height.nA","d15N.permil","C.height.nA","d13C.permil", "dilution",), as.numeric),
         across(IRMS.id, as.character),
         dilution=100-dilution) %>%
  rename("name"="sample.id")

mycolnames <- c("IRMS.ID","Name","Weight.mg","Amplitude","isotope.raw","dilution") #these will be our new column names

nitrogen.df <- epa.df %>% #make a dataframe of the nitrogen data
  select(IRMS.id, name,weight.mg,N.height.nA,d15N.permil, dilution,) %>% #remove the extranious carbon columns
  drop_na()

colnames(nitrogen.df) <- mycolnames #use our column names

nitrogen.df$Species <- factor("N")

#repeat for carbon
carbon.df <- epa.df %>% #make a dataframe of the nitrogen data
  select(IRMS.id, name,weight.mg,C.height.nA,d13C.permil, dilution) %>% #remove the extranious carbon columns
  drop_na()

colnames(carbon.df) <- mycolnames #use our column names

carbon.df$Species <- factor("C")

epa.df <- rbind(nitrogen.df,carbon.df) #create a new dataframe of the nitrogen and carbon data in a long format

epa.df$Facility <- factor("U.S. EPA")

# import elementar linearity data -----------------------------------------

epa_refgas.df <- read.csv("raw/USEPA_refgas.csv")

epa_refgas.df$Name <- "Reference Gas"

epa_refgas.df$weight.mg <- 0

nitrogen.df <- epa_refgas.df %>%
  filter(species=="N") %>%
  mutate(IRMS.ID=paste0(row_number(),"_wg"))

nitrogen.df$isotope.raw <- 1000*((nitrogen.df$X28.29/median(nitrogen.df$X28.29))-1)

carbon.df <- epa_refgas.df %>%
  filter(species=="C") %>%
  mutate(IRMS.ID=paste0(row_number(),"_wg"))

carbon.df$d45C <- 1000*((carbon.df$X45.44/median(carbon.df$X45.44))-1)

carbon.df$d46C <- 1000*((carbon.df$X46.44/median(carbon.df$X46.44))-1)

carbon.df$isotope.raw <- (1.0676*carbon.df$d45C)-((0.0338*1.0010*carbon.df$d46C))/(1-(0.0338*0.0021))

nitrogen.df <- nitrogen.df %>%
  select(c("IRMS.ID", "Name","weight.mg","species","amplitude.nA","isotope.raw"))

carbon.df <- carbon.df %>%
  select(c("IRMS.ID", "Name","weight.mg","species","amplitude.nA","isotope.raw"))

epa_refgas.df <- rbind(nitrogen.df,carbon.df)

colnames(epa_refgas.df) <-c("IRMS.ID", "Name","Weight.mg","Species","Amplitude", "isotope.raw")

epa_refgas.df$Facility <- factor("U.S. EPA")

# bind everything together ------------------------------------------------

iso.df <- bind_rows(csi.df,epa.df)

iso.df <- bind_rows(iso.df,csi_refgas.df)

iso.df <- bind_rows(iso.df,epa_refgas.df)

iso.df <- left_join(iso.df,samples.df)

iso.df <- left_join(iso.df, mg.df)

# calculate 95%CI for raw data --------------------------------------------

CI95.df <- iso.df %>%
  filter(Type!="Reference Gas",
         Name != "Reference Gas") %>%
  group_by(Name, Species, Facility) %>%
  summarize(isotope.raw.95CI = diff(t.test(isotope.raw,conf.level = 0.95)$conf.int)/2) %>%
  group_by(Name, Species) %>%
  summarize(isotope.raw.95CI = mean(isotope.raw.95CI, na.rm=TRUE)) %>%
  ungroup()

iso.df <- left_join(iso.df, CI95.df)

str(iso.df) #check data structure

# estimate mg equivalent for working gas linearity ------------------------

iso.df <- iso.df %>%
  mutate(mg=Weight.mg*(percent.mean/100),
         mg.diluted=ifelse(Species=="C",mg*((100-dilution)/100),mg)
         ) %>%
  group_by(Facility, Species) %>%
  mutate(intercept=coefficients(lm(mg.diluted~Amplitude))[1],
         slope=coefficients(lm(mg.diluted~Amplitude))[2]) %>%
  ungroup() %>%
  mutate(mg.diluted=ifelse(Matrix.1=="Reference Gas",
                   Amplitude*slope+intercept,
                   mg.diluted)) %>%
  select(-c(slope, intercept))


# replace sample names ----------------------------------------------------

iso.df <- iso.df %>%
  mutate(Name = str_replace(Name, "Green Chile Powder", "CSI Chile")) %>%
  mutate(Name = str_replace(Name, "Tuna Muscle", "CSI Tuna")) %>%
  mutate(Name = str_replace(Name, "Casein", "CSI Casein")) %>%
  mutate(Name = str_replace(Name, "Blue Grama", "CSI Blue Grama")) %>%
  mutate_if(is.character,as.factor) #convert character columns to factor

# save data ---------------------------------------------------------------

saveRDS(iso.df,file="Rdata/raw_data.rds")


# summarize the raw data --------------------------------------------------

iso.df <- iso.df%>%
  select(IRMS.ID,Name,Species,Facility,Type,Matrix.1, Matrix.2,isotope.raw,isotope.expected) %>%
  arrange(Species,Facility,Name)

#export the raw data as a csv for external assessment
write.csv(iso.df,"export/raw_data_long.csv",row.names=FALSE)

wide.df <- iso.df %>%
  drop_na(Type) %>%
  pivot_wider(names_from="Species",values_from=c("isotope.raw","isotope.expected"))

#export the raw data as a csv for external assessment
write.csv(wide.df,"export/raw_data_wide.csv",row.names=FALSE)
