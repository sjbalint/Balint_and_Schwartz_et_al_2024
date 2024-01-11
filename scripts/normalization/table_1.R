#clear environment
rm(list = ls())

library(tidyverse) #used for everything

# import data -------------------------------------------------------------

#load raw combined data
iso.df <- readRDS("Rdata/raw_data.rds") %>%
  filter(Type=="Normalization") %>%
  select(Name, Species, isotope.raw.95CI) %>%
  unique() %>%
  pivot_wider(names_from="Species", values_from="isotope.raw.95CI") %>%
  arrange(Name)

#export the  data as a csv
write.csv(iso.df,"export/table1.csv",row.names=FALSE)


#load raw combined data
iso.df <- readRDS("Rdata/raw_data.rds") %>%
  filter(Name!="Reference Gas",
         Type!="Reference Gas") %>%
  select(Name, Type, Species, percent.mean, percent.95CI) %>%
  unique() %>%
  pivot_wider(names_from="Species", values_from=c("percent.mean", "percent.95CI")) %>%
  arrange(Type, Name)

#export the  data as a csv
write.csv(iso.df,"export/table2.csv",row.names=FALSE)
