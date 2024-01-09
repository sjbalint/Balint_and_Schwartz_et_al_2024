#clear environment
rm(list = ls())

# import packages ---------------------------------------------------------

library(tidyverse) #used for everything

# import data -------------------------------------------------------------

normalizations.df <- readRDS("Rdata/normalized_data.rds")
raw_data.df <- readRDS("Rdata/raw_data.rds")

# manipulate data for graphing --------------------------------------------

intermediate.df <- normalizations.df %>%
  filter(Method != "one.point") %>%
  filter(Type=="QC") %>%
  filter(standards != "IAEA600_USGS91") %>%
  mutate(abs.deviation = abs(isotope.deviation)) %>%
  ungroup() %>%
  mutate(deviation_significant = factor(deviation_significant, 
                                        levels=c(TRUE, FALSE),
                                        labels=c("Significant","Nonsignificant")),
         Method.1 = factor(Method,
                           levels=c("one.point","two.point","three.point","four.point"),
                           labels=c("One Point","Two Point","Three Point","Four Point")),
         Method.2 = factor(Method,
                           levels=c("one.point","two.point","three.point","four.point"),
                           labels=c("One","Two","Three","Four")),
         Species.2 = factor(Species,
                            levels=c("N","C"),
                            labels=c("Nitrogen Deviation (‰)",
                                     "Carbon Deviation (‰)")),
         Species.3 = factor(Species,
                            levels=c("N","C"),
                            labels=c("Nitrogen",
                                     "Carbon")),
         Species = factor(Species,
                          levels=c("N","C"),
                          labels=c("Nitrogen\nDeviation\n(‰)",
                                   "Carbon\nDeviation\n(‰)")),
         iso.rangebreaks = cut(isotope.range,
                               breaks=c(0,15,30,45),
                               labels=c("0 to 15","15 to 30","30 to 45"))
  )


# calculate "final" slope -------------------------------------------------

species_list <- raw_data.df %>%
  pull(Species) %>%
  unique()

facility_list <- raw_data.df %>%
  pull(Facility) %>%
  unique()

results.list <- list()

for (species in species_list){
  for (facility in facility_list){
    
    slope <- raw_data.df %>%
      filter(Species==species) %>%
      filter(Facility==facility) %>%
      filter(Type=="Normalization")%>%
      lm(isotope.expected~isotope.raw, data=.)
    
    slope <- slope$coefficients[2]
    
    slope.df <- data.frame("Slope"=slope, "Species"=species, "Facility"=facility)
    
    results.list <- append(results.list,list(slope.df))
  }
}

slope.df <- bind_rows(results.list)

slope.df$Species.3 <- factor(slope.df$Species,
                                    levels=c("N","C"),
                                    labels=c("Nitrogen",
                                             "Carbon"))

# configure graphing ------------------------------------------------------

fig_path <- "figures/v6/"

mywidth <- 6
myheight <- 6

basetheme <- list(
  theme_classic(),
  theme(
    text=element_text(size=12),
    strip.background = element_blank(),
    legend.title=element_blank(),
    #panel.grid.major.x = element_line(linewidth=.1, color="gray"), 
    panel.grid.major.y = element_line(linewidth=.1, color="gray"),
    #axis.title.y = element_text(angle = 0,vjust = 0.5),
    #strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position="top",
    strip.placement = "outside"),
  scale_shape_manual(values=c(21:25)),
  scale_fill_viridis_d(option="cividis", direction=-1
                       #begin = 0.1, end = 0.6
  ),
  scale_color_viridis_d(option="cividis", direction=-1
                        #begin = 0.1, end = 0.6
  ),
  labs(x=NULL,y=NULL)
)


# slope vs range ----------------------------------------------------------

plot.df <- intermediate.df %>%
  filter(Method!="one.point") %>%
  select(Species,Slope,Facility,Method.1, Method.2, Species.3, isotope.range, deviation_significant) %>%
  unique()
#mutate(Method.1=factor(Method.1, levels=c("Four Point","Three Point","Two Point")))

ggplot(plot.df,aes(x=isotope.range, y=Slope, shape=Species.3, fill=Species.3, color=Species.3))+
  basetheme+
  geom_point(alpha=0.5,color="black")+
  geom_hline(color="black", data=slope.df,aes(yintercept=Slope), linetype="dashed")+
  facet_wrap(.~Method.1, nrow=1, strip.position="top")+
  scale_x_continuous(limits=c(0,43), expand = c(0, NA))+
  labs(x="Isotopic Range (‰)")+
  labs(y=expression("Normalization Slope ("*italic(m)*")"))

ggsave(paste0(fig_path,"Fig_7.png"),width=mywidth, height=myheight)
