#clear environment
rm(list = ls())

# import packages ---------------------------------------------------------

library(tidyverse) #used for everything
library(ggpubr) #for regression stats

# import data -------------------------------------------------------------

normalizations.df <- readRDS("Rdata/normalized_data.rds")

# manipulate data for graphing --------------------------------------------

intermediate.df <- normalizations.df %>%
  filter(Method == "one.point") %>%
  filter(Type=="QC") %>%
  mutate(abs.deviation = (isotope.deviation),
         abs.difference = (one.point.mean-isotope.expected),
         Species.3 = factor(Species,
                             levels=c("N","C"),
                             labels=c("Nitrogen",
                                      "Carbon")))

# configure graphing ------------------------------------------------------

fig_path <- "figures/"

mywidth <- 6
myheight <- 6

basetheme <- list(
  theme_classic(),
  theme(
    text=element_text(size=12),
    strip.background = element_blank(),
    legend.title=element_blank(),
    panel.grid.major.x = element_line(linewidth=.1, color="gray"), 
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


# make graph --------------------------------------------------------------

ggplot(intermediate.df,aes(x=abs.difference, y=abs.deviation))+
  basetheme+
  geom_point(aes(shape=Species.3, fill=Species.3), alpha=0.5,color="black")+
  geom_smooth(method="lm", alpha=0.7, se = FALSE, show.legend=FALSE, color="black")+
  stat_cor(p.accuracy = 0.001, show.legend = FALSE)+
  labs(x="Isotopic distance between one-point calibrant and standard (‰)",
       y=expression("Inaccuracy Margin (‰)"))

ggsave(paste0(fig_path,"Fig_S2.png"),width=mywidth, height=myheight)
