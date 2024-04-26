#clear environment
rm(list = ls())

# import packages ---------------------------------------------------------

library(tidyverse) #used for everything

# import data -------------------------------------------------------------

normalizations.df <- readRDS("Rdata/normalized_data.rds")

# manipulate data for graphing --------------------------------------------

intermediate.df <- normalizations.df %>%
  filter(Method != "one.point") %>%
  filter(Type=="QC") %>%
  filter(standards != "IAEA600_USGS91") %>%
  mutate(abs.deviation = abs(isotope.deviation)) %>%
  mutate(deviation_significant = factor(deviation_significant, 
                                        levels=c(TRUE, FALSE),
                                        labels=c("Significant","Nonsignificant")),
         Species.3 = factor(Species,
                            levels=c("N","C"),
                            labels=c("Nitrogen",
                                     "Carbon")))


# set up graphing ---------------------------------------------------------

fig_path <- "figures/"

mywidth <- 6
myheight <- 6

labels <- c(0.001,0.01,0.1,1,10)

theme_set(theme_classic())

basetheme <- list(
  theme(
    strip.background = element_blank(),
    legend.title=element_blank(),
    panel.grid.major.x = element_line(linewidth=.1, color="gray"),
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position="top",
    strip.placement = "outside"),
  scale_shape_manual(values=c(21,4)),
  scale_fill_viridis_d(option="cividis", begin=0.5, aesthetics = c("colour", "fill")),
  labs(x=NULL,y=NULL)
)

# graph -------------------------------------------------------------------

ggplot(intermediate.df, aes(abs.deviation, fill=Species.3))+
  basetheme+
  geom_histogram(alpha=0.5, bins=15, position = "identity", color="black")+
  facet_wrap(~Name, ncol=2, scales="free_y")+
  scale_x_log10(limits=c(0.001,10),breaks=labels,labels=labels)+
  labs(y="Count", x="Deviation from certified value (‰)")

ggsave(paste0(fig_path,"Fig_1.png"), width=mywidth, height=myheight)
