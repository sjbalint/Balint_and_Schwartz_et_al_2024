#clear environment
rm(list = ls())

library(tidyverse) #used for everything

# import data -------------------------------------------------------------

iso.df <- readRDS("Rdata/raw_data.rds")

# configure graphing --------------------------------------r----------------

fig_path <- "figures/"

mywidth <- 6
myheight <- 6

mytheme <- list(
  theme_classic(),
  geom_hline(yintercept=0,linetype="dashed"),
  geom_point(alpha=0.8, size=2),
  theme(
    text=element_text(size=12),
    strip.background = element_blank(),
    legend.title=element_blank(),
    panel.grid.major.x = element_line(size=.1, color="gray"), 
    panel.grid.major.y = element_line(size=.1, color="gray"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position="top",
    strip.placement = "outside"),
  scale_shape_manual(values=c(21:24)),
  scale_fill_viridis_d(option="cividis", begin = 0, end = 1, direction=-1),
  scale_color_viridis_d(option="cividis", begin = 0, end = 1,direction=-1),
  labs(x=NULL,y=NULL),
  lims(y=c(-2,0.5))
)

# graph linearity results -------------------------------------------------

plot.df <- iso.df %>%
  filter(Type=="Linearity") %>%
  select(Name,Matrix.2,Species,Facility,isotope.raw, mg.diluted) %>%
  group_by(Name,Species,Facility) %>%
  mutate(isotope.deviation = isotope.raw-median(isotope.raw)) %>%
  ungroup()

plot.df <- plot.df %>%
  mutate(Species=factor(Species,
                          levels=c("N","C"),
                          labels=c("Nitrogen (mg)",
                                  "Carbon (mg after dilution)")),
         Facility=paste0(Facility,"\nLinearity Effect (‰)"),
         Matrix.2=factor(Matrix.2, levels=c("Reference Gas","High Organic","Plant","Soil")))

ggplot(plot.df,aes(mg.diluted,isotope.deviation,fill=Matrix.2,shape=Matrix.2))+
  mytheme+
  facet_grid(Facility~Species, switch="both", scales="free_x")

ggsave(paste0(fig_path,"Fig_6.png"), width=mywidth, height=myheight)


# plot combustion efficiency ----------------------------------------------

plot.df <- iso.df %>%
  filter(Type=="Linearity") %>%
  #filter(Matrix.2!="Reference Gas") %>%
  select(Name,Weight.mg,mg, Amplitude,Matrix.2,Species,Facility,isotope.raw) %>%
  mutate(relative.amplitude = Amplitude/Weight.mg) %>%
  group_by(Name,Species,Facility) %>%
  mutate(isotope.deviation = isotope.raw-median(isotope.raw),
         amplitude.deviation = relative.amplitude-median(relative.amplitude),
         Matrix.2=factor(Matrix.2, levels=c("Reference Gas","High Organic","Plant","Soil"))) %>%
  ungroup()

plot.df <- plot.df %>%
  mutate(Species=factor(Species,
                        levels=c("N","C"),
                        labels=c("Nitrogen (mg)",
                                 "Carbon (mg)")),
         Facility=paste0(Facility,"\nMajor Amplitude / Weight"))

ggplot(plot.df, aes(mg,amplitude.deviation,fill=Matrix.2,shape=Matrix.2))+
  mytheme+
  facet_grid(Facility~Species, switch="both", scales="free_x")

ggsave(paste0(fig_path,"Fig_S5.png"), width=mywidth, height=myheight)



