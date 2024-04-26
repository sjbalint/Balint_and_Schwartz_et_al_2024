#clear environment
rm(list = ls())

# import packages ---------------------------------------------------------

library(tidyverse) #used for everything
library(rstatix) #for dunn_test
library(rcompanion) #for cld
library(cowplot) #for combining plots
library(scales) #for psuedo_log
library(patchwork) #for combining plots
library(grid) #for adding annotations
library(ggsci) #for colors
library(ggrepel) #for geom_text_repel

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

intermediate.df %>% 
  select(standards, Species, Facility) %>%
  unique() %>%
  count()

intermediate.df %>%
  group_by(deviation_significant) %>%
  count()


# configure graphing ------------------------------------------------------

update_geom_defaults("point", list(shape = 21, fill="grey"))

fig_path <- "figures/"

mywidth <- 6
myheight <- 8

labels <- c(0.001,0.01,0.1,1,10)

theme_set(theme_classic())

basetheme <- list(
  theme(
    strip.background = element_blank(),
    panel.grid.major.y = element_line(linewidth=.1, color="gray"),
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position="top",
    strip.placement = "outside"),
  scale_shape_manual(values=c(21:25)),
  scale_fill_viridis_d(option="cividis", begin=0.3, aesthetics = c("colour", "fill")),
  labs(x=NULL,y=NULL)
)


# create a function to add compact letter display -------------------------

create_cld <- function(plot.df,stat.test,x){
  species_name="Species.3"
  
  results.list <- list()
  
  species.list <- plot.df %>%
    pull(species_name) %>%
    unique()

  for (myspecies in species.list){
    
    stat.temp <- stat.test %>%
      filter(!!rlang::sym(species_name)==myspecies)
    
    stat.temp <- cldList(p.adj~comparison,stat.temp,remove.zero = FALSE)
    
    stat.temp[species_name] <- myspecies
    
    if (myspecies==species.list[2]){
      stat.temp$Letter <- toupper(stat.temp$Letter)
    }
    
    results.list <- append(results.list,list(stat.temp))
  }
  
  cld.df <- bind_rows(results.list)
  
  box.df <- plot.df %>%
    group_by_at(c(species_name,x)) %>%
    summarise(boxplot=list(setNames(boxplot.stats(abs.deviation)$stats,
                                    c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker'))),
              .groups="keep") %>%
    unnest_wider(boxplot) %>%
    ungroup() %>%
    select(all_of(c(x,species_name,"median","upper_hinge","upper_whisker"))) %>%
    mutate(median.label=as.character(sprintf("%.3f",median)),
           Group=gsub(" ", "", !!rlang::sym(x)))
  
  cld.df <- left_join(cld.df,box.df,by = join_by(Group, Species.3)) %>%
    mutate(across(all_of(c(x,species_name)),as.factor))
  
  return (cld.df)
}

make_bar <- function(plot.df, x, percent=TRUE){
  
  percent.df <- plot.df %>%
    group_by_at(c("Species.3", x)) %>%
    count(deviation_significant) %>%
    complete(deviation_significant) %>%
    mutate(sum=sum(n),
           percentage=round(100*n/sum, digits=0),
           percentage=ifelse(is.na(percentage),0, percentage)) %>%
    filter(deviation_significant=="Significant")
  
  if (percent){
    percent.df$label <- paste0(percent.df$percentage,"%")
  } else {
    percent.df$label <- paste0(percent.df$n," / ", percent.df$sum)
  }
  
  p1 <- ggplot(percent.df, aes(x=get(x), y=percentage, fill=Species.3, group=interaction(get(x), Species.3))) +
    geom_text(color="white", data=percent.df,aes(label="X",y=25),
              position = position_dodge(width=-1.6), hjust=0.5)+
    geom_bar(alpha=0.5, width=0.65, stat="identity", position=position_dodge(width=0.75), color="black")+
    basetheme+
    theme(legend.position="none")+
    labs(y="% Signif.")+
    scale_y_continuous(limits=c(0,100), breaks=c(0,50, 100))+
    geom_text(color="black", data=percent.df,aes(label=label, y=percentage),
              size = 3, position = position_dodge(width=0.75), vjust=-1)
    #scale_x_discrete(expand=expansion(mult=c(0.4,0.42)))
}

make_boxplot <- function(plot.df, x, median_vjust=-0.7){
  
  stat.test <- plot.df %>%
    group_by(Species.3) %>%
    dunn_test(as.formula(paste("abs.deviation ~",x)),p.adjust.method="bonferroni") %>%
    mutate(comparison=paste0(group1,"-",group2),
           p.adj=p.adj)
  
  stat.test %>%
    select(Species.3,group1,group2,n1,n2,p.adj) %>%
    print()
  
  cld.df <- create_cld(plot.df,stat.test,x)
  
  cld.df %>%
    select(all_of(c("Group","Species.3",x,"median.label"))) %>%
    print()
  
  percent.df <- plot.df %>%
    group_by_at(c("Species.3", x)) %>%
    count(deviation_significant) %>%
    complete(deviation_significant) %>%
    mutate(sum=sum(n),
           percentage=round(100*n/sum, digits=0),
           percentage=ifelse(is.na(percentage),0, percentage),
           percentage=paste0(percentage, "%")) %>%
    filter(deviation_significant=="Significant")
  
  p2 <- ggplot(plot.df,aes(y=abs.deviation,x=get(x),fill=Species.3, color=Species.3,
                           group=interaction(get(x),Species.3),label=abs.deviation))+
    basetheme+
    #geom_point(alpha=0.3, color="black", aes(shape=deviation_significant), position = position_jitterdodge(jitter.width = 0.8))+
    geom_boxplot(alpha=0.5,
                 outlier.shape=21,
                 color="black")+
    theme(
      legend.title=element_blank(),
      legend.position="top",
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank())+
    scale_y_log10(limits=c(0.001,10), breaks=labels,labels=labels)+
    labs(y=paste("Deviation from certified value (‰)"))+
    geom_text(color="black", data = cld.df, aes(y=median,label=median.label),
              size = 3, position=position_dodge(width=-0.75), vjust = median_vjust)+
    geom_text(color="black", data=cld.df,aes(label=Letter,y=upper_hinge),
              position = position_dodge(width=-1.6), hjust=0.5, vjust=-1)
  
  
  return(p2)
  
}

make_plot <- function(plot.df, x, label=NULL, xlab=NULL, label.position="left", 
                      xaxis=TRUE, yaxis=TRUE, show.legend=TRUE, median_vjust=-0.7,
                      percent.label=TRUE){
  
  remove_y_axis <- function(plot){
    plot <- plot+
      labs(y=NULL)+
      theme(
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
    return(plot)
  }
  
  p1 <- make_boxplot(plot.df, x, median_vjust=median_vjust)
  
  p2 <- make_bar(plot.df, x, percent=percent.label)+
    labs(x=xlab)
  
  if (label.position=="left"){
    p1 <- p1 + annotation_custom(
      grobTree(textGrob(label,
                        gp=gpar(fontsize=10), x=0.1,  y=0.9, just="left")))
  }
  
  if (label.position=="right"){
    p1 <- p1 + annotation_custom(
      grobTree(textGrob(label,
                        gp=gpar(fontsize=10), x=0.9,  y=0.9, just="right")))
  }
  
  if (show.legend==FALSE){
    p1 <- p1 +
      theme(legend.position="none")
  }
  
  if (yaxis==TRUE){
    p1 <- p1 + annotation_logticks(sides = "l")
  } else {
    p1 <- remove_y_axis(p1)
    p2 <- remove_y_axis(p2)
  }
  
  if (xaxis==FALSE){
    p2 <- p2+
      labs(x=NULL)+
      theme(axis.text.x = element_blank())
  }
  
  if (is.null(xlab)){
    p2_height <- 0.2
  } else {
    p2_height <- 0.28
  }
  
  p3 <- plot_grid(p1, p2, ncol=1, align="v", axis="tblr", rel_heights = c(0.9, p2_height))
  
  return(p3)
}

# plot isotopic range of standards ----------------------------------------

plot.df <- raw_data.df %>%
  select(Name, Type, Matrix.2, Species, isotope.expected, isotope.precision) %>%
  unique() %>%
  pivot_wider(names_from="Species", values_from=c("isotope.expected", "isotope.precision")) %>%
  drop_na(Type) %>%
  filter(Matrix.2!="Reference Gas")


ggplot(plot.df,aes(x=isotope.expected_C,y=isotope.expected_N,
                   shape=Type,fill=Matrix.2,label=Name))+
  basetheme+
  theme(panel.grid.major.x=element_line(color="grey", linewidth=0.1),
        legend.position=c(0.2, 0.7))+
  geom_point(size=3,alpha=0.7, color="black")+
  geom_errorbar(aes(xmin=isotope.expected_C-isotope.precision_C,
                    xmax=isotope.expected_C+isotope.precision_C,))+
  geom_errorbar(aes(ymin=isotope.expected_N-isotope.precision_N,
                   ymax=isotope.expected_N+isotope.precision_N,))+
  geom_text_repel(fill="white",point.padding=5, size=2)+
  labs(y=bquote(delta^15*N~'(‰)'),
      x=bquote(delta^13*C~'(‰)'),
      fill="Matrix")

ggsave(paste0(fig_path,"Fig_S1.png"),width=mywidth, height=myheight)


# normalization method comparison -----------------------------------------

plot.df <- intermediate.df %>%
  filter(isotope.range>20)

x <- "Method.2"

label <- "Isotopic Range > 20‰"

p1 <- make_plot(plot.df, x, label, xlab="Number of Reference Materials", show.legend=FALSE)

p1b <- make_plot(plot.df, x, label, xlab="Number of Reference Materials", show.legend=FALSE,
                percent.label=FALSE)

# normalization method comparison in poor cases ---------------------------

plot.df <- intermediate.df %>%
  filter(isotope.range<20)

label <- "Isotopic Range < 20‰"

p2 <- make_plot(plot.df, x, label, label.position="right", 
                xlab="Number of Reference Materials", yaxis=FALSE, show.legend=FALSE)

p2b <- make_plot(plot.df, x, label, label.position="right", percent.label=FALSE,
                xlab="Number of Reference Materials", yaxis=FALSE, show.legend=FALSE)

# combined plot -----------------------------------------------------------

mylegend <- get_legend(make_boxplot(plot.df, x))

p3 <- plot_grid(p1, p2, nrow = 1, align="vh", axis="tblr", rel_widths=c(1,0.9),
                labels="AUTO", label_x=0.12, label_y=0.98) #configure labels for the subplots

plot_grid(mylegend,p3,nrow=2,rel_heights=c(0.1,2))

ggsave(paste0(fig_path,"Fig_2.png"),width=mywidth*2, height=myheight)

p4 <- plot_grid(p1b, p2b, nrow = 1, align="vh", axis="tblr", rel_widths=c(1,0.9),
                labels="AUTO", label_x=0.12, label_y=0.98) #configure labels for the subplots

plot_grid(mylegend,p4,nrow=2,rel_heights=c(0.1,2))

ggsave(paste0(fig_path,"Fig_S3.png"),width=mywidth*2, height=myheight)

# matrix ------------------------------------------------------------------

plot.df <- intermediate.df %>%
  filter(Method=="three.point",
         Extrapolation=="Interpolated",
         isotope.range>20) %>%
  mutate(Matrix.3=ifelse(Matrix.2=="Matched",Matrix.2, "Mixed"))

x="Matrix.3"

label <- "Isotopic Range > 20‰\nThree Point\nInterpolated"

make_plot(plot.df, x, label, median_vjust=-1)

ggsave(paste0(fig_path,"Fig_3.png"),width=mywidth, height=myheight)


# matrix figure for SIs ---------------------------------------------------

plot.df <- intermediate.df %>%
  filter(Method=="three.point",
         Extrapolation=="Interpolated",
         isotope.range>20
         ) %>%
  mutate(
    Matrix.1=gsub(" / ", "\n", Matrix.1),
    Matrix.1=factor(Matrix.1,
    levels=c("Matched","Compound\nMacromolecule","Mixed\nCompound",
             "Mixed\nMacromolecule","Mixed\nPlant","Compound\nPlant")))

x="Matrix.1"

label <- "Isotopic Range > 20‰\nThree Point\nInterpolated"

make_plot(plot.df, x, label,
          xlab="Reference Material Matrix\nQuality Control Matrix")

ggsave(paste0(fig_path,"Fig_S4.svg"),width=mywidth*2, height=myheight)


# extrapolation -----------------------------------------------------------

plot.df <- intermediate.df %>%
  filter(isotope.range>20)%>%
  filter(Method=="three.point") %>%
  filter(Matrix.2=="Matched") %>%
  mutate(Extrapolation=factor(Extrapolation, levels=c("Interpolated", "Extrapolated")))

x="Extrapolation"

label <- "Isotopic Range > 20‰\nThree Point\nMatrix Matched"

make_plot(plot.df, x, label, median_vjust=1.5)

ggsave(paste0(fig_path,"Fig_4.png"),width=mywidth, height=myheight)

# range -------------------------------------------------------------------

plot.df <- intermediate.df %>%
  filter(Method=="three.point") %>%
  filter(Extrapolation=="Interpolated")

x <- "iso.rangebreaks"

label <- "Three Point\nInterpolated"

p1 <- make_plot(plot.df, x, label, xaxis=FALSE, show.legend=FALSE)

plot.df <- intermediate.df %>%
  filter(Method=="two.point") %>%
  filter(Extrapolation=="Interpolated")

label <- "Two Point\nInterpolated"

p2 <- make_plot(plot.df, x, label, label.position="right", 
                xaxis=FALSE, yaxis=FALSE, show.legend=FALSE)

plot.df <- intermediate.df %>%
  filter(Method=="three.point") %>%
  filter(Extrapolation=="Extrapolated")

label <- "Three Point\nExtrapolated"

p3 <- make_plot(plot.df, x, label, label.position="left",
                xlab="Isotopic Range (‰)",show.legend=FALSE, median_vjust=-0.5)

plot.df <- intermediate.df %>%
  filter(Method=="two.point") %>%
  filter(Extrapolation=="Extrapolated")

label <- "Two Point\nExtrapolated"

p4 <- make_plot(plot.df, x, label, yaxis=FALSE, xlab="Isotopic Range (‰)",
                label.position="right", show.legend=FALSE, median_vjust=-0.5)

mylegend <- get_legend(make_boxplot(plot.df, x))

plot <- plot_grid(p1, p2, p3, p4, nrow = 2, align="vh", axis="tblr",
                rel_widths=c(1,0.9, 1, 0.9), rel_heights=c(1,1, 0.9, 0.9),
                labels="AUTO", label_x=0.15, label_y=0.95) #configure labels for the subplots

plot_grid(mylegend,plot, nrow=2, rel_heights=c(0.1,2))

ggsave(paste0(fig_path,"Fig_5.svg"),width=mywidth*1.5, height=myheight*1.5)
