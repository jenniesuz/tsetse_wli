#************************FIGURE 1***********************
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(scales)
library(RColorBrewer)
library(zoo)

dat <- read.csv("countData.csv",header=T) # read in count data

class(dat$spp)
dat$spp <- factor(dat$spp)
levels(dat$spp) <-  c("G. pallidipes" = expression(paste(italic("G.pallidipes")))
                       ,"G. swynnertoni"= expression(paste(italic("G. swynnertoni"))))

var.names<-  c("G. pallidipes" = expression(paste(italic("G.pallidipes")))
                       ,"G. swynnertoni"= expression(paste(italic("G. swynnertoni"))))


tiff("Fig_1.tiff", height = 2.5, width = 4, units = 'in', compression="lzw", res=600)
ggplot(dat, aes(x=trap.dist,y=men_cnt,col=qtr)) + 
  scale_y_continuous(trans='log10',limits=c(0.5, high=301),labels=c("0","10","100"),breaks=c(1,11,101)) +
  labs( y= "Mean tsetse catch per trap per day"
        , x="Distance from protected area boundary (m)" ) + 
  scale_color_manual(values=c("#875777","#eab051","#c0b9ac")) +
  xlim(-4000,10000) +
  geom_vline(xintercept=0) +
  geom_point(mapping = aes(x =trap.dist, y =men_cnt+1),size=0.5,shape=1) +
  scale_shape(solid = FALSE) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=9)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6.5)
        ,legend.key.size = unit(0.4,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=7.5)
        ,legend.position =c(0.93,0.33)
        ,strip.background = element_rect(fill="white",color="white")
        ,legend.title = element_blank()) +
  facet_wrap(~spp, labeller = labeller(.cols=label_parsed,
                                           .multi_line = FALSE))
dev.off()
