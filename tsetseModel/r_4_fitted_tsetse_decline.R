
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(scales)
library(RColorBrewer)

source("r_1_model.R")
source("r_1_paramsForFit.R")

newdat <- read.csv("tsetseCounts.csv")

# model fit
a <- tsetse_run_mod(tsetse_params(alph=(0.5*2)/4
                                  ,l=(1/8)/4
                                  ,md.p= 2.534082e-06
                                            ,mo.val=4.103366e-02
                                          ))
a <- output_func(a)

a <- a[30:60,]

tiff("Fig_2.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
ggplot(newdat, aes(x=dist,y=tot+1)) + 
  geom_point(data=newdat,aes(x=dist,y=tot+1),size=1,alpha=0.2) +
  geom_line(data=a,aes(x =trans, y=vlt/100+1),linetype=1,size=0.5)  +
 scale_y_continuous(trans='log10',limits=c(1, high=501),labels=c("0","10","100","500"),breaks=c(1,11,101,501)) +
  xlim(-5000,10000) +
  labs( y= "Tsetse abundance (trap/day)"
        , x="Distance from protected area boundary (m)" 
        ,col=" "
   ) + 
 # geom_vline(xintercept=0) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position =c(0.9,0.5)
        ,legend.title = element_text(size=6)
  ) 
dev.off()