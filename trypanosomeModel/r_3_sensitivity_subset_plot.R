#********************************************
library(plyr)
library(ggplot2)
#***********T congolense****************
co <- read.csv("sensTcOutput.csv")
tcorp <- read.csv("tcSens.csv")
tcorp <- tcorp[,-c(1,2,12,13,14,15,16,17)]
co1 <- merge(co,tcorp,by="run",all.x=T)  
#***********T brucei****************
br <- read.csv("sensTbOutput.csv")
tbrrp <- read.csv("tbSens.csv")
tbrrp <- tbrrp[,-c(1,2,12,13,14,15,16,17)]
br1 <- merge(br,tbrrp,by="run",all.x=T)  
#*********************************************

co1$spp <- "T. congolense"
br1$spp <- "T. brucei"
dat <- rbind.data.frame(co1,br1)


dat<- dat[!dat$trans < -5000,]
dat$trans <- as.factor(dat$trans)

labs <- c(expression(paste(italic("T. brucei")))
          ,expression(paste(italic("T. congolense")))
)

tiff("Fig_4_trypsInterface.tiff", height =3 , width = 5, units = 'in', compression="lzw", res=400)
ggplot(dat, aes(x=trans,y=perc.host)) +
  geom_boxplot(aes(),outlier.size=0.1) +
  labs( y= "Hosts infected (%)"
        , x="Distance from protected area boundary (m)"
        ,col = "Probability of vector to host transmission") + 
  scale_x_discrete(breaks=c(-5000,-2500,0,2500,5000), labels=c("-5000","-2500","0","2500","5000")) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line") 
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position =c(0.2,0.8)
        ,legend.title =element_text(size=5)
        ,strip.text.x = element_text(size = 8, face = "bold.italic" )
        ,strip.background=element_rect(fill="white",colour="white")
  ) +
  facet_wrap(~spp)
dev.off()
#*********************************************************************

br <- br1[br1$trans == 2000,]
co <- co1[co1$trans == 2000,]

mean(br$perc.host)
mean(co$perc.host)

brMax <- br[max(br$perc.host),]
brMin <- br[br$perc.host %in% min(br$perc.host),]
coMax <- co[max(co$perc.host),]
coMin <- co[min(co$perc.host),]

parmsDat <- rbind.data.frame(brMax,brMin,coMax,coMin)
write.csv(parmsDat,"parmsForTcideITC.csv")
