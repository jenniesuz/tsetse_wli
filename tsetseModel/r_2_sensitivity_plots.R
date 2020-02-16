library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(scales)
library(sensitivity)
library(reshape)

randSnd <- read.csv("tsetseSensitivityResults.csv")

dat <- melt(randSnd,measure.vars=c("larvi" 
                                    ,"emerg"
                                    ,"pden"   
                                    ,"pmort"   
                                    ,"amort"
                                    ,"difcoeff"   
                                    ,"mo.out" ))



dat$variable2 <- factor(dat$variable
                       , levels = c("larvi" 
                                    ,"emerg"
                                    ,"pden"   
                                    ,"pmort"   
                                    ,"amort"
                                    ,"difcoeff"   
                                    ,"mo.out")
                       , ordered = TRUE
                       , labels=c(expression(paste(italic("l")))
                                  ,expression(paste(italic(beta)))
                                  ,expression(paste(italic(delta)))
                                  ,expression(paste(italic(mu["P"])))
                                  ,expression(paste(italic(mu["B"])))
                                  ,expression(paste(italic(alpha)))
                                  ,expression(paste(mu["F"]))
                       )
)



tiff("FigS3.2_sensitivity.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(dat,aes(value,flypop))+
  geom_point(size=0.01,col="black",alpha=0.2)  +
  ylim(0,40) +
  labs( y= "Relative tsetse density (%)"
        , x="Parameter value") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.5)
        ,legend.title = element_text(size=9)
        ,strip.background = element_rect(fill="white",color="white")
  ) +
  facet_wrap(~variable2,scales="free_x",labeller = label_parsed)
dev.off()



#*************PRCC***********
prccNd <- pcc(randSnd[,2:8],randSnd[,9],rank=T,nboot=1000)
prccNdvals <- prccNd$PRCC
names(prccNdvals) <- c("original","bias","st.err","minci","maxci")


tiff("FigS3.3_prcc.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(prccNdvals, aes(x=as.factor(1:length(prccNdvals[,1])),y=original)) +
  geom_point() +
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1) +
  labs( x="Parameter"
        ,y="PRCC"
        , title=" ") +
  scale_x_discrete(breaks=c("1","2","3","4","5","6","7")
                   ,labels=c(expression(paste(italic("l")))
                       ,expression(paste(italic(beta)))
                       ,expression(paste(italic(delta)))
                       ,expression(paste(italic(mu["P"])))
                       ,expression(paste(italic(mu["B"])))
                       ,expression(paste(italic(alpha)))
                       ,expression(paste(mu["F"]))
                   )
  )+
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=9)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.5)
        ,legend.title = element_text(size=9)
        ,strip.background = element_rect(fill="white",color="white")
  )

dev.off()

