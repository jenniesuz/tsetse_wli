library(sensitivity)
library(ggplot2)

inputs <- read.csv("LHSParms.csv")
output <- read.csv("trypanosomeSensitivityResults.csv")
data <- cbind.data.frame(inputs,output)

dat <- melt(data,measure.vars=c("ptrans.ht" 
                                   ,"ptrans.hg"
                                   ,"vinc"   
                                   ,"wane"   
                                   ,"ptrans.vh"
                                   ,"hinc"   
                                   ,"hrec" 
                                   ,"feed"
                               ))



dat$variable2 <- factor(dat$variable
                        , levels = c("ptrans.ht" 
                                     ,"ptrans.hg"
                                     ,"vinc"   
                                     ,"wane"   
                                     ,"ptrans.vh"
                                     ,"hinc"   
                                     ,"hrec" 
                                     ,"feed"
                        )
                        , ordered = TRUE
                        , labels=c(expression(paste(italic(p)[S]))
                                   ,expression(paste(italic(p)[G]))
                                   ,expression(paste(italic(sigma)[v]))
                                   ,expression(paste(italic(gamma)))
                                   ,expression(paste(italic(p)[H]))
                                   ,expression(paste(italic(sigma)[H]))
                                   ,expression(paste(italic(phi)))
                                   ,expression(paste(italic(alpha)))
                        )
)




#**************************Scatter plots**********************************
tiff("FigS5.1_sensitivity.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(dat,aes(value,perc.host))+
  geom_point(size=0.01,col="black",alpha=0.2)  +
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

tiff("FigS5.2_sensitivity.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(dat,aes(value,perc.vector))+
  geom_point(size=0.01,col="black",alpha=0.2)  +
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
#***************************PRCC Plots***********************************************
#************host**********
prccNd <- pcc(inputs[,2:9],output$perc.host,rank=T,nboot=100)
prccNdvals <- prccNd$PRCC
names(prccNdvals) <- c("original","bias","st.err","minci","maxci")


tiff("FigS5.3_hostPRCC.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(prccNdvals, aes(x=as.factor(1:length(prccNdvals[,1])),y=original)) +
  geom_point() +
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1) +
  labs( x="Parameter"
        ,y="PRCC"
        , title=" ") +
  scale_x_discrete(breaks=c("1","2","3","4","5","6","7","8")
                   ,labels=c(expression(paste(italic(p)[S]))
                             ,expression(paste(italic(p)[G]))
                             ,expression(paste(italic(sigma)[v]))
                             ,expression(paste(italic(gamma)))
                             ,expression(paste(italic(p)[H]))
                             ,expression(paste(italic(sigma)[H]))
                             ,expression(paste(italic(phi)))
                             ,expression(paste(italic(alpha)))
                            ) ) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=11)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=11)
        ,legend.key.size = unit(1,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=11)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )
dev.off()

#***********vector************
prccNd <- pcc(inputs[,2:9],output$perc.vector,rank=T,nboot=100)
prccNdvals <- prccNd$PRCC
names(prccNdvals) <- c("original","bias","st.err","minci","maxci")


tiff("FigS5.4_vectorPRCC.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(prccNdvals, aes(x=as.factor(1:length(prccNdvals[,1])),y=original)) +
  geom_point() +
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1) +
  labs( x="Parameter"
        ,y="PRCC"
        , title=" ") +
  scale_x_discrete(breaks=c("1","2","3","4","5","6","7","8")
                   ,labels=c(expression(paste(italic(p)[S]))
                             ,expression(paste(italic(p)[G]))
                             ,expression(paste(italic(sigma)[v]))
                             ,expression(paste(italic(gamma)))
                             ,expression(paste(italic(p)[H]))
                             ,expression(paste(italic(sigma)[H]))
                             ,expression(paste(italic(phi)))
                             ,expression(paste(italic(alpha)))
                   ) ) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=11)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=11)
        ,legend.key.size = unit(1,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=11)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )
dev.off()

#********************************************************************************************

# choose only those outputs with host prevalence between X and vector prevalence between X



outputSubsTb <- data[(data$perc.vector >= 0.09)&(data$perc.vector <= 1.57 ),] # 138
outputSubsTc <- data[(data$perc.vector >= 4.79)&(data$perc.vector <= 5.94 ),] # 150

write.csv(outputSubsTb, "tbSens.csv")
write.csv(outputSubsTc,"tcSens.csv")
