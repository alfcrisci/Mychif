library(ggplot2)
library(dplyr)
library(ggExtra)
library(gridExtra)
library(readxl)
library(cowplot)
library(ggpubr)

#########################################################################################################################à
# require convert of imagemagick

setwd("")
setwd("/home/alf/Scrivania/lav_michyf/repositories/Mychif/co_occurence_plots")
#########################################################################################################################à
blue.plain.8.text <- element_text(face = "plain", color = "blue", angle=-75,size = 8)


fileex="Mycotoxins_co_occurence_fin_arranged.xls"
sheet_names=excel_sheets(fileex)

barley=c("barley DONs NIV","barley DONs ZEN NIV")
maize=c("maize DONs FBs","maize DONs FBs AF")
oat=c("oat DONs T2+HT2","oat DONs NIV","oat DONs T2+HT2 NIV")
wheat=c("wheat DONs NIV","wheat DONs ZEN","wheat DONs ZEN NIV")

wb_vars=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref","Mean_DONs","Mean_ZEN","Mean_NIV")

wheat_sheets=sapply(wheat,function(x) read_xls(fileex,x)) 
barley_sheets=sapply(barley,function(x) read_xls(fileex,x))
barley_sheets[[1]]$Mean_ZEN=0
wheat_sheets[[1]]$Mean_ZEN=0                                    
wheat_sheets[[2]]$Mean_NIV=0   
df_wheat=do.call("rbind",lapply(wheat_sheets,function(x) x[,wb_vars]))
df_barley=do.call("rbind",lapply(barley_sheets,function(x) x[,wb_vars]))

oat_vars=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref","Mean_DONs","Mean_T2+HT2","Mean_NIV")
oat_sheets=sapply(oat,function(x) read_xls(fileex,x))
oat_sheets[[1]]$Mean_NIV=0
oat_sheets[[2]]$`Mean_T2+HT2`=0
df_oat=do.call("rbind",lapply(oat_sheets,function(x) x[,oat_vars]))

mais_vars=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref","Mean_AF","Mean_DONs","Mean_FBs")
maize_sheets=sapply(maize,function(x) read_xls(fileex,x))
maize_sheets[[1]]$Mean_AF=0
df_maize=do.call("rbind",lapply(maize_sheets,function(x) x[,mais_vars]))

################################################################################################################################

df_wheat_gg=reshape2::melt(df_wheat,id=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref"))
df_barley_gg=reshape2::melt(df_barley,id=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref"))
df_oat_gg=reshape2::melt(df_oat,id=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref"))
df_maize_gg=reshape2::melt(df_maize,id=c("plants","mycotoxins","freq_100","ndata_valid","nrecords","Ref"))

#######################################################################################################################################
# wheat

df_wheat_gg$value=as.numeric(df_wheat_gg$value)
df_wheat_gg$Mycotoxin=df_wheat_gg$variable
df_wheat_gg$Co_occurence=df_wheat_gg$mycotoxins
df_wheat_gg$References=gsub(" et al.,",",",gsub(".... et al, ",",",df_wheat_gg$Ref))

id=df_wheat_gg$mycotoxins[as.numeric(sapply(levels(as.factor(df_wheat_gg$References)),function(x) match(x,df_wheat_gg$References)))]

df_wheat_gg$freq_100=df_wheat_gg$freq_100/3
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
p=ggplot(df_wheat_gg) + geom_col(aes(x = References, y = value, fill = Mycotoxin), position = "stack")+ theme(legend.position="bottom")+ scale_x_discrete(labels=id)+theme(axis.text.x =blue.plain.8.text,legend.text=element_text(size=8))+
  labs(title = "Wheat Co-occurence mycotoxins",y="Concentrations",x="")
g <-ggplot(df_wheat_gg, aes(x=References,y=freq_100, fill = Co_occurence)) + geom_bar(stat="identity")+ theme(legend.position="top")+theme(axis.text.x = element_blank(),legend.text=element_text(size=8))+scale_fill_manual(values=my3cols)+labs(y="Data Freq (%)",x="References")
#a=grid.arrange(p, g, ncol = 1,heights=c(5,3))

pdf("wheat.pdf")
plot(p)
dev.off()


png("wheat.png")
plot(p)
dev.off()


#######################################################################################################################################



#######################################################################################################################################
# oat

df_oat_gg$value=as.numeric(df_oat_gg$value)
df_oat_gg$Mycotoxin=df_oat_gg$variable
df_oat_gg$Co_occurence=df_oat_gg$mycotoxins
df_oat_gg$References=gsub(" et al.,",",",gsub(".... et al, ",",",df_oat_gg$Ref))

id=df_oat_gg$mycotoxins[as.numeric(sapply(levels(as.factor(df_oat_gg$References)),function(x) match(x,df_oat_gg$References)))]

df_oat_gg$freq_100=df_oat_gg$freq_100/3
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
p=ggplot(df_oat_gg) + geom_col(aes(x = References, y = value, fill = Mycotoxin), position = "stack")+ theme(legend.position="bottom")+ scale_x_discrete(labels=id)+theme(axis.text.x =blue.plain.8.text,legend.text=element_text(size=8))+
  labs(title = "Oat Co-occurence mycotoxins",y="Concentrations",x="")
g <-ggplot(df_oat_gg, aes(x=References,y=freq_100, fill = Co_occurence)) + geom_bar(stat="identity")+ theme(legend.position="top")+theme(axis.text.x = element_blank(),legend.text=element_text(size=8))+scale_fill_manual(values=my3cols)+labs(y="Data Freq (%)",x="References")
a=grid.arrange(p, g, ncol = 1,heights=c(5,3))

pdf("oat.pdf")
plot(p)
dev.off()
png("oat.png")
plot(p)
dev.off()
#######################################################################################################################################
#######################################################################################################################################
# barley

df_barley_gg$value=as.numeric(df_barley_gg$value)
df_barley_gg$Mycotoxin=df_barley_gg$variable
df_barley_gg$Co_occurence=df_barley_gg$mycotoxins
df_barley_gg$References=gsub(" et al.,",",",gsub(".... et al, ",",",df_barley_gg$Ref))

id=df_barley_gg$mycotoxins[as.numeric(sapply(levels(as.factor(df_barley_gg$References)),function(x) match(x,df_barley_gg$References)))]

df_barley_gg$freq_100=df_barley_gg$freq_100/3
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
p=ggplot(df_barley_gg) + geom_col(aes(x = References, y = value, fill = Mycotoxin), position = "stack")+ theme(legend.position="bottom")+ scale_x_discrete(labels=id)+theme(axis.text.x =blue.plain.8.text,legend.text=element_text(size=8))+
  labs(title = "Barley Co-occurence mycotoxins",y="Concentrations",x="")
g <-ggplot(df_barley_gg, aes(x=References,y=freq_100, fill = Co_occurence)) + geom_bar(stat="identity")+ theme(legend.position="top")+theme(axis.text.x = element_blank(),legend.text=element_text(size=8))+scale_fill_manual(values=my3cols)+labs(y="Data Freq (%)",x="References")
a=grid.arrange(p, g, ncol = 1,heights=c(5,3))

pdf("barley.pdf")
plot(p)
dev.off()
png("barley.png")
plot(p)
dev.off()
#######################################################################################################################################
#######################################################################################################################################
# maize

df_maize_gg$value=as.numeric(df_maize_gg$value)
df_maize_gg$Mycotoxin=df_maize_gg$variable
df_maize_gg$Co_occurence=df_maize_gg$mycotoxins
df_maize_gg$References=gsub(" et al.,",",",gsub(".... et al, ",",",df_maize_gg$Ref))

id=df_maize_gg$mycotoxins[as.numeric(sapply(levels(as.factor(df_maize_gg$References)),function(x) match(x,df_maize_gg$References)))]

df_maize_gg$freq_100=df_maize_gg$freq_100/3
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
p=ggplot(df_maize_gg) + geom_col(aes(x = References, y = value, fill = Mycotoxin), position = "stack")+ theme(legend.position="bottom")+ scale_x_discrete(labels=id)+theme(axis.text.x =blue.plain.8.text,legend.text=element_text(size=8))+
  labs(title = "Maize Co-occurence mycotoxins",y="Concentrations",x="")

g <-ggplot(df_maize_gg, aes(x=References,y=freq_100, fill = Co_occurence)) + geom_bar(stat="identity")+ theme(legend.position="top")+theme(axis.text.x = element_blank(),legend.text=element_text(size=8))+scale_fill_manual(values=my3cols)+labs(y="Data Freq (%)",x="References")
a=grid.arrange(p, g, ncol = 1,heights=c(5,3))

pdf("maize.pdf")
plot(p)
dev.off()
png("maize.png")
plot(p)
dev.off()


#######################################################################################################################################


# @references

# https://deanattali.com/2015/03/29/ggExtra-r-package/
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html#51_arrange_multiple_graphs_on_the_same_page
# http://www.lreding.com/nonstandard_deviations/2017/08/19/cowmarg/
