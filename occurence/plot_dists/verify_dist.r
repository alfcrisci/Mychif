
dists=list()
j=1
idlist=c(NULL)
fw <- try(fitdist(x,"weibull"))
if (!grepl("Error",fw)) {dists[[j]]=fw;idlist=c("weibull");j=j+1}

fg <- try(fitdist(x,"gamma"))
if (!grepl("Error",fg)) {dists[[j]]=fg;idlist=c(idlist,"gamma"); j=j+1 } 

fln <- try(fitdist(x,"lnorm"))
if (!grepl("Error",fln)) {dists[[j]]=fln;idlist=c(idlist,"lognormal");j=j+1 }

flnor <- try(fitdist(x,"norm"))
if (!grepl("Error",flnor)) {dists[[j]]=flnor;idlist=c(idlist,"normal") }

ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot") 
pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
out=ggarrange(ds, qq, cd,pp + rremove("x.text"), ncol = 2, nrow = 2)
options(warn=0)
ggsave(outfilepdf,plot=out,device="pdf")
ggsave(outfilepng,plot=out,device="png")