summary_large=function(x) {data.frame(c(stat.desc(as.numeric(x)),
                                        data.frame(skew=as.numeric(psych::describe(res_pooled[[10]]))[11],
                                                   kurtosis=as.numeric(psych::describe(x))[12])))
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

normFunc <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

