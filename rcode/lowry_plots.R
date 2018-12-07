#==============================================================
# Lowry plots https://www.r-bloggers.com/introducing-the-lowry-plot/
#==============================================================
require(ggplo2)
require(reshape2)

m_xylene_data <- data.frame(
  Parameter = c(
    "BW", "CRE", "DS", "KM", "MPY", "Pba", "Pfaa",
    "Plia", "Prpda", "Pspda", "QCC", "QfaC", "QliC",
    "QPC", "QspdC", "Rurine", "Vfac", "VliC", "Vmax"),
  "Main Effect" = c(
    1.03E-01, 9.91E-02, 9.18E-07, 3.42E-02, 9.27E-3, 2.82E-2, 2.58E-05,
    1.37E-05, 5.73E-4, 2.76E-3, 6.77E-3, 8.67E-05, 1.30E-02,
    1.19E-01, 4.75E-04, 5.25E-01, 2.07E-04, 1.73E-03, 1.08E-03),
  Interaction = c(
    1.49E-02, 1.43E-02, 1.25E-04, 6.84E-03, 3.25E-03, 7.67E-03, 8.34E-05,
    1.17E-04, 2.04E-04, 7.64E-04, 2.84E-03, 8.72E-05, 2.37E-03,
    2.61E-02, 6.68E-04, 4.57E-02, 1.32E-04, 6.96E-04, 6.55E-04
  )
)

#Function to prepare data for lowry plot

fortify_lowry_data <- function(data, 
                               param_var = "Parameter", 
                               main_var = "Main.Effect", 
                               inter_var = "Interaction") {
  #Convert wide to long format
  mdata <- melt(data, id.vars = param_var)
  
  #Order columns by main effect and reorder parameter levels
  o <- order(data[, main_var], decreasing = TRUE)
  data <- data[o, ]
  data[, param_var] <- factor(
    data[, param_var], levels = data[, param_var]
  )
  
  #Force main effect, interaction to be numeric
  data[, main_var] <- as.numeric(data[, main_var])
  data[, inter_var] <- as.numeric(data[, inter_var])
  
  #total effect is main effect + interaction
  data$.total.effect <- rowSums(data[, c(main_var, inter_var)])
  
  #Get cumulative totals for the ribbon
  data$.cumulative.main.effect <- cumsum(data[, main_var])
  data$.cumulative.total.effect <- cumsum(data$.total.effect)
  
  #A quirk of ggplot2 means we need x coords of bars
  data$.numeric.param <- as.numeric(data[, param_var])
  
  #The other upper bound
  #.maximum =  1 - main effects not included
  data$.maximum <- c(1 - rev(cumsum(rev(data[, main_var])))[-1], 1)
  
  data$.valid.ymax <- with(data, 
                          pmin(.maximum, .cumulative.total.effect)
  )
 
  mdata[, param_var] <- factor(
    mdata[, param_var], levels = data[, param_var]
  ) 
  list(data = data, mdata = mdata)
}

#Function to make lowry plot
lowry_plot <- function(data,
                       param_var = "Parameter",
                       main_var = "Main.Effect",
                       inter_var = "Interaction",
                       x_lab = "Parameters",
                       y_lab = "Total Effects (= Main Effects + Interactions)",
                       ribbon_alpha = 0.5,
                       x_text_angle = 25)
{
  #Fortify data and dump contents into plot function environment
  data_list <- fortify_lowry_data(data, param_var, main_var, inter_var)
  list2env(data_list, envir = sys.frame(sys.nframe()))
  
  p <- ggplot(data) +
    geom_bar(aes_string(x = param_var, y = "value", fill = "variable"),
             data = mdata, stat = "identity") +
    #scale_fill_manual(values = c("darkgrey","lightgrey")) +
    geom_ribbon(
      aes(x = .numeric.param, ymin = .cumulative.main.effect, ymax = .valid.ymax),
      data = data,
      alpha = ribbon_alpha) +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(axis.text.x = element_text(angle = x_text_angle, hjust = 1)) +
    scale_fill_grey(start = 0.5,end=0.2) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.direction = "horizontal") +
    guides(fill = guide_legend(reverse=TRUE))
  p
}
