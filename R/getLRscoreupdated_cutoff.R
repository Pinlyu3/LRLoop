#' Calculate an LRscores_LRL-cutoff based on the distribution of a specified column in LRscores_randLRL and create density plots of the specified score column of LRscores_LRL and LRscores_randLRL.
#' @name getLRscoreupdated_cutoff
#' @param myscore_col The name of the column in LRscores_LRL and LRscores_randLRL based on which the cutoff is calculated.
#' @param p_value Specifies the desired fraction of the scores in the column myscore_col in LRscores_randLRL that are above the cutoff.
#' @param LRscores_LRL The LRloops-updated LRscores.
#' @param LRscores_randLRL The random LRloops-updated LRscores.
#' @return
#'  Prints the cutoff value and creates density plots of the specified score column of LRscores_LRL and LRscores_randLRL.
#' @import tidyverse ggplot2
#' @export


# Calculate an LRscores_LRL-cutoff based on the distribution of a specified column in LRscores_randLRL and create density plots of the specified score column of LRscores_LRL and LRscores_randLRL.
# Inputs:
# myscore_col: The name of the column in LRscores_LRL and LRscores_randLRL based on which the cutoff is calculated.
# p_value: Specifies the desired fraction of the scores in the column myscore_col in LRscores_randLRL that are above the cutoff.
# LRscores_LRL: The LRloops-updated LRscores.
# LRscores_randLRL: The random LRloops-updated LRscores.
# Output:
# Prints the cutoff value and creates density plots of the specified score column of LRscores_LRL and LRscores_randLRL.

getLRscoreupdated_cutoff <- function(p_value, LRscores_LRL, LRscores_randLRL, myscore_col) {
  
  myLRscores_LRL = as.numeric(LRscores_LRL[,myscore_col])
  myLRscores_randLRL = as.numeric(LRscores_randLRL[,myscore_col])
  
  Scores = cbind.data.frame(c(myLRscores_LRL, myLRscores_randLRL), c(rep('LRScores_LRL', length(myLRscores_LRL)), 
                                                                           rep('LRScores_randLRL', length(myLRscores_randLRL))))
  colnames(Scores) = c('Score', 'Mode')
  Scores = tibble(Scores)
  
  cutoff = quantile(myLRscores_randLRL, (1 - p_value))
  
  p = Scores%>%
    ggplot(aes(x=Score, fill=Mode)) +
    geom_density(alpha=0.3)+
    labs(x= "Score",
         subtitle="LR-Score distribution",
         caption="caption") +
    geom_vline(xintercept = cutoff, color = 'red', size = 1) +
    annotate(geom = 'text', x = cutoff, y = 0.5, color = 'red', label = sprintf("cutoff = %s (p-value: %s)", round(cutoff,3), p_value), hjust = 0.5)
  
  print(cutoff)
  
  return(p)
}