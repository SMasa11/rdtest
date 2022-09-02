####
#' return chi square test results of indepdendence
#'
#'
#' @param df.check data.frame of data to evaluate
#' @param int.numPartition  integer dimension of partitions
#'
#' @export
####

returnResultsChiIndependence <- function(df.check,int.numPartition)
{
  mat.freqByCell<- matrix(rep(0,int.numPartition*int.numPartition),int.numPartition,int.numPartition)
  for (i in 0:(int.numPartition-1))
  {
    for (j in 0:(int.numPartition-1))
    {
      mat.freqByCell[i+1,j+1]  <- sum((df.check$X <= qnorm((i+1)/int.numPartition)) &
                              (df.check$X > qnorm(i/int.numPartition)) &
                              (df.check$Z <= qnorm((j+1)/int.numPartition)) &
                              (df.check$Z > qnorm(j/int.numPartition)))
    }
  }
  real.statChiSqIndep <- chisq.test(mat.freqByCell)
  return(real.statChiSqIndep)
}
