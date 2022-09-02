library(tictoc)
library(rdtest)
library(gt)
library(ggplot2)
library(dplyr)
library(truncnorm)

test_that("the vignette process for the first loop does not change",{

  int.ns <- 3

  fun.runnerJoint = function(str.spec)
  {
    #tic(str.spec)
    time.start <- Sys.time()
    eval(parse(text = paste0(paste0("list.summary <- returnResultsMCJointTest(",str.spec),")")))
    time.end <- Sys.time()
    #toc()
    return(list.summary)
  }

  # Specification list
  # values for x jumps
  vec.xJump <- c(0) #,0.15,0.3)
  # values for z jumps
  vec.zJump <- c(0) #seq(0,2,length=5)
  # sample sizes
  vec.n <- c(500)
  vec.covZ <- c(0.5)
  vec.dimZ <- c(3)
  # the last one, half, full:
  ## half: for 3, 2:3 jumps.
  vec.frac.jump <- list("'None'") #list("'Last'","'Half'","'All'")

  int.numberSpec <- length(vec.xJump)*length(vec.n)*length(vec.zJump)*length(vec.covZ)*length(vec.dimZ)*length(vec.frac.jump)#*3

  vec.effNzReportSpec <- c(1:int.numberSpec)
  vec.effNxReportSpec <- c(1:int.numberSpec)
  vec.corZReportSpec <- c(1:int.numberSpec)
  vec.dimVecReportSpec <- c(1:int.numberSpec)
  vec.nReportSpec <- c(1:int.numberSpec)
  vec.xjumpReportSpec <- c(1:int.numberSpec)
  vec.zjumpReportSpec <- c(1:int.numberSpec)
  vec.naiveReportResult <- c(1:int.numberSpec)
  vec.bonfeReportResult <- c(1:int.numberSpec)
  vec.jointReportResult <- c(1:int.numberSpec)
  vec.maxTestSpec <- c(1:int.numberSpec)
  vec.avgMaxTestCritical <- c(1:int.numberSpec)
  vec.medianMaxTestCritical <- c(1:int.numberSpec)
  vec.jumpTestSpec <- c(1:int.numberSpec)

  list.speclist <- list()

  counter <- 1
  df.resultJointTest <- data.frame()
  fun.repeaterJoint <- function(n,cov.z,dim,x.jump,z.jump,bool.maxTest,bool.maxTestInv,bool.chisqStd = FALSE,counter,df.resultJointTest,frac.jump,spec) {
    list.sum <- fun.runnerJoint(str.spec=spec)

    real.effNzReportSpec <- round(list.sum$effN.mean,1)
    real.effNxReportSpec <- round(list.sum$effNx.mean,1)
    real.corZReportSpec <- round(list.sum$cor.z,2)
    real.dimVecReportSpec <- dim
    real.nReportSpec <- n
    real.xjumpReportSpec <- x.jump
    real.zjumpReportSpec <- z.jump

    real.naiveReportResult <- round(list.sum$naive,3)
    real.bonfeReportResult <- round(list.sum$bonfe,3)
    real.jointReportResult <- round(list.sum$joint,3)
    real.maxTestSpec <- bool.maxTest
    if (bool.maxTest) {
      real.avgMaxTestCritical <- round(mean(list.sum$maxCritical),3)
      real.medianMaxTestCritical <- round(median(list.sum$maxCritical),3)
    } else {
      real.avgMaxTestCritical <- NA
      real.medianMaxTestCritical <- NA
    }

    list.speclist[counter] <- spec
    #print(spec)

    df.resultJointTest <- rbind(df.resultJointTest,
                                data.frame(vec.effNzReportSpec = real.effNzReportSpec,
                                           vec.effNxReportSpec = real.effNxReportSpec,
                                           vec.corZReportSpec = cov.z,
                                           vec.dimVecReportSpec = real.dimVecReportSpec,
                                           vec.nReportSpec = real.nReportSpec,
                                           vec.xjumpReportSpec = real.xjumpReportSpec,
                                           vec.zjumpReportSpec = real.zjumpReportSpec,
                                           vec.naiveReportResult = real.naiveReportResult,
                                           vec.bonfeReportResult = real.bonfeReportResult,
                                           vec.jointReportResult = real.jointReportResult,
                                           vec.maxTestSpec = real.maxTestSpec,
                                           vec.maxTestInvSpec = bool.maxTestInv,
                                           vec.chisqStd = bool.chisqStd,
                                           vec.avgMaxTestCritical = real.avgMaxTestCritical,
                                           vec.medianMaxTestCritical = real.medianMaxTestCritical,
                                           vec.jumpTestSpec = frac.jump))
    return(df.resultJointTest)
  }

  for (n in vec.n) {
    for (dim in vec.dimZ) {
      if (dim == 1) {
        vec.covZUsed <- c(0)
      } else {
        vec.covZUsed <- vec.covZ
      }
      for (cov.z in vec.covZUsed) {
        for (frac.jump in vec.frac.jump) {
          cat("Running for n=",n,",dim=",dim,",cov.z=",cov.z,",frac=", frac.jump,".\n")
          # bool.maxTest TRUE, bool.maxTestVinv FALSE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",TRUE,",bool.maxTestVinv=",FALSE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = TRUE,bool.maxTestInv=FALSE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec)
                  counter <- counter + 1
              }
            }
          }
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",TRUE,",bool.maxTestVinv=",TRUE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = TRUE,bool.maxTestInv=TRUE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool.maxTest FALSE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",FALSE,",bool.maxTestVinv=",FALSE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = FALSE,bool.maxTestInv=FALSE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool.maxTest FALSE, bool.chisqStd TRUE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",FALSE,",bool.maxTestVinv=",FALSE,",bool.chisqStd=",TRUE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = FALSE,bool.maxTestInv=FALSE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec,bool.chisqStd=TRUE)
                  counter <- counter + 1
              }
            }
          }
        }
      }

    }
  }
  df.coverage <- data.frame(n=df.resultJointTest[,5],cov=df.resultJointTest[,3],dim=df.resultJointTest[,4],p.naive=df.resultJointTest[,8],p.bonferroni=df.resultJointTest[,9],p.joint=df.resultJointTest[,10],isMaxTest=df.resultJointTest[,11],isMaxInv=df.resultJointTest[,12],isChiStd=df.resultJointTest[,13])
  testthat::expect_snapshot_output(df.coverage %>% gt())

  vec.xJump <- c(0.15)
  vec.zJump <- seq(0,2,length=5)
  vec.n <- c(500)
  vec.covZ <- c(0.5)
  vec.dimZ <- c(3)
  vec.frac.jump <- list("'Last'","'Half'","'All'")
  int.numberSpec <- length(vec.xJump)*length(vec.n)*length(vec.zJump)*length(vec.covZ)*length(vec.dimZ)*2*length(vec.frac.jump)

  vec.effNzReportSpec <- c(1:int.numberSpec)
  vec.effNxReportSpec <- c(1:int.numberSpec)
  vec.corZReportSpec <- c(1:int.numberSpec)
  vec.dimVecReportSpec <- c(1:int.numberSpec)
  vec.nReportSpec <- c(1:int.numberSpec)
  vec.xjumpReportSpec <- c(1:int.numberSpec)
  vec.zjumpReportSpec <- c(1:int.numberSpec)
  vec.naiveReportResult <- c(1:int.numberSpec)
  vec.bonfeReportResult <- c(1:int.numberSpec)
  vec.jointReportResult <- c(1:int.numberSpec)
  vec.maxTestSpec <- c(1:int.numberSpec)
  vec.avgMaxTestCritical <- c(1:int.numberSpec)
  vec.medianMaxTestCritical <- c(1:int.numberSpec)
  vec.jumpTestSpec <- c(1:int.numberSpec)

  list.speclist <- list()

  counter <- 1
  df.resultJointTest <- data.frame()
  fun.repeaterJoint <- function(n,cov.z,dim,x.jump,z.jump,bool.maxTest,bool.chisqStd = FALSE,counter,df.resultJointTest,frac.jump,spec) {
    list.sum <- fun.runnerJoint(str.spec=spec)

    real.effNzReportSpec <- round(list.sum$effN.mean,1)
    real.effNxReportSpec <- round(list.sum$effNx.mean,1)
    real.corZReportSpec <- round(list.sum$cor.z,2)
    real.dimVecReportSpec <- dim
    real.nReportSpec <- n
    real.xjumpReportSpec <- x.jump
    real.zjumpReportSpec <- z.jump

    real.naiveReportResult <- round(list.sum$naive,3)
    real.bonfeReportResult <- round(list.sum$bonfe,3)
    real.jointReportResult <- round(list.sum$joint,3)
    real.maxTestSpec <- bool.maxTest
    if (bool.maxTest) {
      real.avgMaxTestCritical <- round(mean(list.sum$maxCritical),3)
      real.medianMaxTestCritical <- round(median(list.sum$maxCritical),3)
    } else {
      real.avgMaxTestCritical <- NA
      real.medianMaxTestCritical <- NA
    }

    list.speclist[counter] <- spec
    #print(spec)

    df.resultJointTest <- rbind(df.resultJointTest,
                                data.frame(vec.effNzReportSpec = real.effNzReportSpec,
                                           vec.effNxReportSpec = real.effNxReportSpec,
                                           vec.corZReportSpec = cov.z,
                                           vec.dimVecReportSpec = real.dimVecReportSpec,
                                           vec.nReportSpec = real.nReportSpec,
                                           vec.xjumpReportSpec = real.xjumpReportSpec,
                                           vec.zjumpReportSpec = real.zjumpReportSpec,
                                           vec.naiveReportResult = real.naiveReportResult,
                                           vec.bonfeReportResult = real.bonfeReportResult,
                                           vec.jointReportResult = real.jointReportResult,
                                           vec.maxTestSpec = real.maxTestSpec,
                                           vec.chisqStd = bool.chisqStd,
                                           vec.avgMaxTestCritical = real.avgMaxTestCritical,
                                           vec.medianMaxTestCritical = real.medianMaxTestCritical,
                                           vec.jumpTestSpec = frac.jump))
    return(df.resultJointTest)
  }

  for (n in vec.n) {
    for (dim in vec.dimZ) {
      for (cov.z in vec.covZ) {
        for (frac.jump in vec.frac.jump) {
          cat("Running for n=",n,",dim=",dim,",cov.z=",cov.z,",frac=", frac.jump,".\n")
          # bool.maxTest TRUE, bool.maxTestVinv FALSE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",TRUE,",bool.maxTestVinv=",FALSE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = TRUE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool.maxTest FALSE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",FALSE,",bool.maxTestVinv=",FALSE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = FALSE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool.maxTest FALSE, bool.chisqStd TRUE
          for (x.jump in vec.xJump) {
            for (z.jump in vec.zJump) {
              if ((x.jump > 0.000001) | (x.jump == 0 & z.jump == 0)) {
                  spec <- paste0("int.ns=",int.ns,",a.2=",0,",x.jump=",x.jump,",jump=",z.jump,",dim=",dim,",n=",n,",cov.z =", cov.z,",frac.jump=",frac.jump,",bool.mutePrint=TRUE,bool.maxTest=",FALSE,",bool.maxTestVinv=",FALSE,",bool.chisqStd=",TRUE)
                  df.resultJointTest <- fun.repeaterJoint(n = n,cov.z = cov.z,dim = dim,x.jump = x.jump,z.jump = z.jump,bool.maxTest = FALSE,counter = counter,df.resultJointTest = df.resultJointTest,frac.jump = frac.jump,spec = spec,bool.chisqStd=TRUE)
                  counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
  testthat::expect_snapshot_output(df.resultJointTest[,c(3,4,6:11)])

  df.resultNull <- df.resultJointTest[df.resultJointTest$vec.xjumpReportSpec == 0,]
  testthat::expect_snapshot_output(df.resultNull[,c(3,4,8:11,14)])

  for (n in vec.n)
  {
    for (dim in unique(df.resultJointTest$vec.dimVecReportSpec))
    {
      for (cor in vec.covZ)
      {
        for (xval in unique(df.resultJointTest$vec.xjumpReportSpec))
        {
          for (jumpTest in list("'Last'","'Half'","'All'")) {
            df.resultEval <- df.resultJointTest[df.resultJointTest$vec.dimVecReportSpec==dim
                                                & df.resultJointTest$vec.corZReportSpec == cor
                                                & df.resultJointTest$vec.nReportSpec== n
                                                & df.resultJointTest$vec.xjumpReportSpec == xval
                                                & df.resultJointTest$vec.jumpTestSpec == jumpTest
                                                & df.resultJointTest$vec.maxTestSpec == TRUE,]
            df.resultEvalMat <- rbind(data.frame(z=df.resultEval$vec.zjumpReportSpec,prob=df.resultEval$vec.jointReportResult,type="joint max"),data.frame(z=df.resultEval$vec.zjumpReportSpec,prob=df.resultEval$vec.bonfeReportResult,type="bonferroni"))

            df.resultEvalMat <- rbind(df.resultEvalMat,data.frame(z=df.resultEval$vec.zjumpReportSpec,prob=df.resultEval$vec.bonfeReportResult,type="bonferroni"))
            df.resultEvalChiStd <- df.resultJointTest[df.resultJointTest$vec.dimVecReportSpec==dim
                                                      & df.resultJointTest$vec.corZReportSpec == cor
                                                      & df.resultJointTest$vec.nReportSpec== n
                                                      & df.resultJointTest$vec.xjumpReportSpec == xval
                                                      & df.resultJointTest$vec.jumpTestSpec == jumpTest
                                                      & df.resultJointTest$vec.maxTestSpec == FALSE
                                                      & df.resultJointTest$vec.chisqStd == TRUE,]
            df.resultEvalMat <- rbind(df.resultEvalMat,data.frame(z=df.resultEvalChiStd$vec.zjumpReportSpec,prob=df.resultEvalChiStd$vec.jointReportResult,type="joint chi Std"))

            p_2 <- ggplot(data = df.resultEvalMat, mapping = aes(x = z, y = prob,color=type)) + geom_point() +
              scale_y_continuous(breaks=c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
              geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
              geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
              geom_hline(yintercept = 1, linetype = "dashed", size = 0.01) + ggtitle(paste0("x jump: ", xval,"; n: ", n,", dim: ",dim,", cor:",cor,", alt Hyp:",jumpTest))
            testthat::expect_snapshot_output(p_2)
          }
        }
      }
    }
  }
})
