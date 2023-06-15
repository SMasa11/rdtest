library(tictoc)
library(rdtest)
library(gt)
library(ggplot2)
library(dplyr)
library(truncnorm)

test_that("the vignette process for the first loop does not change, dim = 3, null coverage",{

  int_ns <- 10

  fun_runnerJoint = function(str_spec)
  {
    #tic(str_spec)
    time_start <- Sys.time()
    eval(parse(
      text = paste0(
        paste0("list_summary <- return_result_MC_joint(",str_spec),")"
        )))
    time_end <- Sys.time()
    #toc()
    return(list_summary)
  }

  # Specification list
  # values for x jumps
  vec_xJump <- c(0) #,0.15,0.3)
  # values for z jumps
  vec_zJump <- c(0) #seq(0,2,length=5)
  # sample sizes
  vec_n <- c(500)
  vec_covZ <- c(0.5)
  vec_dimZ <- c(3)
  # the last one, half, full:
  ## half: for 3, 2:3 jumps.
  vec_frac_jump <- list("'None'")

  int_numberSpec <-
    length(vec_xJump) *
    length(vec_n) *
    length(vec_zJump) *
    length(vec_covZ) *
    length(vec_dimZ) *
    length(vec_frac_jump)

  vec_effNzReportSpec <- c(1:int_numberSpec)
  vec_effNxReportSpec <- c(1:int_numberSpec)
  vec_corZReportSpec <- c(1:int_numberSpec)
  vec_dimVecReportSpec <- c(1:int_numberSpec)
  vec_nReportSpec <- c(1:int_numberSpec)
  vec_xjumpReportSpec <- c(1:int_numberSpec)
  vec_zjumpReportSpec <- c(1:int_numberSpec)
  vec_naiveReportResult <- c(1:int_numberSpec)
  vec_bonfeReportResult <- c(1:int_numberSpec)
  vec_jointReportResult <- c(1:int_numberSpec)
  vec_maxTestSpec <- c(1:int_numberSpec)
  vec_avgMaxTestCritical <- c(1:int_numberSpec)
  vec_medianMaxTestCritical <- c(1:int_numberSpec)
  vec_jumpTestSpec <- c(1:int_numberSpec)

  list_speclist <- list()

  counter <- 1
  df_resultJointTest <- data.frame()
  fun_repeaterJoint <- function(n,
                                cov_z,
                                dim,
                                x_jump,
                                z_jump,
                                bool_max_test,
                                bool_L2_std = FALSE,
                                counter,
                                df_resultJointTest,
                                frac_jump,
                                spec) {
    list_sum <- fun_runnerJoint(str_spec=spec)

    real_effNzReportSpec <- round(list_sum$effN.mean,1)
    real_effNxReportSpec <- round(list_sum$effNx.mean,1)
    real_corZReportSpec <- round(list_sum$cor.z,2)
    real_dimVecReportSpec <- dim
    real_nReportSpec <- n
    real_xjumpReportSpec <- x_jump
    real_zjumpReportSpec <- z_jump

    real_naiveReportResult <- round(list_sum$naive,3)
    real_bonfeReportResult <- round(list_sum$bonfe,3)
    real_jointReportResult <- round(list_sum$joint,3)
    real_maxTestSpec <- bool_max_test
    if (bool_max_test) {
      real_avgMaxTestCritical <- round(mean(list_sum$maxCritical),3)
      real_medianMaxTestCritical <- round(median(list_sum$maxCritical),3)
    } else {
      real_avgMaxTestCritical <- NA
      real_medianMaxTestCritical <- NA
    }

    list_speclist[counter] <- spec
    #print(spec)

    df_resultJointTest <-
      rbind(df_resultJointTest,
            data.frame(vec_effNzReportSpec = real_effNzReportSpec,
                       vec_effNxReportSpec = real_effNxReportSpec,
                       vec_corZReportSpec = cov_z,
                       vec_dimVecReportSpec = real_dimVecReportSpec,
                       vec_nReportSpec = real_nReportSpec,
                       vec_xjumpReportSpec = real_xjumpReportSpec,
                       vec_zjumpReportSpec = real_zjumpReportSpec,
                       vec_naiveReportResult = real_naiveReportResult,
                       vec_bonfeReportResult = real_bonfeReportResult,
                       vec_jointReportResult = real_jointReportResult,
                       vec_maxTestSpec = real_maxTestSpec,
                       vec_L2_std = bool_L2_std,
                       vec_avgMaxTestCritical = real_avgMaxTestCritical,
                       vec_medianMaxTestCritical = real_medianMaxTestCritical,
                       vec_jumpTestSpec = frac_jump))
    return(df_resultJointTest)
  }

  for (n in vec_n) {
    for (dim in vec_dimZ) {
      if (dim == 1) {
        vec_covZUsed <- c(0)
      } else {
        vec_covZUsed <- vec_covZ
      }
      for (cov_z in vec_covZUsed) {
        for (frac_jump in vec_frac_jump) {
          cat("Running for n=",n,
              ",dim=",dim,
              ",cov_z=",cov_z,
              ",frac=", frac_jump,".\n")
          # bool_maxTest TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",TRUE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,
                                      x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = TRUE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",FALSE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,
                                      x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = FALSE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE, bool_L2_std TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",FALSE,
                                 ",bool_L2_std=",TRUE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,
                                      x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = FALSE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec,
                                      bool_L2_std=TRUE)
                  counter <- counter + 1
              }
            }
          }
        }
      }

    }
  }
  df_coverage <- data.frame(n = df_resultJointTest[,5],
                            cov = df_resultJointTest[,3],
                            dim = df_resultJointTest[,4],
                            p.naive = df_resultJointTest[,8],
                            p.bonferroni = df_resultJointTest[,9],
                            p.joint = df_resultJointTest[,10],
                            isMaxTest = df_resultJointTest[,11],
                            isChiStd = df_resultJointTest[,12])
  testthat::expect_equal(df_coverage$p.naive,c(0.2,0.2,0.2))
  testthat::expect_equal(df_coverage$p.bonferroni,c(0,0,0))
  testthat::expect_equal(df_coverage$p.joint,c(0,0,0))
})

test_that("the vignette process for the first loop does not change, dim = 3, alternative power",{
  int_ns <- 10

  fun_runnerJoint = function(str_spec)
  {
    #tic(str_spec)
    time_start <- Sys.time()
    eval(parse(
      text = paste0(
        paste0("list_summary <- return_result_MC_joint(",str_spec),")"
      )))
    time_end <- Sys.time()
    #toc()
    return(list_summary)
  }

  vec_xJump <- c(0.15)
  vec_zJump <- seq(0,2,length=5)
  vec_n <- c(500)
  vec_covZ <- c(0.5)
  vec_dimZ <- c(3)
  vec_frac_jump <- list("'Last'","'Half'","'All'")
  int_numberSpec <-
    length(vec_xJump) *
    length(vec_n) *
    length(vec_zJump) *
    length(vec_covZ) *
    length(vec_dimZ) *
    2 *
    length(vec_frac_jump)

  vec_effNzReportSpec <- c(1:int_numberSpec)
  vec_effNxReportSpec <- c(1:int_numberSpec)
  vec_corZReportSpec <- c(1:int_numberSpec)
  vec_dimVecReportSpec <- c(1:int_numberSpec)
  vec_nReportSpec <- c(1:int_numberSpec)
  vec_xjumpReportSpec <- c(1:int_numberSpec)
  vec_zjumpReportSpec <- c(1:int_numberSpec)
  vec_naiveReportResult <- c(1:int_numberSpec)
  vec_bonfeReportResult <- c(1:int_numberSpec)
  vec_jointReportResult <- c(1:int_numberSpec)
  vec_maxTestSpec <- c(1:int_numberSpec)
  vec_avgMaxTestCritical <- c(1:int_numberSpec)
  vec_medianMaxTestCritical <- c(1:int_numberSpec)
  vec_jumpTestSpec <- c(1:int_numberSpec)

  list_speclist <- list()

  counter <- 1
  df_resultJointTest <- data.frame()
  fun_repeaterJoint <- function(n,
                                cov_z,
                                dim,
                                x_jump,
                                z_jump,
                                bool_max_test,
                                bool_L2_std = FALSE,
                                counter,
                                df_resultJointTest,
                                frac_jump,
                                spec)
    {
    list_sum <- fun_runnerJoint(str_spec=spec)

    real_effNzReportSpec <- round(list_sum$effN.mean,1)
    real_effNxReportSpec <- round(list_sum$effNx.mean,1)
    real_corZReportSpec <- round(list_sum$cor.z,2)
    real_dimVecReportSpec <- dim
    real_nReportSpec <- n
    real_xjumpReportSpec <- x_jump
    real_zjumpReportSpec <- z_jump

    real_naiveReportResult <- round(list_sum$naive,3)
    real_bonfeReportResult <- round(list_sum$bonfe,3)
    real_jointReportResult <- round(list_sum$joint,3)
    real_maxTestSpec <- bool_max_test
    if (bool_max_test) {
      real_avgMaxTestCritical <- round(mean(list_sum$maxCritical),3)
      real_medianMaxTestCritical <- round(median(list_sum$maxCritical),3)
    } else {
      real_avgMaxTestCritical <- NA
      real_medianMaxTestCritical <- NA
    }

    list_speclist[counter] <- spec
    #print(spec)

    df_resultJointTest <-
      rbind(df_resultJointTest,
            data.frame(vec_effNzReportSpec = real_effNzReportSpec,
                       vec_effNxReportSpec = real_effNxReportSpec,
                       vec_corZReportSpec = cov_z,
                       vec_dimVecReportSpec = real_dimVecReportSpec,
                       vec_nReportSpec = real_nReportSpec,
                       vec_xjumpReportSpec = real_xjumpReportSpec,
                       vec_zjumpReportSpec = real_zjumpReportSpec,
                       vec_naiveReportResult = real_naiveReportResult,
                       vec_bonfeReportResult = real_bonfeReportResult,
                       vec_jointReportResult = real_jointReportResult,
                       vec_maxTestSpec = real_maxTestSpec,
                       vec_L2_std = bool_L2_std,
                       vec_avgMaxTestCritical = real_avgMaxTestCritical,
                       vec_medianMaxTestCritical = real_medianMaxTestCritical,
                       vec_jumpTestSpec = frac_jump))
    return(df_resultJointTest)
  }

  for (n in vec_n) {
    for (dim in vec_dimZ) {
      for (cov_z in vec_covZ) {
        for (frac_jump in vec_frac_jump) {
          cat("Running for n=",n,
              ",dim=",dim,
              ",cov_z=",cov_z,
              ",frac=", frac_jump,".\n")
          # bool_maxTest TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",TRUE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = TRUE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",FALSE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,
                                      x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = FALSE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec)
                  counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE, bool_L2_std TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                  spec <- paste0("int_ns=",int_ns,
                                 ",a_2=",0,
                                 ",x_jump=",x_jump,
                                 ",jump=",z_jump,
                                 ",dim=",dim,
                                 ",n=",n,
                                 ",cov_z =", cov_z,
                                 ",frac_jump=",frac_jump,
                                 ",bool_mutePrint=TRUE,bool_max_test=",FALSE,
                                 ",bool_L2_std=",TRUE)
                  df_resultJointTest <-
                    fun_repeaterJoint(n = n,
                                      cov_z = cov_z,
                                      dim = dim,
                                      x_jump = x_jump,
                                      z_jump = z_jump,
                                      bool_max_test = FALSE,
                                      counter = counter,
                                      df_resultJointTest = df_resultJointTest,
                                      frac_jump = frac_jump,
                                      spec = spec,
                                      bool_L2_std = TRUE)
                  counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
  #testthat::expect_snapshot_output(df_resultJointTest[,c(3,4,6:11)])
  testthat::expect_equal(df_resultJointTest[,8], c(0.3,0.5,0.9,1.0,1.0,0.3,0.5,0.9,1.0,1.0,0.3,0.5,0.9,1.0,1.0,0.3,0.5,0.5,0.9,1.0,0.3,0.5,0.5,0.9,1.0,0.3,0.5,0.5,0.9,1.0,0.3,0.2,0.5,0.4,0.5,0.3,0.2,0.5,0.4,0.5,0.3,0.2,0.5,0.4,0.5))
  testthat::expect_equal(df_resultJointTest[,9], c(0.1,0.1,0.5,0.9,0.9,0.1,0.1,0.5,0.9,0.9,0.1,0.1,0.5,0.9,0.9,0.1,0.0,0.3,0.5,0.8,0.1,0.0,0.3,0.5,0.8,0.1,0.0,0.3,0.5,0.8,0.1,0.1,0.1,0.1,0.4,0.1,0.1,0.1,0.1,0.4,0.1,0.1,0.1,0.1,0.4))
  testthat::expect_equal(df_resultJointTest[,10],c(0.1,0.2,0.6,0.9,0.9,0.1,0.2,0.8,1.0,1.0,0.1,0.2,0.4,0.9,0.9,0.1,0.0,0.3,0.6,0.8,0.1,0.1,0.2,0.8,1.0,0.1,0.2,0.2,0.4,0.9,0.1,0.1,0.1,0.3,0.4,0.1,0.0,0.1,0.1,0.3,0.1,0.2,0.1,0.2,0.5))

  list_compare <-
    list(
      case1 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.2,0.6,0.9,0.9,0.1,0.1,0.5,0.9,0.9,0.1,0.1,0.5,0.9,0.9,0.1,0.2,0.4,0.9,0.9),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ),
      case2 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.0,0.3,0.6,0.8,0.1,0.0,0.3,0.5,0.8,0.1,0.0,0.3,0.5,0.8,0.1,0.2,0.2,0.4,0.9),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ),
      case3 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.1,0.1,0.3,0.4,0.1,0.1,0.1,0.1,0.4,0.1,0.1,0.1,0.1,0.4,0.1,0.2,0.1,0.2,0.5),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ))

  count <- 0
  for (n in vec_n)
  {
    for (dim in unique(df_resultJointTest$vec_dimVecReportSpec))
    {
      for (cor in vec_covZ)
      {
        for (xval in unique(df_resultJointTest$vec_xjumpReportSpec))
        {
          for (jumpTest in list("'Last'","'Half'","'All'")) {
            df_resultEval <-
              df_resultJointTest[
                df_resultJointTest$vec_dimVecReportSpec==dim
              & df_resultJointTest$vec_corZReportSpec == cor
              & df_resultJointTest$vec_nReportSpec== n
              & df_resultJointTest$vec_xjumpReportSpec == xval
              & df_resultJointTest$vec_jumpTestSpec == jumpTest
              & df_resultJointTest$vec_maxTestSpec == TRUE,]
            df_resultEvalMat <-
              rbind(data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_jointReportResult,
                               type = "joint max"),
                    data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_bonfeReportResult,
                               type = "bonferroni"))

            df_resultEvalMat <-
              rbind(df_resultEvalMat,
                    data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_bonfeReportResult,
                               type = "bonferroni"))
            df_resultEvalChiStd <-
              df_resultJointTest[
                df_resultJointTest$vec_dimVecReportSpec==dim
              & df_resultJointTest$vec_corZReportSpec == cor
              & df_resultJointTest$vec_nReportSpec== n
              & df_resultJointTest$vec_xjumpReportSpec == xval
              & df_resultJointTest$vec_jumpTestSpec == jumpTest
              & df_resultJointTest$vec_maxTestSpec == FALSE
              & df_resultJointTest$vec_L2_std == TRUE,]
            df_resultEvalMat <-
              rbind(df_resultEvalMat,
                    data.frame(z = df_resultEvalChiStd$vec_zjumpReportSpec,
                               prob = df_resultEvalChiStd$vec_jointReportResult,
                               type = "joint chi Std"))
            count <- count + 1
            testthat::expect_equal(df_resultEvalMat,eval(parse(text=paste0("list_compare$case",count))))
            #
            #
            # p_2 <- ggplot(
            #   data = df_resultEvalMat,
            #   mapping = aes(x = z, y = prob,color=type)
            #   ) +
            #   geom_point() +
            #   scale_y_continuous(
            #     breaks=c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
            #   geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
            #   geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
            #   geom_hline(yintercept = 1, linetype = "dashed", size = 0.01) +
            #   ggtitle(
            #     paste0("x jump: ", xval,
            #            "; n: ", n,
            #            ", dim: ",dim,
            #            ", cor:",cor,
            #            ", alt Hyp:",jumpTest))
            # testthat::expect_snapshot_output(p_2)
          }
        }
      }
    }
  }
})


test_that("the vignette process for the first loop does not change, dim = 1, null coverage",{

  int_ns <- 10

  fun_runnerJoint = function(str_spec)
  {
    #tic(str_spec)
    time_start <- Sys.time()
    eval(parse(
      text = paste0(
        paste0("list_summary <- return_result_MC_joint(",str_spec),")"
      )))
    time_end <- Sys.time()
    #toc()
    return(list_summary)
  }

  # Specification list
  # values for x jumps
  vec_xJump <- c(0) #,0.15,0.3)
  # values for z jumps
  vec_zJump <- c(0) #seq(0,2,length=5)
  # sample sizes
  vec_n <- c(500)
  vec_covZ <- c(0.5)
  vec_dimZ <- c(1)
  # the last one, half, full:
  ## half: for 3, 2:3 jumps.
  vec_frac_jump <- list("'None'")

  int_numberSpec <-
    length(vec_xJump) *
    length(vec_n) *
    length(vec_zJump) *
    length(vec_covZ) *
    length(vec_dimZ) *
    length(vec_frac_jump)

  vec_effNzReportSpec <- c(1:int_numberSpec)
  vec_effNxReportSpec <- c(1:int_numberSpec)
  vec_corZReportSpec <- c(1:int_numberSpec)
  vec_dimVecReportSpec <- c(1:int_numberSpec)
  vec_nReportSpec <- c(1:int_numberSpec)
  vec_xjumpReportSpec <- c(1:int_numberSpec)
  vec_zjumpReportSpec <- c(1:int_numberSpec)
  vec_naiveReportResult <- c(1:int_numberSpec)
  vec_bonfeReportResult <- c(1:int_numberSpec)
  vec_jointReportResult <- c(1:int_numberSpec)
  vec_maxTestSpec <- c(1:int_numberSpec)
  vec_avgMaxTestCritical <- c(1:int_numberSpec)
  vec_medianMaxTestCritical <- c(1:int_numberSpec)
  vec_jumpTestSpec <- c(1:int_numberSpec)

  list_speclist <- list()

  counter <- 1
  df_resultJointTest <- data.frame()
  fun_repeaterJoint <- function(n,
                                cov_z,
                                dim,
                                x_jump,
                                z_jump,
                                bool_max_test,
                                bool_L2_std = FALSE,
                                counter,
                                df_resultJointTest,
                                frac_jump,
                                spec) {
    list_sum <- fun_runnerJoint(str_spec=spec)

    real_effNzReportSpec <- round(list_sum$effN.mean,1)
    real_effNxReportSpec <- round(list_sum$effNx.mean,1)
    real_corZReportSpec <- round(list_sum$cor.z,2)
    real_dimVecReportSpec <- dim
    real_nReportSpec <- n
    real_xjumpReportSpec <- x_jump
    real_zjumpReportSpec <- z_jump

    real_naiveReportResult <- round(list_sum$naive,3)
    real_bonfeReportResult <- round(list_sum$bonfe,3)
    real_jointReportResult <- round(list_sum$joint,3)
    real_maxTestSpec <- bool_max_test
    if (bool_max_test) {
      real_avgMaxTestCritical <- round(mean(list_sum$maxCritical),3)
      real_medianMaxTestCritical <- round(median(list_sum$maxCritical),3)
    } else {
      real_avgMaxTestCritical <- NA
      real_medianMaxTestCritical <- NA
    }

    list_speclist[counter] <- spec
    #print(spec)

    df_resultJointTest <-
      rbind(df_resultJointTest,
            data.frame(vec_effNzReportSpec = real_effNzReportSpec,
                       vec_effNxReportSpec = real_effNxReportSpec,
                       vec_corZReportSpec = cov_z,
                       vec_dimVecReportSpec = real_dimVecReportSpec,
                       vec_nReportSpec = real_nReportSpec,
                       vec_xjumpReportSpec = real_xjumpReportSpec,
                       vec_zjumpReportSpec = real_zjumpReportSpec,
                       vec_naiveReportResult = real_naiveReportResult,
                       vec_bonfeReportResult = real_bonfeReportResult,
                       vec_jointReportResult = real_jointReportResult,
                       vec_maxTestSpec = real_maxTestSpec,
                       vec_L2_std = bool_L2_std,
                       vec_avgMaxTestCritical = real_avgMaxTestCritical,
                       vec_medianMaxTestCritical = real_medianMaxTestCritical,
                       vec_jumpTestSpec = frac_jump))
    return(df_resultJointTest)
  }

  for (n in vec_n) {
    for (dim in vec_dimZ) {
      if (dim == 1) {
        vec_covZUsed <- c(0)
      } else {
        vec_covZUsed <- vec_covZ
      }
      for (cov_z in vec_covZUsed) {
        for (frac_jump in vec_frac_jump) {
          cat("Running for n=",n,
              ",dim=",dim,
              ",cov_z=",cov_z,
              ",frac=", frac_jump,".\n")
          # bool_maxTest TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = TRUE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE, bool_L2_std TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE,
                               ",bool_L2_std=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec,
                                    bool_L2_std=TRUE)
                counter <- counter + 1
              }
            }
          }
        }
      }

    }
  }
  df_coverage <- data.frame(n = df_resultJointTest[,5],
                            cov = df_resultJointTest[,3],
                            dim = df_resultJointTest[,4],
                            p.naive = df_resultJointTest[,8],
                            p.bonferroni = df_resultJointTest[,9],
                            p.joint = df_resultJointTest[,10],
                            isMaxTest = df_resultJointTest[,11],
                            isChiStd = df_resultJointTest[,12])
  testthat::expect_equal(df_coverage$p.naive,c(0,0,0))
  testthat::expect_equal(df_coverage$p.bonferroni,c(0,0,0))
  testthat::expect_equal(df_coverage$p.joint,c(0,0,0))
})

test_that("the vignette process for the first loop does not change, dim = 1, alternative power",{
  int_ns <- 10

  fun_runnerJoint = function(str_spec)
  {
    #tic(str_spec)
    time_start <- Sys.time()
    eval(parse(
      text = paste0(
        paste0("list_summary <- return_result_MC_joint(",str_spec),")"
      )))
    time_end <- Sys.time()
    #toc()
    return(list_summary)
  }

  vec_xJump <- c(0.15)
  vec_zJump <- seq(0,2,length=5)
  vec_n <- c(500)
  vec_covZ <- c(0.5)
  vec_dimZ <- c(1)
  vec_frac_jump <- list("'Last'","'Half'","'All'")
  int_numberSpec <-
    length(vec_xJump) *
    length(vec_n) *
    length(vec_zJump) *
    length(vec_covZ) *
    length(vec_dimZ) *
    2 *
    length(vec_frac_jump)

  vec_effNzReportSpec <- c(1:int_numberSpec)
  vec_effNxReportSpec <- c(1:int_numberSpec)
  vec_corZReportSpec <- c(1:int_numberSpec)
  vec_dimVecReportSpec <- c(1:int_numberSpec)
  vec_nReportSpec <- c(1:int_numberSpec)
  vec_xjumpReportSpec <- c(1:int_numberSpec)
  vec_zjumpReportSpec <- c(1:int_numberSpec)
  vec_naiveReportResult <- c(1:int_numberSpec)
  vec_bonfeReportResult <- c(1:int_numberSpec)
  vec_jointReportResult <- c(1:int_numberSpec)
  vec_maxTestSpec <- c(1:int_numberSpec)
  vec_avgMaxTestCritical <- c(1:int_numberSpec)
  vec_medianMaxTestCritical <- c(1:int_numberSpec)
  vec_jumpTestSpec <- c(1:int_numberSpec)

  list_speclist <- list()

  counter <- 1
  df_resultJointTest <- data.frame()
  fun_repeaterJoint <- function(n,
                                cov_z,
                                dim,
                                x_jump,
                                z_jump,
                                bool_max_test,
                                bool_L2_std = FALSE,
                                counter,
                                df_resultJointTest,
                                frac_jump,
                                spec)
  {
    list_sum <- fun_runnerJoint(str_spec=spec)

    real_effNzReportSpec <- round(list_sum$effN.mean,1)
    real_effNxReportSpec <- round(list_sum$effNx.mean,1)
    real_corZReportSpec <- round(list_sum$cor.z,2)
    real_dimVecReportSpec <- dim
    real_nReportSpec <- n
    real_xjumpReportSpec <- x_jump
    real_zjumpReportSpec <- z_jump

    real_naiveReportResult <- round(list_sum$naive,3)
    real_bonfeReportResult <- round(list_sum$bonfe,3)
    real_jointReportResult <- round(list_sum$joint,3)
    real_maxTestSpec <- bool_max_test
    if (bool_max_test) {
      real_avgMaxTestCritical <- round(mean(list_sum$maxCritical),3)
      real_medianMaxTestCritical <- round(median(list_sum$maxCritical),3)
    } else {
      real_avgMaxTestCritical <- NA
      real_medianMaxTestCritical <- NA
    }

    list_speclist[counter] <- spec
    #print(spec)

    df_resultJointTest <-
      rbind(df_resultJointTest,
            data.frame(vec_effNzReportSpec = real_effNzReportSpec,
                       vec_effNxReportSpec = real_effNxReportSpec,
                       vec_corZReportSpec = cov_z,
                       vec_dimVecReportSpec = real_dimVecReportSpec,
                       vec_nReportSpec = real_nReportSpec,
                       vec_xjumpReportSpec = real_xjumpReportSpec,
                       vec_zjumpReportSpec = real_zjumpReportSpec,
                       vec_naiveReportResult = real_naiveReportResult,
                       vec_bonfeReportResult = real_bonfeReportResult,
                       vec_jointReportResult = real_jointReportResult,
                       vec_maxTestSpec = real_maxTestSpec,
                       vec_L2_std = bool_L2_std,
                       vec_avgMaxTestCritical = real_avgMaxTestCritical,
                       vec_medianMaxTestCritical = real_medianMaxTestCritical,
                       vec_jumpTestSpec = frac_jump))
    return(df_resultJointTest)
  }

  for (n in vec_n) {
    for (dim in vec_dimZ) {
      for (cov_z in vec_covZ) {
        for (frac_jump in vec_frac_jump) {
          cat("Running for n=",n,
              ",dim=",dim,
              ",cov_z=",cov_z,
              ",frac=", frac_jump,".\n")
          # bool_maxTest TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = TRUE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE, bool_L2_std TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE,
                               ",bool_L2_std=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec,
                                    bool_L2_std = TRUE)
                counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
  #testthat::expect_snapshot_output(df_resultJointTest[,c(3,4,6:11)])
  testthat::expect_equal(df_resultJointTest[,8], c(0.3,0.3,0.6,0.9,1.0,0.3,0.3,0.6,0.9,1.0,0.3,0.3,0.6,0.9,1.0,0.3,0.6,1.0,1.0,1.0,0.3,0.6,1.0,1.0,1.0,0.3,0.6,1.0,1.0,1.0,0.3,0.3,0.6,0.9,1.0,0.3,0.3,0.6,0.9,1.0,0.3,0.3,0.6,0.9,1.0))
  testthat::expect_equal(df_resultJointTest[,9], c(0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.4,1.0,1.0,1.0,0.1,0.4,1.0,1.0,1.0,0.1,0.4,1.0,1.0,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0))
  testthat::expect_equal(df_resultJointTest[,10],c(0.1,0.2,0.5,0.7,1.0,0.1,0.2,0.5,0.8,1.0,0.1,0.2,0.5,0.8,1.0,0.1,0.5,1.0,1.0,1.0,0.1,0.5,1.0,1.0,1.0,0.1,0.5,1.0,1.0,1.0,0.1,0.2,0.5,0.7,1.0,0.1,0.2,0.5,0.8,1.0,0.1,0.2,0.5,0.8,1.0))

  list_compare <-
    list(
      case1 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.2,0.5,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.5,0.8,1.0),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ),
      case2 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.5,1.0,1.0,1.0,0.1,0.4,1.0,1.0,1.0,0.1,0.4,1.0,1.0,1.0,0.1,0.5,1.0,1.0,1.0),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ),
      case3 = data.frame(
        z = c(0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0,0.0,0.5,1.0,1.5,2.0),
        prob = c(0.1,0.2,0.5,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.4,0.7,1.0,0.1,0.2,0.5,0.8,1.0),
        type = c("joint max","joint max","joint max","joint max","joint max","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","bonferroni","joint chi Std","joint chi Std","joint chi Std","joint chi Std","joint chi Std")
      ))

  count <- 0
  for (n in vec_n)
  {
    for (dim in unique(df_resultJointTest$vec_dimVecReportSpec))
    {
      for (cor in vec_covZ)
      {
        for (xval in unique(df_resultJointTest$vec_xjumpReportSpec))
        {
          for (jumpTest in list("'Last'","'Half'","'All'")) {
            df_resultEval <-
              df_resultJointTest[
                df_resultJointTest$vec_dimVecReportSpec==dim
                & df_resultJointTest$vec_corZReportSpec == cor
                & df_resultJointTest$vec_nReportSpec== n
                & df_resultJointTest$vec_xjumpReportSpec == xval
                & df_resultJointTest$vec_jumpTestSpec == jumpTest
                & df_resultJointTest$vec_maxTestSpec == TRUE,]
            df_resultEvalMat <-
              rbind(data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_jointReportResult,
                               type = "joint max"),
                    data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_bonfeReportResult,
                               type = "bonferroni"))

            df_resultEvalMat <-
              rbind(df_resultEvalMat,
                    data.frame(z = df_resultEval$vec_zjumpReportSpec,
                               prob = df_resultEval$vec_bonfeReportResult,
                               type = "bonferroni"))
            df_resultEvalChiStd <-
              df_resultJointTest[
                df_resultJointTest$vec_dimVecReportSpec==dim
                & df_resultJointTest$vec_corZReportSpec == cor
                & df_resultJointTest$vec_nReportSpec== n
                & df_resultJointTest$vec_xjumpReportSpec == xval
                & df_resultJointTest$vec_jumpTestSpec == jumpTest
                & df_resultJointTest$vec_maxTestSpec == FALSE
                & df_resultJointTest$vec_L2_std == TRUE,]
            df_resultEvalMat <-
              rbind(df_resultEvalMat,
                    data.frame(z = df_resultEvalChiStd$vec_zjumpReportSpec,
                               prob = df_resultEvalChiStd$vec_jointReportResult,
                               type = "joint chi Std"))
            count <- count + 1
            testthat::expect_equal(df_resultEvalMat,eval(parse(text=paste0("list_compare$case",count))))
            #
            #
            # p_2 <- ggplot(
            #   data = df_resultEvalMat,
            #   mapping = aes(x = z, y = prob,color=type)
            #   ) +
            #   geom_point() +
            #   scale_y_continuous(
            #     breaks=c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
            #   geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
            #   geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
            #   geom_hline(yintercept = 1, linetype = "dashed", size = 0.01) +
            #   ggtitle(
            #     paste0("x jump: ", xval,
            #            "; n: ", n,
            #            ", dim: ",dim,
            #            ", cor:",cor,
            #            ", alt Hyp:",jumpTest))
            # testthat::expect_snapshot_output(p_2)
          }
        }
      }
    }
  }
})

test_that("the vignette process for the first loop does not change, dim = 3, negative corr",{

  int_ns <- 10

  fun_runnerJoint = function(str_spec)
  {
    #tic(str_spec)
    time_start <- Sys.time()
    eval(parse(
      text = paste0(
        paste0("list_summary <- return_result_MC_joint(",str_spec),")"
      )))
    time_end <- Sys.time()
    #toc()
    return(list_summary)
  }

  # Specification list
  # values for x jumps
  vec_xJump <- c(0) #,0.15,0.3)
  # values for z jumps
  vec_zJump <- c(0) #seq(0,2,length=5)
  # sample sizes
  vec_n <- c(500)
  vec_covZ <- c(-0.5)
  vec_dimZ <- c(3)
  # the last one, half, full:
  ## half: for 3, 2:3 jumps.
  vec_frac_jump <- list("'None'")

  int_numberSpec <-
    length(vec_xJump) *
    length(vec_n) *
    length(vec_zJump) *
    length(vec_covZ) *
    length(vec_dimZ) *
    length(vec_frac_jump)

  vec_effNzReportSpec <- c(1:int_numberSpec)
  vec_effNxReportSpec <- c(1:int_numberSpec)
  vec_corZReportSpec <- c(1:int_numberSpec)
  vec_dimVecReportSpec <- c(1:int_numberSpec)
  vec_nReportSpec <- c(1:int_numberSpec)
  vec_xjumpReportSpec <- c(1:int_numberSpec)
  vec_zjumpReportSpec <- c(1:int_numberSpec)
  vec_naiveReportResult <- c(1:int_numberSpec)
  vec_bonfeReportResult <- c(1:int_numberSpec)
  vec_jointReportResult <- c(1:int_numberSpec)
  vec_maxTestSpec <- c(1:int_numberSpec)
  vec_avgMaxTestCritical <- c(1:int_numberSpec)
  vec_medianMaxTestCritical <- c(1:int_numberSpec)
  vec_jumpTestSpec <- c(1:int_numberSpec)

  list_speclist <- list()

  counter <- 1
  df_resultJointTest <- data.frame()
  fun_repeaterJoint <- function(n,
                                cov_z,
                                dim,
                                x_jump,
                                z_jump,
                                bool_max_test,
                                bool_L2_std = FALSE,
                                counter,
                                df_resultJointTest,
                                frac_jump,
                                spec) {
    list_sum <- fun_runnerJoint(str_spec=spec)

    real_effNzReportSpec <- round(list_sum$effN.mean,1)
    real_effNxReportSpec <- round(list_sum$effNx.mean,1)
    real_corZReportSpec <- round(list_sum$cor.z,2)
    real_dimVecReportSpec <- dim
    real_nReportSpec <- n
    real_xjumpReportSpec <- x_jump
    real_zjumpReportSpec <- z_jump

    real_naiveReportResult <- round(list_sum$naive,3)
    real_bonfeReportResult <- round(list_sum$bonfe,3)
    real_jointReportResult <- round(list_sum$joint,3)
    real_maxTestSpec <- bool_max_test
    if (bool_max_test) {
      real_avgMaxTestCritical <- round(mean(list_sum$maxCritical),3)
      real_medianMaxTestCritical <- round(median(list_sum$maxCritical),3)
    } else {
      real_avgMaxTestCritical <- NA
      real_medianMaxTestCritical <- NA
    }

    list_speclist[counter] <- spec
    #print(spec)

    df_resultJointTest <-
      rbind(df_resultJointTest,
            data.frame(vec_effNzReportSpec = real_effNzReportSpec,
                       vec_effNxReportSpec = real_effNxReportSpec,
                       vec_corZReportSpec = cov_z,
                       vec_dimVecReportSpec = real_dimVecReportSpec,
                       vec_nReportSpec = real_nReportSpec,
                       vec_xjumpReportSpec = real_xjumpReportSpec,
                       vec_zjumpReportSpec = real_zjumpReportSpec,
                       vec_naiveReportResult = real_naiveReportResult,
                       vec_bonfeReportResult = real_bonfeReportResult,
                       vec_jointReportResult = real_jointReportResult,
                       vec_maxTestSpec = real_maxTestSpec,
                       vec_L2_std = bool_L2_std,
                       vec_avgMaxTestCritical = real_avgMaxTestCritical,
                       vec_medianMaxTestCritical = real_medianMaxTestCritical,
                       vec_jumpTestSpec = frac_jump))
    return(df_resultJointTest)
  }

  for (n in vec_n) {
    for (dim in vec_dimZ) {
      if (dim == 1) {
        vec_covZUsed <- c(0)
      } else {
        vec_covZUsed <- vec_covZ
      }
      for (cov_z in vec_covZUsed) {
        for (frac_jump in vec_frac_jump) {
          cat("Running for n=",n,
              ",dim=",dim,
              ",cov_z=",cov_z,
              ",frac=", frac_jump,".\n")
          # bool_maxTest TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = TRUE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec)
                counter <- counter + 1
              }
            }
          }
          # bool_maxTest FALSE, bool_L2_std TRUE
          for (x_jump in vec_xJump) {
            for (z_jump in vec_zJump) {
              if ((x_jump > 0.000001) | (x_jump == 0 & z_jump == 0)) {
                spec <- paste0("int_ns=",int_ns,
                               ",a_2=",0,
                               ",x_jump=",x_jump,
                               ",jump=",z_jump,
                               ",dim=",dim,
                               ",n=",n,
                               ",cov_z =", cov_z,
                               ",frac_jump=",frac_jump,
                               ",bool_mutePrint=TRUE,bool_max_test=",FALSE,
                               ",bool_L2_std=",TRUE)
                df_resultJointTest <-
                  fun_repeaterJoint(n = n,
                                    cov_z = cov_z,
                                    dim = dim,
                                    x_jump = x_jump,
                                    z_jump = z_jump,
                                    bool_max_test = FALSE,
                                    counter = counter,
                                    df_resultJointTest = df_resultJointTest,
                                    frac_jump = frac_jump,
                                    spec = spec,
                                    bool_L2_std=TRUE)
                counter <- counter + 1
              }
            }
          }
        }
      }

    }
  }
  df_coverage <- data.frame(n = df_resultJointTest[,5],
                            cov = df_resultJointTest[,3],
                            dim = df_resultJointTest[,4],
                            p.naive = df_resultJointTest[,8],
                            p.bonferroni = df_resultJointTest[,9],
                            p.joint = df_resultJointTest[,10],
                            isMaxTest = df_resultJointTest[,11],
                            isChiStd = df_resultJointTest[,12])
  testthat::expect_equal(df_coverage$p.naive,c(0.1,0.1,0.1))
  testthat::expect_equal(df_coverage$p.bonferroni,c(0,0,0))
  testthat::expect_equal(df_coverage$p.joint,c(0,0,0))
})

test_that("the summary returns the expected display, covariates",{
  set.seed(1)
  N <- 1000
  Z <- data.frame(var1 = rnorm(N),
                  var2 = rnorm(N),
                  var3 = rnorm(N),
                  var4 = rnorm(N),
                  var5 = rnorm(N))
  vec_X <- rnorm(N)
  res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X, bool_joint = FALSE)
  ans <- summary(res_rd)

  tstats <- matrix(NA, 0L, 3L,
                   dimnames = list(NULL,
                                   c("Variable",
                                     "t-stat",
                                     "(Naive) p-value")))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var1",
                `t-stat` = -0.182,
                `(Naive) p-value` = 0.428))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var2",
                `t-stat` = 0.125,
                `(Naive) p-value` = 0.45))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var3",
                `t-stat` = 0.403,
                `(Naive) p-value` = 0.343))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var4",
                `t-stat` = 1.53,
                `(Naive) p-value` = 0.063))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var5",
                `t-stat` = -0.301,
                `(Naive) p-value` = 0.382))
  expect_equal(ans$tstats,tstats)
})

test_that("the summary returns the expected display, stdWald",{
  set.seed(1)
  N <- 1000
  Z <- data.frame(var1 = rnorm(N),
                  var2 = rnorm(N),
                  var3 = rnorm(N),
                  var4 = rnorm(N),
                  var5 = rnorm(N))
  vec_X <- rnorm(N)
  res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X)
  ans <- summary(res_rd)

  tstats <- matrix(NA, 0L, 3L,
                         dimnames = list(NULL,
                                           c("Variable",
                                               "t-stat",
                                               "(Naive) p-value")))
  tstats <-
    rbind(tstats,
          cbind(Variable = "density",
                `t-stat` = -0.303,
                `(Naive) p-value` = 0.381))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var1",
                `t-stat` = -0.182,
                `(Naive) p-value` = 0.428))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var2",
                `t-stat` = 0.125,
                `(Naive) p-value` = 0.45))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var3",
                `t-stat` = 0.403,
                `(Naive) p-value` = 0.343))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var4",
                `t-stat` = 1.53,
                `(Naive) p-value` = 0.063))
  tstats <-
    rbind(tstats,
          cbind(Variable = "var5",
                `t-stat` = -0.301,
                `(Naive) p-value` = 0.382))
  expect_equal(ans$tstats,tstats)

  joint_result <- matrix(NA, 0L, 4L,
                             dimnames = list(NULL,
                                             c("Test Type",
                                               "Test Statistic",
                                               "Critical Value",
                                               "Joint p-value")))

  test_type = "standardized Wald test"
  crit_value <- 12.825
  names(crit_value) <- "95%"
  joint_result <-
    rbind(joint_result,
          cbind(
            `Test Type` = test_type,
            `Test Statistic` = 2.737,
            `Critical Value` = crit_value,
            `Joint p-alue` = 0.834
          ))
  expect_equal(knitr::kable(ans$joint_result,"rst"),knitr::kable(joint_result,"rst"))
})
test_that("the summary returns the expected display, max",{
  set.seed(1)
  N <- 1000
  Z <- data.frame(var1 = rnorm(N),
                  var2 = rnorm(N),
                  var3 = rnorm(N),
                  var4 = rnorm(N),
                  var5 = rnorm(N))
  vec_X <- rnorm(N)
  res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X,bool_max_test = TRUE)
  ans <- summary(res_rd)

  joint_result <- matrix(NA, 0L, 4L,
                         dimnames = list(NULL,
                                         c("Test Type",
                                           "Test Statistic",
                                           "Critical Value",
                                           "Joint p-value")))

  test_type = "Max test"
  crit_value <- 6.794
  names(crit_value) <- "95%"
  joint_result <-
    rbind(joint_result,
          cbind(
            `Test Type` = test_type,
            `Test Statistic` = 2.344,
            `Critical Value` = crit_value,
            `Joint p-alue` = 0.535
          ))
  expect_equal(knitr::kable(ans$joint_result,"rst"),knitr::kable(joint_result,"rst"))
})

test_that("the summary returns the expected display, naive Wald",{
  set.seed(1)
  N <- 1000
  Z <- data.frame(var1 = rnorm(N),
                  var2 = rnorm(N),
                  var3 = rnorm(N),
                  var4 = rnorm(N),
                  var5 = rnorm(N))
  vec_X <- rnorm(N)
  res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X,bool_L2_std = FALSE,bool_max_test = FALSE)
  ans <- summary(res_rd)

  joint_result <- matrix(NA, 0L, 4L,
                         dimnames = list(NULL,
                                         c("Test Type",
                                           "Test Statistic",
                                           "Critical Value",
                                           "Joint p-value")))

  test_type = "(non-standardized) Wald test"
  crit_value <- 12.592
  joint_result <-
    rbind(joint_result,
          cbind(
            `Test Type` = test_type,
            `Test Statistic` = 2.805,
            `Critical Value` = crit_value,
            `Joint p-alue` = 0.833
          ))
  expect_equal(knitr::kable(ans$joint_result,"rst"),knitr::kable(joint_result,"rst"))
})
