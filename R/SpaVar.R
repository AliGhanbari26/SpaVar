
#' @export


setClass("SpaVar", slots=
           list(
             fit.simple   ="stanfit",
             fit.complex ="stanfit",
             fit.simple.summary.X ="matrix",
             fit.complex.summary.X ="matrix",
             data ="list",
             DIC.simple   =  "list",
             DIC.complex  =  "list",
             WAIC.simple  =  "list",
             WAIC.complex =  "list"
           )
)

setMethod("show",
          "SpaVar",
          function(object) {
            cat("DIC complex model: ",object@DIC.complex$DIC, "\n")
            cat("    simple  model: ",object@DIC.simple$DIC, "\n")
            cat("WAIC complex model: ",object@WAIC.complex$WAIC, "\n")
            cat("     simple  model: ",object@WAIC.simple$WAIC, "\n")

            cat("\n")
            cat("Coefficient complex model (central tendancy): ", "\n")
            print( object@fit.complex.summary.X )
            cat("\n")
            cat("Coefficient simple model (central tendancy): ", "\n")
            print( object@fit.simple.summary.X )
            cat("\n")

            cat("      < complex model spatial parameter(s)>","\n")
            cat("Sigma spatial ICAR parameter (central tendancy): ", "\n")
            print(summary(object@fit.complex,"sigma_d")$summary[,c("mean","se_mean")] )
            cat("\n")

            cat("Sigma spatial ICAR parameter (dispersion): ", "\n")
            print(summary(object@fit.complex,"sigma_phi")$summary[,c("mean","se_mean")] )
            cat("\n")
            cat("      < simple model spatial parameter(s)>","\n")
            cat("Sigma spatial ICAR parameter (central tendancy): ", "\n")
            print(summary(object@fit.simple,"sigma_d")$summary[,c("mean","se_mean")] )
            cat("\n")
          }
)

SpaVar = function(formula_mu,
                  formula_var,
                  adjMatMu, adjMatVar ,
                  index_locMu,index_locVar,
                  distribution="normal",
                  SpaStrVar="ICAR", SpaStrMu="ICAR",
                  n.iter=4000, n.warmup=2000,n.core=2,n.chain=2,
                  data        ){

  # stan data ---------------------------------------------------------------

  stan_data = list()
  stan_data$X = model.matrix(formula_mu,data)
  stan_data$Z = model.matrix(formula_var,data)
  stan_data$n       =  nrow(data)
  stan_data$p       =  ncol(stan_data$X)
  stan_data$q       =  ncol(stan_data$Z)

  adjMatVar_one = adjMatVar
  for(i in 1:nrow(adjMatVar)){
    adjMatVar_one[i,] = adjMatVar_one[i,]/sum(adjMatVar_one[i,])
  }
  adjMatMu_one = adjMatMu
  for(i in 1:nrow(adjMatMu)){
    adjMatMu_one[i,] = adjMatMu_one[i,]/sum(adjMatMu_one[i,])
  }
  stan_data$W_dis_one = adjMatMu_one
  stan_data$W_pro_one = adjMatVar_one
  stan_data$n_d      =  nrow(adjMatMu)
  stan_data$n_p     =  nrow(adjMatVar)
  stan_data$index_p  =  index_locVar
  stan_data$index_d  =  index_locMu
  stan_data$nb_count_dis = apply(adjMatMu ,1,sum)
  stan_data$nb_count_pro = apply(adjMatVar,1,sum)
  stan_data$y  = y

  # fit complex model -------------------------------------------------------

  fit.stan.complex = rstan::sampling(compile.complex
                                     ,data = stan_data
                                     ,chain=n.chain
                                     ,core=n.core
                                     ,iter = n.iter
                                     ,warmup  = n.warmup
                                     ,thin = 1
  )

  # fit simple model -------------------------------------------------------

  fit.stan.simple = rstan::sampling(compile.simple
                                    ,data = stan_data
                                    ,chain=n.chain
                                    ,core=n.core
                                    ,iter = n.iter
                                    ,warmup  = n.warmup
                                    ,thin = 1
  )


  # Calc DIC WAIC -----------------------------------------------------------

if(F){
  loop_path = "2_1"
  load(file=paste0("data/compile.simple.RData"))
  load(file=paste0("d:/RData/stan_data_",loop_path,".RData"))
  load(file=paste0("d:/RData/fit_simple_",loop_path,".RData"))
  load(file=paste0("d:/RData/fit_complex_",loop_path,".RData"))
}

  ext.MU  <- t(rstan::extract(fit.stan.complex,"MU")[[1]] );dim(ext.MU)
  ext.PHI <- t(rstan::extract(fit.stan.complex,"PHI")[[1]]);dim(ext.PHI)
  LL.complex = matrix(0,nrow(ext.MU),ncol(ext.MU));dim(LL.complex)
  for(i in 1:nrow(LL.complex)){
    for(j in 1:ncol(LL.complex)){
      LL.complex[i,j] = dnorm(stan_data$y[i]
                              ,  ext.MU[i,j]
                              ,  ext.PHI[i,j]
                              , log=T)

    }
  }
  WAIC.complex = LaplacesDemon::WAIC(LL.complex)
  Dev.complex <- -2*colSums(LL.complex)
  DIC.complex <- list(DIC=mean(Dev.complex) + var(Dev.complex)/2, Dbar=mean(Dev.complex), pV=var(Dev.complex)/2)


  ext.MU  <- t(rstan::extract(fit.stan.simple,"MU")[[1]] );dim(ext.MU)
  ext.PHI <- t(rstan::extract(fit.stan.simple,"PHI")[[1]]);dim(ext.PHI)
  LL.simple = matrix(0,nrow(ext.MU),ncol(ext.MU));dim(LL.simple)
  for(i in 1:nrow(LL.simple)){
    for(j in 1:ncol(LL.simple)){
      LL.simple[i,j] = dnorm(stan_data$y[i]
                             ,  ext.MU[i,j]
                             ,  ext.PHI[i,j]
                             , log=T)

    }
  }
  WAIC.simple = LaplacesDemon::WAIC(LL.simple)
  Dev.simple <- -2*colSums(LL.simple)
  DIC.simple <- list(DIC=mean(Dev.simple) + var(Dev.simple)/2, Dbar=mean(Dev.simple), pV=var(Dev.simple)/2)




  # create output -----------------------------------------------------------

  tmp = summary(object@fit.complex,"beta_x")$summary[,c("mean","se_mean","2.5%","97.5%","n_eff","Rhat")]
  rownames(tmp) = colnames(stan_data$X)
  fit.complex.summary.X =  tmp

  tmp = summary(object@fit.simple,"beta_x")$summary[,c("mean","se_mean","2.5%","97.5%","n_eff","Rhat")]
  rownames(tmp) = colnames(stan_data$X)
  fit.simple.summary.X =  tmp

  output <- new("SpaVar",
                fit.simple    = fit.stan.simple,
                fit.complex   = fit.stan.complex,

                fit.simple.summary.X    = fit.simple.summary.X,
                fit.complex.summary.X   = fit.complex.summary.X,

                data          = stan_data,
                DIC.simple    = DIC.simple,
                DIC.complex   = DIC.complex,
                WAIC.simple   = WAIC.simple,
                WAIC.complex  = WAIC.complex
  )

  output

}

