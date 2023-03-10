\docType{class}
\name{SpaVar-class}
\alias{SpaVar-class}
\title{Fitting two models with and without spatial structure on the dispersion parameter and comparing them.}
\arguments{
\item{formula_mu}{The formula that determines which variables must consider in estimating the central tendency parameter.}

\item{formula_var}{The formula that determines which variables must consider in estimating the dispersion parameter.}

\item{distribution}{The distribution of response variable.}

\item{SpaStrVar}{The spatial structure for dispersion parameter.}

\item{SpaStrMu}{The spatial structure for central tendency parameter.}

\item{adjMatMu}{The adjucent matrix related to central tendency parameter}

\item{adjMatVar}{The adjucent matrix related to dispersen parameter}

\item{index_locMu}{A vector that detrmine each row of data related to which rows of the adjMatMu}

\item{index_locVar}{A vector that detrmine each row of data related to which rows of the adjMatVar}

\item{n.burn}{the number of iterations to use for running model. If not specified, then the default is 4000.}

\item{n.warmup}{the number of samples for warm-up. If not specified, the default is 2000.}

\item{n.core}{The number of cores for running markov chains parallelly.}

\item{n.chain}{The number of markov chains.}


\item{data}{data (dataframe)}
}
\description{
This package aims to emphasize the spatial structure of dispersion parameters in addition to central parameters simultaneously; there are many models and techniques that are used to consider the spatial structure of central tendency parameters in different distributions. In order to consider the dispersion and central tendency parameters independently, two spatial structures are considered separately. This model complexity helps estimate covariates more precisely and has better model selection criteria such as DIC and WAIC. These criteria calculate automatically and are present in the function's result. Moreover, models are implemented in STAN language, and users are able to use STAN objects in the function's output for more exploration.
}
\details{
In the current version, only the ICAR spatial structure and normal distribution are included. The package will be completed more in different aspects, such as supporting more distributions and spatial structures.
}
\value{
An object of S4 class SpaVar.
}
\author{
Ali Ghanbari
}
\examples{


######################################################
####### SIMULATING DATA AND FITTING THE MODEL ########
######################################################

set.seed(26)

# install & load
devtools::install_github("AliGhanbari26/GPRMortality")
library("SpaVar")

# The adjacent matrix for 429 districts of Iran is embedded into the package,
# and it is used as an adjacent matrix for the spatial structure of
# the central tendency parameter
dim(adjMatMu)

# The adjacent matrix for 31 provinces of Iran is embedded into the package,
# and it is used as an adjacent matrix for the spatial structure of
# the dispersion parameter
dim(adjMatVar)

# create a frame of simulated data
data  = expand.grid(Mu = rownames(adjMatMu) , age=1:4  )
n = nrow(data)
data$WI =  rnorm(n)

# create a frame of simulated data
index_locMu  = match(data$Mu,rownames(adjMatMu))
index_locVar = match(substr(data$Mu,1,2),rownames(adjMatVar))

# create X & Z matrix both of them have two covariates ;
# age(categorical) and one continious
X = stats::model.matrix(~1+factor(data$age)+WI,data=data)
Z = X
p = ncol(X) ; q = ncol(Z) ; n_Mu = nrow(adjMatMu) ; n_Var = nrow(adjMatVar)

# simulate the central tendency vector that consists of two parts,
# spatial (ICAR) and non-spatial (covariate) parts
beta_x = rnorm(p,0,1)
D_Mu = diag(apply(adjMatMu,1,sum))
omega_Mu = .5*(   D_Mu -.99*adjMatMu  )
d_Mu =  t(LaplacesDemon::rmvnp(1 , rep(0,n_Mu), omega_Mu ))
MU = X\%*\%beta_x + d_Mu[index_locMu] + rnorm(n,0,.1)

# simulate the dispersion vector that consists of two parts,
# spatial (ICAR) and non-spatial (covariate) parts
beta_z = rnorm(q,0,.5)
PHI_x =  ( Z\%*\%beta_z)+ rnorm(n,0,.05)
D_Var = diag(apply(adjMatVar,1,sum))
omega_Var = .5*(   D_Var -.99*adjMatVar  )
d_Var =  t(LaplacesDemon::rmvnp(1 , rep(0,n_Var), omega_Var ))
PHI_spa = d_Var[index_locVar]
PHI = exp(PHI_x + PHI_spa)

# simulate response
y  = rnorm(n,MU,PHI)

# fit & compare model criterias such as DIC, WAIC, beta (se), beta (bias)

fit = SpaVar(formula_mu  = WI + factor(age),
                  formula_var = WI + factor(age),
                  adjMatMu = adjMatMu,
                  adjMatVar = adjMatVar ,
                  index_locMu = index_locMu,
                  index_locVar = index_locVar,
                  distribution = "normal",
                  SpaStrVar="ICAR", SpaStrMu="ICAR",
                  n.iter=4000,n.warmup=2000,n.core=2,n.chain=2,
                  data        )

fit

}
