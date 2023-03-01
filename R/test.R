
###################################################
####### SIMULATE A DATA AND FIT THE MODEL #########
###################################################

set.seed(26)

# install & load
devtools::install_github("AliGhanbari26/GPRMortality")
library("SpaVar")
library(LaplacesDemon)
library(rstan)
load("data/adjMatMu.RData")
load("data/adjMatVar.RData")

# the adjacent matrix for 429 districts of Iran
# embeded into the package and it is used as adjacent matrix
# for central tendncy parameter spatial structure
dim(adjMatMu)

# the adjacent matrix for 31 provinces of Iran
# embeded into the package and it is used as adjacent matrix
# for dispersion parameter spatial structure
dim(adjMatVar)

# create a frame of simulated data
data  = expand.grid(Mu = rownames(adjMatMu) , age=1:4  )
n = nrow(data)
data$WI =  rnorm(n)

# create a frame of simulated data
index_locMu  = match(data$Mu,rownames(adjMatMu))
index_locVar = match(substr(data$Mu,1,2),rownames(adjMatVar))

# create X & Z matrix
# both of them have two covariates ; age(categorical) and one continious
X = stats::model.matrix(~1+factor(data$age)+WI,data=data)
Z = X
p = ncol(X) ; q = ncol(Z) ; n_Mu = nrow(adjMatMu) ; n_Var = nrow(adjMatVar)

# simulate central tendency vector that
# consist of two part, spatial (ICAR) and non-spatial (covariate) parts
beta_x = rnorm(p,0,1)
D_Mu = diag(apply(adjMatMu,1,sum))
omega_Mu = .5*(   D_Mu -.99*adjMatMu  )
d_Mu =  t(LaplacesDemon::rmvnp(1 , rep(0,n_Mu), omega_Mu ))
MU = X%*%beta_x + d_Mu[index_locMu] + rnorm(n,0,.1)

# simulate dispersion vector that
# consist of two part, spatial (ICAR) and non-spatial (covariate) parts
beta_z = rnorm(q,0,.5)
PHI_x =  ( Z%*%beta_z)+ rnorm(n,0,.05)
D_Var = diag(apply(adjMatVar,1,sum))
omega_Var = .5*(   D_Var -.99*adjMatVar  )
d_Var =  t(LaplacesDemon::rmvnp(1 , rep(0,n_Var), omega_Var ))
PHI_spa = d_Var[index_locVar]
PHI = exp(PHI_x + PHI_spa)

# simulate response
y  = rnorm(n,MU,PHI)

# fit & compare
# model criterias such as DIC, WAIC, beta (se), beta (bias)

formula_mu  = formula("y~WI + factor(age)")
formula_var = formula(" ~WI + factor(age)")

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



