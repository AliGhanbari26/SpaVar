# SpaVar
Fitting two models with and without spatial structure on the dispersion parameter and comparing them.

## Description
This package aims to emphasize the spatial structure of dispersion parameters in addition to central parameters simultaneously; there are many models and techniques that are used to consider the spatial structure of central tendency parameters in different distributions. In order to consider the dispersion and central tendency parameters independently, two spatial structures are considered separately. This model complexity helps estimate covariates more precisely and has better model selection criteria such as DIC and WAIC. These criteria calculate automatically and are present in the function's result. Moreover, models are implemented in STAN language, and users are able to use STAN objects in the function's output for more exploration.

In the current version, only the ICAR spatial structure and normal distribution are included. The package will be completed more in different aspects, such as supporting more distributions and spatial structures.

## Example
We simulate the response variable and two covariates, continuous and discrete. The spatial structure of the ICAR is taken into account for 429 districts. Additionally, the variance of the response variable is affected by the two covariates mentioned above and the ICAR spatial structure (31 provinces). A thin line indicates the district border, and a thick line indicates the province border in the following figure.


![image](https://user-images.githubusercontent.com/30459265/222126706-60659259-f278-4a97-9459-c2c635894529.png)


```r
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
MU = X%*%beta_x + d_Mu[index_locMu] + rnorm(n,0,.1)

# simulate the dispersion vector that consists of two parts,
# spatial (ICAR) and non-spatial (covariate) parts
beta_z = rnorm(q,0,.5)
PHI_x =  ( Z%*%beta_z)+ rnorm(n,0,.05)
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

```
```r
DIC complex model:  10538.18 
    simple  model:  12638.14 
WAIC complex model:  10231.69 
     simple  model:  12540.06 

Coefficient complex model (central tendancy):  
                    mean     se_mean       2.5%      97.5%    n_eff      Rhat
(Intercept)   0.72197492 0.002123780  0.4683816  0.9656499 3523.902 1.0000502
WI            2.20080706 0.001055274  2.0742191  2.3263146 3576.027 0.9999368
factor(age)2  0.40110485 0.002308145  0.1020203  0.7034635 4524.108 0.9998764
factor(age)3  0.05914021 0.002645949 -0.2697902  0.4005724 4221.852 0.9998280
factor(age)4 -1.17831876 0.002210866 -1.4284300 -0.9339829 3229.766 0.9998972

Coefficient simple model (central tendancy):  
                   mean     se_mean       2.5%      97.5%    n_eff      Rhat
(Intercept)   1.3753521 0.009061376  0.6689083  2.0948107 1606.660 0.9996256
WI            2.0244127 0.002947587  1.7531254  2.2944601 2204.706 1.0017068
factor(age)2  0.5109182 0.012509091 -0.4187465  1.4936774 1559.313 1.0007132
factor(age)3 -0.2061040 0.014399381 -1.3672445  0.9065953 1573.793 1.0003234
factor(age)4 -1.6258537 0.009362696 -2.3426071 -0.8953636 1569.403 1.0000438

      < complex model spatial parameter(s)> 
Sigma spatial ICAR parameter (central tendancy):  
       mean     se_mean 
0.844007528 0.002925893 

Sigma spatial ICAR parameter (dispersion):  
       mean     se_mean 
0.879940303 0.001501999 

      < simple model spatial parameter(s)> 
Sigma spatial ICAR parameter (central tendancy):  
      mean    se_mean 
0.54166592 0.07986053 
```

A small amount of DIC, WAIC, and beta demonstrates the superiority of the complex model, which considers spatial structure on variance.
