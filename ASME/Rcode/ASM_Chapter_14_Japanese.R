# 日本語学習用注釈版：元のR・JAGS・NIMBLE・Stan・TMBの実行コードは変更しない。
# 元の英文コメントを参照可能な形で残し、各節・主要処理に日本語解説を付ける。
# Poisson GLMM：群別の正規ランダム切片と傾きを積分または事後推定する。


# ---------------------------------------------------------------
# 日本語注記：Code 対象 the book 'Applied Statistical Modeling 対象 Ecologists'
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -------------------------------------------------------
# Chapter 14  --  Poisson generalized linear mixed model,
#                 or Poisson GLMM
# -------------------------------------------------------

# Last changes: 11 June 2024


# 【節の日本語案内】14.1 はじめに
# 14.1 Introduction
# -----------------
# 日本語注記：(コードなし)
# (no code)

#######################################################################
# --------------------------------------------------------------------#
# 14.2  データ生成 
# --------------------------------------------------------------------#
#######################################################################

# 処理の要点：乱数シードを固定する。
set.seed(14)
nPops <- 16
nYears <- 30
n <- nPops * nYears                  # n = 480
pop <- gl(n = nPops, k = nYears)

orig.year <- rep(1:nYears, nPops)
year <- (orig.year-1)/29             # Squeeze between 0 and 1

# 処理の要点：回帰式から設計行列を作成する。
Xmat <- model.matrix(~ pop * year - 1 - year)
print(Xmat[1:91,], 2)                # Print top 91 rows of 480

# 日本語注記：Choose values 対象 hyperparams and draw Normal random numbers
# Choose values for hyperparams and draw Normal random numbers
mu.alpha <- 3                        # Mean of intercepts
sigma.alpha <- 1                     # SD of intercepts
mu.beta <- -2                        # Mean of slopes
sigma.beta <- 0.6                    # SD of slopes
# 処理の要点：正規分布から乱数を生成する。
alpha <- rnorm(n = nPops, mean = mu.alpha, sd = sigma.alpha)
# 処理の要点：正規分布から乱数を生成する。
beta <- rnorm(n = nPops, mean = mu.beta, sd = sigma.beta)
all.effects <- c(alpha, beta)        # All together

# 日本語注記：真のパラメータ値を保存
# Save true parameter values
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta,
  sigma.alpha = sigma.alpha, sigma.beta = sigma.beta)

lin.pred <- Xmat[,] %*% all.effects  # Value of lin.predictor (eta)
# 処理の要点：ポアソン分布から乱数を生成する。
C <- rpois(n = n, lambda = exp(lin.pred)) # Exponentiate and add Poisson noise
hist(C, col = "grey")                # Inspect what we’ve created (not shown)

# 日本語注記：必要なライブラリを読み込む
# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB); library(lattice);
library(glmmTMB); library(DHARMa); library(abind); library(pracma)

xyplot(C ~ orig.year | pop, ylab = "Red-backed shrike counts", xlab = "Year",
  pch = 16, cex = 1.2, col = rgb(0, 0, 0, 0.4))


#######################################################################
# --------------------------------------------------------------------#
# 14.3 Rの既存関数による尤度解析 
# --------------------------------------------------------------------#
#######################################################################

gtmb.data <- data.frame(C = C, year = year, pop = pop) # bundle data
# 処理の要点：一般化線形混合モデルを既存関数で推定する。
out14.3 <- glmmTMB(C ~ year + (year || pop), data = gtmb.data, family = poisson)
summary(out14.3)                     # Inspect results
sds <- attr(VarCorr(out14.3)$cond$pop, 'stddev')       # Save results
gtmb_est <- c(fixef(out14.3)$cond, sds)

# 日本語注記：表示 推定s of ランダム効果
# Print estimates of random effects
ranef(out14.3) # ran-ef. estimates from zero-mean dist. (not shown)
coef(out14.3) # 'full' ran-ef estimates (not shown)

# 日本語注記：Assess goodness-of-fit of the model (図・出力は省略)
# Assess goodness-of-fit of the model (not shown)
par(mfrow = c(1, 2)) # Distribution of estimates of alpha and beta
hist(coef(out14.3)$cond$pop[,1], main = 'alpha', breaks = 12)
hist(coef(out14.3)$cond$pop[,2], main = 'beta', breaks = 12)
library(DHARMa)
# 処理の要点：シミュレーション残差に基づくモデル診断を実行する。
simOut <- simulateResiduals(out14.3, n = 1000, plot = TRUE)

#######################################################################
# --------------------------------------------------------------------#
# 14.4 JAGSによるベイズ解析 
# --------------------------------------------------------------------#
#######################################################################
# 日本語注記：データをリストにまとめ構造を確認
# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# 日本語注記：JAGSモデルファイルを作成
# Write JAGS model file
cat(file = "model14.4.txt", "
model {
# Priors
for (i in 1:nPops){
  # Models for the sets of random effects alpha and beta
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # Intercepts
  beta[i] ~ dnorm(mu.beta, tau.beta)    # Slopes
}
mu.alpha ~ dnorm(0, 0.0001)             # Hyperparam. for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 0.1, 1)T(0,)        # half-Cauchy prior
mu.beta ~ dnorm(0, 0.0001)              # Hyperparam. for random slopes
tau.beta <- pow(sigma.beta, -2)
sigma.beta ~ dt(0, 0.1, 1)T(0,)         # half-Cauchy prior

# 'Likelihood'
for (i in 1:n){
  C[i] ~ dpois(lambda[i])
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]] * year[i])
  # log(lambda[i]) <- alpha[pop[i]] + beta[pop[i]]* year[i] # same
}
# Posterior predictive simulations (used in DHARMa GoF assessment)
for (i in 1:n) {
  C.new[i] ~ dpois(lambda[i])
}
}
")

# 日本語注記：初期値を生成する関数
# Function to generate starting values
inits <- function(){
# 処理の要点：正規分布から乱数を生成する。
list(mu.alpha = rnorm(1), mu.beta = rnorm(1),
  sigma.alpha = runif(1), sigma.beta = runif(1))
}

# 日本語注記：推定・保存するパラメータ
# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta", "alpha", "beta", "C.new")

# 日本語注記：MCMCの設定
# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# 日本語注記：Call JAGS (ART <1 min), 確認 convergence and summarize posteriors
# Call JAGS (ART <1 min), check convergence and summarize posteriors
# 処理の要点：JAGSによるMCMC推定を行う。
out14.4 <- jags(dataList, inits, params, "model14.4.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out14.4) # not shown
print(out14.4, 3) # partially shown

# 日本語注記：Distribution of ランダム効果 推定s (alpha, beta)
# Distribution of random effects estimates (alpha, beta)
par(mfrow = c(1,2))                  # Very low power with small samples. . .
hist(out14.4$mean$alpha, main = 'alpha', breaks = 12)
hist(out14.4$mean$beta, main = 'beta', breaks = 12)


# 日本語注記：Do quantile residual assessments (図・出力は省略)
# Do quantile residual assessments (not shown)
C.new <- out14.4$sims.list$C.new
sim <- createDHARMa(simulatedResponse = t(C.new), observedResponse = C,
fittedPredictedResponse = apply(C.new, 2, median), integerResponse = T)
plot(sim)

# 日本語注記：Compare 尤度 with Bayesian 推定s and with truth
# Compare likelihood with Bayesian estimates and with truth
jags_est <- unlist(out14.4$mean[1:4])
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est)
print(comp, 4)

#######################################################################
# --------------------------------------------------------------------#
# 14.5 NIMBLEによるベイズ解析 
# --------------------------------------------------------------------#
#######################################################################

library(nimble)

# 日本語注記：データをリストにまとめ構造を確認
# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# 日本語注記：NIMBLEモデルファイルを作成
# Write NIMBLE model file
# 処理の要点：NIMBLEモデルの確率構造を定義する。
model14.5 <- nimbleCode( {
# 日本語注記：事前分布s
  # Priors
  for (i in 1:nPops){    # Priors for random effects
# 処理の要点：正規分布の確率密度を計算する。
    alpha[i] ~ dnorm(mu.alpha, sd = sigma.alpha) # Intercepts
# 処理の要点：正規分布の確率密度を計算する。
    beta[i] ~ dnorm(mu.beta, sd = sigma.beta)    # Slopes
  }
# 処理の要点：正規分布の確率密度を計算する。
  mu.alpha ~ dnorm(0, sd = 100)                  # Same as in JAGS
  sigma.alpha ~ T(dt(0, 0.1, 1), 0,)

# 処理の要点：正規分布の確率密度を計算する。
  mu.beta ~ dnorm(0, sd = 100)
  sigma.beta ~ T(dt(0, 0.1, 1), 0, )

# 日本語注記：'尤度'
  # 'Likelihood'
  for (i in 1:n){
# 処理の要点：ポアソン分布の尤度への寄与を計算する。
    C[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]] * year[i])
  }
} )


# 日本語注記：初期値
# Inits
inits <- function(){
# 処理の要点：正規分布から乱数を生成する。
  list(mu.alpha = rnorm(1), mu.beta = rnorm(1),
       sigma.alpha = runif(1), sigma.beta = runif(1))
}

# 日本語注記：記録するパラメータ: 前と同じ
# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta", "alpha", "beta")

# 日本語注記：MCMCの設定
# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# 日本語注記：Call NIMBLE (ART 50 sec), 確認 convergence, summarize posteriors and 保存 推定s
# Call NIMBLE (ART 50 sec), check convergence, summarize posteriors and save estimates
system.time( out14.5 <- 
# 処理の要点：NIMBLEでMCMCを実行する。
    nimbleMCMC(code = model14.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out14.5)  # not shown
(nsum <- nimble_summary(out14.5, params))    # not shown
nimble_est <- nsum[1:4,1]                    # Save estimates

#######################################################################
# --------------------------------------------------------------------#
# 14.6 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################
# 日本語注記：データをリストにまとめ構造を確認
# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# 日本語注記：Stanモデルファイルを作成
# Write Stan model
cat(file = "model14_6.stan", "
data{
  int n;
  int nPops;
  array[n] int C;
  vector[n] year;
  array[n] int pop;
}

parameters{
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[nPops] alpha;
  vector[nPops] beta;
}

model{
  vector[n] lambda;
  mu_alpha ~ normal(0, 100);
  sigma_alpha ~ cauchy(0, sqrt(10));
  mu_beta ~ normal(0, 100);
  sigma_beta ~ cauchy(0, sqrt(10));
  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
    beta[i] ~ normal(mu_beta, sigma_beta);
  }
  for (i in 1:n){
    lambda[i] = exp(alpha[pop[i]] + beta[pop[i]] * year[i]);
    C[i] ~ poisson(lambda[i]);
  }
}
")

# 日本語注記：HMCの設定
# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# 日本語注記：Call STAN (ART 67/12 sec), assess convergence, show results and 保存 推定s
# Call STAN (ART 67/12 sec), assess convergence, show results and save estimates
system.time(
# 処理の要点：StanによるMCMC推定を行う。
out14.6 <- stan(file = "model14_6.stan", data = dataList,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out14.6) # not shown
print(out14.6, dig = 3) # not shown
stan_est <- summary(out14.6)$summary[1:4,1] # Save estimates

#######################################################################
# --------------------------------------------------------------------#
# 14.7 自作尤度関数による最尤推定 
# --------------------------------------------------------------------#
#######################################################################

# Definition of custom functions f and g (note f is needed in g,
#   and g is needed in definition of NLL function below)
f <- function(alpha_j, beta_j, mu.alpha, sig.alpha, 
  mu.beta, sig.beta, C_j, year_j){
    lambda <- exp(alpha_j + beta_j * year_j)
# 処理の要点：正規分布の確率密度を計算する。
    prod(dpois(C_j, lambda)) * dnorm(alpha_j, mu.alpha, sig.alpha) *
# 処理の要点：正規分布の確率密度を計算する。
      dnorm(beta_j, mu.beta, sig.beta)
}

g <- function(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, year_j){
# 処理の要点：2変数の数値積分を行う。
  integral2(fun = f, xmin = -10, xmax = 10, ymin = -10, ymax = 10,
    vectorized = FALSE,
    mu.alpha = mu.alpha, sig.alpha = sig.alpha, mu.beta = mu.beta,
    sig.beta = sig.beta, C_j = C_j, year_j = year_j)$Q
}

# 日本語注記：Define NLL 対象 Poisson GLMM (using custom functions f and g)
# Define NLL for Poisson GLMM (using custom functions f and g)
NLL <- function(pars, data) {
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sig.alpha <- exp(pars[3])
  sig.beta <- exp(pars[4])
  nll <- 0 # Initialize nll at 0 prior to summation in loop below
  for (j in 1:data$nPops){
    # Subset data to just pop j
    C_j <- data$C[data$pop == j]
    year_j <- data$year[data$pop == j]
    lik <- g(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, year_j)
    nll <- nll - log(lik)
  }
  return(nll)
}

# 日本語注記：負の対数尤度を最小化 to find MLEs, get SEs and CIs and 保存 推定s (ART 122 sec)
# Minimize that NLL to find MLEs, get SEs and CIs and save estimates (ART 122 sec)
inits <- c(mu.alpha = 0, mu.beta = 0, log.sigma.alpha = 0, log.sigma.beta = 0)
system.time(
# 処理の要点：目的関数を数値最適化する。
  out14.7 <- optim(inits, NLL, data = dataList, hessian = TRUE,
    method = 'BFGS', control = list(trace = 1, REPORT = 1)) )
get_MLE(out14.7, 4)
diy_est <- c(out14.7$par[1:2], exp(out14.7$par[3:4])) # Save estimates

#######################################################################
# --------------------------------------------------------------------#
# 14.8 TMBによる尤度解析 
# --------------------------------------------------------------------#
#######################################################################

# 出版社から提供のコードが動かなかったので改変（日本語補注なし）

library(TMB)

# Prepare data
tmbData <- dataList
tmbData$pop <- tmbData$pop - 1
str(tmbData)

# Write C++ model
cat(file = "model14_8.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data
  DATA_VECTOR(C);
  DATA_VECTOR(year);
  DATA_IVECTOR(pop);
  DATA_INTEGER(n);
  DATA_INTEGER(nPops);

  // Parameters
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER(log_sigma_alpha);
  PARAMETER(log_sigma_beta);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type sigma_alpha = exp(log_sigma_alpha);
  Type sigma_beta = exp(log_sigma_beta);

  Type LL = 0.0;

  // Random effects
  for (int i = 0; i < nPops; i++) {
    LL += dnorm(alpha(i), mu_alpha, sigma_alpha, true);
    LL += dnorm(beta(i), mu_beta, sigma_beta, true);
  }

  // Poisson likelihood
  for (int i = 0; i < n; i++) {
    Type lambda = exp(alpha(pop(i)) + beta(pop(i)) * year(i));
    LL += dpois(C(i), lambda, true);
  }

  // Return negative log-likelihood
  return -LL;
}
")

# Compile and load
compile("model14_8.cpp")
dyn.load(dynlib("model14_8"))

# Initial parameters
params <- list(
  mu_alpha = 0,
  mu_beta = 0,
  log_sigma_alpha = 0,
  log_sigma_beta = 0,
  alpha = rep(0, tmbData$nPops),
  beta = rep(0, tmbData$nPops)
)

# Create TMB objective function
out14.8 <- MakeADFun(
  data = tmbData,
  parameters = params,
  random = c("alpha", "beta"),
  DLL = "model14_8",
  silent = TRUE
)

# Optimize
opt <- optim(
  out14.8$par,
  fn = out14.8$fn,
  gr = out14.8$gr,
  method = "BFGS"
)

# Check convergence
print(opt$convergence)
print(opt$message)

# Parameter estimates
print(opt$par)

# Extract estimates in the specified order
tmb_est <- c(
  mu_alpha = unname(opt$par["mu_alpha"]),
  mu_beta = unname(opt$par["mu_beta"]),
  sigma_alpha = exp(unname(opt$par["log_sigma_alpha"])),
  sigma_beta = exp(unname(opt$par["log_sigma_beta"]))
)

print(tmb_est)

#######################################################################
# --------------------------------------------------------------------#
# 14.9 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################

# 日本語注記：推定値を真値と比較
# Compare estimates with truth
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)
