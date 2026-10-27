
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# For more information, see www.kenkellner.org/ASM
# ---------------------------------------------------------------

# ----------------------------------------------------
# 第12章：Overdispersion, zero inflation, and
# ポアソンGLMのオフセット
# ----------------------------------------------------

# 最終更新：11 June 2024

#######################################################################
# --------------------------------------------------------------------#
# 12.2.1 データの生成
# --------------------------------------------------------------------#
#######################################################################


set.seed(122)
nSites <- 50
n <- 2 * nSites
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))
eps <- rnorm(2*nSites, mean = 0, sd = 0.5)             # 正規 変量 効果
lambda.OD <- exp(0.69 +(0.92*(as.numeric(x)-1) + eps) )
lambda.Poisson <- exp(0.69 +(0.92*(as.numeric(x)-1)) ) # For 比較

# 真のパラメータ値を保存する
truth <- c(Intercept = 0.69, arable = 0.92, eps_sd = 0.5)

C_OD <- rpois(n = n, lambda = lambda.OD)               # Counts with OD
C_Poisson <- rpois(n = n, lambda = lambda.Poisson)     # Counts without OD

par(mfrow = c(1,2)) # 図 12.1
boxplot(C_OD ~ x, col = "grey", xlab = "Land-use", main = "With overdispersion",
ylab = "Hare count", las = 1, ylim = c(0, max(C_OD)), frame = FALSE)
boxplot(C_Poisson ~ x, col = "grey", xlab = "Land-use", main = "Without overdispersion",
ylab = "Hare count", las = 1, ylim = c(0, max(C_OD)), frame = FALSE )

# 必要なライブラリを読み込む
library(ASMbook); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
# 12.2.2 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################


glm.fit.no.OD <- glm(C_OD ~ x, family = poisson)
glm.fit.with.OD <- glm(C_OD ~ x, family = quasipoisson)
summary(glm.fit.no.OD)
glm_est <- c(coef(glm.fit.no.OD), NA) # 保存する point 推定値

# Canned function 対象： ポアソン lognormal
# install.packages("PLNmodels")

library(PLNmodels)
PLN.fit <- PLN(matrix(C_OD, ncol = 1) ~ x)
PLN_est <- c(coef(PLN.fit), sqrt(sigma(PLN.fit)))

#######################################################################
# --------------------------------------------------------------------#
# 12.2.3 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# JAGSモデルファイルを書き出す
cat(file = "model12.2.3.txt", "
model {
# 事前分布
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
tau <- pow(sigma, -2)
sigma ~ dt(0, 0.1, 1)T(0, ) # Half-Cauchy

# 尤度
for (i in 1:n) {
  C_OD[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha + beta *x[i] + eps[i]
  eps[i] ~ dnorm(0, tau)
}
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1), sigma = rlnorm(1, 0.1))}

# 推定対象のパラメータ
params <- c("alpha", "beta", "sigma")

# MCMCの設定
na <- 2000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する, summarize and 結果を保存する
out12.2.3 <- jags(dataList, inits, params, "model12.2.3.txt",
n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out12.2.3)                # 図・出力は省略
print(out12.2.3, 3)
jags_est <- as.numeric(out12.2.3$mean[1:3]) # 推定値を保存する

# 尤度法とベイズ法の推定値および真値を比較する
comp <- cbind(truth = truth, glm = glm_est, PLN = PLN_est, JAGS = jags_est)
print(comp, 4)

#######################################################################
# --------------------------------------------------------------------#
# 12.2.4 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


library(nimble)

# Bundle and summarize データ (same as before)
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# NIMBLEモデルコードを定義する
model12.2.4 <- nimbleCode( {
# 事前分布
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
sigma ~ T(dt(0, 0.1, 1), 0, )

# 尤度
for (i in 1:n) {
  C_OD[i] ~ dpois(lambda[i]) 
  log(lambda[i]) <- alpha + beta *x[i] + eps[i]
  eps[i] ~ dnorm(0, sd = sigma)
}
} )

# 初期値
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1), sigma = rlnorm(1, 0.1))}

# 監視するパラメータは前と同じ
params <- c("alpha", "beta", "sigma")

# MCMCの設定
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# NIMBLEを実行する (原著記載の実行時間（ART） 25 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
system.time( out12.2.4 <- 
    nimbleMCMC(code = model12.2.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out12.2.4) # 図・出力は省略
(nsum <- nimble_summary(out12.2.4, params))   # 図・出力は省略
nimble_est <- nsum[1:3,1]                     # 推定値を保存する

#######################################################################
# --------------------------------------------------------------------#
# 12.2.5 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Stanモデルファイルを書き出す
cat(file = "model12_2_5.stan", "
data{
  int n;
  array[n] int C_OD;
  array[n] int x;
}

parameters{
  real alpha;
  real beta;
  real<lower=0> sigma;
  vector[n] eps;
}

model{
  vector[n] lambda;
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  for (i in 1:n){
    eps[i] ~ normal(0, sigma);
    lambda[i] = exp(alpha + beta * x[i] + eps[i]);
    C_OD[i] ~ poisson(lambda[i]);
  }
}
")

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 60/6 sec), 収束を評価する, 結果を出力する table and 推定値を保存する
system.time(
  out12.2.5 <- stan(file = "model12_2_5.stan", data = dataList,
                    warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out12.2.5)          # 図・出力は省略
print(out12.2.5, dig = 3)            # 図・出力は省略
stan_est <- summary(out12.2.5)$summary[1:3,1]

#######################################################################
# --------------------------------------------------------------------#
# 12.2.6 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################


# Bundle データ
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# 負の対数尤度の定義 対象： overdispersed ポアソン GLM with 1 covariate
# Note this is technically a ポアソン GLMM, as in Chapter 14
NLL <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  sigma <- exp(param[3])
  L <- numeric(n)
  for(i in 1:data$n){
    L[i] <- integrate(function(eps){
    lambda <- exp(alpha + beta * data$x[i] + eps)
    dpois(data$C_OD[i], lambda)* dnorm(eps, 0, sigma)},
      lower = -Inf, upper = Inf)$value
  }
  NLL <- -sum(log(L))
  NLL
}

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c('alpha' = log(mean(C_OD)), 'beta(arable)' = 0, 'log.sigma' = 0)
out12.2.6 <- optim(inits, NLL, data = dataList, hessian = TRUE)
get_MLE(out12.2.6, 5)
diy_est <- c(out12.2.6$par[1:2], exp(out12.2.6$par[3]))

#######################################################################
# --------------------------------------------------------------------#
# 12.2.7 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# TMBモデルファイルを書き出す
cat(file = "model12_2_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{


  // 入力データを定義する
  DATA_VECTOR(C_OD);                 // 応答変数
  DATA_VECTOR(x);                    // 共変量
  DATA_INTEGER(n);                   // 観測値数

  // パラメータを定義する
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(eps);

  Type sigma = exp(log_sigma);
  Type LL = 0.0;                     // 対数尤度を0で初期化する
  for (int i= 0; i<n; i++){
    LL += dnorm(eps(i), Type(0.0), sigma, true);
    Type lambda = exp(alpha + beta * x(i) + eps(i));

    // 各観測値の対数尤度を計算し合計に加える
    LL += dpois(C_OD(i), lambda, true);
  }
  return -LL;                        // 負の対数尤度を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model12_2_7.cpp")
dyn.load(dynlib("model12_2_7"))

# Provide dimensions and starting 値 対象： パラメータ (incl. the 変量 効果 'eps')
params <- list(alpha = 0, beta = 0, log_sigma = 0, eps = rep(0, dataList$n))

# TMBの目的関数オブジェクトを作成する
out12.2.7 <- MakeADFun(data = dataList,
                       parameters = params,
                       random = "eps",
                       DLL = "model12_2_7", silent = TRUE)

# TMBオブジェクトを最適化する, print and 結果を保存する
starts <- rep(0, 3)
opt <- optim(starts, fn = out12.2.7$fn, gr = out12.2.7$gr, method = "BFGS")
(tsum <- tmb_summary(out12.2.7)) # 図・出力は省略
tmb_est <- c(tsum[1:2,1], exp(tsum[3,1]))

#######################################################################
# --------------------------------------------------------------------#
# 12.2.8 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################


# 比較する truth with 推定値 から すべての engines and methods
comp <- cbind(truth = truth, glm = glm_est, PLN = PLN_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)

