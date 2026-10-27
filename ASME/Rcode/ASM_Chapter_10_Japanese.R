
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# --------------------------------------------
# 第10章：Linear mixed-effects model
# --------------------------------------------

# 最終更新：11 June 2024


#######################################################################
# --------------------------------------------------------------------#
# 10.2 データの生成
# --------------------------------------------------------------------#
#######################################################################

set.seed(10)
nPops <- 56 # 集団数
nSample <- 10 # Number of vipers in 各 pop
n <- nPops * nSample # データの総数
pop <- gl(n = nPops, k = nSample) # 集団の識別子

# 体長（cm）
orig.length <- runif(n, 45, 70)
mn <- mean(orig.length)
sd <- sd(orig.length)
cat("Mean and sd used to normalise original length:", mn, sd, "\n\n")
lengthN <- (orig.length - mn)/sd # Nは正規化（normalized）を意味する
hist(lengthN, col = "grey") # 表示省略

Xmat <- model.matrix(~ pop * lengthN - 1 - lengthN)
print(Xmat[1:21,], dig = 2) # 出力する top 21 rows (表示省略)

mu.alpha <- 260
sigma.alpha <- 20
mu.beta <- 60
sigma.beta <- 30

alpha <- rnorm(n = nPops, mean = mu.alpha, sd = sigma.alpha)
beta <- rnorm(n = nPops, mean = mu.beta, sd = sigma.beta)
all.pars <- c(alpha, beta) # Put them すべての together

sigma <- 30 # 残差の標準偏差
lin.pred <- Xmat[,] %*% all.pars # 線形予測子の値
eps <- rnorm(n = n, mean = 0, sd = sigma) # 残差
mass <- lin.pred + eps # response = lin.pred + residual
mass <- drop(mass) # 1×nの行列を長さnのベクトルに変換する
hist(mass, col = "grey") # Inspect what we’ve created (表示省略)

# 保存する true 値 対象： comparisons later
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta,
sigma.alpha = sigma.alpha, sigma.beta = sigma.beta, residual.sd = sigma)

library(lattice)
xyplot(mass ~ lengthN | pop, xlab = 'Length', ylab = 'Mass', main = 'Realized mass-length relationships',
  pch = 16, cex = 1.2, col = rgb(0, 0, 0, 0.4))

library(ASMbook); library(jagsUI); library(rstan); library(TMB)


#######################################################################
# --------------------------------------------------------------------#
# 10.3.1 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################


library(glmmTMB)
gtmbData <- data.frame(mass = mass, lengthN = lengthN, pop = pop)
out10.3.ML <- glmmTMB(mass ~ lengthN + (1 | pop), data = gtmbData, REML = FALSE,
                      control = glmmTMBControl(optimizer=optim, 
                                               optArgs=list(method="L-BFGS-B")))  # 最尤法（ML）
out10.3.REML <- glmmTMB(mass ~ lengthN + (1 | pop), data = gtmbData, REML = TRUE) # 制限付き最尤法（REML）
summary(out10.3.ML)
summary(out10.3.REML)

# からの推定値を保存： ML and から REML 当てはめ of 変量-intercepts モデル
gtmbML_est <- c(fixef(out10.3.ML)$cond,
                attr(VarCorr(out10.3.ML)$cond$pop, "stddev"),
                sigma(out10.3.ML))
gtmbREML_est <- c(fixef(out10.3.REML)$cond,
                  attr(VarCorr(out10.3.REML)$cond$pop, "stddev"),
                  sigma(out10.3.REML))

# モデルの適合度を確認する（結果は省略）
hist(coef(out10.3.ML)$cond$pop[,1], main = 'alpha', breaks = 12)
library(DHARMa)
simOut <- simulateResiduals(out10.3.ML, n = 1000, plot = TRUE)

#######################################################################
# --------------------------------------------------------------------#
# 10.3.2 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(mass = mass, pop = as.numeric(pop), lengthN = lengthN, nPops = nPops, n = n) )

# JAGSモデルファイルを書き出す
cat(file = "model10.3.2.txt", "
model {
# 事前分布
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # ランダム切片
}

mu.alpha ~ dnorm(0, 0.000001) # ランダム切片の平均
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 1, 0.0001)I(0,) # ランダム切片の標準偏差

mu.beta ~ dnorm(0, 0.000001) # 共通の傾き
tau.residual <- pow(sigma.residual, -2) # 残差の精度（分散の逆数）
sigma.residual ~ dt(0, 1, 0.0001)I(0,) # 残差の標準偏差

# 尤度
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual) # 観測される確率変数
  mu[i] <- alpha[pop[i]] + mu.beta * lengthN[i] # 期待値
}

# シミュレーション of new データ not conditional on 変量 効果 対象： GoF assessment を用いる DHARM
for (i in 1:nPops){
  alpha.new[i] ~ dnorm(mu.alpha, tau.alpha)
}
for (i in 1:n) {
  mass.new.cond[i] ~ dnorm(mu[i], tau.residual) # 条件付き
  mass.new.un[i] ~ dnorm(mu.new[i], tau.residual) # 無条件
  mu.new[i] <- alpha.new[pop[i]] + mu.beta * lengthN[i]
}
}
")

# 初期値を生成する関数
inits <- function(){list(mu.alpha = rnorm(1, 0, 1),
  mu.beta = rnorm(1, 0, 1), sigma.alpha = rlnorm(1),
  sigma.residual = rlnorm(1)) }

# 推定対象のパラメータ
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.residual", "alpha",
  "mass.new.cond", "mass.new.un")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する, 事後分布を要約する and 推定値を保存する
out10.3.2 <- jags(dataList, inits, params, "model10.3.2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
par(mfrow = c(2, 2)); jagsUI::traceplot(out10.3.2) # 図・出力は省略
print(out10.3.2, 2)
jags_est <- out10.3.2$summary[c(1:4),1]

# 尤度法とベイズ法の推定値および真値を比較する
comp <- cbind(truth = truth[-4], gtmbML = gtmbML_est, JAGS = jags_est)
print(comp, 4)

# Distribution of 変量 効果 推定値 (alpha)
hist(out10.3.2$mean$alpha, main = 'alpha')

# Quantile 残差 assessments (表示省略)
mass.new.un <- out10.3.2$sims.list$mass.new.un # 無条件
mass.new.cond <- out10.3.2$sims.list$mass.new.cond # 条件付き

sim.un <- createDHARMa(simulatedResponse = t(mass.new.un), observedResponse = mass,
  fittedPredictedResponse = apply(mass.new.un, 2, median), integerResponse = FALSE)
plot(sim.un) # Unconditional new データ
sim.cond <- createDHARMa(simulatedResponse = t(mass.new.cond), observedResponse = mass,
  fittedPredictedResponse = apply(mass.new.cond, 2, median), integerResponse = FALSE)
plot(sim.cond) # New データ conditional on estimated 変量 効果

#######################################################################
# --------------------------------------------------------------------#
# 10.3.3 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


library(nimble)

# データ構造を再確認する
str(dataList)

# NIMBLEモデルコードを定義する
model10.3.3 <- nimbleCode( {
# 事前分布
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha)   # ランダム切片
}

mu.alpha ~ dnorm(0, 0.000001)  # ランダム切片の平均
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ T(dt(0, 1, 0.0001), 0, )    # ランダム切片の標準偏差

mu.beta ~ dnorm(0, 0.000001)   # 共通の傾き
tau.residual <- pow(sigma.residual, -2)   # 残差の精度（分散の逆数）
sigma.residual ~ T(dt(0, 1, 0.0001), 0, ) # 残差の標準偏差

# 尤度
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual) # 観測される確率変数
  mu[i] <- alpha[pop[i]] + mu.beta * lengthN[i]  # 期待値
}
} )

# 初期値
inits <- function(){list(mu.alpha = rnorm(1, 0, 1), 
  mu.beta = rnorm(1, 0, 1), sigma.alpha = rlnorm(1), 
  sigma.residual = rlnorm(1)) }

# 監視するパラメータは前と同じ
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.residual", "alpha")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 32 sec), 収束を診断する, 事後分布を要約する and save point 推定値
system.time( out10.3.3 <-
    nimbleMCMC(code = model10.3.3,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out10.3.3)
(nsum <- nimble_summary(out10.3.3, params))   # 図・出力は省略
nimble_est <- nsum[1:4,1]

#######################################################################
# --------------------------------------------------------------------#
# 10.3.4 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# データ構造を再確認する
str(dataList) # 図・出力は省略

# Stanモデルファイルを書き出す
cat(file = "model10_3_4.stan", "
data {
  int n; // 観測値数
  int nPops; // 集団数
  vector[n] mass; // 応答変数
  vector[n] lengthN; // 共変量
  array[n] int pop; // 集団 assignment of 各 obs
}

parameters {
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_residual;
  vector[nPops] alpha;
}

model {
  vector[n] mu; // 期待値

  mu_alpha ~ normal(0, 1000);
  mu_beta ~ normal(0, 1000);
  sigma_alpha ~ cauchy(0, 100);
  sigma_residual ~ cauchy(0, 100);

  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
  }

  // 「尤度」
  for (i in 1:n){
    mu[i] = alpha[pop[i]] + mu_beta * lengthN[i];
    mass[i] ~ normal(mu[i], sigma_residual);
  }
}
")

# 初期値を生成する関数
inits <- function(){list(mu_alpha = rnorm(1, 0, 1),
  mu_beta = rnorm(1, 0, 1), sigma_alpha = rlnorm(1),
  sigma_residual = rlnorm(1)) }

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 50 sec/22 sec) and 結果を保存する
system.time(
out10.3.4 <- stan(file = "model10_3_4.stan", data = dataList, init = inits,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out10.3.4) # 図・出力は省略
print(out10.3.4, dig = 3) # 図・出力は省略
stan_est <- summary(out10.3.4)$summary[1:4,1]

#######################################################################
# --------------------------------------------------------------------#
# 10.3.5 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################


f <- function(alpha.j, mu.alpha, sigma.alpha, mu.beta, sigma.residual,
  mass.j, lengthN.j){
  mu <- alpha.j + mu.beta * lengthN.j
  prod(dnorm(mass.j, mu, sigma.residual)) *
  dnorm(alpha.j, mu.alpha, sigma.alpha)
}

f <- function(alpha.j, mu.alpha, sigma.alpha, mu.beta, sigma.residual,
              mass.j, lengthN.j){
  out <- numeric(length(alpha.j))
  for (i in 1:length(alpha.j)){
    mu <- alpha.j[i] + mu.beta * lengthN.j
    out[i] <- prod(dnorm(mass.j, mu, sigma.residual)) *
    dnorm(alpha.j[i], mu.alpha, sigma.alpha)
  }
  out
}

g <- function(mu.alpha, mu.beta, sigma.alpha, sigma.residual,
              mass.j, lengthN.j){
  integrate(f, lower = 100, upper = 400, mu.alpha = mu.alpha,
            sigma.alpha = sigma.alpha,
            mu.beta = mu.beta, sigma.residual = sigma.residual,
            mass.j = mass.j, lengthN.j = lengthN.j)$value
}

# 負の対数尤度の定義 対象： 変量-intercepts モデル 正規誤差を伴う
NLL <- function(pars, data){
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sigma.alpha <- exp(pars[3])
  sigma.residual <- exp(pars[4])
  LL <- 0 # 負の対数尤度を0で初期化する
  
  for (j in 1:data$nPops){
    mass.j <- data$mass[data$pop == j]
    lengthN.j <- data$lengthN[data$pop == j]
    L <- g(mu.alpha, mu.beta, sigma.alpha, sigma.residual,
           mass.j, lengthN.j) # 尤度 対象： pop j
    LL <- LL + log(L)
  }
  
  return(-LL)
}

# 負の対数尤度を最小化する to 最尤推定値を求める and 標準誤差と信頼区間を求める
inits <- c(mu.alpha = mean(dataList$mass), mu.beta = 0,
           log.sigma.alpha = 1,
           log.sigma.residual = log(sd(dataList$mass)))
out10.3.5 <- optim(inits, NLL, hessian = TRUE, data = dataList)
get_MLE(out10.3.5, 4)

# 保存する MLEs
diy_est <- out10.3.5$par
diy_est[3:4] <- exp(diy_est[3:4]) # Get SD params on ilog scale

#######################################################################
# --------------------------------------------------------------------#
# 10.3.6 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################

# データをリストにまとめ構造を確認する
tmbData <- dataList
tmbData$pop <- tmbData$pop-1 # popの添字を0始まりに変換する
str(tmbData)

# TMBモデルファイルを書き出す
cat(file = "model10_3_6.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(mass); // 応答変数
  DATA_VECTOR(lengthN); // 共変量
  DATA_IVECTOR(pop); // 集団 index
  DATA_INTEGER(nPops); // 集団数
  DATA_INTEGER(n); // 観測値数
  
  // パラメータを定義する
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER(log_sigma_alpha);
  PARAMETER(log_sigma_residual);
  PARAMETER_VECTOR(alpha);

  Type sigma_alpha = exp(log_sigma_alpha);
  Type sigma_residual = exp(log_sigma_residual);
  
  Type LL = 0.0; // 対数尤度を0で初期化する

  // Distribution of 変量 intercepts – add to 対数尤度
  for (int i = 0; i < nPops; i++){ // 添字は1ではなく0から始まる点に注意！
    LL += dnorm(alpha(i), mu_alpha, sigma_alpha, true);
  }
  
  for (int i = 0; i < n; i++){
    Type mu = alpha(pop(i)) + mu_beta * lengthN(i);
    // 各観測値の対数尤度を計算し合計に加える
    LL += dnorm(mass(i), mu, sigma_residual, true);
  }
  
  return -LL; // 負の対数尤度を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model10_3_6.cpp")
dyn.load(dynlib("model10_3_6"))

# パラメータの次元と初期値を指定する
params <- list(mu_alpha = 0, mu_beta = 0,
               log_sigma_alpha = 0, log_sigma_residual = 0,
               alpha = rep(0, tmbData$nPops))

# 作成する TMB object and tell TMB that alpha is 変量
out10.3.6 <- MakeADFun(data = tmbData, parameters = params,
              random = "alpha", # 変量効果として扱うパラメータを指定する
              DLL = "model10_3_6", silent = TRUE)

# TMBオブジェクトを最適化する and 結果を出力する
opt <- optim(out10.3.6$par, fn = out10.3.6$fn, gr = out10.3.6$gr,
             method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out10.3.6)) # 図・出力は省略
# 保存する TMB 推定値
tmb_est <- c(opt$par[1:2], exp(opt$par[3:4]))

#######################################################################
# --------------------------------------------------------------------#
# 10.3.7 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################


# すべての推定結果を真値と比較する
comp <- cbind(truth = truth[-4], gtmbML = gtmbML_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 10.4 Analysis under a 変量-coefficients モデル without correlation between
# 切片 and 傾き
# ---------------------------------------------------------------------------
# （コードなし）

