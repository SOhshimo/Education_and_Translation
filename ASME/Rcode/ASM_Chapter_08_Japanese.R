
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -----------------------------------------------------
# 第8章：Comparisons along two classifications
# in a モデル with two factors
# -----------------------------------------------------

# 最終更新：11 June 2024


#######################################################################
# --------------------------------------------------------------------#
# 8.2 データの生成
# --------------------------------------------------------------------#
#######################################################################


set.seed(8)

# 標本数を設定する
nPops <- 5
nHab <- 3
nSample <- 4
n <- nPops * nHab * nSample

# 作成する factor levels
pop <- gl(n = nPops, k = nHab * nSample, length = n)
hab <- gl(n = nHab, k = nSample, length = n)

# Choose 効果
baseline <- 40                       # 切片
pop.eff <- c(-10, -5, 5, 10)         # 集団 効果
hab.eff <- c(5, 10)                  # Hab 効果
interaction.eff <- c(-2, 3, 0, 4, 4, 0, 3, -2) # Interaction 効果
all.eff <- c(baseline, pop.eff, hab.eff, interaction.eff)
sigma <- 3                           # 残差の標準偏差

eps <- rnorm(n, 0, sigma)            # 残差
Xmat <- as.matrix(model.matrix(~ pop * hab) ) # デザイン行列を作成する
Xmat     # デザイン行列を確認し、その意味を理解する
str(Xmat)                            # 図・出力は省略

wing <- as.vector(Xmat %*% all.eff +  eps)
boxplot(wing ~ hab*pop, col = "grey", xlab = "Habitat-by-Population combination",
  ylab = "Wing length", main = "Simulated data set", las = 1,
  ylim = c(20, 70), frame = FALSE)   # 図示 of generated データ
abline(h = 40)

library(lattice)     # latticeパッケージを読み込む
xyplot(wing ~ hab | pop, ylab = "Wing length", xlab = "Habitat",
  main = "Population-specific relationship between wing and habitat type",
  pch = 16, cex = 2, col = rgb(0,0,0,0.4))
xyplot(wing ~ pop | hab, ylab = "Wing length", xlab = "Population", 
  main = "Habitat-specific relationship between wing and population",    
  pch = 16, cex = 2, col = rgb(0,0,0,0.6)) # 図 8.3

# 後で比較するため真値を保存する
truth <- c(all.eff, sigma)
( true.group.means <- c(baseline, baseline + pop.eff, baseline + hab.eff[1],
  baseline + pop.eff + hab.eff[1] + interaction.eff[1:4], baseline + hab.eff[2],
  baseline + pop.eff + hab.eff[2] + interaction.eff[5:8]) )

# 必要なライブラリ
library(ASMbook); library(jagsUI); library(rstan); library(TMB); library(DHARMa)

#######################################################################
# --------------------------------------------------------------------#
# 8.3 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################

# 当てはめる main-効果 モデル and 推定値を保存する
summary( out83.main <- lm(wing ~ hab + pop) )
lm_est_me <- c(coef(out83.main), sigma = sigma(out83.main))

# Interaction-効果 モデル (means parameterization)
summary(out83.intX <- lm(wing ~ pop:hab-1) ) # 図・出力は省略

# Interaction-効果 モデル (treatment contrast parameterization)
summary(out83.int <- lm(wing ~ pop*hab) )

# GOF
# 分位点残差によって適合度を確認する
# 主効果モデル
plot(out83.main) # 比較のため従来型の残差診断を行う
simOut <- simulateResiduals(out83.main, n = 1000, plot = TRUE) # DHARMa
# 交互作用モデル
plot(out83.int) # 比較のため従来型の残差診断を行う
simOut <- simulateResiduals(out83.int, n = 1000, plot = TRUE) # DHARMa

#######################################################################
# --------------------------------------------------------------------#
# 8.4 An aside: を用いる シミュレーション to assess bias and precision
# of an estimator ... and to understand what a 標準 error is
# --------------------------------------------------------------------#
#######################################################################
# 推定値と真値を比較する
lm_est <- c(coef(out83.int), sigma = sigma(out83.int))
comp <- cbind(truth = truth, lm = lm_est)
print(comp, 4)

simrep <- 10^4 # 反復回数
estimates <- array(dim = c(simrep, length(all.eff))) # データ structure to hold 結果
for(i in 1:simrep) {                 # Run シミュレーション simrep times
  if(i %% 100 == 0) cat(paste('\n iter', i))  # 進捗カウンター
  eps <- rnorm(n, 0, sigma)          # 残差
  y <- as.vector(Xmat %*% all.eff + eps) # Assemble データ
  fit.model <- lm(y ~ pop*hab)       # Break down データ
  estimates[i,] <- fit.model$coef    # Keep 値 of coefficients
}

# 比較する true 値 and 平均 of 推定値
data.frame(params = names(fit.model$coefficients), truth = all.eff,
  mean.of.estimates = round(apply(estimates, 2, mean), 3))

# 比較する theoretical and bootstrapped 標準誤差
SE.lm <- summary(out83.int)$coef[,2]
  data.frame(SE.lm = round(SE.lm, 3), bootstrapped.SE = round(apply(estimates, 2, sd), 3))

# 図示 bootstrapped 標本抽出 distributions of four 推定量 (図 8.4)
par(mfrow = c(2,2), mar = c(5,6,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(estimates[,1], xlab = '', freq = F, main = 'Intercept', col = 'grey', breaks = 40)
hist(estimates[,2], xlab = '', freq = F, main = 'pop2', col = 'grey', breaks = 40)
hist(estimates[,6], xlab = '', freq = F, main = 'hab2', col = 'grey', breaks = 40)
hist(estimates[,8], xlab = '', freq = F, main = ' pop2:hab2', col = 'grey', breaks = 40)

#######################################################################
# --------------------------------------------------------------------#
# 8.5 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# 8.5.1 主効果モデル
# --------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing)))

# JAGSモデルファイルを書き出す
cat(file = "model85.main.txt", "
model {
# 事前分布
alpha ~ dnorm(0, 0.0001)             # 切片
beta.pop[1] <- 0                     # 第1水準の効果をゼロに固定する
for(j in 2:5){                       # 残りの因子水準をループ処理する
  beta.pop[j] ~ dnorm(0, 0.0001)
}
beta.hab[1] <- 0                     # 同様に設定する
for(j in 2:3){                       # 残りの因子水準をループ処理する
  beta.hab[j] ~ dnorm(0, 0.0001)
}
sigma ~ dunif(0, 100)
tau <- pow(sigma, -2)

# 尤度
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau)
  mean[i] <- alpha + beta.pop[pop[i]] + beta.hab[hab[i]]
}
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rnorm(1), sigma = rlnorm(1) )}

# 推定対象のパラメータ
params <- c("alpha", "beta.pop", "beta.hab", "sigma")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する, 事後分布を要約する and 推定値を保存する
out85me <- jags(dataList, inits, params, "model85.main.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out85me) # 図・出力は省略
print(out85me, 3)
jags_est_me <- out85me$summary[c(1, 8:9, 3:6, 10), 1]

# 比較する 尤度 and Bayesian 推定値
comp <- cbind(lm = lm_est_me, JAGS = jags_est_me)
print(comp, 4)


# 8.5.2 交互作用モデル
# ---------------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = length(unique(hab)), nPops = length(unique(pop))) )

# JAGSモデルファイルを書き出す
cat(file = "model85.int.txt", "
model {
# 事前分布
for (i in 1:nPops){
  for(j in 1:nHab) {
    group.mean[i,j] ~ dnorm(0, 0.0001)
  }
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# 尤度
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau)
  mean[i] <- group.mean[pop[i], hab[i]]
}
}
")

# 初期値を生成する関数
inits <- function(){list(sigma = rlnorm(1) )}

# 推定対象のパラメータ
params <- c("group.mean", "sigma")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する, 事後分布を要約する and 推定値を保存する
out85ie <- jags(dataList, inits, params, "model85.int.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out85ie) # 図・出力は省略
print(out85ie, 3)
jags_est_ie <- out85ie$summary[1:16,1]

# Align true vals of すべての params
truthMP <- c(true.group.means, sigma) # Truth 対象： means parameterization

# 尤度法とベイズ法の推定値および真値を比較する
lm_est_ie <- c(coef(out83.intX), sigma = sigma(out83.intX))
comp <- cbind(truth = truthMP, lm = lm_est_ie, JAGS = jags_est_ie)
print(comp, 4)


# 8.5.3 予測値の算出
# ---------------------------
par(mfrow = c(1,1), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
ord <- c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15) # 表示順を指定する
plot(ord, out85ie$mean$group.mean, xlab = "Hab-by-Population", las = 1,
  ylab = "Predicted wing length", cex = 2, ylim = c(20, 70), frame = FALSE,
  pch = 16, col = rgb(0,0,0,0.5))
segments(ord, out85ie$q2.5$group.mean, ord, out85ie$q97.5$group.mean,
  col = rgb(0,0,0,0.5), lwd = 2)
abline(h = 40)

#######################################################################
# --------------------------------------------------------------------#
# 8.6 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# 8.6.1 主効果モデル
# --------------------------

library(nimble)

# Bundle and summarize データ (same as 対象： JAGS)
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing)))

# NIMBLEモデルコードを定義する
model86.main <- nimbleCode( {
# 事前分布
alpha ~ dnorm(0, 0.0001)             # 切片
beta.pop[1] <- 0                     # 第1水準の効果をゼロに固定する
for(j in 2:5){                       # 残りの因子水準をループ処理する
  beta.pop[j] ~ dnorm(0, 0.0001)
}
beta.hab[1] <- 0                     # 同様に設定する
for(j in 2:3){                       # 残りの因子水準をループ処理する
  beta.hab[j] ~ dnorm(0, 0.0001)
}
sigma ~ dunif(0, 100)
tau <- pow(sigma, -2)

# 尤度
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau) 
  mean[i] <- alpha + beta.pop[pop[i]] + beta.hab[hab[i]]
}
} )

# 初期値
inits <- function(){ list(alpha = rnorm(1), sigma = rlnorm(1) )}

# 監視するパラメータは前と同じ
params <- c("alpha", "beta.pop", "beta.hab", "sigma")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 30 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
system.time( out86me <- 
    nimbleMCMC(code = model86.main,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out86me)   # 図・出力は省略
(nsum <- nimble_summary(out86me, params))     # 図・出力は省略
nimble_est_me <- nsum[c(1, 8:9, 3:6, 10), 1]


# 8.6.2 交互作用モデル
# ---------------------------------

# Bundle and summarize データ (same as 対象： JAGS)
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = length(unique(hab)), nPops = length(unique(pop))) )

# NIMBLEモデルコードを定義する
model86.int <- nimbleCode( {
# 事前分布
for (i in 1:nPops){
  for(j in 1:nHab) {
    group.mean[i,j] ~ dnorm(0, 0.0001)
  }
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# 尤度
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau) 
  mean[i] <- group.mean[pop[i], hab[i]]
}
} )

# 初期値
inits <- function(){list(sigma = rlnorm(1) )}

# 監視するパラメータは前と同じ
params <- c("group.mean", "sigma")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 20 sec), 収束を診断する, 事後分布を要約する and 結果を保存する
system.time( out86ie <- 
    nimbleMCMC(code = model86.int,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out86ie)   # 図・出力は省略
(nsum <- nimble_summary(out86ie, params))     # 図・出力は省略
nimble_est_ie <- nsum[1:16, 1]

#######################################################################
# --------------------------------------------------------------------#
# 8.7 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# 8.7.1 主効果モデル
# --------------------------
# Bundle and summarize データ (same as 対象： JAGS)
dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing))


# Write text file with モデル description in Stan
cat(file = "model87_main.stan",      # この行はRコード
"
data {
  int n;
  vector[n] wing;
  int pop[n];
  int hab[n];
}

parameters {
  real alpha;
  vector[5] beta_pop_raw;
  vector[3] beta_hab_raw;
  real<lower=0> sigma;
}

transformed parameters {
  vector[5] beta_pop = beta_pop_raw;
  vector[3] beta_hab = beta_hab_raw;
  beta_pop[1] = 0;
  beta_hab[1] = 0;
}

model {
  vector [n] mu;
  beta_pop_raw ~ normal(0, 100);
  beta_hab_raw ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  for (i in 1:n){
    mu[i] = alpha + beta_pop[pop[i]] + beta_hab[hab[i]];
    wing[i] ~ normal(mu[i], sigma);
  }
}                                    // ここまでがStanコード
" )

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# パラメータ to save
params <- c("alpha", "beta_pop", "beta_hab", "sigma")

# Stanを実行する (原著記載の実行時間（ART） 50/3 sec)
system.time(
  out87me <- stan(file = "model87_main.stan", data = dataList,
    pars = params, chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out87me)            # 図・出力は省略
print(out87me, dig = 2)              # 図・出力は省略
stan_est_me <- summary(out87me)$summary[c(1, 8, 9, 3, 4, 5, 6, 10), 1]


# 8.7.2 交互作用モデル
# ---------------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = nHab, nPops = nPops))

# Write text file with モデル description in Stan
cat(file = "model87_int.stan",       # この行はRコード
"
data {
  int n;
  int nPops;
  int nHab;
  vector[n] wing; int hab[n];
  int pop[n];
}

parameters {
  matrix[nPops, nHab] group_mean;
  real<lower=0> sigma;
}

model {
  to_vector(group_mean) ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  for (i in 1:n) {
    wing[i] ~ normal(group_mean[pop[i], hab[i]], sigma);
  }
}                                    // ここまでがStanコード
" )

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 48/3 sec)
system.time(
  out87ie <- stan(file = "model87_int.stan", data = dataList, chains = nc, iter = ni,
  warmup = nb, thin = nt) )
rstan::traceplot(out87ie)            # 図・出力は省略
print(out87ie, dig = 2)              # 図・出力は省略
# 推定値を保存する
# StanとJAGSでは要約表の並び順が異なるため順序を調整する
stan_est_ie <- summary(out87ie)$summary[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15,16),1]

#######################################################################
# --------------------------------------------------------------------#
# 8.8 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################

# 8.8.1 主効果モデル
# --------------------------

# Design matrix 対象： main-効果 モデル
head( Xmat <- model.matrix(~ hab + pop))

# NLL 対象： main-効果 正規 linear モデル with two factors
NLL <- function(param, y, Xmat) {
  beta <- param[1:7]                 # 切片 and 傾き coefficients
  sigma <- exp(param[8])             # 残差の標準偏差
  mu <- Xmat %*% beta                # Multiply 各 row of Xmat by beta
  LL <- dnorm(y, mu, sigma, log = TRUE) # 対数尤度 対象： 各 obs
  NLL <- -sum(LL)                    # データ全体の負の対数尤度
  return(NLL)
}

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c(mean(wing), rep(0, 6), log(sd(wing)))
names(inits) <- c(names(coef(out83.main)), 'log-sigma')
out88me <- optim(inits, NLL, y = wing, Xmat = Xmat, hessian = TRUE, method = 'BFGS')
get_MLE(out88me, 4)
diy_est_me <- c(out88me$par[-8], exp(out88me$par[8]))


# 8.8.2 交互作用モデル
# ---------------------------------

# Design matrix 対象： interaction-効果 モデル
head( Xmat <- model.matrix(~ pop:hab-1)) # 図・出力は省略

# NLL 対象： interaction-効果 正規 linear モデル with two factors
NLL <- function(param, y, Xmat) {
  beta <- param[1:15]                # Interaction パラメータ
  sigma <- exp(param[16])            # 残差の標準偏差
  mu <- Xmat %*% beta
  LL <- dnorm(y, mu, sigma, log = TRUE) # Log-lik 対象： 各 観測
  NLL <- -sum(LL)                    # データ全体の負の対数尤度
  return(NLL)
}
# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c(rep(mean(wing), 15), log(sd(wing)))
names(inits) <- c(colnames(Xmat), 'log-sigma')
out88.int <- optim(inits, NLL, y = wing, Xmat = Xmat, hessian = TRUE, method = 'BFGS')
get_MLE(out88.int, 4)
diy_est_ie <- c(out88.int$par[-16], exp(out88.int$par[16]))

#######################################################################
# --------------------------------------------------------------------#
# 8.9 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################

# 8.9.1 主効果モデル
# --------------------------
# データをリストにまとめ構造を確認する
Xmat <- model.matrix(~ hab + pop)
str(tmbData <- list(wing = wing, Xmat = Xmat, n = length(wing)))

# TMBモデルファイルを書き出す
cat(file = "model89_main.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(wing);                 // 応答変数
  DATA_MATRIX(Xmat);                 // モデル行列
  DATA_INTEGER(n);                   // 観測数

  // パラメータを定義する
  PARAMETER_VECTOR(beta);            // 切片と傾き
  PARAMETER(log_sigma);              // log(residual standard deviation)

  Type sigma = exp(log_sigma);
  Type LL = 0;                       // 対数尤度の合計を0で初期化する
  
  matrix<Type> mu = Xmat * beta;     // Xmatとbetaの行列積
  
  for (int i=0; i< n; i++){          // 添字は1ではなく0から始まる点に注意
    LL += dnorm(wing(i), mu(i), sigma, true);
  }
  return -LL;                        // 合計対数尤度の負値を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model89_main.cpp")
dyn.load(dynlib("model89_main"))

# パラメータの次元と初期値を指定する
params <- list(beta = rep(0, ncol(Xmat)), log_sigma = 0)

# TMBの目的関数オブジェクトを作成する
out89me <- MakeADFun(data = tmbData, parameters = params,
DLL = "model89_main", silent = TRUE)

# TMBオブジェクトを最適化する, 結果を出力する and 推定値を保存する
opt1 <- optim(out89me$par, fn = out89me$fn, gr = out89me$gr, method="BFGS", hessian = TRUE)
(tsum <- tmb_summary(out89me))
tmb_est <- c(opt1$par[1:7], exp(opt1$par[8]))


# 8.9.2 交互作用モデル
# ---------------------------------

# データをリストにまとめ構造を確認する
str(tmbData <- list(wing = wing,
pop = as.numeric(pop)-1,             # TMBの添字は0から始まる
hab = as.numeric(hab)-1,             # TMBの添字は0から始まる
n= length(wing)))

# TMBモデルファイルを書き出す
cat(file = "model89_int.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(wing);                 // 応答変数
  DATA_IVECTOR(pop);                 // 観測値の集団
  DATA_IVECTOR(hab);                 // 観測値の生息環境区分
  DATA_INTEGER(n);                   // 観測数
  
  // パラメータを定義する
  PARAMETER_MATRIX(group_mean);
  PARAMETER(log_sigma);              // log(residual standard deviation)

  Type sigma = exp(log_sigma);       // Type = match type of function output
  Type LL = 0.0;                     // 対数尤度の合計を0で初期化する

  for (int i= 0; i<n; i++){          // 添字は1ではなく0から始まる点に注意！
    LL += dnorm(wing(i), group_mean(pop(i), hab(i)), sigma, true);
  }
  return -LL;                        // 合計対数尤度の負値を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model89_int.cpp")
dyn.load(dynlib("model89_int"))

# パラメータの次元と初期値を指定する
params <- list(group_mean = matrix(0, nPops, nHab), log_sigma=0)

# TMBの目的関数オブジェクトを作成する
out89ie <- MakeADFun(data = tmbData, parameters = params,
  DLL = "model89_int", silent = TRUE)

# TMBオブジェクトを最適化する, 結果を出力する and 推定値を保存する
starts <- rep(0, nPops*nHab+ 1)
opt2 <- optim(starts, fn = out89ie$fn, gr = out89ie$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out89ie)) # 図・出力は省略
tmb_est_ie <- c(opt2$par[1:15], exp(opt2$par[16]))


# 8.10 パラメータ推定値の比較
# ------------------------------------------

# 比較する interaction-効果 モデル 推定値 with truth
comp <- cbind(truth = truthMP, lm = lm_est_ie, JAGS = jags_est_ie, NIMBLE = nimble_est_ie,
  Stan = stan_est_ie, DIY = diy_est_ie, TMB = tmb_est_ie)
print(comp, 4)

