
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ----------------------------------------------------------------------
# 第7章：Models with a single categoricalcovariate with more
# 2水準を超える場合
# ----------------------------------------------------------------------

# 最終更新：11 June 2024

#######################################################################
# --------------------------------------------------------------------#
# 7.2.1 データの生成
# --------------------------------------------------------------------#
#######################################################################
# シミュレーションデータを生成する
set.seed(72) # 乱数生成器を初期化する
nPops <- 5 # 集団数
nSample <- 10 # 各集団におけるヘビの個体数
pop.means <- c(50, 40, 45, 55, 60) # 集団別のSVL平均
sigma <- 5 # 残差の標準偏差
n <- nPops * nSample # データの総数
eps <- rnorm(n, 0, sigma) # 残差
pop <- factor(rep(1:5, rep(nSample, nPops))) # 集団の識別子
means <- rep(pop.means, rep(nSample, nPops))
X <- as.matrix(model.matrix(~ pop-1)) # デザイン行列を作成する
X # デザイン行列を確認する
y <- as.numeric(X %*% as.matrix(pop.means) + eps)
# %*% は行列積を表す

# 後で比較するため真値を保存する
truth <- c(pop.means, sigma)
names(truth) <- c(paste0("pop", 1:5), "sigma")

# 図7–2を作成する
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ pop, col = "grey", xlab = "Population", ylab = "SVL", main = "", las = 1, frame = FALSE)

# 必要なライブラリを読み込む
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
# 7.2.2 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################
# 基準水準との差で表す通常の処理対比パラメータ化（出力は省略）
summary(out <- lm(y ~ pop))
# 因子の各水準の平均を直接推定するパラメータ化
summary(out72.2 <- lm(y ~ pop - 1))

# 最小二乗推定値を保存する
lm_est <- c(coef(out72.2), sigma = sigma(out72.2))

# 適合度を確認する 通常の残差と分位点残差を用いる
plot(out72.2) # 比較のため従来型の残差診断を行う
simOut <- simulateResiduals(out72.2, n = 1000, plot = TRUE)

#######################################################################
# --------------------------------------------------------------------#
# 7.2.3 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################
# データをリストにまとめ構造を確認する
# JAGSでは因子popを数値型に変換する必要がある
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# JAGSモデルファイルを書き出す
cat(file = "model72.3.txt", "
model {
# 事前分布
for (i in 1:nPops){ # alphaをベクトルとして定義する
  alpha[i] ~ dnorm(0, 0.0001)
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)

# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mean[i], tau)
  mean[i] <- alpha[pop[i]]
}

# 導出量（効果の推定値）
effe2 <- alpha[2] - alpha[1]
effe3 <- alpha[3] - alpha[1]
effe4 <- alpha[4] - alpha[1]
effe5 <- alpha[5] - alpha[1]

# 仮説検定：任意の対比を定義する
test1 <- (effe2 + effe3) - (effe4 + effe5) # 効果2＋効果3が効果4＋効果5に等しいとき0
test2 <- effe5 - 2 * effe4 # effe5がeffe4の2倍に等しいとき0
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rnorm(nPops, mean = mean(y)), sigma = rlnorm(1) )}

# 推定対象のパラメータ
params <- c("alpha", "sigma", "effe2", "effe3", "effe4", "effe5", "test1", "test2")
# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） <1 sec), 収束を診断する and 事後分布を要約する
out72.3 <- jags(dataList, inits, params, "model72.3.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out72.3) # 図・出力は省略
print(out72.3, 3)

# 尤度法とベイズ法の推定値および真値を比較する
jags_est <- out72.3$summary[1:6, 1]
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)

#######################################################################
# --------------------------------------------------------------------#
# 7.2.4 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

library(nimble)

# Bundle and summarize データ (same as before)
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# NIMBLEモデルコードを定義する
model72.4 <- nimbleCode( {
# 事前分布
for (i in 1:nPops){                 # alphaをベクトルとして定義する
  alpha[i] ~ dnorm(0, 0.001)
}
sigma ~ dunif(0, 10)

# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mean[i], sd = sigma) # ここでは標準偏差sigmaを直接指定する
  mean[i] <- alpha[pop[i]]
}

# 導出量
effe2 <- alpha[2] - alpha[1]
effe3 <- alpha[3] - alpha[1]
effe4 <- alpha[4] - alpha[1]
effe5 <- alpha[5] - alpha[1]

# 仮説検定：任意の対比を定義する
test1 <- (effe2+effe3) - (effe4+effe5) # 効果2＋効果3が効果4＋効果5に等しいとき0
test2 <- effe5 - 2 * effe4             # effe5がeffe4の2倍に等しいとき0
}
)

# 初期値
inits <- function(){ list(alpha = rnorm(nPops, mean = mean(y)), sigma = rlnorm(1) )}

# 監視するパラメータは前と同じ
params <- c("alpha", "sigma", "effe2", "effe3", "effe4", "effe5", "test1", "test2")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 24 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
system.time(
  out72.4 <- nimbleMCMC(code = model72.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out72.4)     # 図・出力は省略
(nsum <- nimble_summary(out72.4, params))       # 図・出力は省略
nimble_est <- nsum[1:6,1]

#######################################################################
# --------------------------------------------------------------------#
# 7.2.5 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# Bundle and summarize データ (same as 対象： NIMBLE)
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# Write text file with モデル description in BUGS language
cat(file = "model72_5.stan",
"
data {
  int <lower = 1> n; // Declare すべての データ
  int <lower = 1> nPops;
  vector[n] y;
  array[n] int <lower = 1> pop;
}

parameters { // パラメータの型を定義する
  vector [nPops] alpha; // 各集団の平均値
  real <lower = 0> sigma; // 標準偏差は正の値に制限する
}

model {
  // 事前分布
  alpha ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  // 尤度
  for(i in 1:n) {
    y[i] ~ normal(alpha[pop[i]], sigma);
  }
}

generated quantities {
  // 導出量
  real effe2 = alpha[2] - alpha[1];
  real effe3 = alpha[3] - alpha[1];
  real effe4 = alpha[4] - alpha[1];
  real effe5 = alpha[5] - alpha[1];
  // 仮説検定：任意の対比を定義する
  // 効果2＋効果3が効果4＋効果5に等しいときtest1は0
  real test1 = (effe2 + effe3) - (effe4 + effe5);
  real test2 = effe5 - 2 * effe4;
}
" )

# HMCの設定
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# Stanを実行する (原著記載の実行時間（ART） 46/2 sec)
system.time(
out72.5 <- stan(file = "model72_5.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out72.5) # 図・出力は省略
print(out72.5, dig = 2) # 図・出力は省略
stan_est <- summary(out72.5)$summary[1:6,1]

#######################################################################
# --------------------------------------------------------------------#
# 7.2.6 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################
# 負の対数尤度の定義 （対象： one-factor linear モデル 正規誤差を伴う
NLL <- function(param, y, X) {
  alpha <- param[1:5] # 集団別の平均値
  sigma <- param[6] # 残差の標準偏差
  mu <- X %*% alpha # 各観測値に対応する集団平均を取得する
  LL <- dnorm(y, mu, sigma, log = TRUE) # 各データ点の対数尤度
  NLL <- -sum(LL) # すべてのデータ点の負の対数尤度
  return(NLL)
}

# 平均値によるパラメータ化に対応するデザイン行列を作成する
X <- model.matrix(~ pop - 1)

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c('mu1' = 50, 'mu2' = 50, 'mu3' = 50, 'mu4' = 50, 'mu5' = 50, 'sigma' = 10)
out72.6 <- optim(inits, NLL, y = y, X = X, hessian = TRUE)
get_MLE(out72.6, 4)
diy_est <- out72.6$par

#######################################################################
# --------------------------------------------------------------------#
# 7.2.7 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################
# Bundle and summarize データ (similar as before, except 対象： pop)
str(tmbData <- list(y = y, pop = as.numeric(pop) - 1, n = n))

# TMBモデルファイルを書き出す
cat(file = "model72_7.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(y); // 応答変数
  DATA_IVECTOR(pop); // 集団 index: note IVECTOR, not VECTOR
  DATA_INTEGER(n); // 観測数

  // パラメータを定義する
  PARAMETER_VECTOR(alpha); // 集団別の平均値
  PARAMETER(log_sigma); // log(residual standard deviation)
  Type sigma = exp(log_sigma); // 標準偏差
  Type LL = 0; // 対数尤度の合計を0で初期化する

  // 各データ点の尤度
  for (int i = 0; i < n; i++){ // 添字は1ではなく0から始まる点に注意！
    LL += dnorm(y(i), alpha(pop(i)), sigma, true);
  }

  // Derived 効果 (note index starts at 0!)
  Type effe2 = alpha(1) - alpha(0);
  Type effe3 = alpha(2) - alpha(0);
  Type effe4 = alpha(3) - alpha(0);
  Type effe5 = alpha(4) - alpha(0);
  
  // 任意の対比による仮説検定
  Type test1 = (effe2 + effe3) - (effe4 + effe5);
  Type test2 = effe5 - 2 * effe4;
  
  // Tell TMB to save derived パラメータ 値
  ADREPORT(effe2);
  ADREPORT(effe3);
  ADREPORT(effe4);
  ADREPORT(effe5);
  ADREPORT(test1);
  ADREPORT(test2);
  return -LL; // 合計対数尤度の負値を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model72_7.cpp") # コンパイル時の大量の出力
dyn.load(dynlib("model72_7"))

# パラメータの次元と初期値を指定する
params <- list(alpha = rep(0, nPops), log_sigma = 0)

# TMBの目的関数オブジェクトを作成する
out72.7 <- MakeADFun(data = tmbData, parameters = params,
                     DLL = "model72_7", silent = TRUE)

# TMBオブジェクトを最適化する, 結果を出力する and 推定値を保存する
opt <- optim(out72.7$par, fn = out72.7$fn, gr = out72.7$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out72.7))
tmb_est <- c(opt$par[1:5], exp(opt$par[6]))


# 7.2.8 パラメータ推定値の比較
# ---------------------------------------------

# 推定結果を真値と比較する
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est,
NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 7.3 変量効果モデル
# ------------------------
# （コードなし）


# 7.3.1 データの生成
# -----------------------

# シミュレーションデータを生成する
set.seed(73)
nPops <- 10 # 集団数：5ではなく10集団を設定する
nSample <- 12 # 各集団におけるヘビの個体数
n <- nPops * nSample # データの総数

pop.grand.mean <- 50 # SVLの全体平均
pop.sd <- 3 # 集団平均の全体平均周りの標準偏差
pop.means <- rnorm(n = nPops, mean = pop.grand.mean, sd = pop.sd)
sigma <- 5 # 残差の標準偏差
eps <- rnorm(n, 0, sigma) # 残差を生成する

pop <- factor(rep(1:nPops, rep(nSample, nPops)))
Xmat <- as.matrix(model.matrix(~ pop - 1))
y <- as.numeric(Xmat %*% as.matrix(pop.means) + eps)

# 後で比較するため真値を保存する
truth <- c(pop.grand.mean = pop.grand.mean, pop.sd = pop.sd, residual.sd = sigma)

# 図7–3を作成する
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ pop, col = "grey", xlab = "Population", ylab = "SVL", main = "", las = 1, frame = FALSE)
abline(h = pop.grand.mean)


# 7.3.2 Rの既存関数を用いた尤度解析
# ------------------------------------------------------

library(glmmTMB) # glmmTMBを読み込む

# 当てはめる モデル and inspect 結果
summary(out73.2ML <- glmmTMB(y ~ 1 + 1 | pop, REML = FALSE)) # 最尤法（ML）
summary(out73.2REML <- glmmTMB(y ~ 1 + 1 | pop, REML = TRUE)) # 制限付き最尤法（REML）
ranef(out73.2REML) # Look: THESE are the estimated 変量 効果 !

summary(out73.2ML <- glmmTMB(y ~ 1 + 1 | pop, REML = FALSE)) # 最尤法（ML）
ranef(out73.2REML) # Look: THESE are the estimated 変量 効果 !

# からの推定値を保存： glmmTMB and make interim 比較
gtmb_est <- c(fixef(out73.2ML)$cond, sqrt(VarCorr(out73.2ML)$cond$pop[1]),
sigma(out73.2ML))
gtmb_REMLest <- c(fixef(out73.2REML)$cond, sqrt(VarCorr(out73.2REML)$cond$pop[1]),
sigma(out73.2REML))
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest)
print(comp, 4) # 表示省略, but note REML 推定値 of pop.sd greater

# 分位点残差によって適合度を確認する
simOut <- simulateResiduals(out73.2ML, n = 1000, plot = TRUE)


# 7.3.3 Bayesian による解析： JAGS ... and a quick 固定-変量 比較
# ---------------------------------------------------------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# JAGSモデルファイルを書き出す
cat(file = "model73.3.txt", "
model {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # 集団平均に対する事前分布
  pop.mean[i] ~ dnorm(pop.grand.mean, pop.tau) # 変量効果！
  # Pop. 効果 as derived quantities
  effe[i] <- pop.mean[i] - pop.grand.mean # lme4と同じ形式のランダム効果
}
pop.grand.mean ~ dnorm(0, 0.0001) # SVLの全体平均に対する超事前分布
pop.sd ~ dt(0, 0.01, 1)T(0,) # 集団効果の標準偏差に対する超事前分布
residual.sd ~ dunif(0, 10) # 残差標準偏差の事前分布

pop.tau <- pow(pop.sd, -2)
residual.tau <- pow(residual.sd, -2)

# 「尤度」
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], residual.tau)
  mu[i] <- pop.mean[pop[i]]
}
}
")

# 初期値を生成する関数
inits <- function(){ list(pop.grand.mean = runif(1, 1, 100),
  pop.sd = runif(1), residual.sd = runif(1)) }

# 推定対象のパラメータ
params <- c("pop.grand.mean", "pop.sd", "residual.sd", "pop.mean", "effe")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） 1 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
out73.3 <- jags(dataList, inits, params, "model73.3.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI:: traceplot(out73.3) # 図・出力は省略
print(out73.3, 3)
jags_est <- out73.3$summary[1:3,1]

# 尤度法とベイズ法の推定値および真値を比較する
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest, JAGS = jags_est)
print(comp, 4)

# 確認する shape of 事後分布 distribution 対象： pop.sd (図 7.4)
hist(out73.3$sims.list$pop.sd, col = 'grey', breaks = 30, main = '', freq = F)

# 計算する the 事後分布 mode of pop.sd as a better point 推定値
library(MCMCglmm)
posterior.mode(mcmc(out73.3$sims.list$pop.sd))

# JAGSモデルファイルを書き出す
cat(file = "model73.3fix.txt", "
model {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # 事前分布 対象： 集団 means: no estimated hyperparameters now !
  pop.mean[i] ~ dnorm(0, 0.0001)
}
residual.sd ~ dunif(0, 10) # 残差標準偏差の事前分布
residual.tau <- pow(residual.sd, -2)

# 「尤度」
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], residual.tau)
  mu[i] <- pop.mean[pop[i]]
}
}
")

# 初期値を生成する関数
inits <- function(){ list(pop.mean = runif(nPops, 30, 70), residual.sd = runif(1)) }

# 推定対象のパラメータ
params <- c("pop.mean", "residual.sd")

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する and 事後分布を要約する
out73.3fix <- jags(dataList, inits, params, "model73.3fix.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out73.3fix) # 図・出力は省略
print(out73.3fix, 3) # also 表示省略

# 比較する 固定- and the 変量-効果 集団 means (図 7.5)
off <- 0.15
plot((1:nPops)-off, pop.means, pch = 16, cex = 2, col = 'red', frame = FALSE,
  xlab = 'Population', ylab = 'Mean SVL in population', las = 1, ylim = c(38, 62))
abline(h = out73.3$mean$pop.grand.mean, lty = 2, lwd = 3, col = 'blue')
points(1:nPops, out73.3$mean$pop.mean, pch = 16, cex = 2, col = 'blue')
segments(1:nPops, out73.3$q2.5$pop.mean, 1:nPops, out73.3$q97.5$pop.mean, lwd = 2, col = 'blue')
points((1:nPops) + off, out73.3fix$mean$pop.mean, pch = 16, cex = 2, col = 'black')
segments((1:nPops) + off, out73.3fix$q2.5$pop.mean, (1:nPops) + off,
  out73.3fix$q97.5$pop.mean, lwd = 2, col = 'black')
legend('topleft', pch = 16, cex = 1.5, col = c("red", "blue", "black"), legend = c
("Truth", "Random effects", "Fixed effects"), bty = 'n')


# 7.3.4 NIMBLEによるベイズ解析
# -------------------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# NIMBLEモデルコードを定義する
model73.4 <- nimbleCode( {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # 集団平均に対する事前分布
  pop.mean[i] ~ dnorm(pop.grand.mean, sd = pop.sd)
  # Pop. 効果 as derived quantities
  effe[i] <- pop.mean[i] - pop.grand.mean
}
pop.grand.mean ~ dnorm(0, sd = 1000) # SVLの全体平均に対する超事前分布
pop.sd ~ T(dt(0, 0.01, 1), 0,) # 集団効果の標準偏差に対する超事前分布
residual.sd ~ dunif(0, 10) # 残差標準偏差の事前分布

# 「尤度」
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], sd = residual.sd)
  mu[i] <- pop.mean[pop[i]]
}
}
)

# 初期値
inits <- function(){ list(pop.grand.mean = runif(1, 0, 100), 
  pop.sd = runif(1), residual.sd = runif(1)) }

# 監視するパラメータは前と同じ
params <- c("pop.grand.mean","pop.sd","residual.sd","pop.mean","effe")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 32 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
system.time(
  out73.4 <- nimbleMCMC(code = model73.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(3,3)); coda::traceplot(out73.4)     # 図・出力は省略
(nsum <- nimble_summary(out73.4, params))       # 図・出力は省略
nimble_est <- nsum[1:3,1]


# 7.3.5 Stanによるベイズ解析
# -----------------------------------

# Bundle and summarize データ (same as before)
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# Write text file with モデル description in BUGS language
cat(file = "model73_5.stan",
"data { // 入力データを定義する
  int n;
  int nPops;
  vector[n] y;
  array[n] int pop; // 長さnの整数ベクトル
}

parameters { // パラメータを定義する
  real <lower = 0> pop_grand_mean;
  real <lower = 0> pop_sd;
  real <lower = 0> residual_sd;
  vector <lower = 0> [nPops] pop_mean;
}

model {
  // 各観測値の期待値を格納するベクトル
  vector[n] mu;

  // 事前分布
  pop_grand_mean ~ normal(0, 100);
  pop_sd ~ cauchy(0, 10);
  residual_sd ~ cauchy(0, 10);

  for (i in 1:nPops) {
    pop_mean[i] ~ normal(pop_grand_mean, pop_sd);
  }

  // 「尤度」
  for(i in 1:n) {
    mu[i] = pop_mean[pop[i]];
    y[i] ~ normal(mu[i], residual_sd);
  }
}

generated quantities {
  vector[nPops] effe;
  for (i in 1:nPops){
    effe[i] = pop_mean[i] - pop_grand_mean;
  }
}
" )

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 45/2 sec), 収束を評価する, print and 推定値を保存する
system.time(
out73.5 <- stan(file = "model73_5.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out73.5) # 図・出力は省略
print(out73.5, dig = 2) # 図・出力は省略
stan_est <- summary(out73.5)$summary[1:3,1]


# 7.3.6 自作の最尤推定
# ---------------------------------------------------

# 負の対数尤度の定義 （対象： モデル with one 変量-効果 factor and 正規 errors
NLL <- function(pars) {
  pop.grand.mean <- pars[1]
  pop.sd <- exp(pars[2])
  sigma <- exp(pars[3])
  LL <- 0
  for (i in 1:nPops){
    # Subset データ to just pop i
    ysub <- y[pop == i]
    L <- integrate(function(r){
      tot <- 1 # Starting 値 対象： product over J
      # Iterate over 各 sample j (1–12) in pop i
      for (j in 1:nSample){
        tot = tot * dnorm(ysub[j], pop.grand.mean + r, sigma)
      }
      tot <- tot * dnorm(r, 0, pop.sd)
      tot
    }, lower = -Inf, upper = Inf, subdivisions = 20)$value
    LL <- LL + log(L)
  }
  return(-LL)
}

# 負の対数尤度を最小化する to 最尤推定値を求める, also 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c('pop.grand.mean' = mean(y), 'log(pop.sd)' = 0, 'log(residual.sd)' = log(sd(y)))
out73.6 <- optim(inits, NLL, hessian = TRUE, method = 'BFGS')
get_MLE(out73.6, 4)
exp(out73.6$par[2:3]) # 標準偏差を指数変換で元の尺度に戻す（名称中のlogは無視）
diy_est <- c(out73.6$par[1], exp(out73.6$par[2:3]))


# 7.3.7 TMBを用いた尤度解析
# ------------------------------------

# Bundle and summarize データ; adjust pop so it starts at 0 instead of 1
tmbData <- list(y = y, pop = as.numeric(pop)-1, nPops = nPops, n = n)

# TMBモデルファイルを書き出す
cat(file = "model73_7.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(y); // 応答変数
  DATA_IVECTOR(pop); // 集団 index: note IVECTOR, not VECTOR
  DATA_INTEGER(nPops); // 集団数
  DATA_INTEGER(n); // 観測数
  
  // パラメータを定義する
  PARAMETER(pop_grand_mean);
  PARAMETER(log_pop_sd)
  PARAMETER(log_residual_sd);
  PARAMETER_VECTOR(pop_mean);
  
  Type pop_sd = exp(log_pop_sd);
  Type residual_sd = exp(log_residual_sd);
  
  Type LL = 0; // 対数尤度の合計を0で初期化する
  
  // ランダム切片
  vector <Type> effe(nPops); // 集団効果
  for (int i = 0; i < nPops; i++){
    LL += dnorm(pop_mean(i), pop_grand_mean, pop_sd, true);
    effe(i) = pop_mean(i) - pop_grand_mean;
  }

  // 観測値の尤度
  for (int i = 0; i < n; i++){ // 添字は1ではなく0から始まる点に注意！
    LL += dnorm(y(i), pop_mean(pop(i)), residual_sd, true);
  }
  ADREPORT(effe);
  return -LL; // 合計対数尤度の負値を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model73_7.cpp")
dyn.load(dynlib("model73_7"))

# パラメータの次元と初期値を指定する
params <- list(pop_grand_mean = 0, log_pop_sd = 0, log_residual_sd = 0,
               pop_mean = rep(0, tmbData$nPops))

# TMBの目的関数オブジェクトを作成する
out73.7 <- MakeADFun(data = tmbData,
  parameters = params, random = "pop_mean",
  DLL = "model73_7", silent = TRUE)

# TMBオブジェクトを最適化する and 結果を出力する (including 変量 効果)
opt <- optim(out73.7$par, fn = out73.7$fn, gr = out73.7$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out73.7)) # 図・出力は省略

# からの推定値を保存： TMB
tmb_est <- c(opt$par[1], exp(opt$par[2:3]))


# 7.3.8 パラメータ推定値の比較
# ---------------------------------------------

# 推定結果を真値と比較する and previous 推定値
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 7.4 まとめと展望
# -----------------------
# （コードなし）
