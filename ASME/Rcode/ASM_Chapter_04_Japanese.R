
# ---------------------------------------------------------------
# 書籍『Applied Statistical Modeling for Ecologists』のコード 
# Kéry・Kellner 著、Elsevier、2024年
# ---------------------------------------------------------------

# -----------------------------------------------------------------------
# 第4章 -- 汎用モデル推定エンジン入門 
#                および「平均値のモデル」
# -----------------------------------------------------------------------

# 最終更新：2024年6月11日
#  install.packages("jagsUI")  #  本文にはないので追記
library(jagsUI)

#######################################################################
# --------------------------------------------------------------------#
#   4.2 データ生成
# --------------------------------------------------------------------#
#######################################################################
# オスのハヤブサの体重について、標本サイズの異なる2つのデータを生成
set.seed(4)
y10 <- rnorm(n = 10, mean = 600, sd = 30) # 10羽の標本
y1000 <- rnorm(n = 1000, mean = 600, sd = 30) # 1000羽の標本
# 後で比較するため、データ生成に用いた真のパラメータ値を保存
truth <- c(mean = 600, sd = 30)
# データを描画（図4-2）
xlim = c(450, 750)
par(mfrow = c(1, 2), mar = c(6, 6, 6, 3), cex.lab = 1.5, cex.axis = 1.5)
hist(y10, col = 'grey ', xlim = xlim, main = 'オスのハヤブサ10羽の体重（g）')
hist(y1000, col = 'grey', xlim = xlim, main = 'オスのハヤブサ1000羽の体重（g）')


#######################################################################
# --------------------------------------------------------------------#
#   4.3 Rの既存関数による解析
# --------------------------------------------------------------------#
#######################################################################

summary(out4.3 <- lm(y10 ~ 1)) # 小さいデータセット：表示省略
summary(out4.3 <- lm(y1000 ~ 1)) # 大きいデータセット

# 推定値を保存
lm_est <- c(coef(out4.3), sigma = sigma(out4.3))


#######################################################################
# --------------------------------------------------------------------#
#   4.4 JAGS
# --------------------------------------------------------------------#
#######################################################################
#   4.4.2 JAGSによるモデル推定
# -------------------------------

# library(jagsUI)

# データをリストにまとめ、構造を表示
str(dataList <- list(mass = y1000, n = length(y1000)))

# JAGSのモデルファイルを作成
cat(file = "model4.4.txt", " # この行はRコード
model { # ここからBUGS言語のコード
# 事前分布
pop.mean ~ dunif(0, 5000) # 母平均
precision <- 1 / pop.var # 精度＝1／分散
pop.var <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# 尤度
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
} # BUGS言語のコードはここまで
") # ここから再びRコード

# 初期値を生成する関数
inits <- function(){
  list(pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))
}

# 記録対象のパラメータ
params <- c("pop.mean", "pop.sd", "pop.var")

# MCMCの設定
na <- 1000 # 適応段階の反復回数
ni <- 12000 # 各チェーンの反復回数（バーンインを含む）
nb <- 2000 # バーンインとして破棄する反復回数
nc <- 4 # チェーン数
nt <- 1 # 間引き間隔（nt = 1は間引きなし）

# JAGSを実行（実行時間の目安：1分）、進捗バーを表示
out4.4 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model4.4.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = FALSE)

# JAGSを並列実行（実行時間の目安：1分未満）し、収束を確認
out4.4 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model4.4.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out4.4) # 図4-3を作成

print(out4.4, 3) # 推定されたモデルオブジェクトの要約を表示

names(out4.4)

hist(out4.4$summary[,8]) # 要約表の8列目にあるRhat値
which(out4.4$summary[,8] > 1.1) # 未収束のパラメータを調べる（この例では該当なし）

summary(out4.4$samples) # 表示省略

par(mfrow = c(3,1), mar = c(5,5,4,2)) # 表示省略
pars <- colnames(out4.4$samples[[1]])[1:3]
for (i in pars){
  matplot(as.array(out4.4$samples)[,i,], type = 'l', lty = 1,
          col = c("green", "red", "blue", "orange"),
          xlab = "反復回数", ylab = "値", main = i, frame = FALSE)
}

par(mfrow = c(3, 2), mar = c(5,5,4,2)) # 図4-4
for (i in pars){
  samps <- as.matrix(out4.4$samples[,i])
  hist(samps, col = 'grey', breaks = 100, xlab = i, main = "")
  plot(density(samps), type = 'l', lwd = 2, col = 'gray40', main = "",
       xlab = i, frame = FALSE)
}

sims <- as.matrix(out4.4$samples)
head(sims[,"pop.mean"] < 597) # 多数の論理判定のうち先頭を表示
mean(sims[,"pop.mean"] < 597) # 母平均が597未満である事後確率

# 図4-5：2変量の条件を満たす事後確率を計算し、同時事後分布を描画
test.true <- sims[,"pop.mean"] < 600 & sims[,"pop.sd"] > 30
mean(test.true)
par(mfrow = c(1,1))
plot(sims[,"pop.mean"], sims[,"pop.sd"], pch = 16, col = rgb(0,0,0,0.3),
cex = 0.8, frame = FALSE)
points(sims[test.true,"pop.mean"], sims[test.true,"pop.sd"], pch = 16,
col = rgb(1,0,0,0.6), cex = 0.8)
abline(h = 30, col = 'red')
abline(v = 600, col = 'red')

# 別の方法で2変量の事後分布を表示
pairs(sims[,1:3])

apply(as.matrix(out4.4$samples[,1:2]), 2, summary)

apply(as.matrix(out4.4$samples[,1:2]), 2, sd)

# 最尤推定とベイズ推定の結果を真値と比較
jags_est <- out4.4$summary[1:2,1]
print(cbind(truth = truth, lm = lm_est, JAGS = jags_est), 5)

#######################################################################
# --------------------------------------------------------------------#
#   4.5.2 NIMBLEによるモデル推定
# --------------------------------------------------------------------#
#######################################################################
# NIMBLEを読み込む
# install.packages("nimble")
library(nimble)

# データをリストにまとめる（JAGSと同様）
str(dataList <- list(mass = y1000, n = length(y1000))) # 表示省略

# NIMBLEのモデルコードを定義
model4.5 <- nimbleCode( {
# 事前分布と線形モデル
pop.mean ~ dunif(0, 5000) # 正規分布の精度パラメータによる表現
precision <- 1 / pop.variance # 精度＝1／分散
pop.variance <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# 尤度
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
} )

# JAGSと同じ初期値生成関数を利用できる
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))

# 記録するパラメータ：先ほどと同様
params <- c("pop.mean", "pop.sd", "pop.variance")

# MCMCの設定
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# NIMBLEを実行（実行時間の目安：30秒）
system.time(
out4.5 <- nimbleMCMC(code = model4.5, constants = dataList,
  inits = inits, monitors = params, nburnin = nb, niter = ni,
  thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(out4.5, head) # 各チェーンの先頭の値を表示（表示省略）

par(mfrow = c(1, 3));
  coda::traceplot(out4.5$samples) # 表示省略
out4.5$summary$all.chains # 4チェーンをまとめた事後分布の要約

library(ASMbook)
nsum <- nimble_summary(out4.5$samples, params) # 要約表
round(nsum, 3)
nimble_est <- nsum[1:2,1] # 推定値を保存

# BUGS形式のコードからNIMBLEモデルを作成
rawModel <- nimbleModel(code = model4.5, constants = dataList, inits = inits())

# MCMCアルゴリズムを設定：標準のMCMC設定を作成
# 選択したパラメータを記録対象にする
mcmcConfig <- configureMCMC(rawModel, monitors = params)

# MCMCアルゴリズムの関数を構築
rawMCMC <- buildMCMC(mcmcConfig)

system.time(rawMCMC$run(10)) # 10回サンプリング（約10秒）
as.matrix(rawMCMC$mvSamples) # サンプルを表示（表示省略）

compModel <- compileNimble(rawModel)

compMCMC <- compileNimble(rawMCMC, project = rawModel)

system.time(compMCMC$run(10)) # 表示省略、約0.02秒

system.time(
samples <- runMCMC(compMCMC, niter = ni, nburnin = nb, thin = nt,
  nchains = nc, samplesAsCodaMCMC = TRUE) )

# サンプルを確認
lapply(samples, head) # 表示省略

# 周辺事後分布の要約を作成
(nsum <- nimble_summary(samples)) # 表示省略

# トレースプロット
par(mfrow = c(2,2)); coda::traceplot(samples) # 表示省略

#######################################################################
# --------------------------------------------------------------------#
#   4.6.2 Stanによるモデル推定
# --------------------------------------------------------------------#
#######################################################################
# StanのRパッケージを読み込む
# install.packages("rstan")
library(rstan)

# データをリストにまとめ、構造を表示
str(dataList <- list(n = length(y1000), mass = y1000))

# Stan言語でモデルを記述したテキストファイルを作成
cat(file = "model4_6.stan", # この行まではRコード
"data { // ここからStanコード
  int n; // すべてのデータの型を定義
  vector[n] mass; // ベクトルの次元も定義
}

parameters { // パラメータについても同様に宣言
  real pop_mean;
  real <lower = 0> pop_sd;
}

model {
  // 事前分布
  pop_mean ~ normal(0, 1000);
  pop_sd ~ cauchy(0, 10);

  // 尤度
  for (i in 1:n){
    mass[i] ~ normal(pop_mean, pop_sd);
  }
} // Stanコードはここまで
" )

# HMCの設定
ni <- 1200 ; nb <- 200 ; nc <- 4 ; nt <- 1

# Stanを実行（実行時間の目安：33秒／3秒）
system.time(
out4.6 <- stan(file = "model4_6.stan", data = dataList,
  chains = nc, iter = ni, warmup = nb, thin = nt) )

# 収束を確認
rstan::traceplot(out4.6) # 表示省略

# 事後分布の要約を表示
print(out4.6)

# 推定値を保存
stan_est <- summary(out4.6)$summary[1:2,1]

str(out4.6) # 表示省略

sims <- extract(out4.6)
lapply(sims, head) # 表示省略

sims <- As.mcmc.list(out4.6) # 先頭の「A」が大文字である点に注意
lapply(sims, head) # 表示省略
par(mfrow = c(2,2)); coda::traceplot(sims) # 表示省略

#######################################################################
# --------------------------------------------------------------------#
#   4.7.2 Rで最尤法によりモデルを推定
# --------------------------------------------------------------------#
#######################################################################
# install.packages("ASMbook")
library(ASMbook)

dnorm(y10[1], mean = 600, sd = 30)
dnorm(y10[1], mean = 550, sd = 25)

prod(dnorm(y10, mean = 600, sd = 30))

prod(dnorm(y10, mean = 550, sd = 25))

sum(dnorm(y10, mean = 600, sd = 30, log = TRUE))

# 切片のみの正規線形モデルの対数尤度（LL）を定義
LL <- function(param, y) {
  pop.mean <- param[1] # 第1パラメータは平均（mu）
  pop.sd <- param[2] # 第2パラメータは標準偏差（sigma）
  sum(dnorm(y, mean = pop.mean, sd = pop.sd, log = TRUE))
}

LL(c(600, 30), y10)

starts <- c(pop.mean = 550, pop.sd = 25)

# 平均が一定の正規線形モデルの負の対数尤度（NLL）を定義
NLL <- function(param, y) {
  pop.mean <- param[1] # 第1パラメータは平均（mu）
  pop.sd <- param[2] # 第2パラメータは標準偏差（sigma）
  -sum(dnorm(y, mean = pop.mean,
    sd = pop.sd, log = TRUE)) # マイナス符号に注意
}

out4.7 <- optim(starts, fn = NLL, hessian = TRUE, y = y1000)

out4.7$convergence

out4.7$par
diy_est <- out4.7$par # 推定値を保存

vcov <- solve(out4.7$hessian) # ヘッセ行列の逆行列を計算
ASE <- sqrt(diag(vcov)) # 漸近標準誤差を求める
ASE

get_MLE(out4.7, 5)


# 4.8 Template Model Builderによる最尤推定
# ---------------------------------------------------
# （コードなし）


#   4.8.1 Template Model Builderの紹介
# ----------------------------------------------
# （コードなし）


#   4.8.2 TMBによるモデル推定
# -------------------------------
# install.packages("TMB")
library(TMB)

# データをまとめる
str(dataList <- list(y = y1000, n = length(y1000)))

# TMBのモデルファイルを作成
cat(file = "model4_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type> ::operator() ()
{
  // 入力データを定義
  DATA_INTEGER(n); // 標本サイズ
  DATA_VECTOR(y); // 観測値
  
  // パラメータを定義
  PARAMETER(pop_mean); // 平均
  PARAMETER(log_pop_sd); // log（標準偏差）
  Type pop_sd = exp(log_pop_sd); // Typeは関数の出力型に合わせる

  Type LL = 0.0; // 対数尤度の合計を0で初期化
  for (int i = 0; i < n; i ++ ){ // 添字は1ではなく0から始まる点に注意
    // 観測値の対数尤度を計算
    // trueを指定すると対数密度を返す
    LL += dnorm(y(i), pop_mean, pop_sd, true);
  }
  return -LL; // 対数尤度の合計に負号を付けて返す
}
")

# TMB関数をコンパイルして読み込む
compile("model4_8.cpp") # 大量のコンパイルメッセージが出ても心配不要
dyn.load(dynlib("model4_8"))

# パラメータの初期値などを指定
params <- list(pop_mean = 0, log_pop_sd = 0)

# TMBオブジェクトを作成
out4.8 <- MakeADFun(data = dataList, parameters = params, random = NULL,
  DLL = "model4_8", silent = TRUE)

# TMBオブジェクトを最適化
opt <- optim(out4.8$par, fn = out4.8$fn, gr = out4.8$gr,
  method = "BFGS", hessian = TRUE)

opt$par

summary(sdreport(out4.8))

tsum <- tmb_summary(out4.8) # 表示省略
tmb_est <- tsum[,1] # 結果を保存
tmb_est[2] <- exp(tmb_est[2]) # 標準偏差を対数スケールから元のスケールに変換


# 4.9 推定エンジンの比較とまとめ
# ------------------------------------------------

comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 5)
