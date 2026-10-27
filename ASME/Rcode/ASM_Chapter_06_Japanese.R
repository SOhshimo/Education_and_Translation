
# ---------------------------------------------------------------
# 書籍『Applied Statistical Modeling for Ecologists』のコード 
# Kéry・Kellner著、Elsevier、2024年
# ---------------------------------------------------------------

# -----------------------------------------------------
# 第6章：正規モデルによる2群の比較
# -----------------------------------------------------

# 最終更新：2024年6月11日

#######################################################################
# --------------------------------------------------------------------#
# 6.2 等分散を仮定した2群の比較
# --------------------------------------------------------------------#
#######################################################################
#   6.2.1 データの生成
# -----------------------

# データセットを生成
set.seed(61)
n1 <- 60                             # 雌の個体数
n2 <- 40                             # 雄の個体数
mu1 <- 105                           # 雌の母平均
mu2 <- 77.5                          # 雄の母平均
sigma <- 2.75                        # 2群に共通の母標準偏差
n <- n1+ n2                          # 総標本数
y1 <- rnorm(n1, mu1, sigma)          # 雌のデータ
y2 <- rnorm(n2, mu2, sigma)          # 雄のデータ
y <- c(y1, y2)                       # 2群のデータを結合
x <- rep(c(0,1), c(n1, n2))          # 雄を示す指示変数（雌0、雄1）

# 図を作成（図6.1）
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ x, col = "grey", xlab = "雄（0 = 雌、1 = 雄）", ylab = "翼開長（cm）", las = 1, frame = FALSE)

# 同等のデータを別の表現で生成
set.seed(61)
n <- n1 + n2                         # 総標本数
alpha <- mu1                         # 雌の平均を切片とする
beta <- mu2 - mu1                    # betaは雄の平均と雌の平均との差
E.y <- alpha + beta*x                # 期待値（線形予測子）
y.obs <- rnorm(n = n, mean = E.y, sd = sigma) # ランダムな変動を加える
boxplot(y.obs ~ x, col = "grey", xlab = "雄（0 = 雌、1 = 雄）", ylab = "翼開長（cm）", las = 1, frame = FALSE)

# 比較に用いる真の値を保存
truth <- c(alpha = mu1, beta = mu2 - mu1, sigma = sigma)

# 必要なライブラリを読み込む
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
#   6.2.2 Rの既存関数による尤度に基づく解析
# --------------------------------------------------------------------#
#######################################################################

summary(out6.2.2 <- lm(y ~ x))
lm_est <- c(coef(out6.2.2), sigma(out6.2.2)) # 推定結果を保存

t.test(y ~ x, var.equal = TRUE)

xfac <- factor(rep(c("group1", "group2"), c(n1, n2)))
head(xfac) # 最初の数個の値を確認

summary(lm(y ~ xfac)) # 結果は省略（効果によるパラメータ化／処理対比）
summary(lm(y ~ xfac - 1)) # 結果は省略（群平均によるパラメータ化）

# 残差とモデル適合度の診断（図は省略）
plot(out6.2.2)                       # 従来のピアソン残差による確認
simOut <- simulateResiduals(out6.2.2, n = 1000, plot = TRUE)
plotResiduals(simOut)                # 右側の図と同じ内容
testDispersion(out6.2.2)             # 過分散・過小分散を検定

#######################################################################
# --------------------------------------------------------------------#
#   6.2.3 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################
# データをリストにまとめて構造を確認
str(dataList <- list(y = y, x = x, n = n))

# JAGSのモデルファイルを書き出す
cat(file = "model6.2.3.txt", "
model {
# 事前分布
alpha ~ dnorm(0, 0.0001)             # 切片（雌の平均）
beta ~ dnorm(0, 0.0001)              # 傾き（雄の平均－雌の平均）
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)                 # 共通の標準偏差

# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta * x[i]       # 線形予測子
}

# 派生パラメータ
mu1 <- alpha                         # 雌の平均
mu2 <- alpha + beta                  # 雄の平均
for (i in 1:n){                      # 残差
  residual[i] <- y[i] - mu[i]
}
}")

# 初期値を生成する関数
inits <- function(){list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# 推定するパラメータ
params <- c("alpha", "beta", "sigma", "mu1", "mu2", "residual")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行（参考実行時間：1分未満）し、収束を確認して事後分布を要約
out6.2.3 <- jags(data = dataList, inits = inits, parameters.to.save = params,
model.file = "model6.2.3.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out6.2.3) # 出力は省略
print(out6.2.3, 3)

# 尤度に基づく推定・ベイズ推定・真の値を比較
jags_est <- out6.2.3$summary[1:3,1]
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


# 残差プロット
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(1:100, out6.2.3$mean$residual) # 出力は省略
abline(h = 0)
boxplot(out6.2.3$mean$residual ~ x, col = "grey", xlab = "Male",
  ylab = "翼開長の残差（cm）", las = 1) # 出力は省略
abline(h = 0)

#######################################################################
# --------------------------------------------------------------------#
#   6.2.4 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

library(nimble)

# データをまとめて確認（JAGSと同じ）
dataList <- list(y = y, x = x, n = n)


# NIMBLEのモデルを定義
model6.2.4 <- nimbleCode( {
# 事前分布
alpha ~ dnorm(0, 0.0001)             # 切片（雌の平均）
beta ~ dnorm(0, 0.0001)              # 傾き（雄の平均－雌の平均）
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)                 # 共通の標準偏差

# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta*x[i]         # 線形予測子
}

# 派生量
mu1 <- alpha                         # 雌の平均
mu2 <- alpha + beta                  # 雄の平均
for (i in 1:n){                      # 残差
  residual[i] <- y[i] - mu[i]
}
} )

# 初期値を生成する関数
inits <- function(){list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# 推定するパラメータ
params <- c("alpha", "beta", "sigma", "mu1", "mu2", "residual")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行（参考実行時間：30秒）し、収束確認・事後分布の要約・推定値の保存を行う
system.time(  
  out6.2.4 <- nimbleMCMC(code = model6.2.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out6.2.4) # 出力は省略
(nsum <- nimble_summary(out6.2.4, params))   # 出力は省略
nimble_est <- nsum[1:3,1]                    # 推定値を保存

#######################################################################
# --------------------------------------------------------------------#
#   6.2.5 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをまとめて確認（前と同じ）
dataList <- list(y = y, x = x, n = n)

# モデル記述ファイルを書き出す（原文ではBUGS言語と表記）
# バージョン1：応答変数に対してループを使用
cat(file = "model6_2_5.stan",        # この行はRコード
"data {                              // ここからStanコード
  int<lower=0> n;                    // 全データの型・形式を指定
  vector[n] y;                       // ベクトルの次元も指定
  vector[n] x;                       //
}

parameters {                         // 全パラメータの型・形式を指定
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters{
  vector[n] mu;
  for (i in 1:n){
    mu[i] = alpha + beta * x[i];     // 線形予測子を計算
  }
}

model {
  // 事前分布
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  // 尤度
  for (i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
}

generated quantities {
  vector[n] residuals;
  for (i in 1:n){
    residuals[i] = y[i] - mu[i];
  }
}                                    // Stanコードはここまで
")


# モデル記述ファイルを書き出す（原文ではBUGS言語と表記）
# バージョン2：ベクトル化
cat(file = "model6_2_5.stan",        # この行はRコード
"data {                              // ここからStanコード
  int<lower=0> n;                    // 全データの型・形式を指定
  vector[n] y;                       // ベクトルの次元も指定
  vector[n] x;                       //
}

parameters {                         // 全パラメータの型・形式を指定
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters{
  vector[n] mu;
  mu = alpha + beta * x;             // 線形予測子を計算
}

model {
  // 事前分布
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  // 尤度
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[n] residuals;
  residuals = y - mu;
}                                    // Stanコードはここまで
")

# HMCの設定
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行（原文の参考実行時間：57/4秒）し、収束確認・事後分布の要約・推定値の保存を行う
system.time(
  out6.2.5 <- stan(file = "model6_2_5.stan", data = dataList,
    chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out6.2.5)           # 出力は省略

print(out6.2.5, dig = 2)             # 出力は省略
stan_est <- summary(out6.2.5)$summary[1:3,1]

#######################################################################
# --------------------------------------------------------------------#
#   6.2.6 自作関数による最尤推定
# --------------------------------------------------------------------#
#######################################################################


NLL <- function(params, y, x) {
  alpha <- params[1]
  beta <- params[2]
  sigma <- exp(params[3])            # 対数スケールから元のスケールへ変換
  n <- length(y)                     # データ数
  mu <- LL <- numeric(n)             # 空のベクトルを用意
  # このループはベクトル化すると高速化できる
  for (i in 1:n){
    mu[i] <- alpha + beta * x[i]
    LL[i] <- dnorm(y[i], mean = mu[i], sd = sigma, log = TRUE)
  }
  -sum(LL)
}

# 負の対数尤度を最小化し、最尤推定値と標準誤差を求める
inits <- c('alpha' = 50, 'beta' = 10, 'log_sigma' = 0)
out6.2.6 <- optim(inits, NLL, y = y, x = x, hessian=TRUE)
diy_est <- c(out6.2.6$par[1:2], exp(out6.2.6$par[3])) # 推定値を保存
get_MLE(out6.2.6, 5)

#######################################################################
# --------------------------------------------------------------------#
#   6.2.7 TMBによる尤度に基づく解析
# --------------------------------------------------------------------#
#######################################################################

# データをまとめて確認（前と同じ）
str(dataList <- list(y = y, x = x, n = n))   # 出力は省略

# TMBモデルファイルを書き出す
cat(file = "model6_2_7.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //入力データを定義
  DATA_VECTOR(y);                    // 応答変数
  DATA_VECTOR(x);                    // 共変量
  DATA_INTEGER(n);                   // 観測数

  //パラメータを定義
  PARAMETER(alpha);                  // 切片
  PARAMETER(beta);                   // 傾き
  PARAMETER(log_sigma);              // 残差標準偏差の対数

  Type sigma = exp(log_sigma);       // Typeは関数の出力型に合わせる
  Type LL = 0.0;                     // 対数尤度の合計を0で初期化
  Type mu;

  for (int i= 0; i<n; i++){          // 添字は1ではなく0から始まることに注意
    mu = alpha + beta * x(i);
    // 各観測の対数尤度を計算し、合計に加える
    LL += dnorm(y(i), mu, sigma, true); // 観測iの対数尤度を加算
  }

  // 派生パラメータ
  Type mu1 = alpha;
  Type mu2 = alpha + beta;
  ADREPORT(mu1);                     // mu1とmu2を出力に含める
  ADREPORT(mu2);
  return -LL;                        // 対数尤度の合計に負号を付けて返す
}
")

# TMB関数をコンパイルして読み込む
compile("model6_2_7.cpp")
dyn.load(dynlib("model6_2_7"))

# パラメータの初期値とデータ構造を用意
params <- list(alpha = 0, beta = 0, log_sigma = 0)

# TMBオブジェクトを作成
out6.2.7 <- MakeADFun(data = dataList, parameters = params,
DLL = "model6_2_7", silent = TRUE)

# TMBを最適化し、結果を表示して推定値を保存
opt <- optim(out6.2.7$par, fn = out6.2.7$fn, gr = out6.2.7$gr, method = "BFGS", hessian = TRUE)
tmb_est <- c(opt$par[1:2], exp(opt$par[3]))  # 推定値を保存
(tsum <- tmb_summary(out6.2.7))              # 出力を確認（表示省略）


#   6.2.8 パラメータ推定値の比較
# ---------------------------------------------

# 真の値と各手法の推定結果を比較
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est, NIMBLE = nimble_est,
  Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 6.3 不等分散を仮定した2群の比較
# -----------------------------------------------
# コードなし


#   6.3.1 データの生成
# -----------------------

set.seed(63)
# データセットを生成
n1 <- 60                             # 雌の個体数
n2 <- 40                             # 雄の個体数
mu1 <- 105                           # 雌の母平均
mu2 <- 77.5                          # 雄の母平均
sigma1 <- 3                          # 雌の母標準偏差
sigma2 <- 2.5                        # 雄の母標準偏差
n <- n1+ n2                          # 総標本数
y1 <- rnorm(n1, mu1, sigma1)         # 雌のデータ
y2 <- rnorm(n2, mu2, sigma2)         # 雄のデータ
y <- c(y1, y2)                       # 2群のデータを結合
x <- rep(c(0,1), c(n1, n2))          # 雄を示す指示変数

truth <- c(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)

# 図を作成（図6.2）
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ x, col = "grey", xlab = "雄（0 = 雌、1 = 雄）", ylab = "翼開長（cm）", las = 1, frame = FALSE)


#   6.3.2 Rの既存関数による頻度主義的解析
# -------------------------------------------------------

# ウェルチのt検定
(out6.3.2 <- t.test(y ~ x))
tt_est <- c(out6.3.2$est, NA, NA) # 平均の推定値を保存

# 残差とモデルの診断
summary(fm <- lm(y ~ as.factor(x)-1) ) # 6.2.2のモデルを当てはめる
simOut <- simulateResiduals(fm, n = 1000, plot = TRUE)


#   6.3.3 JAGSによるベイズ解析
# -----------------------------------

# データをリストにまとめて構造を確認
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# JAGSのモデルファイルを書き出す
cat(file = "model6.3.3.txt", "
model {
# 事前分布
mu1 ~ dnorm(0, 0.0001)
mu2 ~ dnorm(0, 0.0001)
tau1 <- pow(sigma1, -2)
sigma1 ~ dunif(0, 1000)
tau2 <- pow(sigma2, -2)
sigma2 ~ dunif(0, 1000)

# 尤度
for (i in 1:n1) {                    # 第1標本（雌）
  y1[i] ~ dnorm(mu1, tau1)
}

for (i in 1:n2) {                    # 第2標本（雄）
  y2[i] ~ dnorm(mu2, tau2)
}

# 派生量
delta <- mu2 - mu1
}
")

# 初期値を生成する関数
inits <- function(){ list(mu1 = rnorm(1), mu2 = rnorm(1), sigma1 = rlnorm(1), sigma2 = rlnorm(1))}

# 推定するパラメータ
params <- c("mu1", "mu2", "delta", "sigma1", "sigma2")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行（参考実行時間：1分未満）し、収束確認・事後分布の要約・推定値の保存を行う
out6.3.3 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model6.3.3.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out6.3.3) # 出力は省略
jags_est <- out6.3.3$summary[c(1, 2, 4, 5), 1] # 推定値を保存
print(out6.3.3, 3)

samps <- as.matrix(out6.3.3$samples)
quantile(samps[,'sigma2'] - samps[,'sigma1'], c(0.025,0.975))


#   6.3.4 NIMBLEによるベイズ解析
# -------------------------------------

# データをリストにまとめて構造を確認
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# NIMBLEのモデルを定義
model63.4 <- nimbleCode( {
# 事前分布
mu1 ~ dnorm(0, 0.0001)
mu2 ~ dnorm(0, 0.0001)
tau1 <- pow(sigma1, -2)
sigma1 ~ dunif(0, 1000)
tau2 <- pow(sigma2, -2)
sigma2 ~ dunif(0, 1000)

# 尤度
for (i in 1:n1) {                    # 第1標本（雌）
  y1[i] ~ dnorm(mu1, tau1) 
}

for (i in 1:n2) {                    # 第2標本（雄）
  y2[i] ~ dnorm(mu2, tau2) 
}

# 派生量
delta <- mu2 - mu1
} )


# 初期値
inits <- function(){ list(mu1 = rnorm(1), mu2 = rnorm(1), sigma1 = rlnorm(1), sigma2 = rlnorm(1))}

# 監視するパラメータ：前と同じ
params <- c("mu1", "mu2", "delta", "sigma1", "sigma2")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行（参考実行時間：30秒）し、収束確認・事後分布の要約・推定値の保存を行う
system.time(  
  out6.3.4 <- nimbleMCMC(code = model63.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out6.3.4)     # 出力は省略
(nsum <- nimble_summary(out6.3.4, params))       # 出力は省略
nimble_est <- nsum[c(1,2,4,5),1] # 推定値を保存


#   6.3.5 Stanによるベイズ解析
# -----------------------------------

# データをリストにまとめて構造を確認
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# モデル記述ファイルを書き出す（原文ではBUGS言語と表記）
cat(file = "model6_3_5.stan",        # この行はRコード
"data {                              //ここからStanコード
  int n1;                            //全データの型・形式を指定
  int n2;
  vector[n1] y1;
  vector[n2] y2;
}

parameters {                         //全パラメータの型・形式を指定
  real mu1;
  real mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

model {                              //ここに事前分布と尤度を記述
  //事前分布
  mu1 ~ normal(0, 100);
  mu2 ~ normal(0, 100);
  sigma1 ~ cauchy(0, 10);
  sigma2 ~ cauchy(0, 10);
  //尤度
  y1 ~ normal(mu1, sigma1);
  y2 ~ normal(mu2, sigma2);
}

generated quantities{
  real delta = mu2 - mu1;
}                                    //Stanコードはここまで
" )

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行（原文の参考実行時間：55/1秒）
system.time(
  out6.3.5 <- stan(file = "model6_3_5.stan", data = dataList,
    chains = nc, iter = ni, warmup = nb, thin = nt) )

rstan::traceplot(out6.3.5) # 出力は省略
stan_est <- summary(out6.3.5)$summary[1:4, 1] # 推定値を保存
print(out6.3.5, dig = 2)             # 出力は省略


#   6.3.6 自作関数による最尤推定
# ---------------------------

NLL <- function(param, y1, y2){
  mu1 <- param[1]
  mu2 <- param[2]
  sigma1 <- exp(param[3])
  sigma2 <- exp(param[4])
  LL1 <- dnorm(y1, mu1, sigma1, log = TRUE) # 雌
  LL2 <- dnorm(y2, mu2, sigma2, log = TRUE) # 雄
  LL <- sum(LL1) + sum(LL2)                 # 合計
  -LL                                       # 負の対数尤度
}

# 負の対数尤度を最小化して最尤推定値・標準誤差・信頼区間を求め、推定値を保存
inits <- c('mu1' = 50, 'mu2' = 50, 'log_sigma1' = 1, 'log_sigma2' = 1)
out6.3.6 <- optim(inits, NLL, y1 = y1, y2 = y2, hessian = TRUE)
diy_est <- c(out6.3.6$par[1:2], exp(out6.3.6$par[3:4]))
get_MLE(out6.3.6, 5)


#   6.3.7 TMBによる尤度に基づく解析
# ------------------------------------

# データをまとめて確認（前と同じ）
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# TMBモデルファイルを書き出す
cat(file = "model6_3_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //入力データを定義
  DATA_VECTOR(y1);                   //第1群の応答変数
  DATA_VECTOR(y2);                   //第2群の応答変数
  DATA_INTEGER(n1);                  //第1群の観測数
  DATA_INTEGER(n2);                  //第2群の観測数
  
  //パラメータを定義
  PARAMETER(mu1);                    //第1群の切片
  PARAMETER(mu2);                    //第2群の切片
  PARAMETER(log_sigma1);             //第1群の残差標準偏差の対数
  PARAMETER(log_sigma2);             //第2群の残差標準偏差の対数

  Type sigma1 = exp(log_sigma1);
  Type sigma2 = exp(log_sigma2);
  Type LL = 0.0;                     //対数尤度の合計を0で初期化
  
  for (int i= 0; i<n1; i++){         //添字は1ではなく0から始まることに注意
    LL += dnorm(y1(i), mu1, sigma1, true); //観測iの対数尤度を加算
  }

  for (int i= 0; i<n2; i++){         //添字は1ではなく0から始まることに注意
    LL += dnorm(y2(i), mu2, sigma2, true); //観測iの対数尤度を加算
  }
  
  Type delta = mu2 - mu1;
  ADREPORT(delta);

  return -LL;                        //対数尤度の合計に負号を付けて返す
}
")

# TMB関数をコンパイルして読み込む
compile("model6_3_7.cpp")
dyn.load(dynlib("model6_3_7"))

# パラメータの初期値とデータ構造を用意
params <- list(mu1 = 0, mu2 = 0, log_sigma1 = 0, log_sigma2 = 0)

# TMBオブジェクトを作成
out6.3.7 <- MakeADFun(data = dataList, parameters = params,
  DLL = "model6_3_7", silent = TRUE)

# TMBを最適化し、結果を表示して推定値を保存
opt <- optim(out6.3.7$par, fn = out6.3.7$fn, gr = out6.3.7$gr,
method = "BFGS", hessian = TRUE)
tmb_est <- c(opt$par[1:2], exp(opt$par[3:4])) # 推定値を保存
(tsum <- tmb_summary(out6.3.7))


#   6.3.8 パラメータ推定値の比較
# ---------------------------------------------

# 真の値と各手法の推定結果を比較
comp <- cbind(truth = truth, tt = tt_est, JAGS = jags_est, NIMBLE = nimble_est,
 Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)


# 6.4 まとめと分散のモデル化に関する考察
# ------------------------------------------------------
# コードなし
