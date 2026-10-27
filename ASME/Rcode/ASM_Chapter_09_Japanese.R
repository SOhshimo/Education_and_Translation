
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ----------------------------------------------------------------
# 第9章：General linear model for a normal response with
# continuous andcategorical explanatory variables
# ----------------------------------------------------------------

# 最終更新：11 June 2024

#######################################################################
# --------------------------------------------------------------------#
# 9.2 データの生成
# --------------------------------------------------------------------#
#######################################################################

set.seed(9)
nPops <- 3
nSample <- 10
n <- nPops * nSample # データの総数
x <- rep(1:nPops, rep(nSample, nPops)) # 集団の識別子
pop <- factor(x, labels = c("Pyrenees", "Massif Central", "Jura"))
length <- runif(n, 45, 70) # 成体の体長は45 cm未満になることが少ない
lengthC <- length-mean(length) # 中心化した体長を使う

Xmat <- model.matrix(~ pop * lengthC)
print(Xmat, dig = 2) # 表示は省略するが、構造を確認して理解すること
beta.vec <- c(80, -30, -20, 6, -3, -4)

sigma <- 10 # Choose 残差 SD
lin.pred <- Xmat[,] %*% beta.vec # 線形予測子の値
eps <- rnorm(n = n, mean = 0, sd = sigma) # 残差
mass <- as.numeric(lin.pred + eps) # response = lin.pred + residual
hist(mass) # Inspect what we’ve created (表示省略)

par(mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
matplot(cbind(length[1:10], length[11:20], length[21:30]), cbind(mass[1:10],
mass[11:20], mass[21:30]), ylim = c(0, max(mass)), ylab = "Body mass (g)",
xlab = "Body length (cm)", col = c("Red","Green","Blue"), pch = c("P", "M", "J"),
las = 1, cex = 1.6, cex.lab = 1.5, frame = FALSE) # 図 9–2

# 後の比較用に設定値を保存する
truth <- c(beta.vec, sigma)

# 必要なライブラリを読み込む
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
# 9.3 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################
summary(out9.3 <- lm(mass ~ pop * lengthC))

# 最小二乗推定値を保存する
lm_est <- c(coef(out9.3), sigma = summary(out9.3)$sigma)

# 分位点残差によって適合度を確認する
plot(out9.3) # 比較のため従来型の残差診断を行う
simOut <- simulateResiduals(out9.3, n = 1000, plot = TRUE)

par(mfrow = c(1,2))
# 図示 効果 of length, holding pop constant (left 図示)
# make sequence of 100 length 値
length.seq <- seq(min(length), max(length), length.out = 100)
# Center the length sequence を用いる the original 平均
lengthC.seq <- length.seq - mean(length)
# make newdata データ.frame, predict and 計算する 95% CI around prediction
newdata <- data.frame(lengthC = lengthC.seq, pop = "Pyrenees")
pr.length <- predict(out9.3, newdata = newdata, se.fit = TRUE)
LCL <- pr.length$fit - 1.96 * pr.length$se.fit # Wald法による95%信頼区間の下限
UCL <- pr.length$fit + 1.96 * pr.length$se.fit # Wald法による95%信頼区間の上限

# 描く the 図示 (図 9-3 left)
plot(length[pop == "Pyrenees"], mass[pop == "Pyrenees"], pch = 16,
ylab = "Mass", xlab = "Length", main = "Effect of length on mass for Pyrenees",
frame = FALSE, ylim = c(0, 180), cex = 1.5, col = rgb(0,0,0,0.5))
lines(length.seq, pr.length$fit, col = 'red', lwd = 3)
polygon(c(length.seq, rev(length.seq)), c(LCL, rev(UCL)), col = rgb(1,0,0,0.2), border = NA)

# 図示 効果 of pop, holding length constant at its median (right 図示)
newdata <- data.frame(lengthC = median(lengthC), pop = levels(pop))
pr.pop <- predict(out9.3, newdata = newdata, se.fit = TRUE)
LCL <- pr.pop$fit - 1.96 * pr.pop$se.fit
UCL <- pr.pop$fit + 1.96 * pr.pop$se.fit
plot(as.numeric(pop), mass, pch = 16, xaxt = 'n', xlab = "Pop",
ylab = "Mass", main = "Effect of pop on mass at median length",
frame = FALSE, cex = 1.5, col = rgb(0,0,0,0.4)) # 描く the 図示 (図 9.3 right)
axis(1, at = 1:3, labels = levels(pop))
points(1:3, pr.pop$fit, pch = 16, col = "red", cex = 1.5)
segments(1:3, LCL, 1:3, UCL, lwd = 2, col = "red")
par(mfrow = c(1,1))

#######################################################################
# --------------------------------------------------------------------#
# 9.4 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC,
    nPops = nPops, n = n) )

# JAGSモデルファイルを書き出す
cat(file = "model9.4.txt", "
model {
# 事前分布
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, 0.0001) # 切片
  beta[i] ~ dnorm(0, 0.0001) # 傾き
}
sigma ~ dunif(0, 100) # 残差の標準偏差
tau <- pow(sigma, -2)

# 尤度
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthC[i]
}

# 定義する 効果 relative to baseline level
a.effe2 <- alpha[2] - alpha[1] # 切片の差：中央高地とピレネー
a.effe3 <- alpha[3] - alpha[1] # 切片の差：ジュラとピレネー
b.effe2 <- beta[2] - beta[1] # 傾きの差：中央高地とピレネー
b.effe3 <- beta[3] - beta[1] # 傾きの差：ジュラとピレネー
# 任意の対比
test1 <- beta[3] - beta[2] # 傾きの差：ジュラと中央高地
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rnorm(nPops, 0, 2), beta = rnorm(nPops, 1, 1), sigma = runif(1))}

# 推定対象のパラメータ
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） 2 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
out9.4 <- jags(dataList, inits, params, "model9.4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out9.4) # 図・出力は省略
print(out9.4, 3)
jags_est <- out9.4$summary[c(1, 8, 9, 4, 10, 11, 7), 1]

# 尤度法とベイズ法の推定値および真値を比較する
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)

#######################################################################
# --------------------------------------------------------------------#
# 9.5 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


library(nimble)

# データをリストにまとめ構造を確認する
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC, nPops = nPops, n = n) )

# NIMBLEモデルコードを定義する
model9.5 <- nimbleCode( {
# 事前分布
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, 0.0001)     # 切片
  beta[i] ~ dnorm(0, 0.0001)      # 傾き
}
sigma ~ dunif(0, 100)             # 残差の標準偏差
tau <- pow(sigma, -2)

# 尤度
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthC[i]
}

# 定義する 効果 relative to baseline level
a.effe2 <- alpha[2] - alpha[1]    # 切片の差：中央高地とピレネー
a.effe3 <- alpha[3] - alpha[1]    # 切片の差：ジュラとピレネー
b.effe2 <- beta[2] - beta[1]      # 傾きの差：中央高地とピレネー
b.effe3 <- beta[3] - beta[1]      # 傾きの差：ジュラとピレネー

# 任意の対比
test1 <- beta[3] - beta[2]        # 傾きの差：ジュラと中央高地
} )

# 初期値
inits <- function(){ list(alpha = rnorm(nPops, 0, 2), beta = rnorm(nPops, 1, 1), sigma = runif(1))}

# 監視するパラメータは前と同じ
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 25 sec), 収束を診断する, 事後分布を要約する and 推定値を保存する
system.time( out9.5 <- 
    nimbleMCMC(code = model9.5,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out9.5)
(nsum <- nimble_summary(out9.5, params))
nimble_est <- nsum[c(1, 8, 9, 4, 10, 11, 7), 1]

#######################################################################
# --------------------------------------------------------------------#
# 9.6 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC,
  nPops = nPops, n = n) )

# Stanモデルファイルを書き出す
cat(file = "model9_6.stan", "
data {
  int n; // 標本数
  int nPops; // 集団数
  vector[n] mass; // 応答変数
  vector[n] lengthC; // 共変量
  array[n] int pop; // 集団 of 各 観測
}

parameters {
  vector[nPops] alpha; // intercepts
  vector[nPops] beta; // slopes
  real <lower = 0> sigma; // 残差 標準偏差
}

model {
  vector[n] mu; // expected 値 of 観測値
  
  // 事前分布
  for (i in 1:nPops){
    alpha[i] ~ normal(0, 100);
    beta[i] ~ normal(0, 100);
  }
  sigma ~ uniform(0, 100);

  // 尤度
  for (i in 1:n){
    mu[i] = alpha[pop[i]] + beta[pop[i]] * lengthC[i];
    mass[i] ~ normal(mu[i], sigma);
  }
}

generated quantities {
  real a_effe2;
  real a_effe3;
  real b_effe2;
  real b_effe3;
  real test1;
  a_effe2 = alpha[2] - alpha[1]; // 切片の差：中央高地とピレネー
  a_effe3 = alpha[3] - alpha[1]; // 切片の差：ジュラとピレネー
  b_effe2 = beta[2] - beta[1]; // 傾きの差：中央高地とピレネー
  b_effe3 = beta[3] - beta[1]; // 傾きの差：ジュラとピレネー
  test1 = beta[3] - beta[2]; // 傾きの差：ジュラと中央高地
}
")

# 推定対象のパラメータ
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# HMCの設定
ni <- 1000 ; nb <- 500 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 40/2 sec), print output and 推定値を保存する
system.time(
out9.6 <- stan(file = "model9_6.stan", data = dataList,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out9.6) # 図・出力は省略
print(out9.6, dig = 3) # 図・出力は省略
stan_est <- summary(out9.6)$summary[c(1, 8, 9, 4, 10, 11, 7), 1]

#######################################################################
# --------------------------------------------------------------------#
# 9.7 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################


NLL <- function(param, y, Xmat) {
  beta <- param[1:6] # 切片/slopes matching columns of Xmat
  sigma <- exp(param[7]) # 残差の標準偏差
  mu <- Xmat %*% beta # Matrix-multiply beta and design matrix
  LL <- dnorm(y, mu, sigma, log = TRUE) # 観測値ごとの対数尤度
  -sum(LL) # データ全体の負の対数尤度
}

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と95%信頼区間を求める and 推定値を保存する
inits <- c(mean(mass), rep(0, 5), log(sd(mass)))
names(inits) <- c(colnames(Xmat), 'log-sigma')
out9.7 <- optim(inits, NLL, y = mass, Xmat = Xmat, hessian = TRUE, method = "BFGS")
get_MLE(out9.7, 4)
diy_est <- c(out9.7$par[1:6], exp(out9.7$par[7]))

#######################################################################
# --------------------------------------------------------------------#
# 9.8 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################

# データをリストにまとめ構造を確認する
tmbData <- dataList
tmbData$pop <- tmbData$pop-1
str(tmbData)

# TMBモデルファイルを書き出す
cat(file = "model9_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  // 入力データを定義する
  DATA_INTEGER(n); // 観測数
  DATA_VECTOR(mass); // 応答変数
  DATA_VECTOR(lengthC); // 共変量
  DATA_IVECTOR(pop); // 集団 index

  // パラメータを定義する
  PARAMETER_VECTOR(alpha); // 切片
  PARAMETER_VECTOR(beta); // 傾き
  PARAMETER(log_sigma); // 残差 sd on log scale

  Type sigma = exp(log_sigma); // 残差の標準偏差
  
  vector <Type> mu(n); // Expected 値 of 観測値
  
  Type LL = 0; // 対数尤度を0で初期化する

  // 観測値ごとに繰り返す
  for (int i = 0; i < n; i++){ // 添字は1ではなく0から始まる点に注意
    mu(i) = alpha(pop(i)) + beta(pop(i)) * lengthC(i);
    // 計算する 対数尤度 of 観測 and add to 合計 loglik
    LL += dnorm(mass(i), mu(i), sigma, true);
  }

  // 導出パラメータを計算する
  Type a_effe2 = alpha(1) - alpha(0); // 切片の差：中央高地とピレネー
  Type a_effe3 = alpha(2) - alpha(0); // 切片の差：ジュラとピレネー
  Type b_effe2 = beta(1) - beta(0); // 傾きの差：中央高地とピレネー
  Type b_effe3 = beta(2) - beta(0); // 傾きの差：ジュラとピレネー
  Type test1 = beta(2) - beta(1); // 傾きの差：ジュラと中央高地

  // 導出パラメータをTMBから報告させる
  ADREPORT(a_effe2);
  ADREPORT(a_effe3);
  ADREPORT(b_effe2);
  ADREPORT(b_effe3);
  ADREPORT(test1);
  
  return -LL; // 合計対数尤度の負値を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model9_8.cpp")
dyn.load(dynlib("model9_8"))

# パラメータの次元と初期値を指定する
# Starting 値 need to be semi-informative (use same ones as DIY)
params <- list(alpha = c(mean(mass), 0, 0), beta = rep(0,3),
log_sigma = log(sd(mass)))

# TMBの目的関数オブジェクトを作成する
out9.8 <- MakeADFun(data = tmbData, parameters = params,
  DLL = "model9_8", silent = TRUE)

# TMBオブジェクトを最適化する and print and 結果を保存する
opt <- optim(out9.8$par, fn = out9.8$fn, gr = out9.8$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out9.8))
tmb_est <- c(tsum[c(1, 8, 9, 4, 10, 11), 1], exp(tsum[7, 1]))

#######################################################################
# --------------------------------------------------------------------#
# 9.9 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################


# 推定結果を真値と比較する and previous 推定値
comp <- cbind(cbind(truth = truth, lm = lm_est, JAGS = jags_est,
NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est))
print(comp, 3)

