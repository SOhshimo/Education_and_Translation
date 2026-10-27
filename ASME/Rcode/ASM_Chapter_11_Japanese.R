
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ------------------------------------------------------
# 第11章：Introduction to the generalized linear
# model (GLM): comparing two groups in a
# ポアソン回帰
# ------------------------------------------------------

# 最終更新：11 June 2024


#######################################################################
# --------------------------------------------------------------------#
# 11.4 データの生成
# --------------------------------------------------------------------#
#######################################################################

set.seed(11)
nSites <- 30
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))
n <- 2 * nSites

lambda <- exp(0.69 + 0.92*(as.numeric(x)-1))           # x has levels 1 and 2, not 0 and 1

# 保存する true 値 of パラメータ 対象： later 比較
truth <- c(intercept = 0.69, arable = 0.92)

y <- rpois(n = n, lambda = lambda)                     # ポアソン乱数による変動を加える
y[c(1, 10, 35)] <- NA                                  # 作成する some 観測値 NA
y                                                      # 出力する the counts (表示省略)

aggregate(y, by = list(x), FUN = mean, na.rm = TRUE)   # Observed means
boxplot(y ~ x, col = "grey", xlab = "Land-use",
  ylab = "Hare count", las = 1, frame = FALSE)         # 図 11–2

# Required packages
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
# 11.5 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################


out11.5 <- glm(y ~ x, family = poisson)                # 当てはめる the モデル
summary(out11.5)                                       # Two-group 比較
anova(out11.5, test = "Chisq")                         # 尤度 ratio test (LRT)

# 推定値を保存する
glm_est <- coef(out11.5)

# 適合度を確認する 通常の残差と分位点残差を用いる
plot(out11.5)                                          # 比較のため従来型の残差診断を行う
simOut <- simulateResiduals(out11.5, n = 1000, plot = TRUE)

#######################################################################
# --------------------------------------------------------------------#
# 11.6 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# JAGSモデルファイルを書き出す
cat(file = "model11.6.txt", "
model {
# 事前分布
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# 尤度
for (i in 1:n) {
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha + beta *x[i]
  # lambda[i]) <- exp(alpha + beta *x[i])              # Same
}

# 適合度の評価
for (i in 1:n) {
  Presi[i] <- (y[i]-lambda[i]) / sqrt(lambda[i])       # ピアソン残差
  y.new[i] ~ dpois(lambda[i])                          # 複製データを生成する
  Presi.new[i] <- (y.new[i]-lambda[i]) / sqrt(lambda[i]) # Pearson resi
  D[i] <- pow(Presi[i], 2)
  D.new[i] <- pow(Presi.new[i], 2)
}
# 不一致度の指標を合計する
fit <- sum(D[])
fit.new <- sum(D.new[])
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1))}

# 推定対象のパラメータ
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit.new", "y")

# MCMCの設定
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する and 事後分布を要約する
out11.6 <- jags(dataList, inits, params, "model11.6.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out11.6) # 表示省略 (yet)
print(out11.6, 3) # 図・出力は省略


# 11.6.1 Assessment of モデル adequacy
# -------------------------------------

# 図 11-3
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(out11.6$mean$Presi, ylab = "Pearson residual", las = 1, pch = 16,
  cex = 1.5, col = rgb(0,0,0,0.5), frame = FALSE)
abline(h = 0)
plot(out11.6$sims.list$fit, out11.6$sims.list$fit.new, main = "",
  xlab = "Discrepancy actual data set", ylab = "Discrepancy perfect data sets",
  pch = 16, cex = 1.5, col = rgb(0,0,0,0.2), frame = FALSE)
abline(0, 1, lwd = 2, col = "black")

# Bayesian p-値
mean(out11.6$sims.list$fit.new > out11.6$sims.list$fit)

# 図 11.4
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(out11.6$sims.list$beta, col = "grey", las = 1,
  xlab = "Coefficient", main = "Effect of 'arable'",
  breaks = 50, freq = FALSE, xlim= c(0,2.2))
abline(v=0, col='red', lty=2, lwd = 3)
hist(exp(out11.6$sims.list$alpha), main= "Comparison of land covers",
  col = "red", xlab = "Expected hare count", xlim = c(0,7),
  breaks = 20, freq = FALSE)
hist(exp(out11.6$sims.list$alpha + out11.6$sims.list$beta),
  col = "blue", breaks = 20, freq = FALSE, add = TRUE)
  legend('topright', legend=c('Grassland', 'Arable'), 
  fill = c("red", "blue"), cex = 1.5)

#######################################################################
# --------------------------------------------------------------------#
# 11.6.2 Inference under the モデル
# --------------------------------------------------------------------#
#######################################################################


print(out11.6, 2)

# 尤度法とベイズ法の推定値および真値を比較する
jags_est <- unlist(out11.6$mean[1:2])
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est)
print(comp, 4)


# 11.7 NIMBLEによるベイズ解析
# ----------------------------------

library(nimble)

# Bundle and summarize データ (same as 対象： JAGS)
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# NIMBLEモデルコードを定義する
model11.7 <- nimbleCode( {
# 事前分布
alpha ~ dnorm(0, sd = 100)
beta ~ dnorm(0, sd = 100)

# 尤度
for (i in 1:n) {
  y[i] ~ dpois(lambda[i]) 
  log(lambda[i]) <- alpha + beta * x[i]
}

# 適合度の評価
for (i in 1:n) {
  Presi[i] <- (y[i]-lambda[i]) / sqrt(lambda[i])   # ピアソン残差
  y.new[i] ~ dpois(lambda[i])                      # 複製データを生成する
  Presi.new[i] <- (y.new[i]-lambda[i]) / sqrt(lambda[i]) # Pearson resi
  D[i] <- pow(Presi[i], 2)
  D.new[i] <- pow(Presi.new[i], 2)
}

# 不一致度の指標を合計する
 fit <- sum(D[1:20])
 fit.new <- sum(D.new[1:20])
} )

# 初期値
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1))}

# 監視するパラメータは前と同じ
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit.new")

# MCMCの設定
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行する (原著記載の実行時間（ART） 28 sec), 収束を診断する and 事後分布を要約する
system.time( out11.7 <- 
    nimbleMCMC(code = model11.7,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out11.7)
(nsum <- nimble_summary(out11.7, params))

# からの推定値を保存： NIMBLE
nimble_est <- nsum[1:2,1]

#######################################################################
# --------------------------------------------------------------------#
# 11.8 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
not_na <- which(!is.na(y))
x_stan <- as.numeric(x) - 1
x_stan <- x_stan[not_na]
str(dataList <- list(y = y[not_na], x = x_stan, n = length(x_stan)))

# Stanモデルファイルを書き出す
cat(file = "model11_8.stan", "
data {
  int n;
  array[n] int y;
  vector[n] x;
}

parameters {
  real alpha;
  real beta;
}

transformed parameters {
  vector[n] lambda;
  for (i in 1:n) {
    lambda[i] = exp(alpha + beta * x[i]);
  }
}

model {

  // 事前分布
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  
  // 尤度
  for (i in 1:n) {
    y[i] ~ poisson(lambda[i]);
  }
}

generated quantities {
  vector[n] Presi;
  vector[n] D;
  array[n] int y_new;
  vector[n] Presi_new;
  vector[n] D_new;
  real fit;
  real fit_new;
  for (i in 1:n) {
    Presi[i] = (y[i] - lambda[i])/sqrt(lambda[i]);
    D[i] = pow(Presi[i], 2);
    y_new[i] = poisson_rng(lambda[i]);
    Presi_new[i] = (y_new[i] - lambda[i])/sqrt(lambda[i]);
    D_new[i] = pow(Presi_new[i], 2);
  }
  fit = sum(D);
  fit_new = sum(D_new);
}
")

# 監視するパラメータは前と同じ
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit_new")

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 32/3 sec), 収束を評価する, print 推定値 and save point 推定値
system.time(
  out11.8 <- stan(file = "model11_8.stan", data = dataList, pars = params,
                  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out11.8) # 図・出力は省略
print(out11.8, dig = 3) # 図・出力は省略
stan_est <- summary(out11.8)$summary[1:2,1]

#######################################################################
# --------------------------------------------------------------------#
# 11.9 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################



# 負の対数尤度の定義 （対象： ポアソン GLM with 1 covariate
NLL <- function(param, y, x) {
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x) # 計算する lambda 対象： 各 datum
  LL <- dpois(y, lambda, log = TRUE) # 対数尤度 対象： 各 datum
  -sum(LL, na.rm = TRUE) # 全観測値の負の対数尤度
}

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
inits <- c('Intercept' = 0, 'beta(arable)' = 0)
out11.9 <- optim(inits, NLL, y = y, x = as.numeric(x)-1, hessian = TRUE)
get_MLE(out11.9, 5)
diy_est <- out11.9$par

#######################################################################
# --------------------------------------------------------------------#
# 11.10 TMBを用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################



# データをリストにまとめ構造を確認する
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# TMBモデルファイルを書き出す
cat(file = "model11_10.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(y); // 応答変数
  DATA_VECTOR(x); // 共変量
  DATA_INTEGER(n); // 観測値数

  // パラメータを定義する
  PARAMETER(alpha);
  PARAMETER(beta);
  
  Type LL = 0.0; // 対数尤度を0で初期化する
  Type lambda; // Initialize lambda

  for (int i= 0; i<n; i++){
    if(R_IsNA(asDouble(y(i)))) continue; // Here accommodate NAs
    lambda = exp(alpha + beta * x(i)); // Expected count

    // 各観測値の対数尤度を計算し合計に加える
    LL += dpois(y(i), lambda, true);
  }

  return -LL; // 負の対数尤度を返す
}
")


# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model11_10.cpp")
dyn.load(dynlib("model11_10"))

# パラメータの次元と初期値を指定する
params <- list(alpha = 0, beta = 0)

# TMBの目的関数オブジェクトを作成する
out11.10 <- MakeADFun(data = dataList,
                      parameters = params,
                      DLL = "model11_10", silent = TRUE)

# TMBオブジェクトを最適化する and 結果を出力する
opt <- optim(out11.10$par, fn = out11.10$fn, gr = out11.10$gr, method = "BFGS")
(tsum <- tmb_summary(out11.10)) # 図・出力は省略
tmb_est <- c(opt$par[1:2]) # 推定値を保存する

#######################################################################
# --------------------------------------------------------------------#
# 11.11 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################


# 推定値と真値を比較する
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)
