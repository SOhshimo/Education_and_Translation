
# ---------------------------------------------------------------
# 『Applied Statistical Modeling for Ecologists』のコード 
# 著者：Kéry & Kellner、Elsevier、2024年
# ---------------------------------------------------------------

# ------------------------------------------
# 第5章 ― 正規線形回帰
# ------------------------------------------

# 最終更新：2024年6月11日

#######################################################################
# --------------------------------------------------------------------#
#   5.2 データ生成
# --------------------------------------------------------------------#
#######################################################################

set.seed(5)
n <- 16 # 年数
a <- 40 # 切片
b <- -0.5 # 傾き
sigma2 <- 25 # 残差分散
# 後で比較するため真の値を保存
truth <- c(alpha = a, beta = b, sigma = sqrt(sigma2))

# 類似した図は図5-3（左）を参照
x <- 1:16 # 共変量である年の値
eps <- rnorm(n, mean = 0, sd = sqrt(sigma2)) # 残差
y <- a + b*x + eps # データセットを作成
plot((x + 1989), y, xlab = "年", las = 1,
ylab = "占有率 (%)", cex = 1.5,
pch = 16, col = rgb(0,0,0,0.6), frame = FALSE) # 図示・出力は本文に掲載なし

# 必要なパッケージを読み込む
library(ASMbook); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
#   5.3 Rの既存関数による解析
# --------------------------------------------------------------------#
#######################################################################

#   5.3.1 モデルの当てはめ
# -------------------------

# モデルを当てはめ、lmの推定値を保存
summary(out5.3 <- lm(y ~ x))
lm_est <- c(coef(out5.3), sigma = sigma(out5.3))


#   5.3.2 従来の残差と分位点残差を使った適合度評価
# ---------------------------------------------------------------------------

# 従来の残差（＝生残差）を表示
print(residuals(out5.3), 4)

# 当てはめた線形回帰モデルの妥当性を視覚的に評価
hist(residuals(out5.3), breaks = 20) # 図示・出力は本文に掲載なし
plot(out5.3) # 図示・出力は本文に掲載なし

# 分位点残差を計算し、適合度を評価
library(DHARMa)
simOut <- simulateResiduals(out5.3, n = 1000)
residuals(simOut) # スケーリングされた分位点残差を表示

# 正規分布に変換した分位点残差
print(residuals(simOut, quantileFunction = qnorm), 4)

par(mfrow = c(1,2))
plotQQunif(simOut)
plotResiduals(simOut, rank = FALSE) # 上のコメントを参照
par(mfrow = c(1,1))
testDispersion(simOut) # 過分散または過小分散を検定
# （本文に掲載なし。なお、この検定が本来
# 意味を持つのは非正規GLMの場合である）


#   5.3.3 予測値の計算
# -----------------------------

xpred <- seq(1, 16, length.out = 100) # xの100個の値について予測
newdata <- list(x = xpred)
pred <- predict(out5.3, newdata = newdata, se.fit = TRUE) # 標準誤差も取得
LCL <- pred$fit - 1.96 * pred$se.fit # 95%信頼区間の下限
UCL <- pred$fit + 1.96 * pred$se.fit # 95%信頼区間の上限
pred.table <- cbind(x = newdata$x, "point_est" =
pred$fit, "se" = pred$se.fit, "LCL" = LCL, "UCL" = UCL) # xと予測値をまとめた表
print(pred.table) # 一部だけを表示。
# 図5.4も参照

#######################################################################
# --------------------------------------------------------------------#
#   5.4 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################
#   5.4.1 モデルの当てはめ
# -------------------------

# データをリストにまとめ、構造を確認
str(dataList <- list(y = y, x = x, n = n))

# JAGSモデルファイルを書き出す
cat(file = "model5.4.txt", "
model {
# 事前分布
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
sigma ~ dunif(0, 100) # より狭い範囲も検討できる
tau <- pow(sigma, -2)
# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau) # 応答変数yは正規分布に従う
  mu[i] <- alpha + beta * x[i] # 線形予測子
}
# 平方和に基づく乖離量を用いてモデルの適合度を評価
for (i in 1:n) {
  residual[i] <- y[i]-mu[i] # 観測データの生残差
  predicted[i] <- mu[i] # 予測値（当てはめ値）
  sq[i] <- residual[i]^2 # 観測データの残差平方
  # 事後予測チェック
  # データの事後予測分布に基づく
  # 複製データを生成し、その適合度統計量を計算
  y.new[i] ~ dnorm(mu[i], tau) # MCMCの各反復で新たなデータセットを1つ生成：
  # これがyの事後予測分布となる
  sq.new[i] <- (y.new[i]-predicted[i])^2 # 新たなデータの残差平方
}
fit <- sum(sq[]) # 実際のデータセットの残差平方和
fit.new <- sum(sq.new[]) # 新たなデータセットの残差平方和
}
")

# 初期値を生成する関数
inits <- function(){ list(alpha = rnorm(1),
beta = rnorm(1), sigma = rlnorm(1))}
# sigmaの初期値には対数正規分布を使用
# 0未満の値を避けるため
# 推定するパラメータ
params <- c("alpha","beta", "sigma", "fit", "fit.new", "residual", "predicted", "y.new")
# MCMCの設定
na <- 1000 ; ni <- 6000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# JAGSを実行（実行時間1分未満）、収束を確認し事後分布を要約
out5.4 <- jags(dataList, inits, params, "model5.4.txt", n.iter = ni, n.burnin = nb,
n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out5.4) # 図示・出力は本文に掲載なし
print(out5.4, 2)

# JAGSの推定値を保存
jags_est <- out5.4$summary[1:3, 1]

# lm()とJAGSの結果を比較
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


#   5.4.2 ベイズ解析における適合度評価
# -------------------------------------------------------

# 図5-3（左）
par(mfrow = c(1, 2), mar = c(5,5,2,3))
plot(out5.4$mean$predicted, out5.4$mean$residual, main = "", las = 1,
xlab = "予測値", ylab = "残差", frame = FALSE,
pch = 16, col = rgb(0,0,0,0.5), cex = 1.5)
abline(h = 0)

# 図5-3（右）
lim <- c(0, 4300)
plot(out5.4$sims.list$fit, out5.4$sims.list$fit.new, main = "", las = 1,
xlab = "実測データの残差平方和", ylab = "新しいデータの残差平方和", xlim = lim,
ylim = lim, frame = FALSE, pch = 16, col = rgb(0,0,0,0.3), cex = 1.5)
abline(0, 1)
mean(out5.4$sims.list$fit.new > out5.4$sims.list$fit) # ベイズp値

# JAGSで当てはめたモデルにDHARMaを適用し、分位点残差を計算する方法
y.new <- out5.4$sims.list$y.new
sim <- createDHARMa(simulatedResponse = t(y.new), observedResponse = y,
fittedPredictedResponse = apply(y.new, 2, median))
par(mfrow = c(1,2)) # Graphical results as in Fig. 5.2 (図示・出力は本文に掲載なし)
plotQQunif(sim)
plotResiduals(sim, rank = FALSE)


#   5.4.3 予測値の計算
# -------------------------------------

# ベイズモデルによる予測値を取得
sims <- out5.4$sims.list # JAGSで得た事後分布のサンプルを取得
predictions <- array(dim = c(length(xpred), length(sims$alpha)))
for(i in 1:length(xpred)){
predictions[i,] <- sims$alpha + sims$beta * xpred[i]
}
BPE <- rowMeans(predictions) # ベイズ点推定値（事後平均）
LPB <- apply(predictions, 1, quantile, probs = 0.025) # 下限
UPB <- apply(predictions, 1, quantile, probs = 0.975) # 上限

# 両モデルの予測値と95%不確実性区間を描画
# Fig. 5.4
plot(1990:2005, y, xlab = "年", las = 1, ylab = "占有率 (%)", pch = 16,
ylim = c(25, 50), col = rgb(0,0,0,0.5), frame = FALSE, cex = 1.5)
# lm()による頻度論的推論の信頼区間（5.3.3節）
lines(xpred + 1989, pred.table[,2], lwd = 2, col = 'blue')
polygon(c(newdata$x + 1989, rev(newdata$x + 1989)), c(pred.table[,4],
rev(pred.table[,5])), col = rgb(0, 0, 1, 0.1), border = NA)
# ベイズ推論による信用区間（ここで計算）
lines(xpred + 1989, BPE, lwd = 2, col = 'red')
polygon(c(xpred+ 1989, rev(xpred+ 1989)), c(LPB, rev(UPB)), col = rgb(1, 0, 0, 0.2), border = NA)
legend('bottomleft', legend = c("最尤推定", "事後分布に基づく推論"),
cex = 1.2, bty = 'n', lty = 1, col = c('blue', 'red'), lwd = 3)


#   5.4.4 信頼区間と信用区間の解釈
# ------------------------------------------------------------------------

par(mfrow = c(1, 1), mar = c(5,5,2,3)) # Fig. 5–5
hist(out5.4$sims.list$beta, main = "", col = "grey", xlab = "傾向の推定値",
xlim = c(-2, 1), breaks = 100, freq = FALSE)
abline(v = 0, col = "black", lwd = 2, lty = 3)
mean(out5.4$sims.list$beta < 0) # 個体群減少の確率


#   5.4.5 予測区間についてはどうか？
# --------------------------------------------

# 細かい刻みでyの事後予測分布を取得
sims <- out5.4$sims.list # すべての事後サンプルを取得
str(sims) # MCMC出力の形式を再確認
xpred <- seq(1, 16, length.out = 1000) # xの1000個の値について予測
y_tilde <- array(NA, dim = c(1000, 20000)) # 事後予測分布
for(k in 1:1000){ # xpredの1000個の値すべてについて繰り返す
mu_k <- sims$alpha + sims$beta * xpred[k]
y_tilde[k,] <- rnorm(20000, mean = mu_k, sd = sims$sigma)
}
# y_tildeの95%予測区間を取得
LPI <- apply(y_tilde, 1, quantile, probs = 0.025) # 下限
UPI <- apply(y_tilde, 1, quantile, probs = 0.975) # 上限
# データの95%予測区間を描画（図5.6）
plot(1990:2005, y, xlab = "年", las = 1, ylab = "占有率 (%)",
pch = 16, ylim = c(18, 58), col = rgb(0,0,0,0.5), frame = FALSE,cex = 1.5)
polygon(c(xpred + 1989, rev(xpred + 1989)), c(LPI, rev(UPI)),
col = rgb(1, 0, 0, 0.2), border = NA)

#######################################################################
# --------------------------------------------------------------------#
#   5.5 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################
library(nimble)

# データをリストにまとめ、構造を確認
str(dataList <- list(y = y, x = x, n = n)) # 図示・出力は本文に掲載なし

# NIMBLEモデルを記述
model5.5 <- nimbleCode( {
# 事前分布
alpha ~ dnorm(0, sd = 100)
beta ~ dnorm(0, sd = 100)
sigma ~ dunif(0, 100)

# 尤度
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], sd = sigma) 
  mu[i] <- alpha + beta * x[i]
}

# 平方和に基づく乖離量を用いてモデルの適合度を評価
for (i in 1:16) {
  residual[i] <- y[i]-mu[i]     # 残差 for observed data
  predicted[i] <- mu[i]         # Predicted values
  sq[i] <- residual[i]^2        # 観測データの残差平方

# 複製データを生成し、その適合度統計量を計算
  y.new[i] ~ dnorm(mu[i], sd = sigma) # MCMCの各反復で新しいデータセットを1つ生成
  sq.new[i] <- (y.new[i]-predicted[i])^2 # 新しいデータの残差平方
}
fit <- sum(sq[1:16])     # 実際のデータセットの残差平方和
fit.new <- sum(sq.new[1:16])# 新たなデータセットの残差平方和
} )

# 初期値
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# 記録するパラメータ：前と同じ
params <- c("alpha", "beta", "sigma", "fit", "fit.new",  "residual", "predicted")

# MCMCの設定
#  返されるサンプル数はfloor((niter-nburnin)/thin)となる。
ni <- 6000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# NIMBLEを実行（実行時間約25秒）、収束を確認し事後分布を要約
system.time(  
  out5.5 <- nimbleMCMC(code = model5.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out5.5)     # 図示・出力は本文に掲載なし
(nsum <- nimble_summary(out5.5, params))       # 図示・出力は本文に掲載なし

# ベイズp値
mean(as.matrix(out5.5[,"fit.new"]) > as.matrix(out5.5[,"fit"]))

# NIMBLEの推定値を保存
nimble_est <- nsum[1:3,1]

#######################################################################
# --------------------------------------------------------------------#
#  5.6 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################

# データをリストにまとめ、構造を確認（前と同じ）
str(dataList <- list(y = y, x = x, n = n)) # 図示・出力は本文に掲載なし

# Stan言語で記述したモデルをテキストファイルに保存
cat(file = "model5_6.stan", # この行はRコード
"data { // ここからStanコードが始まる
  int <lower = 0> n; // すべてのデータの形式を定義
  vector[n] y; //ベクトルの次元も指定
  vector[n] x; //
}
parameters { // すべてのパラメータの形式を定義
  real alpha;
  real beta;
  real <lower = 0> sigma; // sigma（標準偏差）は負にならない
}
transformed parameters {
  vector[n] mu;
  for (i in 1:n){
    mu[i] = alpha + beta * x[i]; // 線形予測子を計算
  }
}
model {
  // 事前分布
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  // 尤度（ベクトル化すると高速化できる）
  for (i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
}
generated quantities {
  vector[n] residuals;
  vector[n] sq;
  vector[n] sq_new;
  vector[n] y_new;
  real fit;
  real fit_new;
  for (i in 1:n){
    residuals[i] = y[i] - mu[i];
    sq[i] = residuals[i]^2;
    y_new[i] = normal_rng(mu[i], sigma);
    sq_new[i] = (y_new[i] - mu[i])^2;
  }
  fit = sum(sq);
  fit_new = sum(sq_new);
} // ここでStanコードが終わる
" )

# HMCの設定
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# Stanを実行（実行時間30/3秒）、収束を評価し推定値を保存
system.time(
out5.6 <- stan(file = "model5_6.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt))
rstan::traceplot(out5.6) # 図示・出力は本文に掲載なし
print(out5.6, dig = 2) # 図示・出力は本文に掲載なし
stan_est <- summary(out5.6)$summary[1:3,1]

#######################################################################
# --------------------------------------------------------------------#
#  5.7 最尤推定を自作する
# --------------------------------------------------------------------#
#######################################################################


NLL <- function(params, y, x) {
alpha <- params[1]
beta <- params[2]
sigma <- exp(params[3]) # 対数尺度から元の尺度に変換
n <- length(y) # データ点数
mu <- LL <- numeric(n) # 空のベクトル
# このループはベクトル化すると高速化できる
for (i in 1:n){
mu[i] <- alpha + beta * x[i]
LL[i] <- dnorm(y[i], mean = mu[i], sd = sigma, log = TRUE)
}
-sum(LL)
}

# 負の対数尤度を最小化して最尤推定値と標準誤差を求める
inits <- c('alpha' = 0, 'beta' = 0, 'log.sigma' = 0)
out5.7 <- optim(inits, NLL, y = y, x = x, hessian = TRUE, method = 'BFGS')
get_MLE(out5.7, 4)
# 自作の最尤推定による推定値を保存
diy_est <- c(out5.7$par[1:2], exp(out5.7$par[3]))

#######################################################################
# --------------------------------------------------------------------#
# 5.8 TMBによる尤度解析
# --------------------------------------------------------------------#
#######################################################################

# データをリストにまとめ、構造を確認（前と同じ）
str(tmbData <- list(y = y, x = x, n = n)) # 図示・出力は本文に掲載なし

# TMBモデルファイルを書き出す
cat(file = "model5_8.cpp",
"#include <TMB.hpp> // (1) まず定型的なコードを記述

template <class Type>
Type objective_function <Type>::operator() ()
{
  // (2) 入力データを定義
  DATA_VECTOR(y); // 応答変数
  DATA_VECTOR(x); // 共変量
  DATA_INTEGER(n); // 観測数

  // (3) パラメータを定義
  PARAMETER(alpha); // 切片
  PARAMETER(beta); // 傾き
  PARAMETER(log_sigma); // 残差標準偏差の対数
  Type sigma = exp(log_sigma); // Typeは関数出力の型に合わせる

  // (4) 対数尤度を計算
  Type LL = 0.0; // 対数尤度の合計を0に初期化
  Type mu; //期待値を宣言

  for (int i = 0; i < n; i++){ // 添字は1でなく0から始まる点に注意！
    mu = alpha + beta * x(i);
    LL += dnorm(y(i), mu, sigma, true);
    //LL = LL + dnorm(...)と同じ意味
  }

  return -LL; // 合計対数尤度の負値を返す
}
")

# TMBの関数をコンパイルして読み込む
compile("model5_8.cpp")
dyn.load(dynlib("model5_8"))

# パラメータの次元・名前・初期値を指定
params <- list(alpha = 0, beta = 0, log_sigma = 0)

# TMBオブジェクトを作成
out5.8 <- MakeADFun(data = tmbData, parameters = params,
DLL = "model5_8", silent = TRUE)

# TMBオブジェクトを最適化して結果を表示
starts <- rep(0, length(unlist(params)))
opt <- optim(starts, fn = out5.8$fn, gr = out5.8$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out5.8))

# TMBの推定値を保存
tmb_est <- c(opt$par[1:2], exp(opt$par[3]))


# 5.9 パラメータ推定値の比較
# -------------------------------------

# すべての推定エンジンの結果を真の値と比較
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est, NIMBLE = nimble_est,
Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 5.10 まとめと展望
# ------------------------
# （コードなし）
