
# ---------------------------------------------------------------
# 書籍付属コード： 'Applied Statistical Modeling 対象： Ecologists'
# 著者：Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# 第13章：Poisson GLM with continuous and categorical
# explanatory variables
# ---------------------------------------------------------------

# 最終更新：11 June 2024



#######################################################################
# --------------------------------------------------------------------#
# 13.2 データの生成
# --------------------------------------------------------------------#
#######################################################################


set.seed(13)
nPops <- 3 # Three populations
nSample <- 100 # . . .with 100 individuals 各
n <- nPops * nSample # Total sample size
x <- rep(1:nPops, rep(nSample, nPops)) # 集団 indicator
pop <- factor(x, labels = c("Pyrenees", "Massif Central", "Jura"))
orig.length <- runif(n, 4.5, 7.0) # Wing length (cm)
wing.length <- orig.length - mean(orig.length) # Center 値

Xmat <- model.matrix(~ pop * wing.length)
head(Xmat, 10) # Look at first 10 rows of design matrix

# 保存する truth 対象： comparisons
truth <- beta.vec <- c(-2, 1, 2, 4, -2, -5)

lin.pred <- Xmat[,] %*% beta.vec # 線形予測子の値
lambda <- exp(lin.pred) # ポアソン 平均: expected count
load <- rpois(n = n, lambda = lambda) # ポアソン乱数による変動を加える

# Inspect what we’ve created (図 13.2)
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.axis = 1.5, cex.lab = 1.5)
hist(load, col = "grey", breaks = 30, xlab = "Parasite load", main = "", las = 1)
plot(wing.length, load, pch = rep(c("P", "M", "J"), each = nSample), las = 1,
col = rep(c("Red", "Green", "Blue"), each = nSample), ylab = "Parasite load",
xlab = "Wing length", cex = 1.2, frame = FALSE)

# 必要なライブラリを読み込む
library(ASMbook); library(jagsUI); library(rstan); library(TMB)

#######################################################################
# --------------------------------------------------------------------#
# 13.3 Rの既存関数を用いた尤度解析
# --------------------------------------------------------------------#
#######################################################################


# 当てはめる モデル を用いる least-squares and 推定値を保存する
summary(out13.3 <- glm(load ~ pop * wing.length, family = poisson))
glm_est <- coef(out13.3)

# モデル assessment by checks of residuals (表示省略)
plot(out13.3) # 確認する of traditional (Pearson) residuals
library(DHARMa) # 計算する quantile residuals based on シミュレーション
simOut <- simulateResiduals(out13.3, n = 1000, plot = TRUE) # 確認する

#######################################################################
# --------------------------------------------------------------------#
# 13.4 JAGSによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################



# 13.4.1 Fitting the モデル
# --------------------------

# データをリストにまとめ構造を確認する
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops,
wing.length = wing.length, n = n) )

# JAGSモデルファイルを書き出す
cat(file = "model13.4.txt", "
model {
# 事前分布
for (i in 1:nPops){ # Loop over 3 populations
  alpha[i] ~ dnorm(0, 0.0001) # 切片
  beta[i] ~ dnorm(0, 0.0001) # 傾き
}
# 尤度
for (i in 1:n) { # Loop over すべての 300 データ points
  load[i] ~ dpois(lambda[i]) # The response variable (C above)
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]]* wing.length[i])
} # Note nested indexing: alpha[pop[i]] and beta[pop[i]]
# 導出量
# Recover 効果 relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1] # 切片の差：中央高地とピレネー
a.effe3 <- alpha[3] - alpha[1] # 切片の差：ジュラとピレネー
b.effe2 <- beta[2] - beta[1] # 傾きの差：中央高地とピレネー
b.effe3 <- beta[3] - beta[1] # 傾きの差：ジュラとピレネー
# Custom test
test1 <- beta[3] - beta[2] # 傾きの差：ジュラと中央高地
}
")

# 初期値を生成する関数
inits <- function(){list(alpha = rlnorm(nPops, 3, 1), beta = rlnorm(nPops, 2, 1))}

# 推定対象のパラメータ
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")
# MCMCの設定
na <- 2000; ni <- 25000; nb <- 5000; nc <- 4; nt <- 5
# JAGSを実行する (原著記載の実行時間（ART） <1 min), 収束を診断する and 事後分布を要約する
out13.4 <- jags(dataList, inits, params, "model13.4.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out13.4) # 図・出力は省略
print(out13.4, 3)

# 比較する glm() estimes with Bayesian 推定値 and with truth
jags_est <- unlist(out13.4$mean)[c(1, 7, 8, 4, 9, 10)]
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est)
print(comp, 4)


# 13.4.2 予測値の算出
# ----------------------------

# 作成する a vector with 100 wing lengths
orig.wlength <- sort(orig.length)
wlength <- orig.wlength - mean(orig.length)
# 作成する matrices to contain prediction 対象： 各 winglength and MCMC iteration
nsamp <- out13.4$mcmc.info$n.samples # Get size of 事後分布 sample
mite.load.Pyr <- mite.load.MC <- mite.load.Ju <- array(dim = c(nsamp, 300))
# Fill in these vectors: this is clumsy, but it works
alpha <- out13.4$sims.list$alpha # 事後分布 of alpha
beta <- out13.4$sims.list$beta # 事後分布 of beta
for (i in 1:300){
  mite.load.Pyr[,i] <- exp(alpha[,1] + beta[,1] * wlength[i])
  mite.load.MC[,i] <- exp(alpha[,2] + beta[,2] * wlength[i])
  mite.load.Ju[,i] <- exp(alpha[,3] + beta[,3] * wlength[i])
}

# 計算する 95% Bayesian credible intervals around predictions
LCB.Pyr <- apply(mite.load.Pyr, 2, quantile, prob = 0.025)
UCB.Pyr <- apply(mite.load.Pyr, 2, quantile, prob = 0.975)
LCB.MC <- apply(mite.load.MC, 2, quantile, prob = 0.025)
UCB.MC <- apply(mite.load.MC, 2, quantile, prob = 0.975)
LCB.Ju <- apply(mite.load.Ju, 2, quantile, prob = 0.025)
UCB.Ju <- apply(mite.load.Ju, 2, quantile, prob = 0.975)

# 計算する 事後分布 means
mean.rel <- cbind(exp(out13.4$mean$alpha[1] + out13.4$mean$beta[1] * wlength),
exp(out13.4$mean$alpha[2] + out13.4$mean$beta[2] * wlength),
exp(out13.4$mean$alpha[3] + out13.4$mean$beta[3] * wlength))
covar <- cbind(orig.wlength, orig.wlength, orig.wlength)

# Plot (Fig. 13.3)
par(mar = c(6, 6, 5, 3), cex.lab = 1.5, cex.axis = 1.5)
matplot(orig.wlength, mean.rel, col = c("red", "green", "blue"), type = "l",
lty = 1, lwd = 2, las = 1, ylab = "Expected mite load",
xlab = "Wing length (cm)", frame = FALSE, ylim = c(0, 25))
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.Pyr, rev(UCB.Pyr)),
col = rgb(1,0,0,0.2), border = NA)
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.MC, rev(UCB.MC)),
col = rgb(0,1,0,0.2), border = NA)
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.Ju, rev(UCB.Ju)),
col = rgb(0,0,1,0.2), border = NA)
matplot(orig.wlength, mean.rel, col = c("red", "green", "blue"),
type = "l", lty = 1, lwd = 3, add = TRUE)

#######################################################################
# --------------------------------------------------------------------#
# 13.5 NIMBLEによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


library(nimble)

# データをリストにまとめ構造を確認する
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops, 
                     wing.length = wing.length, n = n) )

# NIMBLEモデルコードを定義する
model13.5 <- nimbleCode( {
# 事前分布
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, sd = 100)  # 切片
  beta[i] ~ dnorm(0, sd = 100)   # 傾き
}

# 尤度
for (i in 1:n) {
  load[i] ~ dpois(lambda[i])     # 応答変数
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]]* wing.length[i])
}   # Note nested indexing: alpha[pop[i]] and beta[pop[i]]

# 導出量
# Recover 効果 relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1]   # 切片の差：中央高地とピレネー
a.effe3 <- alpha[3] - alpha[1]   # 切片の差：ジュラとピレネー
b.effe2 <- beta[2] - beta[1]     # 傾きの差：中央高地とピレネー
b.effe3 <- beta[3] - beta[1]     # 傾きの差：ジュラとピレネー

# 任意の対比
test1 <- beta[3] - beta[2]       # 傾きの差：ジュラと中央高地
} )

# 初期値
inits <- function(){list(alpha=rlnorm(nPops, 3, 1), beta = rlnorm(nPops, 2, 1))}

# 監視するパラメータは前と同じ
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMCの設定
ni <- 25000  ;  nb <- 5000  ; nc <- 4  ; nt <- 5

# NIMBLEを実行する (原著記載の実行時間（ART） 40 sec), 収束を診断する and 事後分布を要約する
system.time( out13.5 <- 
    nimbleMCMC(code = model13.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out13.5)
(nsum <- nimble_summary(out13.5, params))
nimble_est <- nsum[c(1,7,8,4,9,10),1]    # 保存する param 推定値


#######################################################################
# --------------------------------------------------------------------#
# 13.6 Stanによるベイズ解析
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops,
wing_length = wing.length, n = n) )

# Stanモデルファイルを書き出す
cat(file = "model13_6.stan", "
data {
  int n;
  int nPops;
  array[n] int load;
  vector[n] wing_length;
  array[n] int pop;
}

parameters {
  vector[nPops] alpha;
  vector[nPops] beta;
}

model {
  vector[n] lambda;

  for (i in 1:nPops){
    alpha[i] ~ normal(0, 100);
    beta[i] ~ normal(0, 100);
  }

  for (i in 1:n) {
    lambda[i] = exp(alpha[pop[i]] + beta[pop[i]] * wing_length[i]);
    load[i] ~ poisson(lambda[i]);
  }
}

generated quantities {
  real a_effe2 = alpha[2] - alpha[1];
  real a_effe3 = alpha[3] - alpha[1];
  real b_effe2 = beta[2] - beta[1];
  real b_effe3 = beta[3] - beta[1];
  real test1 = beta[3] - beta[2];
}
")

# HMCの設定
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Stanを実行する (原著記載の実行時間（ART） 60/10 sec), 収束を評価する and 結果を出力する table
system.time(
out13.6 <- stan(file = "model13_6.stan", data = dataList,
warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out13.6) # 図・出力は省略
print(out13.6, dig = 3) # 図・出力は省略
stan_est <- summary(out13.6)$summary[c(1,7,8,4,9,10),1] # 推定値を保存する

#######################################################################
# --------------------------------------------------------------------#
# 13.7 自作の最尤推定
# --------------------------------------------------------------------#
#######################################################################


head(Xmat, 10) # Look at first 10 rows of design matrix

# 定義する NLL 対象： general ポアソン回帰
NLL <- function(beta, y, Xmat) {
  lambda <- exp(Xmat %*% beta) # Multiply design matrix by beta
  LL <- dpois(y, lambda, log = TRUE) # 対数尤度 of 各 obs
  NLL <- -sum(LL) # NLL 対象： すべての 観測値 in データ set
  return(NLL)
}

# 負の対数尤度を最小化する to 最尤推定値を求める, 標準誤差と信頼区間を求める and 推定値を保存する
# Need to use method = "BFGS" here
inits <- rep(0, 6)
names(inits) <- colnames(Xmat)
out13.7 <- optim(inits, NLL, y = load, Xmat = Xmat, hessian = TRUE, method = "BFGS")
get_MLE(out13.7, 4)
diy_est <- out13.7$par

#######################################################################
# --------------------------------------------------------------------#
# 13.8 尤度 による解析： Template モデル Builder
# --------------------------------------------------------------------#
#######################################################################


# データをリストにまとめ構造を確認する
str(tmbData <- list(load = load, pop = as.numeric(pop) - 1,
nPops = nPops, wing_length = wing.length, n = n) ) # 図・出力は省略

# TMBモデルファイルを書き出す
cat(file = "model13_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type> ::operator() ()
{
  // 入力データを定義する
  DATA_VECTOR(load); // 応答変数
  DATA_VECTOR(wing_length); // 共変量
  DATA_IVECTOR(pop); // 集団 index (of integers; IVECTOR)
  DATA_INTEGER(n); // 観測値数

  // パラメータを定義する
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type LL = 0; // 対数尤度を0で初期化する
  for (int i = 0; i < n; i ++ ){
    Type lambda = exp(alpha(pop(i)) + beta(pop(i)) * wing_length(i));
    // 各観測値の対数尤度を計算し合計に加える
    LL += dpois(load(i), lambda, true);
  }

  // 導出パラメータ
  Type a_effe2 = alpha(1) - alpha(0);
  Type a_effe3 = alpha(2) - alpha(0);
  Type b_effe2 = beta(1) - beta(0);
  Type b_effe3 = beta(2) - beta(0);
  Type test1 = beta(2) - beta(1);
  ADREPORT(a_effe2);
  ADREPORT(a_effe3);
  ADREPORT(b_effe2);
  ADREPORT(b_effe3);
  ADREPORT(test1);
  return -LL; // 負の対数尤度を返す
}
")

# TMBモデルをコンパイルして動的ライブラリを読み込む
compile("model13_8.cpp")
dyn.load(dynlib("model13_8"))

# パラメータの次元と初期値を指定する
params <- list(alpha = rep(0, tmbData$nPops), beta = rep(0, tmbData$nPops))

# TMBの目的関数オブジェクトを作成する
out13.8 <- MakeADFun(data = tmbData,
                     parameters = params,
                     DLL = "model13_8", silent = TRUE)

# TMBオブジェクトを最適化する and 結果を出力する
opt <- optim(out13.8$par, fn = out13.8$fn, gr = out13.8$gr, method = "BFGS")
(tsum <- tmb_summary(out13.8))
tmb_est <- tsum[c(1, 7, 8, 4, 9, 10), 1] # 保存する パラメータ 推定値

#######################################################################
# --------------------------------------------------------------------#
# 13.9 パラメータ推定値の比較
# --------------------------------------------------------------------#
#######################################################################

# 推定結果を真値と比較する and previous 推定値
comp <- cbind(truth = beta.vec, glm = glm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)
