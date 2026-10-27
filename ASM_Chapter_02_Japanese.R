
# ---------------------------------------------------------------
# 書籍 'Applied Statistical Modeling for Ecologists' のコード 
# Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------


# ----------------------------------------------------
# 第2章 -- 統計的推論への導入
# ----------------------------------------------------


# このファイルには、この章で使用するすべてのコードが含まれている。
# 書籍本文に掲載されているコードに加えて、いくつかの拡張コードも含まれる。
# それらの結果は本文中で説明・議論されているが、
# コード自体は書籍には掲載されていないか、一部のみ掲載されている。
# これらの追加コードは2.6節に含まれている。


# 最終更新: 2024年7月8日


# 2.1 統計的推論の基礎としての確率
# ------------------------------------------------------
# （コードなし）


# 2.2 確率変数と確率分布による記述
# -----------------------------------------------------------------------
# （コードなし）


#   2.2.1 離散型確率変数と連続型確率変数
# ------------------------------------------------
# （コードなし）


#   2.2.2 確率分布による確率変数の記述
# --------------------------------------------------------------------

# 離散確率分布の例としてPoisson分布を扱う
# Poisson分布のヘルプを表示する
?dpois # '?Poisson' でもよい

# y = 0 におけるPoisson分布の確率質量関数（PMF）を計算する
dpois(0, lambda = 2)                 # 0が得られる確率
dpois(0, lambda = 2, log = TRUE)     # 同じ値を対数スケールで計算

# dpois() が省略表現であることを確認するため、同じ確率を式から直接計算する
lam <- 2; y <- 0
(lam^y)*exp(-lam) / factorial(y)     # 0が得られる確率

# 累積分布関数（CDF）を計算する
ppois(3, lambda = 2)                 # 0, 1, 2, 3 のいずれかが得られる確率
sum(dpois(0:3, 2))                   # 各確率を直接足し合わせても同じ
plot(ppois(0:9, 2), type = 'h', lend = 'butt', lwd = 20, ylim = c(0, 1), 
  xlab = 'Yの値', ylab = 'P(Y <= y)', main = "Yの累積分布関数")  # 書籍には非掲載

# 分位点関数：CDF(y)が指定した値に達するyの値を求める
qpois(0.85, lambda = 2)

# 乱数生成関数（RNG）
set.seed(2016)                       # 同じ乱数を再現するため乱数シードを固定
rpois(n = 10, lambda = 2)            # Poisson(2)から乱数を10個生成

# lambda = 2 のPoisson分布の確率質量関数を図示する（図2.2）
par(mar = c(6,8,5,4), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
plot(0:10, dpois(0:10, lambda = 2), xlab = '確率変数の実現値',
  ylab = '確率', main = 'Poisson(2) の確率質量関数', type = 'h', lend = 'butt',
  lwd = 20, col = 'gray30', frame = FALSE, las = 1)


# 連続確率分布の例として正規分布（Gaussian分布）を扱う
?dnorm # '?Normal' でもよい：正規分布のヘルプを表示

# y = 650 における正規分布の確率密度関数（PDF）を計算する
dnorm(650, mean = 600, sd = 30)
dnorm(650, mean = 600, sd = 30, log = TRUE) # 同じ値を対数スケールで計算

# dnorm() が省略表現であることを確認するため、同じ確率密度を式から直接計算する
mu <- 600; sig <- 30; y <- 650
1/(sig*sqrt(2*pi)) * exp(-(y - mu)^2 / (2 * sig^2))

# 累積分布関数（CDF）を計算する
pnorm(600, mean = 600, sd = 30)      # -Infから600までの累積確率

# 次は同じ累積確率を、PDFを直接積分して求める
integrate(dnorm, lower = -Inf, upper = 600, mean = 600, sd = 30)
plot(pnorm(seq(400, 800, by = 1), 600, 30), type = 'l', lwd = 5, ylim = c(0, 1),
xlab = 'Yの値', ylab = '累積確率', main = "Yの累積分布関数", frame = FALSE) # 書籍には非掲載

# 分位点関数：CDF(y)が指定した値に達するyの値を求める
qnorm(0.95, mean = 600, sd = 30)

# 乱数生成関数（RNG）
set.seed(2016)
rnorm(n = 10, mean = 600, sd = 30)   # Normal(600, 30)から乱数を10個生成


# 区間分けした正規確率変数について各区間の確率を計算する
limits <- seq(500, 700, by = 5)
midpts <- seq(500 + 5/2, 700-5/2, by = 5)
cumProb <- pnorm(limits, mean = 600, sd = 30)
probs <- diff(cumProb)

# Normal(600, 30) の確率密度関数（PDF）を図示する
par(mfrow = c(1, 2), mar = c(6,6,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5) # Fig. 2.3
curve(dnorm(x, mean = 600, sd = 30), 500, 700, xlab = '体重 (g)',
  ylab = '確率密度', main = 'Normal(600, 30) の確率密度関数', type = 'l',
  lwd = 3, col = 'gray30', frame = FALSE, las = 1)
# 5 g区分した体重の確率変数について各区間の確率を図示する
plot(midpts, probs, xlab = '5 g区分した体重 (g)', ylab = '確率',
  main = '区分化した正規確率変数の確率質量関数', type = 'h', lend = 'butt', lwd = 5,
  col = 'gray30', frame = FALSE, las = 1)


#   2.2.3 確率変数の位置とばらつき、および「パラメータをモデル化する」とは何か
# ------------------------------------------------
# （コードなし）


#   2.2.4 確率変数と確率分布についての短いまとめ
# ------------------------------------------------
# （コードなし）


# 2.3 統計モデルとその利用
# ---------------------------------------

# パラメータ値によってBeta分布の形がどのように変わるかを見る
par(mfrow = c(2,2)) # Fig. 2.5
curve(dbeta(x,10,10), 0, 1, lwd = 3, col = 1, main = "beta(10, 10)")
curve(dbeta(x,0.1,0.1), 0, 1, lwd = 3, col = 1, main = "beta(0.1, 0.1)")
curve(dbeta(x,1,100), 0, 1, lwd = 3, col = 1, main = "beta(1, 100)")
curve(dbeta(x,100,1), 0, 1, lwd = 3, col = 1, main = "beta(100, 1)")


# 2.4 尤度関数
# ---------------------------
# （コードなし）


# 2.5 最尤法による古典的推論と1パラメータモデルへの適用
# ---------------------------------------------------------------------------------------------
# （コードなし）


#   2.5.1 最尤法への導入
# ------------------------------------------

# 期待個体数 lam をもつカウントデータを2組作成する
set.seed(2)
lam <- 2
y1 <- rpois(10, lambda = lam)        # 小標本データ（n = 10）
y2 <- rpois(100, lambda = lam)       # 大標本データ（n = 100）
table(y1) ; table(y2)

# Rで尤度関数を定義する
L <- function(lambda, y) {
  Li <- dpois(y, lambda)             # 各データ点iの尤度への寄与
  L <- prod(Li)                      # データ全体の尤度は各寄与の積
  return(L)
}

# Rで対数尤度関数を定義する
LL <- function(lambda, y) {
  LLi <- dpois(y, lambda, log = TRUE) # 各データ点iの対数尤度への寄与
  LL <- sum(LLi)                      # データ全体の対数尤度は各寄与の和
  return(LL)
}

# Rで負の対数尤度関数を定義する
NLL <- function(lambda, y) {
  LL <- dpois(y, lambda, log =TRUE)   # 各データ点iの対数尤度への寄与
  NLL <- -sum(LL)                     # データ全体の負の対数尤度は和にマイナスを付けたもの
  return(NLL)
}

# lambda の広い候補範囲について3つの関数を評価する
possible.lam <- seq(0.001, 4, by = 0.001)
like1 <- loglike1 <- negloglike1 <- numeric(length(possible.lam))
for(i in 1:length(like1)){
  like1[i] <- L(possible.lam[i], y1)         # Likelihood
  loglike1[i] <- LL(possible.lam[i], y1)     # log-Likelihood
  negloglike1[i] <- NLL(possible.lam[i], y1) # negative log-likelihood
}

# 尤度・対数尤度・負の対数尤度のプロファイルを図示する（図2.6）
par(mfrow = c(1, 3), mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(possible.lam, like1, xlab = 'lambda', ylab = 'L', main = '尤度', frame = FALSE)
abline(v = possible.lam[which(like1 == max(like1))], col = 'gray', lwd = 3, lty = 3)
plot(possible.lam, loglike1, xlab = 'lambda', ylab = 'LL', main = '対数尤度', frame = FALSE)
abline(v = possible.lam[which(loglike1 == max(loglike1))], col = 'gray', lwd = 3, lty = 3)
plot(possible.lam, negloglike1, xlab = 'lambda', ylab = '-LL', main = '負の対数尤度', frame = FALSE)
abline(v = possible.lam[which(negloglike1 == min(negloglike1))], col = 'gray', lwd = 3, lty = 3)


# MLEを「手作業」で求める
possible.lam[which(like1 == max(like1))]             # Maximum likelihood
possible.lam[which(loglike1 == max(loglike1))]       # Maximum log-likelihood
possible.lam[which(negloglike1 == min(negloglike1))] # Minimum negative log-likelihood


# Rの optim() による関数最小化を使ってMLEを求める
inits <- c('lambda' = 1)
out <- optim(inits, NLL, y = y1)     # y1についてlambdaを変化させて関数を最適化
out
(MLE <- out$par) # MLEを取り出す


# 大標本についてL、LL、NLLを評価する
like2 <- loglike2 <- negloglike2 <- numeric(length(possible.lam))
for(i in 1:length(like2)){
  like2[i] <- L(possible.lam[i], y2)         # Likelihood
  loglike2[i] <- LL(possible.lam[i], y2)     # log-Likelihood
  negloglike2[i] <- NLL(possible.lam[i], y2) # negative log-likelihood
}

# 2つのデータセットについて最大化された対数尤度の値を計算する
(lmax1 <- max(loglike1, na.rm = TRUE))
(lmax2 <- max(loglike2, na.rm = TRUE))

# 2つの標本についてMLE近傍の対数尤度を図示する（図2.7）
par(mfrow = c(1, 1))
ylim <- c(-5, 0) # y軸の範囲を設定
plot(possible.lam, loglike1-max(loglike1), xlab = 'lambda', ylab = '基準化した対数尤度',
   type = 'l', frame = FALSE, col = 'gray', lwd = 5, ylim = ylim, xlim = c(1, 3))
lines(possible.lam, loglike2-max(loglike2), col = 'gray40', lwd = 5)
abline(v = possible.lam[which(loglike1 == max(loglike1))], lwd = 3, lty = 3, 
  col = 'gray') # 小標本データのMLE
abline(v = possible.lam[which(loglike2 == max(loglike2))], lwd = 3, lty = 3, 
  col = 'grey40') # 大標本データのMLE
legend('bottomright', col = c('gray', 'gray40'), lwd = 4, bty = 'n', 
  legend = c("n = 10", "n = 100"), cex = 1.5)


#   2.5.2 尤度比検定の逆転による信頼区間
# -----------------------------------------------------

# 信頼区間の方法1：プロファイル尤度またはLRTに基づく95%信頼区間
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––
(lmax <- max(loglike1, na.rm = TRUE)) # 最大化された対数尤度
(llim <- lmax - 1.92)                 # 最大値から3.84/2だけ低い値

# MLE周辺に範囲を絞ってプロファイルを図示する（図2.8）
par(mar = c(6, 6, 6, 3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
ooo <- which(abs(loglike1-lmax) < 3.5)
xlim <- c(0.9 * possible.lam[min(ooo)], 1.1*possible.lam[max(ooo)])
ylim <- c(-20, -16)
plot(possible.lam[ooo], loglike1[ooo], type = 'l', lwd = 3,
xlab = 'lambda', ylab = '対数尤度', frame = FALSE, xlim = xlim,
ylim = ylim, main = 'MLE近傍の対数尤度曲線\n（LRTに基づく95%信頼区間）')
# グラフに、最大対数尤度からカイ二乗値/2だけ下がった位置を示す
idx <- which(loglike1 == max(loglike1)) # 最大値を含む位置
arrows(possible.lam[idx], llim, possible.lam[idx], lmax, length = 0.25, angle = 30, code = 3, lwd = 2)
abline(h = llim, lty = 2, lwd = 3)
text(2.4, mean(c(llim, lmax)), "1.92", cex = 2, srt = 90)
# プロファイルの左側と右側で信頼区間の2つの限界を別々に求める
l.left <- loglike1[1:idx]                      # 左側の枝
l.right <- loglike1[(idx+ 1):length(loglike1)] # 右側の枝

# llimと左枝・右枝の値を比較する
idx.left <- which(abs(l.left - llim) == min(abs(l.left - llim), na.rm = TRUE))
idx.right <- which(abs(l.right - llim) == min(abs(l.right - llim), na.rm = TRUE))

# プロファイル尤度による信頼限界を計算する
(LCL <- possible.lam[1:idx][idx.left])
(UCL <- possible.lam[(idx+1):length(loglike1)][idx.right])

# 信頼限界を図に追加する
segments(possible.lam[idx.left], -21, possible.lam[idx.left], llim, lwd = 2, lty = 3)
segments(possible.lam[idx+idx.right], -21, possible.lam[idx+idx.right], llim, lwd = 2, lty = 3)
arrows(possible.lam[idx.left], -20, possible.lam[idx+idx.right], -20, length = 0.25,
angle = 30, code = 3, lwd = 2)
text(2.4, -19.5, "95% CI", cex = 2, srt = 0)

# Rに用意されているPoisson GLMのプロファイル信頼区間と比較する
exp(confint(glm(y1 ~ 1, family = 'poisson')))


#   2.5.3 MLEの近似正規性に基づく標準誤差と信頼区間
# -------------------------------------------------------------------------------------------

# 標準誤差と信頼区間の方法2：Hessianの逆行列を用いる
# –––––––––––––––––––––––––––––––––––––––––––––––––
inits <- c(lambda = 1)
# 小標本データ（n = 10）
(out1 <- optim(inits, NLL, y = y1, hessian= TRUE))
(MLE1 <- out1$par)                   # Grab MLEs
(VC1 <- solve(out1$hessian))         # 分散共分散行列を得る
(VC1 <- 1/out1$hessian)              # 1パラメータモデルでは同じ
(ASE1 <- sqrt(diag(VC1)))            # 漸近標準誤差を取り出す

# 大標本データ（n = 100）
(out2 <- optim(inits, NLL, y = y2, hessian=TRUE))
(MLE2 <- out2$par)                   # Grab MLEs
(VC2 <- solve(out2$hessian))         # 分散共分散行列を得る
(VC2 <- 1/out2$hessian)              # 1パラメータモデルでは同じ
(ASE2 <- sqrt(diag(VC2)))            # 漸近標準誤差を取り出す

# 真値、MLE、漸近標準誤差（ASE）を小標本・大標本で比較する
print(cbind('truth' = lam, MLE1, ASE1, MLE2, ASE2), 5)

# 小標本についてWald型信頼区間を自分で計算する
print(CI1 <- c(MLE1-1.96 * ASE1, MLE1 + 1.96 * ASE1), 4)

# 小標本についてLRTに基づく信頼区間をRで求める
print(exp(confint(glm(y1 ~ 1, family = 'poisson'))), 4)

# 大標本についてWald型信頼区間を自分で計算する
print(CI2 <- c(MLE2-1.96 * ASE2, MLE2 + 1.96 * ASE2), 4)

# 大標本についてLRTに基づく信頼区間をRで求める
print(exp(confint(glm(y2 ~ 1, family = 'poisson'))), 4)


#   2.5.4 ブートストラップによる標準誤差と信頼区間
# ------------------------------------------------------------------------------

# 標準誤差と信頼区間の方法3.1：パラメトリック・ブートストラップ
# ––––––––––––––––––––––––––––––––––––––––––––––––
simrep <- 10000                      # ブートストラップ標本数：十分大きな値
estiPB <- array(NA, dim = c(simrep, 2)) # 推定値を保存する配列
colnames(estiPB) <- c('小標本データ', '大標本データ')

for(i in 1:simrep){
  if(i %% 500 == 0) cat(paste("iter", i, "\n")) # 進捗表示
  yb1 <- rpois(10, lambda = MLE1)    # 小標本データを新たに生成
  fm1 <- summary(glm(yb1~1, family = 'poisson')) # 同じモデルを再適合
  estiPB[i,1] <- exp(fm1$coef[1])    # 自然尺度で推定値を保存
  yb2 <- rpois(100, lambda = MLE2)   # 大標本データを新たに生成
  fm2 <- summary(glm(yb2~1, family = 'poisson')) # 同じモデルを再適合
  estiPB[i,2] <- exp(fm2$coef[1])    # 推定値を保存
}

# パラメトリック・ブートストラップによるMLEの標本分布（図2.9）
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)
xlim <- range(estiPB[,1])
hist(estiPB[,1], col = 'gray', main = 'Bootstrapped sampling distribution of MLE 
  of lambda (small data set; parametric bootstrap)', xlim = xlim)
hist(estiPB[,2], col = 'gray', main = 'Bootstrapped sampling distribution of MLE
  of lambda (large data set; parametric bootstrap)', xlim = xlim)

# パラメトリック・ブートストラップによる標準誤差
se1.pb <- sd(estiPB[,1])             # Small data set
se2.pb <- sd(estiPB[,2])             # Large data set

# パラメトリック・ブートストラップによる信頼区間
CI1.pb <- quantile(estiPB[,1], c(0.025, 0.975)) # Small data set
CI2.pb <- quantile(estiPB[,2], c(0.025, 0.975)) # Large data set


# 標準誤差と信頼区間の方法3.2：ノンパラメトリック・ブートストラップ
# –––––––––––––––––––––––––––––––––––––––––––––––––––
simrep <- 10000                          # ブートストラップ標本数：十分大きな値
estiNPB <- array(NA, dim = c(simrep, 2)) # 推定値を保存する配列
colnames(estiNPB) <- c('小標本データ', '大標本データ')
for(i in 1:simrep){
  if(i %% 500 == 0) cat(paste("iter", i, "\n")) # 進捗表示
  yb1 <- sample(y1, 10, replace = TRUE)  # 小標本データを復元抽出
  fm1 <- summary(glm(yb1~1, family = 'poisson')) # モデルを再適合
  estiNPB[i,1] <- exp(fm1$coef[1])       # 推定値を保存
  yb2 <- sample(y2, 100, replace = TRUE) # 大標本データを復元抽出
  fm2 <- summary(glm(yb2~1, family = 'poisson')) # モデルを再適合
  estiNPB[i,2] <- exp(fm2$coef[1])       # 推定値を保存
}

# ノンパラメトリック・ブートストラップによるMLEの標本分布
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6) # 書籍には非掲載
xlim <- range(estiNPB[,1])
hist(estiNPB[,1], col = 'gray', main = 'Bootstrapped sampling distribution of MLE
  of lambda (small data set; non-parametric bootstrap)', xlim = xlim)
hist(estiNPB[,2], col = 'gray', main = 'Bootstrapped sampling distribution of MLE
  of lambda (large data set; non-parametric bootstrap)', xlim = xlim)
  
# ノンパラメトリック・ブートストラップによる標準誤差
se1.npb <- sd(estiNPB[,1])                        # Small data set
se2.npb <- sd(estiNPB[,2])                        # Large data set

# ノンパラメトリック・ブートストラップによる信頼区間
CI1.npb <- quantile(estiNPB[,1], c(0.025, 0.975)) # Small data set
CI2.npb <- quantile(estiNPB[,2], c(0.025, 0.975)) # Large data set


# 標準誤差
# –––––––––––––––
# 小標本：漸近正規性と2種類のブートストラップを比較する
print(cbind(ASE1, 'PB-SE' = se1.pb, 'NPB-SE' = se1.npb), 4)

# 大標本：漸近正規性と2種類のブートストラップを比較する
print(cbind(ASE2, 'PB-SE' = se2.pb, 'NPB-SE' = se2.npb), 4)

# 95%信頼区間
# ––––––––––––––––––––––––
# LRTに基づく信頼区間を再計算する
CIp1 <- exp(confint(glm(y1 ~ 1, family = 'poisson')))
CIp2 <- exp(confint(glm(y2 ~ 1, family = 'poisson')))

# 小標本の信頼区間：LRT、正規近似、ブートストラップを比較する
print(cbind(CIp1, CI1, CI1.pb, CI1.npb), 4)

# 大標本の信頼区間：LRT、正規近似、ブートストラップを比較する
print(cbind(CIp2, CI2, CI2.pb, CI2.npb), 4)


#   2.5.5 最尤推定についての短いまとめ
# ------------------------------------------------------
# （コードなし）


# 2.6 2パラメータモデルにおける最尤推定
# --------------------------------------------------------

# Swiss Bee-eaterのデータ（1990–2020年の国内つがい数）
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1990–2020年に国内で確認されたつがい数
y <- c(0, 2, 7, 5, 1, 2, 4, 8, 10, 11, 16, 11, 10, 13, 19, 31, 20, 26, 19, 21, 34,
35, 43, 53, 66, 61, 72, 120, 102, 159, 199)
year <- 1990:2020                    # 年の範囲を定義
x <- (year-1989)                     # スケーリングしたが中心化していない年共変量
x <- x-16                            # 年共変量をさらに中心化

# bee-eaterデータを図示する（図2.11）
par(mfrow = c(1, 2), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(y, xlab = '個体数 (y)', ylab = '頻度', breaks = 20,
  col = 'gray', main = '個体数の頻度分布')
plot(year, y, xlab = '年 (x)', ylab = '個体数 (y)', frame = FALSE, cex = 1.5,
  pch = 16, col = 'gray20', main = 'yとxの関係')
fm <- glm(y ~ x, family = 'poisson') # Poisson GLMの当てはめ線を追加
lines(year, predict(fm, type = 'response'), lwd = 3, col = 'red', lty = 3)


# 過分散を簡単に確認する
# ------------------------------
summary(fmX <- glm(y ~ x, family = 'quasipoisson'))
(c_hat <- summary(fmX)$disp)
# （quasipoisson族の分散パラメータは2.375374）
pchisq(c_hat, df = 31 - 3)
# [1] 4.219326e-11


# 過分散の検定としてパラメトリック・ブートストラップを行う
simrep <- 10000
c_hat_boot <- numeric(simrep)

for(i in 1:simrep){
  yb <- rpois(length(y), lambda=exp(predict(fmX))) # 新しいbee-eaterデータを生成
  fmb <- summary(glm(yb~x, family = 'quasipoisson')) # モデルを再適合してc_hatを推定
  c_hat_boot[i] <- fmb$disp          # 推定値を保存 of OD factor c_hat
}

# 過分散係数 c_hat のパラメトリック・ブートストラップ標本分布
par(mfrow = c(1, 1), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)
hist(c_hat_boot, col = 'grey', main = 'c_hat のブートストラップ標本分布', breaks = 60)
abline(v = c_hat, col = 'red', lwd = 3)
# c_hat のブートストラップ分布における百分位を検定のp値として用いる
mean(c_hat_boot < c_hat)
#[1] 0.9999               # ここでも過分散を強く示す
# -----------------------------------------------------------------

# RからMLE（IWLSによる解）を得る
fm <- glm(y ~ x, family = 'poisson')
summary(fm)

#  Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 2.886821   0.052659   54.82   <2e-16 ***
# x           0.147230   0.004864   30.27   <2e-16 ***


# 対数尤度関数を定義する（Poisson対数線形モデル／Poisson GLMの対数尤度）
LL <- function(param, y, x) {
  alpha <- param[1]        # 切片
  beta <- param[2]         # 傾き
  lambda <- exp(alpha + beta * x)   # 線形予測子
  LLi <- dpois(y, lambda, log=TRUE) # 各観測iの対数尤度への寄与
  LL <- sum(LLi)           # yの全観測についての対数尤度
  return(LL)
}


# Poisson対数線形モデル／Poisson GLMの負の対数尤度（NLL）を定義する
NLL <- function(param, y, x) {
  alpha <- param[1]        # 切片
  beta <- param[2]         # 傾き
  lambda <- exp(alpha + beta * x)   # 線形予測子
  LLi <- dpois(y, lambda, log=TRUE) # 各観測iの対数尤度への寄与
  LL <- -sum(LLi)          # yの全観測についての負の対数尤度（マイナス符号に注意）
  return(LL)
}

# パラメータの組について尤度（厳密には負の対数尤度）を評価する
# 値が小さいほど適合がよい
NLL(c(0, 0), y, x)
NLL(c(0, 0.5), y, x)
NLL(c(1, 0), y, x)


# NLLが最小となるパラメータ値を大域的なグリッドサーチで探索する
# NLLを評価するパラメータ値を設定する
nside <- 1000   # グリッドの分解能、したがって近似精度を決める
try.alpha <- seq(2.5, 3.5, length.out = nside)     # 切片用
try.beta <- seq(0.1, 0.2, length.out = nside)      # 傾き用
nll.grid <- array(NA, dim = c(nside, nside))

# グリッド上のすべてのパラメータ組についてNLLを評価する
for(i in 1:nside){
  for(j in 1:nside){  
     nll.grid[i,j] <- NLL(c(try.alpha[i], try.beta[j]), y, x)
  }
}

# 最小値を求める
(best <- which(nll.grid == min(nll.grid)))
(best.point <- cbind(rep(try.alpha, nside), rep(try.beta, each=nside))[best,])

#> (best <- which(nll.grid == min(nll.grid)))
#[1] 472387
#> (best.point <- cbind(rep(try.alpha, nside), rep(try.beta, each=nside))[best,])
#[1] 2.8863864 0.1472472


# NLLを最小化し、Hessianの逆行列による二次近似から標準誤差を求める
inits <- c(alpha=0, beta=0)
# out <- optim(inits, NLL, y=y, x=x, hessian=TRUE, method = 'Nelder-Mead')
out <- optim(inits, NLL, y=y, x=x, hessian=TRUE, method = 'BFGS')
# out <- optim(inits, NLL, y=y, x=x, hessian=TRUE, method = 'CG')
# out <- optim(inits, NLL, y=y, x=x, hessian=TRUE, method = 'SANN')
MLE <- out$par                 # Grab MLEs
VC <- solve(out$hessian)       # 分散共分散行列を得る
ASE <- sqrt(diag(VC))          # 漸近標準誤差を取り出す
print(cbind(MLE, ASE), 4)      # Print MLEs and SE's

getMLE(out, 5)           # short, but with CIs in addition

# Nelder-Mead
#         MLE      ASE
#alpha 2.8867 0.052661
#beta  0.1472 0.004864

# BFGS
#         MLE      ASE
#alpha 2.8868 0.052659
#beta  0.1472 0.004864

# CG
#         MLE     ASE
#alpha 0.8537 0.13197
#beta  0.3043 0.01046

# SANN
#         MLE      ASE
#alpha 2.8935 0.052454
#beta  0.1467 0.004848

# Nelder-Mead法とBFGS法はいずれも解を見つけるが、結果は完全には同一ではない

#          MLE       ASE LCL.95  UCL.95
#alpha 2.88683 0.0526585 2.7836 2.99004
#beta  0.14723 0.0048635 0.1377 0.15676


# 図を作成する
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100), axes = FALSE,
  xlab = "切片 (alpha)", ylab = "傾き (beta)")
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 2, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.05))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.02))
box()
title(main = "負の対数尤度曲面", font.main = 1)
abline(v = MLE[1], lwd = 3, lty = 2, col = 'blue')    # NLLの数値最小化から得た値
abline(h = MLE[2], lwd = 3, lty = 2, col = 'blue')
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # グリッドサーチから得た値


# 同じ図に95%および99%信頼領域を追加する
# ------------------------------------------------------------------
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100), axes = FALSE,
  xlab = "切片 (alpha)", ylab = "傾き (beta)")
#contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 2, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.05))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.02))
box()
title(main = "負の対数尤度曲面", font.main = 1)
abline(v = MLE[1], lwd = 3, lty = 2, col = 'blue')    # NLLの数値最小化から得た値
abline(h = MLE[2], lwd = 3, lty = 2, col = 'blue')
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # グリッドサーチから得た値

out$value  # 関数の最小値
(x95 <- qchisq(0.95, 1))     # 自由度1のカイ二乗分布の95%点
(x99 <- qchisq(0.99, 1))     # 自由度1のカイ二乗分布の99%点
levs <- c(out$value + x95/2, out$value + x99/2)  # 95%および99%信頼領域を定義する水準
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 2, labcex = 1.5,
    levels = levs, labels = c("95%", "99%") )
# 計算上はうまくいくが、この図では見栄えがあまりよくない


# より局所的なグリッドサーチを行い、MLE周辺を拡大する
# NLLを評価するパラメータ値を設定する
# この設定が最も見やすい図になる
nside <- 1000   # グリッドの分解能、したがって近似精度を決める
try.alpha <- seq(2.8, 2.96, length.out = nside)    # 切片
try.beta <- seq(0.135, 0.16, length.out = nside)   # 傾き
# やや見にくい設定
#try.alpha <- seq(2.75, 3, length.out = nside)     # 切片
#try.beta <- seq(0.135, 0.16, length.out = nside)  # 傾き
nll.grid <- array(NA, dim = c(nside, nside))

# グリッド上のすべてのパラメータ組についてNLLを評価する
for(i in 1:nside){
  for(j in 1:nside){  
     nll.grid[i,j] <- NLL(c(try.alpha[i], try.beta[j]), y, x)
  }
}

# 最小値を求める
(best <- which(nll.grid == min(nll.grid)))
(best.point <- cbind(rep(try.alpha, nside), rep(try.beta, each=nside))[best,])
# [1] 2.8868068 0.1472372

# 図を作成する
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100), axes = FALSE,
  xlab = "切片 (alpha)", ylab = "傾き (beta)")
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 1.5, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.05))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.005))
box()
title(main = "負の対数尤度曲面", font.main = 1)
abline(v = MLE[1], lwd = 3, lty = 2, col = 'blue')    # NLLの数値最小化から得た値
abline(h = MLE[2], lwd = 3, lty = 2, col = 'blue')
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # グリッドサーチから得た値


# 同じ図に95%および99%信頼領域を追加する
# ------------------------------------------------------------------
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100), axes = FALSE,
  xlab = "切片 (alpha)", ylab = "傾き (beta)")
#contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 2, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.05))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.02))
box()
title(main = "負の対数尤度曲面", font.main = 1)
abline(v = MLE[1], lwd = 3, lty = 2, col = 'blue')    # NLLの数値最小化から得た値
abline(h = MLE[2], lwd = 3, lty = 2, col = 'blue')
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # グリッドサーチから得た値

out$value  # 関数の最小値
(x95 <- qchisq(0.95, 1))     # 自由度1のカイ二乗分布の95%点
(x99 <- qchisq(0.99, 1))     # 自由度1のカイ二乗分布の99%点
levs <- c(out$value + x95/2, out$value + x99/2)  # 95%および99%信頼領域を定義する水準
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 2, labcex = 1.5,
    levels = levs, labels = c("95%", "99%") )
# 改善したが、この信頼領域を示すにはより広いグリッド範囲の方がよい


# さらに拡大し、グリッドサーチの各軸の分解能を1000にする
# NLLを評価するパラメータ値を設定する
nside <- 1000   # グリッドの分解能、したがって近似精度を決める
try.alpha <- seq(2.85, 2.93, length.out = nside)     # 切片
try.beta <- seq(0.145, 0.15, length.out = nside)     # 傾き
nll.grid <- array(NA, dim = c(nside, nside))

# グリッド上のすべてのパラメータ組についてNLLを評価する
for(i in 1:nside){
  for(j in 1:nside){  
     nll.grid[i,j] <- NLL(c(try.alpha[i], try.beta[j]), y, x)
  }
}

# 最小値を求める
(best <- which(nll.grid == min(nll.grid)))
(best.point <- cbind(rep(try.alpha, nside), rep(try.beta, each=nside))[best,])
# [1] 2.8868368 0.1472272

# 図を作成する
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100), axes = FALSE,
  xlab = "切片 (alpha)", ylab = "傾き (beta)")
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 1.5, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.01))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.001))
box()
title(main = "負の対数尤度曲面", font.main = 1)
abline(v = MLE[1], lwd = 3, lty = 2, col = 'blue')    # NLLの数値最小化から得た値
abline(h = MLE[2], lwd = 3, lty = 2, col = 'blue')
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # グリッドサーチから得た値
# この図はやや見にくい


# パラメトリック・ブートストラップによる標準誤差と信頼区間
# -----------------------------------
simrep <- 10000        # ブートストラップ標本数：十分大きな値
estiPB <- array(NA, dim = c(simrep, 2)) # 推定値を保存する配列
colnames(estiPB) <- c('alpha', 'beta')
# 各年についてlambdaを計算する
lam <- exp(MLE[1] + MLE[2] * x)

for(i in 1:simrep){
  yb <- rpois(n, lambda=lam)         # 新しいbee-eaterデータを生成
  fmb <- summary(glm(yb~x, family = 'poisson')) # 同じモデルを再適合
  estiPB[i,] <- fmb$coef[,1]    # 自然尺度で推定値を保存
}

# パラメトリック・ブートストラップによるMLEの標本分布（図2–7）
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)
hist(estiPB[,1], col = 'grey', main = 'MLEのブートストラップ標本分布\n切片 alpha（パラメトリック）', breaks = 60)
hist(estiPB[,2], col = 'grey', main = 'MLEのブートストラップ標本分布\n傾き beta（パラメトリック）', breaks = 50)

# パラメトリック・ブートストラップによる標準誤差
se.alpha.pb <- sd(estiPB[,1])        # 切片 (alpha)
se.beta.pb <- sd(estiPB[,2])         # 傾き (beta)

# パラメトリック・ブートストラップによる信頼区間
CI.alpha.pb <- quantile(estiPB[,1], c(0.025, 0.975)) # 切片 (alpha)
CI.beta.pb <- quantile(estiPB[,2], c(0.025, 0.975)) # 傾き (beta)


# ノンパラメトリック・ブートストラップによる標準誤差と信頼区間
# ---------------------------------------
simrep <- 10000        # ブートストラップ標本数：十分大きな値
estiNPB <- array(NA, dim = c(simrep, 2)) # 推定値を保存する配列
colnames(estiNPB) <- c('alpha', 'beta')

for(i in 1:simrep){
  idx <- sample(1:n, n, replace = TRUE)    # 各データ組の添字を復元抽出
  yb <- y[idx]
  xb <- x[idx]  
  fmb <- summary(glm(yb~xb, family = 'poisson')) # 同じモデルを再適合
  estiNPB[i,] <- fmb$coef[,1]    # 推定値を保存
}

# ノンパラメトリック・ブートストラップによるMLEの標本分布
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)
hist(estiNPB[,1], col = 'grey', main = 'MLEのブートストラップ標本分布\n切片 alpha（ノンパラメトリック）', breaks = 60)
hist(estiNPB[,2], col = 'grey', main = 'MLEのブートストラップ標本分布\n傾き beta（ノンパラメトリック）', breaks = 50)

# ノンパラメトリック・ブートストラップによる標準誤差
se.alpha.npb <- sd(estiNPB[,1])        # 切片 (alpha)
se.beta.npb <- sd(estiNPB[,2])         # 傾き (beta)

# ノンパラメトリック・ブートストラップによる信頼区間
CI.alpha.npb <- quantile(estiNPB[,1], c(0.025, 0.975)) # 切片 (alpha)
CI.beta.npb <- quantile(estiNPB[,2], c(0.025, 0.975)) # 傾き (beta)

# プロファイル尤度による区間を求める
PCI <- confint(fm)

# 自作MLE、Wald標準誤差・信頼区間、2種類のブートストラップの標準誤差・信頼区間をまとめる
SE.PB <- c(se.alpha.pb, se.beta.pb)
SE.NPB <- c(se.alpha.npb, se.beta.npb)
Wald.LCL <- c(MLE - 1.96 * ASE)
Wald.UCL <- c(MLE + 1.96 * ASE)
PB.LCL <- c(CI.alpha.pb[1], CI.beta.pb[1])
PB.UCL <- c(CI.alpha.pb[2], CI.beta.pb[2])
NPB.LCL <- c(CI.alpha.npb[1], CI.beta.npb[1])
NPB.UCL <- c(CI.alpha.npb[2], CI.beta.npb[2])

sumtab <- cbind(MLE, ASE, SE.PB, SE.NPB, PCI, Wald.LCL, Wald.UCL, PB.LCL, PB.UCL, NPB.LCL, NPB.UCL)
print(sumtab, 3)

#> print(sumtab, 3)
#        MLE     ASE   SE.PB  SE.NPB 2.5 % 97.5 %
#alpha 2.887 0.05266 0.05339 0.07608 2.781  2.988
#beta  0.147 0.00486 0.00492 0.00913 0.138  0.157

#      Wald.LCL Wald.UCL PB.LCL PB.UCL NPB.LCL NPB.UCL
#alpha    2.784    2.990  2.777  2.985   2.729   3.025
#beta     0.138    0.157  0.138  0.157   0.126   0.163




# 2.7 事後分布を用いたベイズ推論
# ----------------------------------------------------
# （コードなし）


#   2.7.1 不確実性または不完全な知識の一般的尺度としての確率
# --------------------------------------------------------------------------------
# （コードなし）

#   2.7.2 推論の基本的な対象としての事後分布
# ---------------------------------------------------------------------
# （コードなし）

#   2.7.3 事前分布
# -------------------------
# （コードなし）


# 2.8 Markov chain Monte Carlo（MCMC）によるベイズ計算
# -----------------------------------------------------
# （コードなし）

#   2.8.1 Monte Carlo積分
# ------------------------------------------------

# Monte Carlo積分の考え方を図示する（図2.14）
set.seed(2)
par(mfrow = c(2,2), mar = c(3,5,4,2))

# 分布から10個サンプリングする
out1 <- rnorm(10, mean = 5, sd = 2)
hist(out1, breaks = 50, col = "gray", main = "正規分布には見えない！")
mean(out1) ; sd(out1)

# 分布から100個サンプリングする
out2 <- rnorm(100, mean = 5, sd = 2)
hist(out2, breaks = 50, col = "gray", main = 'まあ、そうかもしれない…')
mean(out2) ; sd(out2)

# 分布から1000個サンプリングする
out3 <- rnorm(1000, mean = 5, sd = 2)
hist(out3, breaks = 50, col = "gray", main = 'なるほど！')
mean(out3) ; sd(out3)

# 分布から1,000,000個サンプリングする
out4 <- rnorm(10^6, mean = 5, sd = 2)
hist(out4, breaks = 100, col = "gray", main = 'おお！')
mean(out4) ; sd(out4)


#   2.8.2 Swiss bee-eaterデータに対するMarkov chain Monte Carlo
# ------------------------------------------------

# Poisson GLMの対数尤度を計算するR関数を定義する
LL <- function(param, y, x) {
  alpha <- param[1]                  # 切片
  beta <- param[2]                   # 傾き
  lambda <- exp(alpha + beta * x)    # Log-linear model
  LLi <- dpois(y, lambda, log=TRUE)  # LL contribution of each datum i
  LL <- sum(LLi)                     # LL for all observations in the data vector y
  return(LL)
}

# モデルの未正規化対数事後密度を計算するR関数を定義する
log.posterior <- function(alpha, beta, y, x){
  loglike <- LL(c(alpha, beta), y, x)
  logprior.alpha <- dnorm(alpha, mean = 0, sd = 100, log = TRUE)
  logprior.beta <- dnorm(beta, mean = 0, sd = 100, log = TRUE)
  return(loglike + logprior.alpha + logprior.beta)
}

# アルゴリズムを実行する反復回数を設定する
niter <- 10000

# 得られた事後分布サンプルを保存するRオブジェクトを作成する
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# 2つのパラメータのchainを初期化する
alpha <- rnorm(1)                    # alphaの初期値
beta <- rnorm(1)                     # betaの初期値

# 2つのパラメータについて採択の有無を保存するRオブジェクト
acc <- matrix(0, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# 初期値における現在の対数事後密度を計算する
logpost.curr <- log.posterior(alpha, beta, y, x)

# アルゴリズムのチューニングパラメータを設定する
sd.tune <- c(0.1, 0.01)              # 1つ目がalpha、2つ目がbeta


# MCMCアルゴリズムを実行する
for(t in 1:niter){
  if(t %% 1000 == 0)                 # 進捗表示
    cat("iter", t, "\n")

  # 最初に対数線形モデルの切片（alpha）を更新する
  # ––––––––––––––––––––––––––––––––––––––––––
  # alphaの候補値を提案する
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # チューニングパラメータに注意
  # 提案された新しいalphaと
  # 現在のbetaを用いてデータyと共変量xに対する対数事後密度を評価する
  logpost.cand <- log.posterior(alpha.cand, beta, y, x)
  # Metropolis採択比 r を計算する
  r <- exp(logpost.cand - logpost.curr)
  # 条件（u < r）を満たせば候補alphaを採択する
  if(runif(1) < r){
    alpha <- alpha.cand
    logpost.curr <- logpost.cand
    acc[t,1] <- 1 # 候補alphaが採択されたかを示す指標
  }

  # 次に対数線形モデルの傾き（beta）を更新する
  # ––––––––––––––––––––––––––––––––––––––
  beta.cand <- rnorm(1, beta, sd.tune[2]) # ここでもチューニングパラメータに注意

  # 提案された新しいbetaと
  # 現在のalphaを用いてデータyと共変量xに対する対数事後密度を評価する
  logpost.cand <- log.posterior(alpha, beta.cand, y, x)
  # Metropolis採択比を計算する
  r <- exp(logpost.cand - logpost.curr)
  # 条件（u < r）を満たせば候補値を採択する
  if(runif(1) < r){
    beta <- beta.cand
    logpost.curr <- logpost.cand
    acc[t,2] <- 1 # 候補betaが採択されたかを示す指標
  }
  out.mcmc[t,] <- c(alpha, beta) # 各反復のサンプルを保存
  # 注：提案値が採択されなかった場合、このステップでは
  # 以前の値をそのままコピーしてoutオブジェクトに保存する
}


# alphaとbetaそれぞれについて全体の採択率を計算する
(acc.ratio <- apply(acc, 2, mean)) # 25～40%程度が目安

# trace plotを用いてchainの収束状況を確認する
# 最初からすべてのMonte Carloサンプルを図示する（図2.15）
par(mfrow = c(2, 1))
plot(1:niter, out.mcmc[,1], main = paste('alpha（採択率 = ', round(acc.ratio[1], 2),')'),
  xlab = "MCMC反復", ylab = '事後分布からのサンプル', type = 'l')
plot(1:niter, out.mcmc[,2], main = paste('beta（採択率 = ', round(acc.ratio[2], 2),')'),
  xlab = "MCMC反復", ylab = '事後分布からのサンプル', type = 'l')


# この図に基づいてburn-inの長さを決める
nb <- 500                          # 最初のこの数のサンプルを破棄
# burn-in後の採択率を計算する
(acc.ratio.pb <- apply(acc[nb:niter,], 2, mean))

# burn-in後のサンプルだけでtrace plotを再度描く（図2.16）
par(mfrow = c(2, 1))
plot(nb:niter, out.mcmc[nb:niter,1], main = paste('alpha（burn-in後の採択率 = ',
  round(acc.ratio.pb[1], 2),')'), xlab = "反復回数", ylab = '事後分布からのサンプル', type = 'l')
plot(nb:niter, out.mcmc[nb:niter,2], main = paste('beta（burn-in後の採択率 = ',
  round(acc.ratio.pb[2], 2),')'), xlab = "反復回数", ylab = '事後分布からのサンプル', type = 'l')

# burn-in後のサンプルから事後平均、標準偏差、95%信用区間を求める
post.mn <- apply(out.mcmc[nb:niter,], 2, mean) # 事後平均
post.sd <- apply(out.mcmc[nb:niter,], 2, sd)   # 事後標準偏差
CRI <- apply(out.mcmc[nb:niter,], 2,
function(x) quantile(x, c(0.025, 0.975)))      # 95%信用区間

# burn-in後のchainから周辺事後分布の要約を作成する
tmp <- cbind(post.mn, post.sd, t(CRI))
print(tmp, 4)

# ベイズ推論と自作MLEによる推論を比較する
library(ASMbook)
get_MLE(out, 4)

# alphaとbetaの同時事後分布を図示する（図2.17）
par(mfrow = c(1, 1), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(out.mcmc[nb:niter,1], out.mcmc[nb:niter,2], xlab = 'alpha', ylab = 'beta',
  frame = FALSE, cex = 1.2, pch = 16, col = rgb(0,0,0,0.1))

# alphaとbetaの周辺事後分布を図示する（図2.18）
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
hist(out.mcmc[nb:niter,1], xlab = 'alpha', col = 'gray', main = '')
hist(out.mcmc[nb:niter,2], xlab = 'beta', col = 'gray', main = '')


#   2.8.3 Markov chain Monte Carloの実習
#     demoMCMC()関数を用いたbee-eaterデータの解析
# -----------------------------------------------
# RStudioユーザーへの注意：この節では通常のRコンソールを使う方がよい

# デフォルトより短い反復回数で実行する
# （図示用のtrue.valsには上で得たMLEを利用できる）
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(0, 0),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.1, 0.1),
  niter = 2000, nburn = 1000)) # mixingが悪く、自己相関が強い

# 採択率を上げるためチューニングパラメータを小さくする
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(0,0),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.1, 0.01),
  niter = 2000, nburn = 1000))

# ステップ幅（= tuning.params）を大きくする
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(0, 0),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.5, 0.5),
  niter = 2000, nburn = 1000)) # Manhattan状のtrace plotとなり、効率が悪い

# 1000ステップ後でも収束していない
# 極端に離れた初期値を与えて挙動を試す
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(40, 40),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.1, 0.1),
  niter = 5000, nburn = 2000) )

# さらに極端な初期値を与えて試す
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(-40, 40),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.1, 0.1),
  niter = 5000, nburn = 2000) ) # より長いburn-inが必要

# さらに試す
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(40, 0),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.1, 0.1),
  niter = 5000, nburn = 2000) ) # 特徴的な探索軌跡

# より大きなステップ幅と短いchainで同じことを試す
str(tmp <- demoMCMC(y = y, x = x, true.vals = MLE, inits = c(40, 0),
  prior.sd.alpha = 100, prior.sd.beta = 100, tuning.params = c(0.5, 0.2),
  niter = 3000, nburn = 2000) ) # 特徴的な探索軌跡


# 2.9 では、ベイズ派と頻度論派のどちらになるべきか？
# -----------------------------------------------------
# （コードなし）


# 2.10 まとめと展望
# ------------------------
# （コードなし）
