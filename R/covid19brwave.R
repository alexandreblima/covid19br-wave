# Times Series Analysis
# Total number of deaths per month - Brazil 
# Jan/2015 - jul/2020

rm(list = ls()) # remove todos os objetos 

library(readxl)
library(ggplot2)
library(gvlma)
library(tseries)
library(UsingR)
library(waveslim)
library(longmemo)
library(segmented)
library(changepoint)
library(rmatio)
library(splus2R)
library(rio)

dados <- read_xls("obitosBRjan2015jul2020.xls", sheet = "serie1", col_names = TRUE, col_types=c("numeric", "text", "text"))
dados1 <- read_xls("obitosBRjan2015jul2020.xls", sheet = "serie2", col_names = TRUE, col_types=c("numeric", "text", "text"))
x <- dados[1] # 2015/jan to 2020/jul
y <- dados1[1] # 2015/mar to 2020/jul
x1 <- unlist(x, use.names = FALSE)
y1 <- unlist(y, use.names = FALSE)
N = length(x1) # 67 samples
M = length(y1) # 65 samples
t <- 1:N
t1 <- 1:(M-1)
fit <- lm(x1 ~ t) # linear regression
residuos <- fit$residuals
results <-summary(fit)
print(results)

# quartz() # open graphics in MacOs, if you will
# Plot time series and sample ACF
# par(mfrow=c(2,1))
plot(t, x1, type = "o", xlab = "time: 2015/jan to 2020/jul", ylab = "deaths per month") 
abline(fit)
acfx1 <- acf(x,100)   
plot(acfx1, main = "")

# quartz()
# Analysis of the residuals
plot(fit) 

# plot residuals time series
plot(t, residuos, type = "o", xlab = "time: 2015/jan to 2020/jul", ylab = "residuals") 
acf(residuos,100)    
# QQ-plot
qqnorm(residuos)

# Tests the null of normality for residuals using the Jarque-Bera test 
test1 <- jarque.bera.test(residuos)
print(test1)

# Tests the null of normality for ydiff using the Shapiro-Wilk test 
test2 <- shapiro.test(residuos)
print(test2)

# detrend time series y1: take first difference of time series from 2015/mar to 2020/jul
ydiff <- diff(y1) # 64 samples (= 2^6 samples)
mu <- mean(ydiff)
ydiff <- ydiff - mu
ydiff1 <- as.data.frame(ydiff)
export(ydiff1, "ydiff.mat")

# quartz()
# Plot detrended time series and SACF
#par(mfrow=c(2,1))
plot(t1, ydiff, type = "o", xlab = "time: 2015/apr to 2020/jul", ylab = "Detrended times series") 
acf_ydiff <- acf(ydiff,100)  
plot(acf_ydiff, main = "")

# quartz()
# QQ-plot
qqnorm(ydiff)

# Histogram and normal PDF
ydiff.hist <- hist(ydiff, nclass = 10, freq = FALSE, main = "Histogram")
x3 <- seq(min(ydiff), max(ydiff), length.out=100)
y3 <- with(galton, dnorm(x3, mean(ydiff), sd(ydiff)))
lines(x3, y3, col = "red", main = "")

# Tests the null of normality for ydiff using the Jarque-Bera test 
test3 <- jarque.bera.test(ydiff)
print(test3)

# Tests the null of normality for ydiff using the Shapiro-Wilks test 
test4 <- shapiro.test(ydiff)
print(test4)

# Statistical model selection: AIC criterion (Akaike)
ydiff.ar = ar(ydiff, order.max = 30) # obtain AICs

# quartz() 
# plot AICs
plot(0:30, ydiff.ar$aic, type="l")

ydiff.ar # selected model: AR(11)

# Spectral Analysis

#quartz() 
# Plot periodogram (Daniell Method) and the estimated PSD of the AR(11) model
spec.pgram(ydiff,log="no",demean=F, spans=c(2,2), taper = 0.1, main = "Periodogram - Daniell Method") 
# Plot estimated PSD of the AR(11) model
spec.ar(ydiff,log="no", main = "AR(11) model estimated PSD")

# Stationarity testing
# Computes the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test 
# for the null hypothesis that the time series is level stationary
# Is the series I(0)? 
result1.kpss <- kpss.test(x1, null="Trend")
print(result1.kpss)
result2.kpss <- kpss.test(ydiff, null="Level")
print(result2.kpss)
# The smaller the p-value, the stronger the evidence that you should reject the null hypothesis.

# quartz() # open graphics
# Identifying heteroskedasticity in the residuals of the linear model
plot(fit) # identifying heteroskedasticity graphically

# Assessment of the linear models assumptions (level of significance = 0.05)
gvlma(fit) 

####################### Wavelet Analysis ############################################

# Testing for a change in the slope
cp.davies <- davies.test(fit,~t)
print(cp.davies)
#testing for one changepoint
#use the simple null fit
pscore.test(fit,~t) #compare with davies.test(o,~z)..

# wavelet-based approach for detecting changes in second order structure of time series
# pacote "changepoint" do R de Rebecca Killick

cp.ydiff = cpt.var(ydiff,penalty="SIC",method="AMOC",class=TRUE)
print(cp.ydiff)

# MRA: use Haar or db1
## Haar MRA
ydiff.haar <- mra(ydiff, "haar", 4, "dwt")
names(ydiff.haar) <- c("D1", "D2", "D3", "D4", "S4")

## MRA using Daubechies la8 wavelet -- least asymmetric wavelet
## The scaling and wavelet filters have an approximated linear phase
## with permit us to better align the wavelet and details with the signal
ydiff.haar <- mra(ydiff, "la8", 4, "dwt")
names(ydiff.haar) <- c("D1", "D2", "D3", "D4", "S4")

## plot MRA of ydiff
par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
plot(t1, ydiff, type = "o", xlab = "2015/jan to 2020/jul", ylab = "first difference")
for(i in 1:5)
  plot.ts(ydiff.haar[[i]], axes=FALSE, ylab=names(ydiff.haar)[i])
axis(side=1, at=seq(0,368,by=23), 
     labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))

wf <- "d4"
J <- 5
demusd.modwt <- modwt(ydiff, wf, n.levels = J)
demusd.modwt.bw <- brick.wall(demusd.modwt, wf)
jpyusd.modwt <- modwt(ydiff, wf, J)
jpyusd.modwt.bw <- brick.wall(jpyusd.modwt, wf)
returns.modwt.cov <- wave.covariance(demusd.modwt.bw, jpyusd.modwt.bw)
par(mfrow=c(1,1), las=0, mar=c(5,4,4,2)+.1)
matplot(2^(0:(J-1)), returns.modwt.cov[-(J+1),], type="b", log="x",
        pch="*LU", xaxt="n", lty=1, col=c(1,4,4), xlab="Wavelet Scale", 
        ylab="Wavelet Covariance")
axis(side=1, at=2^(0:7))
abline(h=0)

graphics.off() 