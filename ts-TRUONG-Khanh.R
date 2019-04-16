# ======================================================================= #
#                       Time Series -  December 2018
#                               Khanh TRUONG
# ======================================================================= #

# This project contains analysis on two data sets:
# SNCF Traffic and UK Stock Index FTSE

### =======================================================================
### SNCF Traffic ==========================================================
### =======================================================================
rm(list=ls())
library(tseries)
library(TSPred)
library(forecast)

# Prepare data ============================================================
Y <- scan('trafic.dat')
Y <- ts(Y,start=1963,frequency=12)
X <- window(Y,end=1979+11/12)
X.test <- window(Y,start=1980)


# Fit models ==============================================================
ts.plot(X, ylab='') # not stationary: there is trend and seasonality

diffX <- diff(X)
par(mfrow=c(3,1))
ts.plot(diffX, main='diffX') # there is still seasonality
abline(h=0, col='red', lty=2) 
acf(diffX, main='ACF') # acf pops up at lag 12
pacf(diffX, main='PACF')

diff12X <- diff(diff(X),lag=12)
par(mfrow=c(3,1))
ts.plot(diff12X, main='diff12X') # stationary
abline(h=0, col='red', lty=2) # there is still seasonality.
acf(diff12X, main='ACF') # acf vanish after lag 13 
pacf(diff12X, main='PACF') # pacf decays exponentially (less important)

# => apply MA(13) to diff12X
# => equivalently: apply SARIMA (0,1,1) (0,1,1) 12 to X
(sarima <- arima(X,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12)))


# Fit model automatically by function 'auto.arima' ========================
(sarima.auto <- auto.arima(X, trace=TRUE, ic='bic')) 
# choose smallest BIC model
summary(sarima.auto)
confint(sarima.auto, level=0.95) # CI of sar1 contains 0

t.statistic <- sarima.auto$coef['sar1']/ # estimate coefficient
  sqrt(sarima.auto$var.coef['sar1','sar1']) # standard error
pt(abs(t.statistic), df=1, lower.tail=FALSE)*2 # p-value
# insignificant different from 0, so come to the same model in sarima above


### Study residuals =======================================================
residuals <- sarima$residuals

# residuals are stationary? -----------------------------------------------
par(mfrow=c(3,1))
ts.plot(residuals, main='Residuals')
abline(h=0, col='red',lty=2) # visually stationary

# residuals have no auto-correlation? -------------------------------------
acf(residuals, main='ACF') # no visual correlation
pacf(residuals, main='PACF') # no visual partial correlation
Box.test(residuals, lag=12, type='Ljung-Box') 
# H0: no serial correlation upto n lags
# cannot reject the null hypothesis

# volatility of residuals is stable? --------------------------------------
Box.test(residuals^2, lag=12, type='Ljung-Box')
# H0: no arch effect
# cannot reject the null hypothesis

# distribution of residuals? ----------------------------------------------
par(mfrow=c(1,1))
hist(residuals, prob=TRUE, main=NULL, xlab=NULL, col='grey')
curve(dnorm(x, mean(residuals), sd(residuals)), lty=2, add=TRUE) # gaussian
jarque.bera.test(residuals)
# H0: residuals follow Gaussian distribution
# reject the null hypothesis


### Forecast ==============================================================
(forecast <- forecast(sarima, h=12))
# will have error if import library(TSA) before,
# if so, restart R and run from the beginning.

plot(forecast)
abline(v=c(1977,1981), lty=c(2,2)) # line to zoom-in
plot(forecast, xlim=c(1977,1981), ylim=c(2500,4500)) # zoom-in

par(mfrow=c(1,1))
plotarimapred(X.test, sarima, xlim=c(1977,1981), range.percent=0.05)
legend("topleft", legend=c("Actual","Forecast"), col=c("black","blue"),
       lwd=2, cex=0.7, lty=c(2,1)) # legend
abline(v=c(1979,1981), lty=c(2,2)) # line to zoom-in

par(mfrow=c(1,1))
plotarimapred(X.test, sarima, xlim=c(1979,1981), range.percent=0.05)
legend("topleft", legend=c("Actual","Forecast"), col=c("black","blue"),
       lwd=2, cex=0.7, lty=c(2,1)) # legend

accuracy(forecast, X.test) # some accuracy metrics


### Compare to method 'decompose' =========================================
plot(decompose(X))
# plot(decompose(X, type='multiplicative')) # method multiplicative
names(decompose(X))
random <- decompose(X)$random # residuals of decompose method

# residuals are stationary? -----------------------------------------------
par(mfrow=c(2,1))
ts.plot(random, main='Decompose')
abline(h=0, col='red',lty=2) # visually stationary
ts.plot(residuals, main='SARIMA')
abline(h=0, col='red',lty=2) # visually stationary

# residuals have no auto-correlation? -------------------------------------
random <- na.omit(random)

par(mfrow=c(2,1))
acf(random, main='Decompose') # decay exponentially (with seasonality)
acf(residuals, main='SARIMA') # no visual correlation

par(mfrow=c(2,1))
pacf(random, main='Decompose') # partial correlated upto lag 12
pacf(residuals, main='SARIMA') # no visual correlation

Box.test(random, lag=12, type='Ljung-Box') 
# H0: no serial correlation upto n lags
# reject the null hypothesis

# volatility of residuals is stable? --------------------------------------
Box.test(random^2, lag=12, type='Ljung-Box')
# H0: no arch effect
# reject the null hypothesis

# distribution of residuals? ----------------------------------------------
par(mfrow=c(1,1))
hist(random, prob=TRUE, main=NULL, xlab='Residuals', col='grey',
     ylim=c(0,0.003))
curve(dnorm(x, mean(random), sd(random)), lty=2, add=TRUE) # gaussian
jarque.bera.test(random)
# H0: residuals follow Gaussian distribution
# reject the null hypothesis

# good prediction? --------------------------------------------------------
T <- time(X)
T2 <- T^2  
reg <- lm(X~T+T2) # regression traffic to T and T2
summary(reg)

T.test <- data.frame(T=time(X.test), T2=time(X.test)^2) # new data

trend.predict <- predict(reg,T.test) # predict trend based on regression
season.predict <- window(decompose(X)$seasonal,start=1979) # seasonality
X.predict <- trend.predict + season.predict # prediction based on decompose

# actual values
plot(X.test, ylim=c(2500,4500), type='l', lty=2, lwd=2, col='black',
     xlab=NULL, ylab='Traffic')
# decompose prediction
par(new=TRUE)
plot(X.predict, ylim=c(2500,4500), type='l', lty=1, lwd=2, col='red',
     ylab=NULL, xlab=NULL, xaxt='n', yaxt='n') 
# SARIMA prediction
par(new=TRUE)
plot(forecast$mean, ylim=c(2500,4500), type='l', lty=1, lwd=2, col='blue',
     ylab=NULL, xlab=NULL, xaxt='n', yaxt='n')
# legend
legend("topleft",legend=c("Actual","Decompose Forecast","SARIMA Forecast"),
       col=c("black","red","blue"), lwd=2, cex=0.7, lty=c(2,1,1))
# => 'decompose' is not that bad, especially at first and last quarters



### =======================================================================
### ARCH and GARCH ========================================================
### =======================================================================
rm(list=ls())
library(TSA)
library(rugarch)

# Simulate ARCH process ===================================================
set.seed(1234)
arch <- garch.sim(n=1000,alpha = c(0.05,0.3))

# process is stationary? --------------------------------------------------
par(mfrow=c(3,1))
ts.plot(arch, main='ARCH(1)')
abline(h=0, col='red', lty=2) # visually stationary(!)

# process have no auto-correlation? ---------------------------------------
acf(arch, main='ACF') # no correlation
pacf(arch, main='PACF') # no partial correlation
Box.test(arch, lag=20, type='Ljung-Box') 
# H0: no serial correlation upto n lags
# cannot reject the null hypothesis

# squared series follows AR(1)? -------------------------------------------
par(mfrow=c(1,1))
pacf(arch^2, main='') # partial correlation at lag 1

# volatility of process is stable? ----------------------------------------
Box.test(arch^2, lag=20, type='Ljung-Box') 
# H0: no arch effect
# reject the null hypothesis

# distribution of process? ------------------------------------------------
par(mfrow=c(1,1))
hist(arch, prob=TRUE, main=NULL, xlab='process', col='grey')
curve(dnorm(x, mean(arch), sd(arch)), lty=2, add=TRUE) # gaussian
jarque.bera.test(arch)
# H0: process follow Gaussian distribution
# reject the null hypothesis


# Simulate GARCH process ==================================================
set.seed(1234)
garch <- garch.sim(n=1000, alpha = c(0.05,0.3), beta=c(0.1,0.15))

# process is stationary? --------------------------------------------------
par(mfrow=c(3,1))
ts.plot(garch, main='GARCH(1,2)')
abline(h=0, col='red', lty=2) # visually stationary(!)

# process have no auto-correlation? ---------------------------------------
acf(garch, main='ACF') # some appoach to the bounds of CI 
pacf(garch, main='PACF') # some appoach to the bounds of CI
Box.test(garch, lag=20, type='Ljung-Box') 
# H0: no serial correlation upto n lags
# cannot reject the null hypothesis

# squared series follows ARMA? --------------------------------------------
par(mfrow=c(2,1))
acf(garch^2, main='') # no partial correlation
pacf(garch^2, main='') # no partial correlation

# volatility of process is stable? ----------------------------------------
Box.test(garch^2, lag=20, type='Ljung-Box') 
# H0: no garch effect
# reject the null hypothesis

# distribution of process? ------------------------------------------------
par(mfrow=c(1,1))
hist(garch, prob=TRUE, main=NULL, xlab='process', col='grey')
curve(dnorm(x, mean(garch), sd(garch)), lty=2, add=TRUE) # gaussian
jarque.bera.test(garch)
# H0: process follow Gaussian distribution
# reject the null hypothesis


# Validate SARIMA model ===================================================
(sarima <- arima(Y,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12)))



### =======================================================================
### UK Stock Index FTSE ===================================================
### =======================================================================
FTSE <- EuStockMarkets[,'FTSE']
tsp(FTSE)
FTSE.log.return <- diff(log(FTSE))*100

par(mfrow=c(2,1))
ts.plot(FTSE, ylab='', main='FTSE Index') # there is trend
ts.plot(FTSE.log.return, ylab='', main='FTSE Log Return') # no trend

# test <- window(FTSE.log.return,start=1998+0.61) # last 10 observations
# train <- window(FTSE.log.return,end=1998+0.61) # remaining

test <- tail(FTSE.log.return, 10) # last 10 observations
train <- FTSE.log.return[!(FTSE.log.return %in% test)] # remaining

par(mfrow=c(3,1))
ts.plot(train, main='FTSE') 
acf(train, main='ACF') # there is still seasonality (5 lags)
pacf(train, main='PACF')

# diff.train <- diff(train) # differentiate order 1
# par(mfrow=c(3,1))
# ts.plot(diff.train, main='FTSE') 
# acf(diff.train, main='ACF') # there is still seasonality (5 lags)
# pacf(diff.train, main='PACF')
# 
# diff.train <- diff(train, lag=5)
# par(mfrow=c(3,1))
# ts.plot(diff.train, main='FTSE')
# acf(diff.train, main='ACF')
# pacf(diff.train, main='PACF')


# Fit models ==============================================================
(arma <- auto.arima(train, ic='bic'))
Box.test(arma$residuals^2, lag=10, type='Ljung-Box')
# H0: no arch effect
# reject the null hypothesis => use GARCH model

ug.spec <- ugarchspec(mean.model=list(armaOrder=c(0,1))) # specify garch
ug.fit <- ugarchfit(ug.spec, train) # fit the model garch

par(mfrow=c(1,1))
plot(ug.fit, which=3) # wrong x-axis, should be from 1992 to 1998
legend("topright", legend=c("Conditional SD","Absolute Return"),
       col=c("blue","grey"), lwd=2, cex=0.8, lty=c(1,1))

names(ug.fit@fit)
ug.fit@fit$coef
ug.res2 <- ug.fit@fit$residuals^2 # squared residuals
ug.var <- ug.fit@fit$var # conditional variance

par(mfrow=c(1,1))
plot(ug.res2, type='l', ylab='', xlab='Days', ylim=c(0,10))
lines(ug.var, col='green')
legend("topright", legend=c("Squared residual","Conditional variance"),
       col=c("black","green"), lwd=2, cex=0.8, lty=c(1,1))
# typical garch type pattern:
# periods with large volatility (green line),
# (squared) residuals tends to larger (black line)


# Forecast ================================================================
(ug.fore <- ugarchforecast(ug.fit, n.ahead=10)) # forecasted

# Value Forecast ----------------------------------------------------------
# actual values
plot(test, type='l', lty=1, lwd=2, col='black',
     xlab=NULL, ylab=NULL, ylim=c(-3,3), xaxt='n', yaxt='n')
# GARCH prediction
par(new=TRUE)
plot(ug.fore@forecast$seriesFor, type='l', lty=2, lwd=2, col='blue',
     ylab='FTSE', xlab='Days in Future', ylim=c(-3,3)) 
# legend
legend("topleft",legend=c("Actual","GARCH Prediction"),
       col=c("black","blue"), lwd=2, cex=0.7, lty=c(1,2))

# Volatility Forecast -----------------------------------------------------
ug.f <- ug.fore@forecast$sigmaFor # standard error of forecasted values
plot(ug.f, type = "l", ylab='',
     xlab='Days in Future') # standard error decreasing

# gets the last 20 observations of squared residuals
ug.res2.t <- c(tail(ug.res2,20),rep(NA,10))
# gets the last 20 observations of conditional variance
ug.var.t <- c(tail(ug.var,20),rep(NA,10))  
# variance of forecasted values
ug.f <- c(rep(NA,20),(ug.f)^2)

plot(ug.res2.t, type="l", xlab='Days in Future', ylab='')
lines(ug.var.t, col="blue", lwd=2)
lines(ug.f, col="red", lty=2, lwd=2)
legend("topright", legend=c("Squared residuals","Conditional variance",
  "Forecast variance"), col=c("black","blue","red"), lwd=2, cex=0.8,
  lty=c(1,1,2))

# ======================================================================= #