rm(list=ls())

library(ggplot2)
library(forecast)

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)

load("SeaIce.Rdata")

##################
# Exercise for all
##################

gp = autoplot(y.ts) + ggtitle("Baffin Bay Sea Ice Extent") + xlab("Year") + ylab("Sea Ice Extent (km^2)")
gp

# Exploratory analysis

seaice.monthly = aggregate(c(y.ts), list(month = cycle(y.ts)), mean) 
seaice.annual = aggregate(c(y.ts), list(year = floor(time(y.ts))), mean) 

seaice.monthly = ts(seaice.monthly[,2], start = 1)
autoplot(seaice.monthly) + ggtitle("Baffin Bay Sea Ice Extent Intra-annual trend") + xlab("Month") + ylab("Sea Ice Extent(km^2)")

seaice.annual = ts(seaice.annual[,2], start = seaice.annual[2,1], end = seaice.annual[42,1])
autoplot(seaice.annual) + ggtitle("Baffin Bay Sea Ice Extent Inter-annual trend") + ylab("Sea Ice Extent (km^2)") + xlab("Year")

# Detrending seasonality and trend

year.ts = floor(time(y.ts))

X1 = fourier(y.ts, K = 1)
X2 = fourier(y.ts, K = 2)
X3 = fourier(y.ts, K = 3)

mod.h1 = tslm(y.ts ~ X1 + year.ts)
summary(mod.h1)
mod.h2 = tslm(y.ts ~ X2 + year.ts)
summary(mod.h2)
mod.h3 = tslm(y.ts ~ X3 + year.ts)
summary(mod.h3)

# Plot the fitted values

df = data.frame(X = time(y.ts), Y = mod.h2$fitted.values)
gp + geom_line(data = df, aes(X,Y), color = "red", linetype = "dashed") + ggtitle("Baffin Bay Sea Ice Extent Fitted Values (Seasonal Harmonic and Linear Trend)")

# Store and plot the residuals

res.ts = y.ts - mod.h2$fitted.values
autoplot(res.ts) + ggtitle("Residuals from Seasonal and Linear Trend") + xlab("Year") + ylab("Residual")

# Fit a non-seasonal ARMA model to the residuals
# (the information left over after detrending)

seaice.res.fit = auto.arima(res.ts, d = 0, seasonal = FALSE, max.p = 4, max.q = 4, stepwise = FALSE)
summary(seaice.res.fit)
aic = NULL
for (p in 0:3)
{
  for (q in 0:3)
  {
    dummy = arima(res.ts, order = c(p, 0, q), include.mean = FALSE)
    aic = rbind(aic,cbind(p, q, dummy$aic))
  }
}
aic

ggAcf(seaice.res.fit$residuals) + ggtitle("ACF of ARMA Model Residuals")

# Forecast for 2020-2021

horizon = 22
forc = forecast(seaice.res.fit, h = horizon, level = c(80, 95))

# Add back the trend

year.ts = c(kronecker(2019,rep(1,10)), kronecker(2020,rep(1,12)))
trend.predict = forecast(mod.h2, data.frame(fourier(y.ts, K = 2, h = horizon), year.ts))

forc$x = forc$x + mod.h2$fitted.values
forc$mean = forc$mean + trend.predict$mean
forc$lower = forc$lower + trend.predict$mean
forc$upper = forc$upper + trend.predict$mean

# Plot final results
autoplot(forc) + xlab("Year") + ylab("Sea Ice Extent (km^2)") + ggtitle("Forecast of Sea Ice Extent")
autoplot(forc, include = 100) + xlab("Year") + ylab("Sea Ice Extent (km^2)") + ggtitle("Forecast of Sea Ice Extent")
autoplot(forc, include = 24) + xlab("Year") + ylab("Sea Ice Extent (km^2)") + ggtitle("Forecast of Sea Ice Extent")

################################
# Exercise for graduate students
################################

set.seed(1)

len = 1000
nsim = 1000

results = NULL

# Test 1 - ARMA(1,0) where phi = 0.8

successes = 0
for (n in 1:nsim)
{
  test = arima.sim(model = list(ar = 0.8), n = len)
  fit = arima(test[1:(nsim-1)], order = c(1,0,0), include.mean = FALSE)
  check = forecast(fit, h = 1, level = 95)
  check.low = check$lower[1]
  check.high = check$upper[1]
  if (check.low <= test[len] & test[len] <= check.high)
    successes = successes + 1
}
results = rbind(results, cbind(1, 0, successes/nsim))

# Test 2 - ARMA(0,1) where theta = 0.8

successes = 0
for (n in 1:nsim)
{
  test = arima.sim(model = list(ma = 0.8), n = len)
  fit = arima(test[1:(nsim-1)], order = c(0,0,1), include.mean = FALSE)
  check = forecast(fit, h = 1, level = 95)
  check.low = check$lower[1]
  check.high = check$upper[1]
  if (check.low <= test[len] & test[len] <= check.high)
    successes = successes + 1
}
results = rbind(results, cbind(0, 1, successes/nsim))

# Test 3 - ARMA(1,1) where phi = 0.6 and theta = -0.2

successes = 0
for (n in 1:nsim)
{
  test = arima.sim(model = list(ar = 0.6, ma = -0.2), n = len)
  fit = arima(test[1:(nsim-1)], order = c(1,0,1), include.mean = FALSE)
  check = forecast(fit, h = 1, level = 95)
  check.low = check$lower[1]
  check.high = check$upper[1]
  if (check.low <= test[len] & test[len] <= check.high)
    successes = successes + 1
}
results = rbind(results, cbind(1, 1, successes/nsim))

# Test 4 - ARMA(2,2) where phi1 = -0.4, phi2 = -0.6, theta1 = 0.8, and theta2 = -0.2

successes = 0
for (n in 1:nsim)
{
  test = arima.sim(model = list(ar = c(-0.4, -0.6), ma = c(0.8, -0.2)), n = len)
  fit = arima(test[1:(nsim-1)], order = c(2,0,2), include.mean = FALSE)
  check = forecast(fit, h = 1, level = 95)
  check.low = check$lower[1]
  check.high = check$upper[1]
  if (check.low <= test[len] & test[len] <= check.high)
    successes = successes + 1
}
results = rbind(results, cbind(2, 2, successes/nsim))

results
























