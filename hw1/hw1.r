'''
1.17 (Brockwell et al., 2016, p. 36)

Load the dataset deaths in R using the read.table function. Plot the data. Also create a histogram of
the data using the R function hist. Plot the sample autocorrelation function using the acf function. The
presence of a strong seasonal component with period 12 is evident in the graph of the data and in the sample
autocorrelation function.
'''

data = read.table("C:/prosjekt/Time_Series_stat211/data/deaths.txt", skip = 9, header=FALSE)
names(data) = c("month","year","deaths")
data$time = as.Date(paste(data$year,data$month,"15",sep="-"))
attach(data)
plot(time,deaths/1000, type="l", xlab="", ylab = "deaths (thousands)")

acf(deaths ,lag.max = 24) #computes (and by default plots) estimates of the autocovariance or autocorrelation function

hist(deaths, breaks = 15)

'''
Problem 1.6
[BD, Exercise 1.18, page 37]
We are still studying the dataset deaths. In this exercise, you are supposed to reproduce the figures 1-24 and 1-25 in [BD, pp. 27-28] . 
In 1.17, we found a period of length 12. 
Fit a seasonal component using the procedure described in section 1.5.2.1 on page 26.
You may use the following functions or write your own:
'''
# Function for calculating a moving average when d is even
ma <- function(x, n=12){filter(x,c(.5,rep(1,n-1),.5)/n, sides=2)}
# Function for finding the seasonal component
seasonal.component <- function(x){ # x-  deaths
	# First step: detrending
	detrended <- x - ma(x)
	# Second step: Calculating sesonal component from detrended data
	wt <- rowMeans(matrix(detrended[!is.na(detrended)], nrow=12,byrow=FALSE))
	st<-(wt - mean(wt))[c(7:12,1:6)] #seasonal component
	return(st)
}


#sol
library(itsmr)
detrended <- deaths - seasonal.component(deaths) #season(deaths, 12) #removing season

'''
Plot the deseasonalized data (as in figure 1-24). Fit a quadratic trend (polynomial of
order two) to the deseasonalized data and add the curve to the plot you just created. The
trend should be mb = 9952 − 71.82t + 0.8260t2 for 1 ≤ t ≤ 72. This can be done using the
following code:
'''

M <- poly(1:72, degree=2, raw=TRUE)
trend <- lm(detrended ~ M) # Re-estimating trend of the detrended data
'''
Plot the sample autocorrelation function of Yt. Forecast the data for the next 24 months
without allowing for this dependence, based on the assumption that the estimated seasonal
and trend components are true values and that Yt is a white noise sequence with zero
mean. Calculate sb72+k for k = 1, . . . , 24 and do the forecasting by
Xbt = mb 72+k + sb72+k, k = 1, . . . , 24.
'''

#solution
plot(data$time, detrended/1000, type="b", xlab="",ylab="(thousands)") #plotting deseasonalized data

lines(data$time, trend(detrended, 2)/1000, col=2)

plot(data$time,season(deaths,12)/1000, type="b",xlab="",ylab="(thousands)") # Plotting seasonal component
abline(h=0) # adding horizontal line at zero

##########
'''
Plot the original data with the forecasts appended. Later we shall see how to improve on
these forecasts by taking into account the dependence in the series Yt. Hint: To calculate
mb 72+k the following code may be useful:
'''
M <- poly(72 + 1:24, 2, raw=TRUE)
yhat <- predict(trend, newdata= M)

#######
p.time <- seq(as.Date("1979-01-15"),by="month",length.out=24)

#plot(time,deaths/1000,xlim=range(time,p.time),type="b",
plot(time,detrended/1000,xlim=range(time,p.time),type="b",
	#ylim = range(deaths, yhat)/1000,
	ylim = range(trend(detrended, 2), yhat)/1000,
	ylab = "(thousands)", xlab ="")

lines(data$time, trend(detrended, 2)/1000, col=2)	
lines(p.time,yhat/1000, type="b", col = "blue")





'''
References
Peter J Brockwell and Richard A Davis. Introduction to time series and forecasting;
3rd ed. Springer texts in statistics. Springer, Cham, 2016. doi: 10.1007/
'''
