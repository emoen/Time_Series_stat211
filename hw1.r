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