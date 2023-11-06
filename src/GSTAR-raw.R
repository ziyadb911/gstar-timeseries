library(gstar)
library(xts)

#-----import dataset----#
df <- read.csv("dataset.csv", header=TRUE)
df <- df[,c("DateTime","Zone1","Zone2","Zone3")]

#-----adjust date time format----#
df$DateTime <- as.POSIXct(df$DateTime,format='%m/%d/%Y %H:%M')
df$Date <- strftime(df$DateTime,format='%Y-%m-%d')
df$Time <- strftime(df$DateTime,format='%H:%M')

df <- df[df$Time == "00:00",c("DateTime","Zone1","Zone2","Zone3")]
df


#-----Use data with xts object----#
x = xts(df[, -1], order.by = df$DateTime)
s <- round(nrow(x) * 0.8) 

## split into training and testing (80:20)
x_train <- x[1:s, ]
x_test <- x[-c(1:s), ]

weight = matrix(c(0, 1, 1, # create the uniform weight.
                  1, 0, 1,
                  1, 1, 0), ncol = 3, nrow = 3)

weight = weight/(ncol(x) - 1) #the sum of weight is equal to 1 every row.
fit <- gstar(x_train, weight = weight, p = 1, d = 0, est = "OLS")

summary(fit)

performance(fit)
performance(fit, x_test) ## to check the performance with testing data
predict(fit, n = 10) #forecast 10 data ahead

plot(fit)
plot(fit, n_predict = 10) #plot with 10 forecasting data
plot(fit, testing = x_test)
