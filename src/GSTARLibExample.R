library(gstar)
library(xts)

data("LocationCPI")

#-----Use data with xts object----#

x = xts(LocationCPI[, -1], order.by = as.Date(LocationCPI[, 1]))
s <- round(nrow(x) * 0.8) 

## split into training and testing (80:20)
x_train <- x[1:s, ]
x_test <- x[-c(1:s), ]

weight = matrix(c(0, 1, 1, 1, # create the uniform weight.
                  1, 0, 1, 1,
                  1, 1, 0, 1,
                  1, 1, 1, 0), ncol = 4, nrow = 4)

weight = weight/(ncol(x) - 1) #the sum of weight is equal to 1 every row.
fit <- gstar(x_train, weight = weight, p = 1, d = 0, est = "OLS")

summary(fit)

performance(fit)
performance(fit, x_test) ## to check the performance with testing data
predict(fit, n = 10) #forecast 10 data ahead

plot(fit)
plot(fit, n_predict = 10) #plot with 10 forecasting data
plot(fit, testing = x_test)

#---- Use dataframe or matrix---#
x2 <- LocationCPI
x2$Date <- NULL # remove the date column
data(Loc)

dst <- as.matrix(dist(Loc[, -1], diag = TRUE, upper = TRUE))
dst1 <- matrix(0, nrow = nrow(dst), ncol = ncol(dst))

for(i in 1:nrow(dst)) {
  for(j in 1:ncol(dst)){
    if(j == i) next
    dst1[i, j] <- sum(dst[i, -j])/sum(dst[i,])
  }
}

weight_inverse_distance <- matrix(0, nrow = nrow(dst), ncol = ncol(dst))

for(i in 1:nrow(dst)) {
  for(j in 1:ncol(dst)){
    if(j == i) next
    weight_inverse_distance[i, j] <- sum(dst1[i, j])/sum(dst1[i,])
  }
}

fit_inverse_distance <- gstar(x2, weight = weight_inverse_distance, p = 2, d = 1, est = "OLS")

summary(fit_inverse_distance)
performance(fit_inverse_distance)
predict(fit_inverse_distance)
plot(fit_inverse_distance)