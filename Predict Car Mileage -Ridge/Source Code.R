## Part 1a
require(reshape2)
require(ggplot2)
require()
X_train <- read.csv("X_train.csv", header = FALSE)
Y_train <- read.csv("y_train.csv", header = FALSE)
X <- as.matrix(X_train)
y <- as.matrix(Y_train)
## taking the SVD
SVD <- svd(X) 
## setting values of lambda
lambda = seq(0, 5000, 1)
## estimating w_rr
w_rr = matrix(0, length(lambda), ncol(X))
for(i in 1:length(lambda))
{
  w_rr[i, ] = solve(t(X)%*%X + lambda[i]*diag(ncol(X)))%*%t(X)%*%y
}
df_ridge <- data.frame(cbind(lambda, w_rr))
colnames(df_ridge) <- c("Lambda", "w1", "w2", "w3", "w4", "w5", "w6", "w7")
## estimating df_lambda
df_lambda = matrix(0, length(lambda), 1)
for(i in 1:length(lambda))
{
  df_lambda[i, ] = sum(diag(solve(lambda[i]*(1/(SVD$d)^2)*diag(ncol(X)) + diag(ncol(X)))))
}
## combining lambda, df_lambda, and w_rr together into one dataframe
df_ridge2 <- data.frame(cbind(lambda, df_lambda, w_rr))
colnames(df_ridge2) <- c("Lambda", "df_lambda", "w1", "w2", "w3", "w4", "w5", "w6", "w7")
df2 <- melt(df_ridge2, id=c("Lambda","df_lambda"))
colnames(df2) <- c("Lambda", "df_lambda", "Variable", "w_rr")
## plotting the graph
ggplot(data=df2, aes(x=df_lambda, y=w_rr, colour=Variable)) + geom_line()

## Part 1c
X_test <- read.csv("X_test.csv", header = FALSE)
y_test <- read.csv("y_test.csv", header = FALSE)
X_test <- as.matrix(X_test)
y_test <- as.matrix(y_test)
w_rr <- as.matrix(w_rr)
## setting values for lambda
lambda_c = c(0:50)
## getting predictions for y
predicted_y = X_test%*%t(w_rr[c(0:51), ])
colnames(predicted_y) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                           "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                           "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43",
                           "44", "45", "46", "47", "48", "49", "50")
## finding RMSE for all values of lambda in test data
RMSE = matrix(0, 51, 1)
for(i in 1:51)
{
  RMSE[i, ] = sqrt(mean((y_test - (X_test%*%w_rr[i, ])) ^ 2))
}
RMSE = data.frame(cbind(RMSE, lambda_c))
## plotting the graph
gg <- ggplot(data=RMSE, aes(x=lambda_c, y = V1)) + geom_line() 
gg + labs(x = "lambda", y = "RMSE")

## Part 1d
## setting the values for lambda
lambda_d = c(0:500)
## RMSE for p = 1
RMSE = matrix(0, 501, 1)
for(i in 1:501)
{
  RMSE[i, ] = sqrt(mean((y_test - (X_test%*%w_rr[i, ])) ^ 2))
}
## finding RMSE for p = 2
X_sq = X^2
X_sq <- cbind(X, X_sq)
X_sq <- X_sq[ , -7]
w_rrsq = matrix(0, length(lambda_d), ncol(X_sq))
for(i in 1:length(lambda_d))
{
  w_rrsq[i, ] = solve(t(X_sq)%*%X_sq + lambda[i]*diag(ncol(X_sq)))%*%t(X_sq)%*%y
}
X_testsq = X_test^2
X_testsq <- cbind(X_test, X_testsq)
X_testsq <- X_testsq[ , -7]
RMSE_sq = matrix(0, 501, 1)
for(i in 1:501)
{
  RMSE_sq[i, ] = sqrt(mean((y_test - (X_testsq%*%w_rrsq[i, ])) ^ 2))
}
## finding RMSE for p = 3
X_cube = X^3
X_cube <- cbind(X_sq, X_cube)
X_cube <- X_cube[ , -13]
w_rrcube = matrix(0, length(lambda_d), ncol(X_cube))
for(i in 1:length(lambda_d))
{
  w_rrcube[i, ] = solve(t(X_cube)%*%X_cube + lambda[i]*diag(ncol(X_cube)))%*%t(X_cube)%*%y
}
X_testcube = X_test^3
X_testcube <- cbind(X_testsq, X_testcube)
X_testcube <- X_testcube[ , -13]
RMSE_cube = matrix(0, 501, 1)
for(i in 1:501)
{
  RMSE_cube[i, ] = sqrt(mean((y_test - (X_testcube%*%w_rrcube[i, ])) ^ 2))
}
## plotting RMSE for p = 1, 2, 3 against lambda
RMSE_all = data.frame(cbind(RMSE, RMSE_sq, RMSE_cube, lambda_d))
gg <- ggplot(data=RMSE_all, aes(x=lambda_d, y = value, color = variable)) +
  geom_line(aes(y = V1, col = "p = 1")) +
  geom_line(aes(y = V2, col = "p = 2")) +
  geom_line(aes(y = V3, col = "p = 3"))
gg + labs(x = "lambda", y = "RMSE")
gg
## finding minimum RMSE for all values of p and lambda
which(RMSE_all[,1:3] == min(RMSE_all[, 1:3]), arr.ind = TRUE)