# Part 2a
X_train <- read.csv("X_train.csv", header = FALSE)
Y_train <- read.csv("y_train.csv", header = FALSE)
r_num <- c(1:4508)
X_train2 <- cbind(r_num, X_train)
Y_train2 <- cbind(r_num, Y_train)
XY_train <- merge(X_train2, Y_train2, by = "r_num")
# Subsetting to Y = 1 and Y = 0 to obtain class conditional parameters below
XY_train_Y1 <- XY_train[XY_train$V1.y == 1, ]
XY_train_Y0 <- XY_train[XY_train$V1.y == 0, ]
# Finding pi
pi <- sum(XY_train$V1.y)/nrow(XY_train)
# Finding theta for Y = 1, using Bernoulli distribution for first 54 dimensions and Pareto for last 3
theta_Y1_54 <- t(colSums(XY_train_Y1[, -c(1, 59, 58, 57, 56)])/nrow(XY_train_Y1))
theta_Y1_3 <- t(nrow(XY_train_Y1)/colSums(log(XY_train_Y1[, -c(1:55, 59)])))
# Finding theta for Y = 0
theta_Y0_54 <- t(colSums(XY_train_Y0[, -c(1, 59, 58, 57, 56)])/nrow(XY_train_Y0))
theta_Y0_3 <- t(nrow(XY_train_Y0)/colSums(log(XY_train_Y0[, -c(1:55, 59)])))


# Calculating objective function 
X_test <- read.csv("X_test.csv", header = FALSE)
Y_test <- read.csv("y_test.csv", header = FALSE)
r_num_test <- c(1:93)
X_test2 <- cbind(r_num_test, X_test)
Y_test2 <- cbind(r_num_test, Y_test)
XY_test <- merge(X_test2,Y_test2, by = "r_num_test")[,-1]
test_thetaY1 <- cbind(X_test, theta_Y1_54)
test_thetaY1 <- cbind(test_thetaY1, theta_Y1_3)
test_thetaY0 <- cbind(X_test, theta_Y0_54)
test_thetaY0  <- cbind(test_thetaY0, theta_Y0_3)
# Getting odds with parameter for Y = 1
predict_Y1 <- matrix(0, nrow = 93, ncol = 57)
for(i in 1:54){
predict_Y1[,i] <- ((test_thetaY1[, i+57])^(test_thetaY1[, i]))*(1-test_thetaY1[, i+57])^(1-test_thetaY1[, i])
}
for(i in 55:57){
predict_Y1[,i] <- test_thetaY1[, i+57]*(test_thetaY1[,i])^(-(test_thetaY1[, i+57]+1))
}
library(tidyverse)
odds_Y1 <- as.data.frame(predict_Y1) %>% 
  mutate(Prod = Reduce(`*`, .)*pi)
odds_Y1 <- odds_Y1[,58]
# Getting odds with parameter for Y = 0
predict_Y0 <- matrix(0, nrow = 93, ncol = 57)
for(i in 1:54){
  predict_Y0[,i] <- ((test_thetaY0[, i+57])^(test_thetaY0[, i]))*(1-test_thetaY0[, i+57])^(1-test_thetaY0[, i])
}
for(i in 55:57){
  predict_Y0[,i] <- test_thetaY1[, i+57]*(test_thetaY0[,i])^(-(test_thetaY0[, i+57]+1))
}
library(tidyverse)
odds_Y0 <- as.data.frame(predict_Y0) %>% 
  mutate(Prod = Reduce(`*`, .)*(1-pi))
odds_Y0 <- odds_Y0[,58]
# Generating model predictions
log_odds<- log(odds_Y1/odds_Y0)
model_predict <- as.data.frame(ifelse(log_odds>=0,1,0))
colnames(model_predict)[1] <- "Y_hat"
y_table <- table(model_predict$Y_hat, Y_test$V1)
# horizontal is Y_test, Vertical is model_predict such that model_predict predicts 57 0 and 36 1 of which 53 0 correct etc
# Prediction accuracy is 86/93 = 0.925


# Part 2b
library(ggplot2)
library(reshape)
df_stemplot <- as.data.frame(cbind(c(1:54), t(theta_Y0_54), t(theta_Y1_54)))
colnames(df_stemplot)[1] <- c("x")
colnames(df_stemplot)[2] <- c("Class Y=0")
colnames(df_stemplot)[3] <- c("Class Y=1")
df_stemplot <- melt(df_stemplot, id = c("x"))
p <- ggplot(aes(group = variable, color = variable, shape = variable), data=df_stemplot) +
  geom_hline(aes(yintercept=0)) +
  geom_segment(aes(x,value,xend=x,yend=value-value)) + 
  geom_point(aes(x,value),size=3) 
p
# 16 and 52 both have theta for class y = 1 higher than for class y = 0, implying it's more correlated to spam 


# Part 2c (Euclidean distance)
library(DMwR)
XY_train <- XY_train[ ,-1]
XY_knn <- rbind(XY_train, XY_test)
train_knn <- as.data.frame(XY_knn)[1:4508,]
test_knn <- as.data.frame(XY_knn)[4509:4601,]
knn_accuracy_table <- matrix(0, nrow = 20, ncol = 1)
for(i in 1:20){
  fit_knn <- kNN(V1.y ~ .,train_knn,test_knn,norm=FALSE,k=i)
  knn_table <- table(test_knn$V1.y,fit_knn)
  knn_accuracy_table[i,] <- sum(diag(knn_table)/93)
}
k <- as.data.frame(c(1:20))
knn_accuracy_table <- as.data.frame(knn_accuracy_table)
knn_accuracy_table2 <- cbind(k, knn_accuracy_table)
colnames(knn_accuracy_table2)[1] <- "k"
colnames(knn_accuracy_table2)[2] <- "Prediction_Accuracy"
p1 <- ggplot(knn_accuracy_table2, mapping = aes(x=k,y= Prediction_Accuracy)) + geom_line()
p1


# Part 2c (Absolute distance 2)
library(kknn)
XY_knn[,"V1.y"] <- as.factor(XY_knn[,"V1.y"])
train_knn <- as.data.frame(XY_knn)[1:4508,]
test_knn <- as.data.frame(XY_knn)[4509:4601,]
knn_accuracy_table <- matrix(0, nrow = 20, ncol = 1)
for(i in 1:20){
  fit_knn <- kknn(V1.y ~ ., train_knn, test_knn, k=i, distance = 1, kernel="rectangular", scale=FALSE)
  knn_table <- table(test_knn[,"V1.y"],fit_knn$fitted.values)
  knn_accuracy_table[i,] <- sum(diag(knn_table)/93)
}
k <- as.data.frame(c(1:20))
knn_accuracy_table <- as.data.frame(knn_accuracy_table)
knn_accuracy_table2 <- cbind(k, knn_accuracy_table)
colnames(knn_accuracy_table2)[1] <- "k"
colnames(knn_accuracy_table2)[2] <- "Prediction_Accuracy"
p1 <- ggplot(knn_accuracy_table2, mapping = aes(x=k,y= Prediction_Accuracy)) + geom_line()
p1

# Part 2d
Y_train_d <- Y_train
Y_train_d[Y_train_d==0] <- -1
Y_train_d <- as.matrix(Y_train_d)
Y_train_d58 <- matrix(Y_train_d, nrow = 4508, ncol = 58)
Y_test_d <- Y_test
Y_test_d[Y_test_d==0] <- -1
Y_test_d <- as.matrix(Y_test_d)
X_train_d <- X_train
X_train_d$newcolumn <- 1
X_train_d <- as.matrix(X_train_d)
X_test_d <- X_test
X_test_d$newcolumn <- 1
X_test_d <- as.matrix(X_test_d)
w <- matrix(0, nrow = 58, ncol = 1)
L <- matrix(0, nrow = 10, ncol = 1)
for(i in 1:2){
  sigmoid = (exp(Y_train_d*(X_train_d%*%w)))/(1+exp(Y_train_d*(X_train_d%*%w)))
  if(is.infinite(sigmoid))
  {
    if(sigmoid > 0)
    {
      sigmoid = exp(10^2)/(1+exp(10^2))
    }
    if(sigmoid < 0)
      sigmoid = exp(10^-5)/(1+exp(10^-5))
  }
  sigmoid_inv = 1-sigmoid
  sigmoid_58 = matrix(sigmoid_inv, nrow = 4508, ncol = 58)
  w = w + (1/((10^5)*((i+1)^0.5)))*(colSums((sigmoid_58)*Y_train_d58*X_train_d))
  value <- colSums(log(sigmoid))
  if(is.infinite(value))
     {
       value = log(10^(-2))
  }
          L[i,] <- value
    print(i)
    flush.console()
}


# Part 2e
w_e <- matrix(0, nrow = 58, ncol = 1)
L_e <- matrix(0, nrow = 10, ncol = 1)
for(i in 1:10){
  sigmoid = (exp(Y_train_d*(X_train_d%*%w_e)))/(1+exp(Y_train_d*(X_train_d%*%w_e)))
  sigmoid_inv = 1-sigmoid
  sigmoid_58 = matrix(sigmoid_inv, nrow = 4508, ncol = 58)
  sigmoid2 = (exp((X_train_d%*%w_e)))/(1+exp((X_train_d%*%w_e)))
  sigmoid582 = matrix(sigmoid2, nrow=4508, ncol=58)
  sigmoid_inv2 = 1-sigmoid2
  sigmoid_582 = matrix(sigmoid_inv, nrow = 4508, ncol = 58)
  w_e = w_e - (1/(((i+1)^0.5)))*(solve((-colSums(sigmoid582*(sigmoid_582)*t(X_train_d)%*%(X_train_d)))))%*%(colSums((sigmoid_58)*Y_train_d58*X_train_d))
  
  L[i,] <- colSums(log(sigmoid))
  print(i)
  flush.console()
}




