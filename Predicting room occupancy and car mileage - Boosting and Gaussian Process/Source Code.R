# Problem 1a
# Gaussian Kernel
kernel_gauss = function(b, x1,x2) {
  kernel_out <- matrix(rep(0, nrow(x1)*nrow(x2)), nrow=nrow(x1))
  for (i in 1:nrow(kernel_out)) {
    for (j in 1:ncol(kernel_out)) {
      kernel_out[i,j] <- exp(-(1/b)*((x1[i,1]-x2[j,1])^2+(x1[i,2]-x2[j,2])^2+(x1[i,3]-x2[j,3])^2
                                      +(x1[i,4]-x2[j,4])^2+(x1[i,5]-x2[j,5])^2+(x1[i,6]-x2[j,6])^2
                                      +(x1[i,7]-x2[j,7])^2))
    }
  }
  return(kernel_out)
}
# Code to implement the Gaussian Process
gauss_solve = function(x_train, y_train, x_pred, b, sigmasq) {
  solution = list()
  # compute training covariance matrix (used to get relationships in training)
  k_xx = kernel_gauss(b, x_train,x_train)
  # compute covariance between training and testing (used to predict weights into new data)
  k_xp_x = kernel_gauss(b, x_pred,x_train)
  # compute covariance between testing (used to estimate covariance of predictions in new data)
  k_xp_xp = kernel_gauss(b, x_pred,x_pred)
  # Mean and covariance functions
  Vinv = solve( sigmasq * diag(1, ncol(k_xx))+ k_xx)
  # compute the estimate and variance for the prediction 
  solution[["miu"]] = k_xp_x %*% Vinv %*% y_train
  solution[["var"]] = sigmasq + k_xp_xp - k_xp_x %*% Vinv %*% t(k_xp_x)
  return( solution )
}


# Problem 1b
X_train_gp <- read.csv("X_train_gp.csv", header = FALSE)
Y_train_gp <- read.csv("y_train_gp.csv", header = FALSE)
X_train_gp <- as.matrix(X_train_gp)
y_train_gp <- as.matrix(Y_train_gp)
X_test_gp <- read.csv("X_test_gp.csv", header = FALSE)
y_test_gp <- read.csv("y_test_gp.csv", header = FALSE)
X_test_gp <- as.matrix(X_test_gp)
y_test_gp <- as.matrix(y_test_gp)
RMSE = matrix(0, 6, 10)
b <- c(5, 7, 9, 11, 13, 15)
sigmasq <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (i in b){
  for (j in sigmasq){
    RMSE[which(b==i),which(sigmasq==j)] = sqrt(mean((y_test_gp - as.matrix(unlist(gauss_solve(X_train_gp, y_train_gp, X_test_gp,
                                                    i, j)["miu"]))) ^ 2))
  }
}
RMSE

# Problem 1c
min(RMSE)
which(RMSE == min(RMSE), arr.ind = TRUE)
## From homework 1, the lowest RMSE in part c was 2.633644 at lambda = 0 (least squares > ridge as optimal regulariazation parameter = 0)
## lowest RMSE in part d was 2.192574 at p = 2 and lambda = 23


# Problem 1d
# Gaussian Kernel for only 4th dimension
kernel_gauss = function(b, x1,x2) {
  kernel_out <- matrix(rep(0, nrow(x1)*nrow(x2)), nrow=nrow(x1))
  for (i in 1:nrow(kernel_out)) {
    for (j in 1:ncol(kernel_out)) {
      kernel_out[i,j] <- exp(-(1/b)*((x1[i,4]-x2[j,4])^2))
    }
  }
  return(kernel_out)
}
# Code to implement the Gaussian Process for only 4th dimension
gauss_solve = function(x_train, y_train, x_pred, b, sigmasq) {
  solution = list()
  # compute training covariance matrix (used to get relationships in training)
  k_xx = kernel_gauss(b, x_train,x_train)
  # compute covariance between training and testing (used to predict weights into new data)
  k_xp_x = kernel_gauss(b, x_pred,x_train)
  # compute covariance between testing (used to estimate covariance of predictions in new data)
  k_xp_xp = kernel_gauss(b, x_pred,x_pred)
  # Mean and covariance functions
  Vinv = solve( sigmasq * diag(1, ncol(k_xx))+ k_xx)
  # compute the estimate and variance for the prediction (note, does not depend on y) [eq 2.19 or 2.23]
  solution[["miu"]] = k_xp_x %*% Vinv %*% y_train
  solution[["var"]] = sigmasq + k_xp_xp - k_xp_x %*% Vinv %*% t(k_xp_x)
  return( solution )
}
# Plotting scatter plot of y against x[4]
library(ggplot2)
XY_train_merged <- cbind(X_train_gp[,4], y_train_gp)
colnames(XY_train_merged)[1] <- "x4"
colnames(XY_train_merged)[2] <- "y"
gg <- ggplot(data=as.data.frame(XY_train_merged), aes(x=x4, y = y)) + geom_point()
gg
# Computing predictive mean using gaussian
miu_4 <- as.matrix(unlist(gauss_solve(X_train_gp, y_train_gp, X_train_gp,
                             5, 2)["miu"]))
Xmiu4_train_merged <- cbind(X_train_gp[,4], miu_4)
Xmiu4_train_merged <- as.data.frame(Xmiu4_train_merged)
colnames(Xmiu4_train_merged)[1] <- "x4"
colnames(Xmiu4_train_merged)[2] <- "miu4"
gg2 <- gg + geom_smooth(data=Xmiu4_train_merged, aes(x=x4, y = miu_4), color = 'red')
gg2

# Problem 2a
library(dplyr)
library(ggplot2)
Xtrain<-read.csv("X_train_boost.csv", header=FALSE)
ytrain<-read.csv("y_train_boost.csv", header=FALSE)
Xtest<-read.csv("X_test_boost.csv", header=FALSE)
ytest<-read.csv("y_test_boost.csv", header=FALSE)
## Adding a dimension of 1s
V0<-rep(1,nrow(Xtrain))
Xtrain<-cbind(V0, Xtrain)
V0<-rep(1,nrow(Xtest))
Xtest<-cbind(V0, Xtest)
rm(V0)
n<-nrow(Xtrain)
n
weightStart<-rep(1/n, nrow(Xtrain))
str(weightStart)
id<-rownames(Xtrain)
Xtrain_1<-cbind(id=id, Xtrain, weightStart)
Xtrain_1<-as.matrix(Xtrain_1)
class(Xtrain_1)
id<-rownames(Xtest)
Xtest_1<-cbind(id=id, Xtest)
XtrainMat<-as.matrix(Xtrain)
ytrainMat<-as.matrix(ytrain)
XtestMat<-as.matrix(Xtest)
ytestMat<-as.matrix(ytest)
## Function to generate epsilon
findError<-function(w, wt)
{
  yhat<-sign(XtrainMat%*%w)
  temp<-as.data.frame(cbind(ytrainMat, yhat, wt))
  colnames(temp)<-c("y","yhat", "wt")
  size<-nrow(temp)
  temp1<-temp
  temp1$difference=temp1$yhat-temp1$y
  temp1<-subset(temp1, difference!=0)
  epsilon<-sum(temp1$wt)
  output<-list(epsilon, temp)
  return (output)
}
## Function to generate prediction
prediction<-function(X,w,Alpha)
{
  Alpha<-as.vector(Alpha)
  X<-as.matrix(X)
  w<-as.matrix(w)
  ifelse(ncol(w)>1, w<-t(w), w<-w)
  A<-sign(X%*%w) 
  B<-Alpha*t(A)
  C<-t(apply(B,2, sum))
  C<-t(apply(C,2, sign))
  dim(C)
  predictvalues<-C[1,]
  return (predictvalues)
}
## Function to generate prediction error
perform<-function(y, yhat)
{
  data<-cbind(y, yhat)
  data$error<-data[,1]-data[,2]
  data$error[data$error!=0]<-1
  pred_error<-mean(data$error)
  pred_error
  
}
## ADABOOST Algorithm
T<-1500
trainError<-rep(0,T)
testError<-rep(0,T)
W<-matrix(0, nrow=T, ncol=6)
ALPHA<-rep(0,T)
EPSILON<-rep(0,T)
weight<-matrix(0,nrow=nrow(Xtrain), ncol=T)
wt<-weightStart
count <- matrix(0, nrow = nrow(Xtrain), ncol = T)
for (u in 1:T)
{
  weight[,u]<-wt  
  vec<-Xtrain_1[,"id"]
  sample_vec<-sample(vec, size=n, replace=TRUE, prob=wt)
  sample_vec<-as.numeric(sample_vec) 
  sample_vec2<- as.data.frame(sample_vec)
  count[,u] <- sample_vec2$sample_vec
  XtrainMat<-as.matrix(Xtrain)
  Xbs<-XtrainMat[sample_vec, ]
  yid<-rownames(ytrain)
  yid<-as.numeric(yid)
  ytrainNew<-cbind(ytrain, yid)
  ytrainNew<-ytrainNew[sample_vec, ]
  ytrainNew<-ytrainNew$V1
  ybs<-as.matrix(ytrainNew)
  rm(ytrainNew)
  w<-(solve(t(Xbs)%*%Xbs))%*%(t(Xbs)%*%ybs)
  W[u, ]<-t(w)  
  result<-findError(w, wt)
  epsilon<-as.numeric(result[1])
  epsilon
  temp<-as.data.frame(result[2])
  ## If epsilon>0.5, change the sign 
  if (epsilon>0.5)
  {
    w<- -(w)
    W[u,]<-w
    result<-findError(w, wt)
    epsilon<-as.numeric(result[1])
    temp<-as.data.frame(result[2])
  }
  alpha<-0.5*(log((1-epsilon)/epsilon))
  ALPHA[u]<-alpha
  ALPHA
  EPSILON[u]<-epsilon
  ## Updating weights for next iteration
  temp<-temp%>%
    mutate(alpha=alpha)%>%
    mutate (wt_new=wt*exp(yhat*y*-(alpha)))
  constant<-sum(temp$wt_new)
  temp$wt_new<-(temp$wt_new)/constant
  wt<-temp$wt_new           
  rm(temp)  
  print(paste("After round: ", u))
  print(paste("EPSILON: ", EPSILON[u]))
  print(paste("ALPHA: ", ALPHA[u]))
} 
## Computing training and testing error for each iteration
for (iteration in 1:T)
{
  ALPHAtemp<-data.frame(ALPHA)
  index<-as.numeric(rownames(ALPHAtemp))
  ALPHAtemp<-ALPHA[1:iteration]
  ALPHAtemp<-as.numeric(ALPHAtemp)
  ALPHAtemp
  Wtemp<-as.matrix(W[1:iteration, ])
  Wtemp
  yhatTrain<-prediction(Xtrain, Wtemp, ALPHAtemp)
  pred_error_train<-perform(ytrain, yhatTrain)
  trainError[iteration]<-pred_error_train
  yhatTest<-prediction(Xtest, Wtemp, ALPHAtemp)
  pred_error_test<-perform(ytest, yhatTest)
  testError[iteration]<-pred_error_test
}
## Plotting the errors
compareError<-data.frame(trainError, testError)
compareError$id<-as.numeric(rownames(compareError))
str(compareError)
colnames(compareError)<-c("trainError", "testError", "t")
gg <- ggplot (compareError, aes(t))+
  geom_line(aes(y=trainError, color="trainError"))+
  geom_line(aes(y=testError, color="testError"))

# Problem 2b
upperbound <- matrix(0, nrow = T, ncol = 1)
D <- rep(0, T)
for (i in 1:T){
  D[i] <- (0.5 - EPSILON[i,1])^2
  upperbound[i,] <- exp(-2*sum(D))
}
upperbound <- as.data.frame(upperbound)
upperbound$t <- c(1:1500)
gg2 <- ggplot (upperbound, aes(t))+
  geom_line(aes(y=V1, color="upperbound"))

# Problem 2c
hist <- table(count)
hist <- as.data.frame(hist)
gg3 <- ggplot (hist, aes(count))+
  geom_col(aes(y=Freq))

# Problem 2d
ALPHA <- as.data.frame(ALPHA)
ALPHA$t <- c(1:1500)
gg4 <- ggplot (ALPHA, aes(t))+
  geom_line(aes(y=ALPHA, color="ALPHA"))

EPSILON <- as.data.frame(EPSILON)
EPSILON$t <- c(1:1500)
gg5 <- ggplot (as.data.frame(EPSILON), aes(t)) +
  geom_line(aes(y=EPSILON, color="EPSILON"))
