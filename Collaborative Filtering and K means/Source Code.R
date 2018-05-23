# Problem 1a
# Generating the 500 observations
mu1<-as.vector(c(0,0))
sigma1<-matrix(c(1,0,0,1), nrow=2, ncol=2)
mu2<-as.vector(c(3,0))
sigma2<-matrix(c(1,0,0,1), nrow=2, ncol=2)
mu3<-as.vector(c(0,3))
sigma3<-matrix(c(1,0,0,1), nrow=2, ncol=2)
pi<-c(0.2,0.5,0.3)
library(mvtnorm)
set.seed(1)
N=500
dists = runif(N)
data = matrix(0, nrow=N, ncol=2)
for(i in 1:N){
  if(dists[i]<pi[1]){
    data[i,] = rmvnorm(1, mean=mu1, sigma = sigma1)
  } else if(dists[i]<pi[2]){
    data[i,] = rmvnorm(1, mean=mu2, sigma = sigma2)
  } else {
    data[i,] = rmvnorm(1, mean=mu3, sigma = sigma3)
  }
}
# Implementing the K-means algorithm
K_means <- function (k) {
  X <- data
  c <- rep(0, N)
  mu<-matrix(runif(2*k), ncol=k) # random initialization of mu
  T <- 20
  obj <- c()
  for (t in 1:T){
    for (j in 1:N){
      dist_m <- rbind(t(mu), X[j,])
      distance <- as.matrix(dist(dist_m))
      distance <- distance[k+1, 1:k]
      c[j] <- which.min(distance) # choosing the class that minimizes euclidean distance for each j
    }
    counts <- table(c)
    indic <- matrix(0, nrow=N, ncol=k)
    for (i in 1:N)
    {
      for (j in 1:k)
      {
        if(c[i] == j)
        {
          indic[i, j] = 1
        }
      }
    }
    mu2 <- t(X)%*%indic
    for (j in 1:k)
      {
      n <- counts[j]
      mu2[,j] <- (1/n)*mu2[,j] # calculating mu for each k
    }
    for(j in 1:k){
      mu[ ,j] <- mu2[ , j] # updating mu before next iteration
    }
    distance <- rep(0, k)
    data1 <- t(X)
    for (i in 1:N)
    {
      vec1<-as.matrix(data1[,i]) # extract the ith observation from the data matrix
      for (j in 1:k)
      {
        if (c[i] == j)
        {
          vec2<-as.matrix(vec1-mu[,j]) # the difference between the ith observation and the jth cluster mean
          squared_distance<-as.numeric(t(vec2)%*%vec2)
          distance[j]<-distance[j]+squared_distance
        }  
      }
    }
    obj[t] <- sum(distance) # calculating the objective function for each iteration
    print(t)
  }
  output<-list(obj, c)
  return(output)
}
values <- matrix(0, nrow = 20, ncol = 4)
for (k in 2:5) # running function for k=2,3,4,5 and plotting objective function against iteration
{
  values[,k-1]<-K_means(k)[[1]]
}
colnames(values)<-c("2", "3", "4", "5")
iterations <-c(1:20)
values<-cbind(values, iterations)
library(ggplot2)
library(reshape)
values_long <- melt(values, id="iterations")
colnames(values_long)[1] <- "iteration"
colnames(values_long)[2] <- "k"
values_long <- values_long[values_long$k!="iterations",]
p <-ggplot(data=values_long,
              aes(x=iteration, y=value, colour=k)) +
  geom_line()+
  ylab("Objective Function Value")
p


# Problem 1b
# Clustering when K = 3
output<-K_means(3)
class<-output[2]
class_data<-data.frame(data)
class_data<-cbind(class_data,class)
colnames(class_data)<-c("x", "y", "cluster")
p2<-ggplot(class_data, aes(x,y))+
  geom_point(aes(color=factor(cluster)))+
  labs(xlab="x")+
  labs(ylab="y")+
  labs(color = "Cluster Number") +
  labs(title ="Clustering when K = 3")
p2 
# Clustering when K = 5
output<-K_means(5)
class<-output[2]
class_data<-data.frame(data)
class_data<-cbind(class_data,class)
colnames(class_data)<-c("x", "y", "cluster")
p3<-ggplot(class_data, aes(x,y))+
  geom_point(aes(color=factor(cluster)))+
  labs(xlab="x")+
  labs(ylab="y")+
  labs(color = "Cluster Number") +
  labs(title ="Clustering when K = 5")
p3 

# Problem 2a
library(tidyr)
ratings_train<-read.csv("ratings.csv", header = FALSE)
colnames(ratings_train)[1] <- "user_id"
colnames(ratings_train)[2] <- "movie_id"
colnames(ratings_train)[3] <- "rating"
ratings_train<-spread(ratings_train,movie_id, rating)[,-1]
ratings_test<-read.csv("ratings_test.csv", header = FALSE)
colnames(ratings_test)[1] <- "user_id"
colnames(ratings_test)[2] <- "movie_id"
colnames(ratings_test)[3] <- "rating"
ratings_test<-spread(ratings_test,movie_id, rating)[,-1]
library(mvtnorm)
lambda<-1
d<-10
sigmasqr<-0.25
runs<-10
iterations<-100
RMSE<-rep(NA,runs)
# initializing two lists to store U and V matrices for each run
listU<-list()
listV<-list()
obj<-as.data.frame((matrix(nrow=iterations, ncol=runs)))
N1<-nrow(ratings_train)
N2<-ncol(ratings_train)
# creating empty matrices for U and V
U<-matrix(0,nrow=N1, ncol=d)
V<-matrix(0,nrow=d, ncol=N2)
for (t in 1: runs){
  mu<-as.vector(rep(0,d))
  sigma<-lambda*diag(d)
  f<-function(vector)
  {
    vector=rmvnorm(1, mean=mu, sigma = sigma) # initializing U and V
  }
  U<-t(apply(U,1,f))
  V<-apply(V,2,f)
  mat_train<-as.matrix(ratings_train)
  obj2<-c()
  for (num in 1:iterations){
    for (user in 1:nrow(U)){ # running inner loop for U
      given_u<-as.matrix(mat_train[user,])
      given_u<-t(given_u)
      index_movies_rated<-which(!is.na(given_u))
      Vtemp<-V[,index_movies_rated]
      A<-Vtemp%*%t(Vtemp)
      B<-lambda*sigmasqr*diag(d)
      C<-solve(A+B)
      M<-as.matrix(t(given_u[index_movies_rated]))
      D<-Vtemp%*%t(M)
      E<-C%*%D
      U[user,]<-t(E)
    }
    for (movie in 1:ncol(V)){ # running inner loop for V
      given_v<-as.matrix(mat_train[,movie])
      index_users_whorated<-which(!is.na(given_v))
      Utemp<-as.matrix(U[index_users_whorated,])
      if (length(index_users_whorated)>1){
        Utemp<-t(Utemp)
      }
      W<-Utemp%*%t(Utemp)
      X<-lambda*sigmasqr*diag(d)
      Y<-solve(W+X)
      M2<-as.matrix(t(given_v[index_users_whorated]))
      Z<-(Utemp)%*%t(M2)
      z<-Y%*%Z
      V[,movie]<-z
    }
    predicted_ratings<-U%*%V # getting predicted ratings
    missing_count<-(sum(is.na(mat_train)))
    total<-nrow(mat_train)*ncol(mat_train)
    num1<-total-missing_count
    errorMat<-predicted_ratings-mat_train
    errorMat<-errorMat^2
    first<-(sum(errorMat, na.rm=TRUE))/(2*sigmasqr)
    missing_count_new<-(sum(is.na(errorMat)))
    missing_count_new==missing_count
    U1<-as.matrix(c(t(U)))
    second<-(t(U1)%*%U1)*(lambda/2)
    V1<-as.matrix(c((V)))
    third<-(t(V1)%*%V1)*(lambda/2)
    obj2[num] <- -(first+second+third)
  }
  listU[[t]]<-U
  listV[[t]]<-V
  obj[,t] <- obj2
  mat_test<-as.matrix(ratings_test)
  predicted_ratings<-U%*%V
  missing_count<-(sum(is.na(mat_test)))
  predicted_data<-predicted_ratings[1:nrow(mat_test),1:ncol(mat_test)]
  error_mat<-(predicted_data-mat_test)^2
  total<-nrow(mat_test)*ncol(mat_test)
  num<-total-missing_count
  MSE<-(sum(error_mat, na.rm=TRUE))/num
  rmse<-sqrt(MSE)
  RMSE[t]<-rmse
  print(t)
}
iterations2<-c(1:100)
obj_ans<-cbind(obj, iterations2)
colnames(obj_ans)<-c("1","2", "3", "4", "5", "6","7","8","9","10", "iterations")
obj_ans<-obj_ans[-1,] # to get iterations 2 to 10
library(reshape2)
library(ggplot2)
results_long <- melt(obj_ans, id="iterations")
colnames(results_long)[2] <- "Run_Number"
colnames(results_long)[3] <- "objectivevalue"
p <-ggplot(data=results_long,
           aes(x=iterations, y=objectivevalue, colour=Run_Number)) +
  geom_line()+
  ylab("Log Joint Likelihood")
p
final_values<-obj[100,]
final_values<-as.data.frame(t(final_values))
RMSE <- sqrt(RMSE)
results_1<-cbind(final_values, RMSE)
colnames(results_1)[1]<-"final_obj_function_values"
rownames(results_1)<-c("1","2","3","4","5","6","7","8","9","10")      
results_final<-results_1[order(results_1[,1],decreasing=TRUE),]


# Problem 2b
movies_data<-readLines("movies.txt", n =-1)
movies_data<-data.frame(movies_data)
Vbest<-as.matrix(unlist(listV[[6]])) # since the 6th run had the highest obj, extract its data
closest_movies<-as.matrix(dist(t(Vbest), method="euclidean", p=2))
interested_movies <-c(50,485,182) # indexes of the 3 movies in question
combined_closest_movies<-closest_movies[interested_movies,]
StarWarsNN<-data.frame(combined_closest_movies[1,])
StarWarsNN <- cbind(movie_id = rownames(StarWarsNN), StarWarsNN)
colnames(StarWarsNN)[2]<-"Euclidean_distance"
StarWarsNN<-StarWarsNN[order(StarWarsNN$Euclidean_distance),]
StarWarsNN<-StarWarsNN[1:11,]
StarWarsNN$movie_id<-as.numeric(levels(StarWarsNN$movie_id))[StarWarsNN$movie_id]
StarWarsNNids<-StarWarsNN$movie_id
StarWarsNNnames<-movies_data[StarWarsNNids,]
StarWarsNNfinal<-cbind(as.character(StarWarsNNnames), StarWarsNN)
colnames(StarWarsNNfinal)[1]<-"movie_name"
StarWarsNNfinal<-StarWarsNNfinal[c(2,1,3)]

FairLadyNN<-data.frame(combined_closest_movies[2,])
FairLadyNN <- cbind(movie_id = rownames(FairLadyNN), FairLadyNN)
colnames(FairLadyNN)[2]<-"Euclidean_distance"
FairLadyNN<-FairLadyNN[order(FairLadyNN$Euclidean_distance),]
FairLadyNN<-FairLadyNN[1:11,]
FairLadyNN$movie_id<-as.numeric(levels(FairLadyNN$movie_id))[FairLadyNN$movie_id]
FairLadyNNids<-FairLadyNN$movie_id
FairLadyNNnames<-movies_data[FairLadyNNids,]
FairLadyNNfinal<-cbind(as.character(FairLadyNNnames), FairLadyNN)
colnames(FairLadyNNfinal)[1]<-"movie_name"
FairLadyNNfinal<-FairLadyNNfinal[c(2,1,3)]

GoodFellasNN<-data.frame(combined_closest_movies[3,])
GoodFellasNN <- cbind(movie_id = rownames(GoodFellasNN), GoodFellasNN)
colnames(GoodFellasNN)[2]<-"Euclidean_distance"
GoodFellasNN<-GoodFellasNN[order(GoodFellasNN$Euclidean_distance),]
GoodFellasNN<-GoodFellasNN[1:11,]
GoodFellasNN$movie_id<-as.numeric(levels(GoodFellasNN$movie_id))[GoodFellasNN$movie_id]
GoodFellasNNids<-GoodFellasNN$movie_id
GoodFellasNNnames<-movies_data[GoodFellasNNids,]
GoodFellasNNfinal<-cbind(as.character(GoodFellasNNnames), GoodFellasNN)
colnames(GoodFellasNNfinal)[1]<-"movie_name"
GoodFellasNNfinal<-GoodFellasNNfinal[c(2,1,3)]
