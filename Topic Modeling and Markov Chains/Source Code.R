# Problem 1a
d <- read.csv("CFB2017_scores.csv")
teamnames<-readLines("TeamNames.txt")
# unnormalized matrix of zeros
Mhat<-matrix(0, nrow=763, ncol=763)
colnames(d)<-c("teamA_id", "teamA_points", "teamB_id", "teamB_points")
# processing games and updating matrix per the question
for (a in 1:nrow(d)) {
  i <- d$teamA_id[a]
  j <- d$teamB_id[a]
  points_i <- d$teamA_points[a]
  points_j <- d$teamB_points[a]
  if (points_i > points_j){
    Mhat[i,i] <- Mhat[i,i] + 1 + (points_i/(points_i + points_j))
    Mhat[j,j] <- Mhat[j,j] + 0 + (points_j/(points_i + points_j))
    Mhat[i,j] <- Mhat[i,j] + 0 + (points_j/(points_i + points_j))
    Mhat[j,i] <- Mhat[j,i] + 1 + (points_i/(points_i + points_j))
  }
  if (points_i < points_j){
    Mhat[i,i] <- Mhat[i,i] + 0 + (points_i/(points_i + points_j))
    Mhat[j,j] <- Mhat[j,j] + 1 + (points_j/(points_i + points_j))
    Mhat[i,j] <- Mhat[i,j] + 1 + (points_j/(points_i + points_j))
    Mhat[j,i] <- Mhat[j,i] + 0 + (points_i/(points_i + points_j))
  }
}
M<-Mhat
for (j in 1:nrow(Mhat))
{
  total <- sum(Mhat[j,])
  M[j,] <- M[j,]/total
}

# Initializing and updating w0
w0<-rep(1/763,763)
w0<-as.matrix(w0)
stationary <- function(iterations, w0, M){
for (i in 1:iterations){
  w0 <- w0%*%M
}
  w0<-w0/sum(w0)
  return(w0)
}
# for t = 10
w10 <- stationary(10,t(w0),M)
Index <- order(w10, decreasing = TRUE)[1:25]
top25 <- w10[Index]
Teams<-as.data.frame(teamnames)
colnames(Teams)<-"team"
top25_teams <- as.data.frame(Teams[Index,])
result10 <- cbind(Index, top25, top25_teams)
result10$Rank <- rownames(result10)
result10 <- result10[,c(4,3,2)]
colnames(result10)[2] <- "Team"
colnames(result10)[3] <- "wt"
# for t = 100
w100 <- stationary(100,t(w0),M)
Index <- order(w100, decreasing = TRUE)[1:25]
top25 <- w100[Index]
Teams<-as.data.frame(teamnames)
colnames(Teams)<-"team"
top25_teams <- as.data.frame(Teams[Index,])
result100 <- cbind(Index, top25, top25_teams)
result100$Rank <- rownames(result100)
result100 <- result100[,c(4,3,2)]
colnames(result100)[2] <- "Team"
colnames(result100)[3] <- "wt"
# for t = 1000
w1000 <- stationary(1000,t(w0),M)
Index <- order(w1000, decreasing = TRUE)[1:25]
top25 <- w1000[Index]
Teams<-as.data.frame(teamnames)
colnames(Teams)<-"team"
top25_teams <- as.data.frame(Teams[Index,])
result1000 <- cbind(Index, top25, top25_teams)
result1000$Rank <- rownames(result1000)
result1000 <- result1000[,c(4,3,2)]
colnames(result1000)[2] <- "Team"
colnames(result1000)[3] <- "wt"
# for t = 10000
w10000 <- stationary(10000,t(w0),M)
Index <- order(w10000, decreasing = TRUE)[1:25]
top25 <- w10000[Index]
Teams<-as.data.frame(teamnames)
colnames(Teams)<-"team"
top25_teams <- as.data.frame(Teams[Index,])
result10000 <- cbind(Index, top25, top25_teams)
result10000$Rank <- rownames(result10000)
result10000 <- result10000[,c(4,3,2)]
colnames(result10000)[2] <- "Team"
colnames(result10000)[3] <- "wt"


# Problem 1b
# getting the first eigenvector and eigenvalue of M transposed
eigenM <- eigen(t(M))
eigenvector1 <- as.matrix(eigenM$vectors[,1])
# geting w infinity
w_infinity <- t(eigenvector1/sum(eigenvector1))
# getting distance measure to plot
iterations <- 10000
distance <- rep(NA,iterations)
state_vec <- t(w0)
for (i in 1:iterations){
  state_vec <- state_vec%*%M
  total_dist <- sum(abs(state_vec-w_infinity))
  distance[i] <- total_dist
}
library(ggplot2)
df <- as.data.frame(distance)
df$t <- as.numeric(rownames(df))
g <- ggplot(df, aes(x=t, y=distance, group=1))+geom_line()

# Problem 2a
# Reading in the data and creating matrix X
X<-matrix(0,nrow=3012, ncol=8447)
nyt_d <- readLines("nyt_data.txt")
for (i in 1:length(nyt_d)){
  doc_no <- nyt_d[i] 
  word_count <- strsplit(doc_no,",")[[1]] # gives the index and count
  for (j in 1:length(word_count)){
    word <- strsplit(word_count[j], ":")[[1]] # separating index and count
    idx <- as.numeric(word[1]) 
    cnt <- as.numeric(word[2])
    X[idx,i] <- cnt # filling X with counts
  }
}
# initialize W and H
K <- 25
N <- nrow(X)
M<-ncol(X)
NK_val <- N*K
KM_val <- K*M
# initializing values from Uniform (1,2) distribution
W <-matrix(runif(NK_val,min=1, max=2), nrow=N, ncol=K) 
H <- matrix(runif(KM_val,min=1, max=2), nrow=K, ncol=M)
# running NMF with divegence penalty
iterations <- 100
divergence <- rep(NA,iterations)
for (t in 1:iterations){
  Purple <- X/((W%*%H+(10^-16)))
  Wt <- t(W)
  Wt <- t(apply(Wt, 1, function(x) x/sum(x)))
  H <- H*(Wt%*%(Purple)) #update H
  Purple <- X/((W%*%H+(10^-16)))
  Ht <- t(H)
  Ht <- (apply(Ht, 2, function(x) x/sum(x)))
  W <- W*(Purple%*%Ht) #update W
  divergence[t] <- -(sum(X*log(W%*%H+(10^-16))-(W%*%H+10^-16)))
  print(t)
}
# plot of objective function against iteration
library(ggplot2)
df <- data.frame(divergence)
df$iterations <- as.numeric(rownames(df))
ggplot(df, aes(x=iterations, y=divergence, group = 1)) + geom_line()


# Problem 2b
# normalizing column of W
W <- (apply(W, 2, function(x) x/sum(x)))
nyt_vocab <- readLines("nyt_vocab.dat")
top_words <- data.frame(matrix(0,nrow=10))
# finding 10 words with highest weight
for (topic in 1:ncol(W))
{
  vec <- W[,topic]
  index <- order(vec, decreasing = TRUE)
  top10 <- index[1:10]
  words <- nyt_vocab[top10]
  top_words[,paste0("Words in topic ",topic)] <- words
  wordWeights <- vec[top10]
  top_words[,paste0("Weights ",topic)] <- wordWeights
  
}
top_words <- top_words[,-1]

