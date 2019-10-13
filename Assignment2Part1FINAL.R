library(dplyr)
library(tidyverse)
library(ggplot2)
library(devtools)
library(modelr)

##PART1 - AssIGNMENT2 

######1.2

ridge <- function(data, dep, indep, lambda) {
  y <- as.matrix(data[,dep])
  x <- as.matrix(data[,indep])
  x <- cbind(1,x)
  Ix <- diag(x)
  beta_ridge <- as.vector(c(solve((crossprod(x,x)) + (lambda*Ix)) %*% (crossprod(x,y))))
  
  names(beta_ridge) <- colnames(x)
  return_obj <- list(beta = beta_ridge, dep = dep , indep = indep)
  
  class(return_obj) <- "ridge"
  return(return_obj)
  
  
  
}


ridge_model<-ridge(data = prostate, 
                dep = "lpsa",
                indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                lambda = 1)


#####1.3

pred <- function(data, indep, beta_hat){
  x <- as.matrix(data[,indep]) 
  beta <- as.matrix(beta_hat)
  preds <- (x %*% beta)
  return(preds)
}



prediction <- pred(data= prostate, indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                   beta_hat = ridgemod$beta[c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45")])

#####1.4
### first model wITHOUT intercept

cv4 <- function(data, dep, indep, lambda) {
  A <- matrix(0, nrow = 5, ncol = 1)
  for (i in 1:5){
    train <- data %>% filter(group != i)
    test <- data %>% filter(group == i)
    y <- as.matrix(train[,dep])
    x <- as.matrix(train[,indep])
    y_test <- as.matrix(test[,dep])
    x_test <- as.matrix(test[,indep])
    x <- cbind(1,x)
    Ix <- diag(x)
    beta_ridge <- as.vector(c(solve((crossprod(x,x)) + (lambda*Ix)) %*% (crossprod(x,y))))
    beta_ridge2 <-beta_ridge[ -c(1) ]
    y_hat <- x_test %*% (as.matrix(beta_ridge2))
    group_mse <- mean((as.vector(y_hat)-as.vector(y_test))^2)
    A[i,] <- group_mse
    
  }
  mse<-sum(A)/5
  return(mse)
}

MSE<-cv4(data = prostate, 
         dep = "lpsa",
         indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
         lambda = 1)
### second model WITH intercept

cv4intercept <- function(data, dep, indep, lambda) {
  A <- matrix(0, nrow = 5, ncol = 1)
  for (i in 1:5){
    train <- data %>% filter(group != i)
    test <- data %>% filter(group == i)
    y <- as.matrix(train[,dep])
    x <- as.matrix(train[,indep])
    y_test <- as.matrix(test[,dep])
    x_test <- as.matrix(test[,indep])
    x <- cbind(1,x)
    x_test <- cbind(1,x_test)
    Ix <- diag(x)
    beta_ridge <- as.vector(c(solve((crossprod(x,x)) + (lambda*Ix)) %*% (crossprod(x,y))))
    y_hat <- x_test %*% (as.matrix(beta_ridge))
    group_mse <- mean((as.vector(y_hat)-as.vector(y_test))^2)
    A[i,] <- group_mse
    
  }
  mse<-sum(A)/5
  return(mse)
}

MSEintercept<-cv4intercept(data = prostate, 
         dep = "lpsa",
         indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
         lambda = 1)
#### 1.5 and 1.6
#### first model and ggplot WITHOUT intercept
cv5 <- function(data, dep, indep, lambda) {
  B <- matrix(0, nrow = length(lambda), ncol = 2)
  for(j in lambda){
    A <- matrix(0, nrow = 5, ncol = 1)
    for (i in 1:5){
      train <- data %>% filter(group != i)
      test <- data %>% filter(group == i)
      y <- as.matrix(train[,dep])
      x <- as.matrix(train[,indep])
      y_test <- as.matrix(test[,dep])
      x_test <- as.matrix(test[,indep])
      x <- cbind(1,x)
      Ix <- diag(x)
      beta_ridge <- as.vector(c(solve((crossprod(x,x)) + (j*Ix)) %*% (crossprod(x,y))))
      beta_ridge2 <-beta_ridge[ -c(1) ]
      y_hat <- x_test %*% (as.matrix(beta_ridge2))
      group_mse <- mean((as.vector(y_hat)-as.vector(y_test))^2)
      A[i,] <- group_mse
      
    }
    mse<-sum(A)/5
    B[j+1,2]<-mse
    B[j+1,1]<-j
    returnobject<-as.data.frame(B)
    names(returnobject)<- c("lambda","mse")
  }
  return(returnobject)
}


MSE_lambda<-cv5(data = prostate, 
               dep = "lpsa",
               indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
               lambda = 0:50)

ggplot(MSE_lambda, aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point()
##### second model and ggplot WITH intercept

cv5intercept <- function(data, dep, indep, lambda) {
  B <- matrix(0, nrow = length(lambda), ncol = 2)
  for(j in lambda){
    A <- matrix(0, nrow = 5, ncol = 1)
    for (i in 1:5){
      train <- data %>% filter(group != i)
      test <- data %>% filter(group == i)
      y <- as.matrix(train[,dep])
      x <- as.matrix(train[,indep])
      y_test <- as.matrix(test[,dep])
      x_test <- as.matrix(test[,indep])
      x <- cbind(1,x)
      x_test <- cbind(1,x_test)
      Ix <- diag(x)
      beta_ridge <- as.vector(c(solve((crossprod(x,x)) + (j*Ix)) %*% (crossprod(x,y))))
      y_hat <- x_test %*% (as.matrix(beta_ridge))
      group_mse <- mean((as.vector(y_hat)-as.vector(y_test))^2)
      A[i,] <- group_mse
      
    }
    mse<-sum(A)/5
    B[j+1,2]<-mse
    B[j+1,1]<-j
    returnobject<-as.data.frame(B)
    names(returnobject)<- c("lambda","mse")
  }
  return(returnobject)
}


MSE_lambdaintercept<-cv5intercept(data = prostate, 
                dep = "lpsa",
                indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                lambda = 0:50)

ggplot(MSE_lambdaintercept, aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point()

