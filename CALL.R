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
  
  return(beta_ridge)
  
  
  
}


ridge_model<-ridge(data = prostate, 
                   dep = "lpsa",
                   indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                   lambda = 1)


#####1.3

pred <- function(data, indep, beta_hat){
  x <- as.matrix(data[,indep]) 
  beta <- as.matrix(beta_hat)
  x <- cbind(1,x)
  preds <- (x %*% beta)
  return(preds)
}



prediction <- pred(data= prostate, indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                   beta_hat = ridge_model)

#### 1.4

cv4_Aintercept <- function(data, dep, indep, lambda) {
  k <- max(data$group)
  A <- matrix(0, nrow = k, ncol = 1)
  for (i in 1:k){
    train <- data %>% filter(group != i)
    test <- data %>% filter(group == i)
    y_test <- as.matrix(test[,dep])
    beta_hat <- ridge(data=train, dep, indep, lambda  )
    y_hat <- pred(data = test, indep , beta_hat )
    group_mse <- mean((as.vector(y_hat)-as.vector(y_test))^2)
    A[i,] <- group_mse
    
  }
  mse<-sum(A)/5
  return(mse)
}

MSEintercept_B<-cv4_Aintercept(data = prostate, 
                           dep = "lpsa",
                           indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                           lambda = 1)

########### 1.5

cv5intercept_B <- function(data, dep, indep, lambda) {
  B <- matrix(0, nrow = length(lambda), ncol = 2)
  for(j in lambda){
    mse<-cv4_Aintercept(data, dep, indep, lambda = j)
    B[j+1,2]<-mse
    B[j+1,1]<-j
    returnobject<-as.data.frame(B)
    names(returnobject)<- c("lambda","mse")
  }
  return(returnobject)
}


MSE_lambdaintercept_B<-cv5intercept_B(data = prostate, 
                                      dep = "lpsa",
                                      indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                                      lambda = 0:50)

#####1.6
ggplot(MSE_lambdaintercept_B, aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point()

