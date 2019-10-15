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
  
  
  beta_ridge <- c((solve(crossprod(x,x) + lambda*diag(ncol(x)))) %*% crossprod(x,y))
  
  names(beta_ridge) <- colnames(x)
  
  return(beta_ridge)
  
  
  
}


ridge_model<-ridge(data = prostate, 
                   dep = "lpsa",
                   indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                   lambda = 1)
ridge_model


#####1.3

pred <- function(data, indep, beta_hat){
  x <- as.matrix(data[,indep]) 
  beta <- as.matrix(beta_hat)
  x <- cbind(1,x)
  preds <- x %*% beta
  return(preds)
}



prediction <- pred(data= prostate, indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                   beta_hat = ridge_model)
prediction

#### 1.4

cv1 <- function(data, dep, indep, lambda) {
  k <- max(data$group)
  A <- matrix(0, nrow = 1, ncol = k)
  for (i in 1:k){
    train <- data %>% filter(group != i)
    test <- data %>% filter(group == i)
    y_test <- as.matrix(test[,dep])
    beta_hat <- ridge(data=train, dep, indep, lambda  )
    y_hat <- pred(data = test, indep , beta_hat )
    group_mse <- mean(((y_test)-mean(y_hat))^2)
    A[,i] <- group_mse
    
  }
  mse<-sum(A)/k
  return(mse)
}

MSE<-cv1(data = prostate, 
         dep = "lpsa",
         indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
         lambda = 1)

########### 1.5

cv2 <- function(data, dep, indep, lambda) {
  B <- matrix(0, nrow = length(lambda), ncol = 2)
  for(j in lambda){
    mse<-cv1(data, dep, indep, lambda = j)
    B[j+1,2]<-mse
    B[j+1,1]<-j
    returnobject<-as.data.frame(B)
    names(returnobject)<- c("lambda","mse")
  }
  return(returnobject)
}


MSE_lambda<-cv2(data = prostate, 
                dep = "lpsa",
                indep = c("lcavol", "lweight", "age", "lbph" , "svi", "lcp", "gleason" , "pgg45"),
                lambda = 0:50)

#####1.6
ggplot(MSE_lambda, aes(x = lambda, y = mse)) +
  geom_line() +
  geom_point() +
  labs( title = "MSE as a variable dependent on lambda") +
  scale_y_continuous(name= "Mean Squared Error") +
  scale_x_continuous(name = "Value of lambda" , breaks = c(0,5,10,15,20,25,30,35,40,45,50))
