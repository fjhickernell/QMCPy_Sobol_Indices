modelE <- function(X){
  
  # distribution transformations
  X[,1] <- (1000-100)*X[,1]+100
  X[,2] <- exp((log(10^(-2))-log(10^(-3)))*X[,2]+log(10^-(3)))
  X[,3] <- exp((log(10^(-5))-log(10^(-6)))*X[,3]+log(10^-(6)))
  X[,4] <- exp((log(10^(-1))-log(10^(-3)))*X[,4]+log(10^-(3)))
  X[,5] <- (500-100)*X[,5]+100
  X[,6] <- (5-1)*X[,6]+1
  X[,7] <- (30-3)*X[,7]+3
  X[,8] <- exp((log(10^(-1))-log(10^(-2)))*X[,8]+log(10^-(2)))
  X[,9] <- (200-50)*X[,9]+50
  X[,10] <- (5-1)*X[,10]+1
  X[,11] <- (30-3)*X[,11]+3
  X[,12] <- exp((log(10^(7))-log(10^(5)))*X[,12]+log(10^-(5)))
  
  # call to level E model
  path = "/Users/Laurent/Desktop/levelEmodel/src"
  write.table(x=nrow(X),file = paste(path,"/try.txt",sep=""),col.names = FALSE,row.names = FALSE)
  write.table(x=X,file = paste(path,"/try.txt",sep=""),append=TRUE,col.names = FALSE,row.names = FALSE)
  system(paste(paste("cd", path,sep=" "),"; ./levele.out",sep=""),intern=TRUE)
  
  # retrieve the outputs
  y <- read.table(file = paste(path,"/LevelE.txt",sep=""))
  y <- as.matrix(y)
  
  return(y)
}

library(boot)
# designs of experiments
n <- 14350
d <- 12
X1 <- data.frame(matrix(runif(d * n), nrow = n))
X2 <- data.frame(matrix(runif(d * n), nrow = n))

# sensitivity analysis : sobolSalt
x <- sobolSalt(model = NULL, X1, X2, scheme="A", nboot = 100)
y_modelE <- modelE(X = x$X)

AS1 <- x
tell(x = AS1,y=y_modelE[,1])

AS26 <- x
tell(x = AS26,y=y_modelE[,26])


# sensitivity analysis : sobolowen

library(boot)
n <- 14350
d <- 12
X1 <- data.frame(matrix(runif(d * n), nrow = n))
X2 <- data.frame(matrix(runif(d * n), nrow = n))
X3 <- data.frame(matrix(runif(d * n), nrow = n))

# sensitivity analysis

## Not run: 
SA_owen <- sobolowen(model = NULL, X1, X2, X3, nboot = 0)
y_modelE_owen <- modelE(X = SA_owen$X)
