## ------------------------------------------------------------------------
l <- list(a= matrix(1:4,2,2),b = 0.4)
sapply(l,as.vector) %>% do.call(c,.)

## ------------------------------------------------------------------------
vpack(l)

## ------------------------------------------------------------------------
vunpack <- gen.vunpack(l)
vunpack(1:5)
vpack(l) %>% vunpack

## ------------------------------------------------------------------------
x <- seq(0,1,l=20)
y <- 3*x - 2 +rnorm(length(x))
cfun <- function(a,b) sum((y-(a*x+b))^2)

## ------------------------------------------------------------------------
guess <- list(a=0,b=1)
up <- gen.vunpack(guess) #Unpacking function
cfun.lifted <- lift_up(cfun,up) #Now cfun accepts a vector argument
opt <- optim(vpack(guess),cfun.lifted)
up(opt$par)

## ------------------------------------------------------------------------
opt <- vpoptim(guess,cfun)
opt$par

## ------------------------------------------------------------------------
X <- cov(USArrests)
cfun <- function(A,sigma2)
{
    lrank <- tcrossprod(A,A) + sigma2*diag(nrow(A)) 
    sum((X-lrank)^2)
}
guess <- list(A=matrix(0,4,2),sigma2=1)
opt <- vpoptim(guess,cfun)
opt$par

