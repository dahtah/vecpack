vecpack is an experimental package for R, meant to help you stuff various things into a single vector (which is easy), and get them back in the right shape (which isn't). 

It lets you interface with optim in a more natural way: in the following example we compute a low-rank + diagonal approximation to a matrix. The most natural way of writing the cost function involves two arguments: a matrix and a scalar. Using vecpack you can send that cost function directly to optim. 

```{r}
devtools::install_github("dahtah/vecpack")
library(vecpack)
## We'll compute a rank-2 + diagonal approximation to the matrix X
X <- cor(USArrests ) 
##A cost function over two parameters: a is a scalar and B is a matrix
cfun <- function(a,B)
    {
        lowrank <- a*diag(4) + tcrossprod(B)
        sum((lowrank-X)^2)
    }
guess <- list(a=0,B=matrix(0,4,2))
#vpoptim wraps optim
res <- vpoptim(guess,cfun)
#Results are returned as a list 
with(res$par,tcrossprod(B)+a*diag(4)) %>% round(2)
```

vecpack is extensible and it should be easy enough to roll out similar interfaces for other optimisers/MCMC samplers. Have a look at the vignette for more information. 
