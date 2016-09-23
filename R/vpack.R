#' vecpack: a package for packing things into a vector
#'
#' vecpack is an experimental package for R, meant to help you stuff various things into a single vector (which is easy), and get them back in the right shape (which isn't). 
#' It lets you interface with optim in a more natural way: in the following example we compute a low-rank + diagonal approximation to a matrix. The most natural way of writing the cost function involves two arguments: a matrix and a scalar. Using vecpack you can send that cost function directly to optim. 
#' @docType package
#' @examples
#' #Packing things into a vector is easy
#' list(A=matrix(0,2,2),b=1:3) %>% sapply(as.vector) %>% do.call("c",.)
#' #Going the other way isn't: many objects correspond to the same vector
#' #gen.vunpack generates an "unpacking" function from an example
#' l <- list(A=matrix(0,2,2),b=1:3)
#' vunpack <- gen.vunpack(l)
#' vpack(l)
#' vpack(l) %>% vunpack
#'
#' #The most common use for vpack is for interfacing with optimisation packages
#' #In this example we compute a rank-2 + diagonal approximation 
#' #to the matrix X
#' X <- cor(USArrests ) 
#' #A cost function over two parameters: a is a scalar and B is a matrix
#' cfun <- function(a,B)
#'    {
#'        lowrank <- a*diag(4) + tcrossprod(B)
#'        sum((lowrank-X)^2)
#'    }
#' guess <- list(a=0,B=matrix(0,4,2))
#' #vpoptim wraps optim
#' res <- vpoptim(guess,cfun)
#' @name vecpack
#' @author Simon Barthelme
NULL

#' @importFrom purrr map map_dbl map_lgl map_df map2
#' @import magrittr
#' @importFrom stats optim
#' @export "%>%"
NULL

#' Pipe operator from magrittr package
#'
#' The pipe operator lets you string function calls together. See help from magrittr package for more.
#' 
#' @name %>%
#' @rdname pipe
#' @param lhs data
#' @param rhs a function
NULL


utils::globalVariables(c(".", "%>%"))


##' Pack a list of objects into a vector
##'
##' Turns a list of objects into a concatenated vector. Objects must have an "as.vector" method.
##' If given a single object, vpack just runs as.vector. 
##' @param l a list of objects or a single object
##' @return a vector
##' @author Simon Barthelme
##' @seealso gen.vunpack
##' @examples
##' list(a=1:3,b=matrix(1:4,2,2)) %>% vpack
##' matrix(1:4,2,2) %>% vpack
##' @export
vpack <- function(l)
{
    if (is.list(l))
    {
        map(l,as.vector) %>% do.call(c,.)
    }
    else
    {
        as.vector(l)
    }
}
##' S3 generic for converting an object back into a vector
##'
##' When converting an object into a vector part of the meta-data is lost (for example, when converting a matrix to a vector, the dimensions of the original matrix are lost). These generic functions return a "reshaping function" appropriate for turning a vector back into the original object. They're used internally by vecpack, see Vignette.
##' @param x an object
##' @param ... ignored
##' @export 
reshapefun <- function(x,...)    UseMethod("reshapefun", x)

##' @export 
reshapefun.numeric <- function(x,...) function(k) k
##' @export 
reshapefun.integer <- function(x,...) function(k) k
##' @export 
reshapefun.vector <- function(x,...) function(k) k
##' @export 
reshapefun.matrix <- function(x,...) function(k) matrix(k,nrow(x),ncol(x))
##' @export 
reshapefun.cimg <- function(x,...) function(k) imager::as.cimg(k,dim=dim(x))
##' @export 
reshapefun.array <- function(x,...) function(k) array(k,dim(x))




##' Change function to accept a packed vector as argument
##'
##' Similar to the family of lift functions in purrr. An example works better than an explanation, so please see below. 
##' 
##' @param f a function 
##' @param up an unpacking method
##' @return a function 
##' @examples
##' f <- function(a,b) a + b #A function of two arguments
##' l <- list(a= matrix(1:4,2,2),b = 5)
##' up <- gen.vunpack(l)
##' f.lifted <- lift_up(f,up)
##' f(l$a,l$b)
##' f.lifted(vpack(l))
##' @author Simon Barthelme
##' @export
lift_up <- function(f,up) {
    if (!attr(up,"direct"))
    {
        f <- purrr::lift(f)
    }
    . %>% up %>% f
}

wrap <- function(f)
{
    function(x)
    {
        out <- vpack(x) %>% f
        gen.vunpack(x)(out)
    }
}

##' Generate an unpacking method for vectors 
##'
##' @param l a list
##' @return an unpacking method
##' @author Simon Barthelme
##' @examples
##' #A list of values that can all be coerced into a vector
##' l <- list(a= matrix(1:4,2,2),b = 5)
##' #Pack elements into a vector
##' vpack(l)
##' #For the inverse operation (getting the list back from our packed vector)
##' #we create an unpacking method
##' unpack <- gen.vunpack(l)
##' vpack(l) %>% unpack
##' @export
gen.vunpack <- function(l)
    {
        if (is.list(l))
        {
            ss <- map_dbl(l,length)
            n <- length(l)
            nm <- names(l)
            cs <- cumsum(ss)
            unvec <- function(vec)
            {
                map2(c(1,1+cs[-n]),cs,function(a,b) vec[a:b])
            }
            reshape <- map(l,reshapefun)
            vunpack <- function(vec)
            {
                l <- unvec(vec)
                names(l) <- nm
                map2(l,reshape,function(l,f) f(l))
            }
            attr(vunpack,"direct") <- FALSE
        }
        else
        {
            vunpack <- reshapefun(l)
            attr(vunpack,"direct") <- TRUE

        }
        vunpack
    }

##' An interface to optim for optimisation over multiple arguments
##'
##' vpoptim uses vpack to let you call optim on cost functions involving multiple arguments. See examples for usage.
##' @param guess list of initial parameters
##' @param cfun a cost function (the names of the arguments should be the *same ones* as in the "guess" argument)
##' @param ... passed on to optim
##' @examples
##' #A cost function in two arguments:
##' cost <- function(a,b)
##' {
##'   (3*a-b+2)^2 
##' }
##' res <- vpoptim(list(a=1,b=0),cost)
##' res$par
##' ##Won't work: the names in the list and the cost function need to be the same
##' try(res <- vpoptim(list(a=1,v=0),cost))
##' ##Use optim for low-rank + diagonal approximation
##' X <- cov(USArrests)
##' cfun <- function(A,sigma2)
##' {
##'     lrank <- tcrossprod(A,A) + sigma2*diag(nrow(A)) 
##'     sum((X-lrank)^2)
##' }
##' guess <- list(A=matrix(0,4,2),sigma2=1)
##' opt <- vpoptim(guess,cfun)
##' opt$par
##' #Extra parameters are passed on to optim:
##' opt <- vpoptim(guess,cfun,method="BFGS")
##' @export
vpoptim <- function(guess,cfun,...)
{
    up <- gen.vunpack(guess) 
    cfun.lifted <- lift_up(cfun,up) 
    opt <- optim(vpack(guess),cfun.lifted,...)
    opt$par <- up(opt$par)
    opt
}


vpoptimx <- function(guess,cfun,...)
{
    up <- gen.vunpack(guess) 
    cfun.lifted <- lift_up(cfun,up) 
    opt <- optimx(vpack(guess),cfun.lifted,...)
    n <- length(vpack(guess))
    opt$par <- apply(opt,1,function(v) up(v[1:n]))
    opt
}


##' An interface to numDeriv::grad for numerical gradients over multiple arguments
##'
##' Same as vpoptim, pretty much. See examples
##'
##' @param fun function to differentiate (the names of the arguments should be the *same ones* as in the "guess" argument)
##' @param x compute gradient at x
##' @param ... passed on to grad
##' @examples
##' #gradient of a function of two arguments
##' vpgrad(function(x,y) x+y,list(x=1,y=3))
##' #gradient of a matrix function, returned as a matrix
##' vpgrad(det,diag(2))
##' @export
vpgrad <- function(fun,x,...)
    {
        up <- gen.vunpack(x) 
        fun.lifted <- lift_up(fun,up) 
        gr <- numDeriv::grad(fun.lifted,vpack(x),...)
        up(gr)
    }
