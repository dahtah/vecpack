#' @importFrom purrr map map_dbl map_lgl map_df map2
#' @importFrom magrittr "%>%"
NULL


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
##' @param guess: list of initial parameters
##' @param fun: a cost function (the names of the arguments should be the *same ones* as in the "guess" argument)
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
##' @export
vpoptim <- function(guess,cfun,...)
{
    up <- gen.vunpack(guess) 
    cfun.lifted <- lift_up(cfun,up) 
    opt <- optim(vpack(guess),cfun.lifted,...)
    opt$par <- up(opt$par)
    opt
}
