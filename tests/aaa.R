library(magic)
n <- 10

#first test magic() for magicness, standardness, and normality for magic(3)..magic(n):
stopifnot(is.magic(magic(3:n)))
stopifnot(is.standard(magic(3:n)))
stopifnot(is.normal(magic(3:n)))

#now test some of the specific algorithms:
stopifnot(is.magic(strachey(1:n)))
stopifnot(is.magic(lozenge (1:n)))
stopifnot(is.magic(magic.4n(1:n)))

stopifnot(sapply(1:n,function(i){is.square.palindromic(circulant(i))}))


#now test for magic.2np1() giving a generalized circulant:
stopifnot(sapply(1:n,function(i){is.circulant(magic.2np1(i)%%(2*i+1),c(2,-1))}))


## Now test that is.diagonally.correct() in fact extracts the correct elements,
## using a function that returns just the first element:
test.corners <- function(i)
{
  a <- magic(i)
  identical(a[c(1,i),c(1,i)],is.diagonally.correct(a,FUN=function(x){x[1]},g=TRUE)$diag.sums)
}

stopifnot(all(sapply(3:n,test.corners)))


## Now check that the first eigenvalue of a magic square is indeed equal to its magic constant:
f <- function(i){minmax(c(eigen(magic(i),FALSE,TRUE)$values[1],magic.constant(i)))}
stopifnot(sapply(3:n,f))

## Now check that the sum of eigenvalues 2,...,n of a magic square is zero:
f <- function(i){minmax(c(1,1+sum(eigen(magic(i),FALSE,TRUE)$values[2:i])))}
stopifnot(sapply(3:n,f))


## Check hudson() for 6n+1 and 6n-1:
stopifnot(sapply(c(6*(1:n)+1,6*(1:n)-1),function(i){is.magic(hudson(i))}))

## Check magichypercube.4n() for a range of dimensions and orders:
stopifnot(apply(expand.grid(m=1:2,d=2:4),1,function(x){is.magichypercube(magichypercube.4n(x[1],x[2]))}))

## Check magiccube.2np1():
stopifnot(sapply(1:n,function(i){is.magichypercube(magiccube.2np1( i))}))

## Sundry tests for transf;

## is transf(a,0) == a?
stopifnot(sapply(3:n , function(i){a <- magic(i);identical(a,transf(a,0))}))

## is transf(a,X) equal (ie eq()) to "a" for different X?
stopifnot(sapply(3:n , function(i){a <- magic(i);eq(a,transf(a,i%%8  ))}))
stopifnot(sapply(3:n , function(i){a <- magic(i);eq(a,transf(a,i%%8+1))}))


data(magiccubes)
stopifnot(unlist(lapply(magiccubes,is.magichypercube)))

data(Ollerenshaw)
stopifnot(is.mostperfect(Ollerenshaw))

data(cube2)
stopifnot(is.magichypercube(cube2))

data(hendricks)
stopifnot(is.perfect(hendricks))

data(perfectcube5)
stopifnot(is.perfect(perfectcube5))

data(perfectcube6)
stopifnot(is.perfect(perfectcube6))

data(Frankenstein)
stopifnot(is.perfect(Frankenstein))

#stopifnot(apply(magic.8(),3,is.magic))

## Now check magic.product() works:
f <- function(x){is.magic(magic.product(x[1],x[2]))}
stopifnot(apply(expand.grid(3:5,3:5),1,f))


## Now check some identities for adiag():
a <- matrix(1:6,2,3)
a2 <- matrix(1,2,2)
a3 <- matrix(1,3,3)

x <- 0
dim(x) <- rep(1,7)

stopifnot(identical(dim(adiag(x,x,x)),rep(3:3,7)))
stopifnot(identical(adiag(a,t(a)),t(adiag(t(a),a))))
stopifnot(identical(adiag(1,1,1,1,1),diag(5)))
stopifnot(identical(adiag(a2,a2),kronecker(diag(2),a2)))
stopifnot(identical(adiag(a3,a3,a3),kronecker(diag(3),a3)))
stopifnot(identical(adiag(matrix(1,0,5),matrix(1,5,0),pad=1:5), kronecker(t(rep(1,5)),1:5)))


#now some tests for is.circulant():
a <- array(0,rep(2,10))
a[1] <- a[1024] <- 1
stopifnot(is.circulant(a))

#"break" a by changing just one (randomly chosen) element:
a[1,1,1,1,2,1,2,1,1,1] <- 1
stopifnot(!is.circulant(a))
