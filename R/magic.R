"adiag" <-
function (..., pad = 0, do.dimnames = TRUE) 
{
    args <- list(...)
    if (length(args) == 1) {
        return(args[[1]])
    }
    if (length(args) > 2) {
        jj <- do.call("Recall", c(args[-1], list(pad = pad)))
        return(do.call("Recall", c(list(args[[1]]), list(jj), 
            list(pad = pad))))
    }
    a <- args[[1]]
    b <- args[[2]]
    if (is.null(b)) {
        return(a)
    }
    if (is.null(dim(a)) & is.null(dim(b))) {
        dim(a) <- rep(1, 2)
        dim(b) <- rep(1, 2)
    }
    if (is.null(dim(a)) & length(a) == 1) {
        dim(a) <- rep(1, length(dim(b)))
    }
    if (is.null(dim(b)) & length(b) == 1) {
        dim(b) <- rep(1, length(dim(a)))
    }
    if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
        stop("a and b must have identical number of dimensions")
    }
    seq.new <- function(i) {
        if (i == 0) {
            return(NULL)
        }
        else {
            return(1:i)
        }
    }
    s <- array(pad, dim.a + dim.b)
    s <- do.call("[<-", c(list(s), lapply(dim.a, seq.new), list(a)))
    ind <- lapply(seq(dim.b), function(i) seq.new(dim.b[[i]]) + 
        dim.a[[i]])
    out <- do.call("[<-", c(list(s), ind, list(b)))
    n.a <- dimnames(a)
    n.b <- dimnames(b)
    if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
        dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
        names(dimnames(out)) <- names(n.a)
    }
    return(out)
}
"allsubhypercubes" <-
function (a) 
{
    if (!minmax(dim(a))) {
        stop("only cubes of equal dimensions allowed")
    }
    n <- dim(a)[1]
    d <- length(dim(a))
    tri <- c("", "i", "n-i+1")
    q <- expand.grid(sapply(1:d, function(x) {
        tri
    }, simplify = FALSE))
    jj <- apply(apply(q, 2, paste), 1, paste, collapse = ",")
    wanted <- grep("i.*i", jj)
    jj <- jj[wanted]
    number.of.subhypercubes <- length(jj)
    f <- function(i, a, string) {
        n <- dim(a)[1]
        execute.string <- paste("jj <- a[", string, "]", collapse = "")
        eval(parse(text = execute.string))
        d <- log(length(jj))/log(n)
        return(array(jj, rep(n, d)))
    }
    dummy <- function(p) {
        x <- sapply(1:n, f, a = a, string = jj[p], simplify = FALSE)
        along.dim <- 1 + sum(dim(x[[1]]) > 1)
        return(do.call("abind", c(x, along = along.dim)))
    }
    out <- lapply(1:number.of.subhypercubes, dummy)
    names(out) <- jj
    return(out)
}
"allsums" <-
function (m, FUN = sum) 
{
    n <- nrow(m)
    rowsums <- apply(m, 1, FUN)
    colsums <- apply(m, 2, FUN)
    f1 <- function(i) {
        FUN(diag.off(m, i, nw.se = TRUE))
    }
    f2 <- function(i) {
        FUN(diag.off(m, i, nw.se = FALSE))
    }
    majors <- sapply(0:(n - 1), f1)
    minors <- sapply(0:(n - 1), f2)
    return(list(rowsums = rowsums, colsums = colsums, majors = majors, 
        minors = minors))
}
"ashift" <-
function (a, v) 
{
    if (is.vector(a)) {
        return(shift(a, v))
    }
    f <- function(i) {
        shift(1:(dim(a)[i]), v[i])
    }
    do.call("[", c(list(a), sapply(1:length(dim(a)), f, simplify = FALSE)))
}
"as.standard" <-
function (a) 
{
    d.a <- dim(a)
    if (max(d.a) <= 1) {
        return(a)
    }
    d <- length(d.a)
    corners <- as.matrix(do.call("expand.grid", lapply(1:d, function(i) c(1, 
        dim(a)[i]))))
    pos.of.min <- corners[which.min(a[corners]), ]
    d.a[pos.of.min > 1] <- -d.a[pos.of.min > 1]
    f <- function(n) {
        if (n > 0) {
            1:n
        }
        else {
            (-n):1
        }
    }
    a <- do.call("[", c(list(a), lapply(d.a, f)))
    return(aperm(a, order(-a[1 + diag(nrow = d)])))
}
"circulant" <-
function (n,vec=1:n)
{
   a <- matrix(0,n,n)
   out <- process(1-row(a)+col(a),n)
   out[] <- vec[out]
   return(out)
}
"diag.off" <-
function (a, offset = 0, nw.se = TRUE) 
{
    n <- dim(a)[1]
    if (nw.se == TRUE) {
        indices <- cbind(1:n, 1:n)
    }
    else {
        indices <- cbind(1:n, n:1)
    }
    jj <- process(sweep(indices, 2, c(0, offset), "+"), n)
    return(a[jj])
}
"dimension" <-
function (a, i) 
{
    p <- 1:prod(dim(a))
    n <- length(dim(a))
    d <- dim(a)[i]
    permute <- c(i, (1:n)[-i])
    a <- aperm(a, permute)
    a[] <- p
    permute[permute] <- 1:n
    return(force.integer(aperm(process(a, d), permute)))
}
"eq" <-
function (m1, m2) 
{
    all(as.standard(m1) == as.standard(m2))
}
"%eq%" <-
function (m1, m2) 
{
    eq(m1, m2)
}
"force.integer" <-
function (x) 
{
    out <- as.integer(x)
    attributes(out) <- attributes(x)
    return(out)
}
"ge" <-
function (m1, m2) 
{
    return(eq(m1, m2) || gt(m1, m2))
}
"%ge%" <-
function (m1, m2) 
{
    return(ge(m1, m2))
}
"gt" <-
function (m1, m2) 
{
    jj <- t(as.standard(m1) - as.standard(m2))
    return(ne(m1, m2) && jj[min(which(jj != 0))] > 0)
}
"%gt%" <-
function (m1, m2) 
{
    return(gt(m1, m2))
}
"hudson" <-
function (n = NULL, a = NULL, b = NULL) 
{
    if (is.null(n)) {
        n <- length(a)
    }
    if (is.null(a)) {
        a <- c(n - 1, 0:(n - 2))
    }
    if (is.null(b)) {
        b <- c(2:(n - 1), n, 1)
    }
    perm <- c(n - 1, n, 1:(n - 2))
    f <- function(i) {
        recurse(start = a, perm = perm, i)
    }
    g <- function(i) {
        recurse(start = b, perm = perm, i)
    }
    jj <- 0:(n - 1)
    aa <- t(sapply(jj, f))
    bb <- t(sapply(-jj, g))
    return(n * aa + bb)
}
"is.2x2.correct" <-
function (m, give.answers = FALSE) 
{
    window <- c(2, 2)
    sums <- subsums(m, window)
    answer <- minmax(sums)
    if (give.answers == FALSE) {
        return(answer)
    }
    else {
        return(list(answer = answer, tbt.sums = sums))
    }
}
"is.associative" <-
function (m) 
{
    is.magic(m) & minmax(c(m + rev(m), nrow(m)^2 + 1))
}
"is.square.palindromic" <-
function (m, base=10, give.answers=FALSE)
{
    n <- nrow(m)
    S <- function(i){ashift(diag(n),c(i,0))}
    P <- function(k){
      out <- matrix(0,n,n)
      out[k,k] <- 1
      return(force.integer(out))
    }
    f.maj <- function(i){
      is.persymmetric(m %*% S(i) %*% t(m))
    }

    f.min <- function(i){
      is.persymmetric(t(m) %*% S(i) %*% m)
    }

    row.sufficient <- is.persymmetric(t(m) %*% m)
    col.sufficient <- is.persymmetric(m %*% t(m))
    major.diag.sufficient <- all(sapply(1:nrow(m),f.maj))
    minor.diag.sufficient <- all(sapply(1:nrow(m),f.min))

    
    sufficient <- row.sufficient & col.sufficient & major.diag.sufficient & minor.diag.sufficient
      
    b <- base^(0:(n-1))
    R <- diag(n)[n:1,]

    is.necessary <- function(mat,tol=1e-8){
      as.vector(abs(
      (crossprod(b, R %*% mat %*% R) %*% b)/
        (crossprod(b, mat) %*% b)-1)<tol)
    }
    
    row.necessary <- is.necessary(crossprod(m,m))
    col.necessary <- is.necessary(m %*% t(m))

    f1 <- function(i) {
      diag.off(m, i, nw.se = TRUE)
    }
    f2 <- function(i) {
      diag.off(m, i, nw.se = FALSE)
    }
    m.tilde.major <- sapply(0:(n-1),f1)
    m.tilde.minor <- sapply(0:(n-1),f2)

    major.diag.necessary <-  is.necessary(crossprod(m.tilde.major %*% t(m.tilde.major)))
    minor.diag.necessary <-  is.necessary(crossprod(m.tilde.minor %*% t(m.tilde.minor)))

    necessary=row.necessary & col.necessary & major.diag.necessary & minor.diag.necessary

    if(give.answers){
      return(list(
                  necessary = necessary,
                  sufficient = sufficient,
                  row.necessary = row.necessary,
                  col.necessary = col.necessary,
                  major.diag.necessary = major.diag.necessary,
                  minor.diag.necessary = minor.diag.necessary,
                  row.sufficient = row.sufficient,
                  col.sufficient = col.sufficient,
                  major.diag.sufficient = major.diag.sufficient,
                  minor.diag.sufficient = minor.diag.sufficient
                  ))
    } else {
      return( sufficient )
    }
}  
"is.bree.correct" <-
function (m, give.answers = FALSE) 
{
    diag.dist <- nrow(m)%/%2
    offsets <- matrix(c(0, diag.dist), 2, 2)
    diag.sums <- subsums(m, offsets)
    answer <- minmax(diag.sums)
    if (give.answers == FALSE) {
        return(answer)
    }
    else {
        return(list(answer = answer, diag.sums = diag.sums))
    }
}
"is.centrosymmetric" <-
function (m)
{
    all(m==rev(m))
}
"is.circulant" <-
function(m)
{
  n <- nrow(m)
  S <- ashift(diag(n),c(1,0))
  return(all(force.integer(m%*%S)==force.integer(S%*%m)))
}
    
  
"is.diagonally.correct" <-
function (a, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    if (!minmax(dim(a))) {
        stop("only cubes of equal dimensions allowed")
    }
    n <- dim(a)[1]
    d <- length(dim(a))
    b <- c(1, -1)
    f <- function(dir) {
        (dir > 0) * (1:n) + (dir < 0) * (n:1)
    }
    g <- function(jj) {
        FUN(a[sapply(jj, f)])
    }
    ans <- do.call("expand.grid", lapply(1:d, function(...) {
        b
    }))
    diag.sums <- apply(ans, 1, g)
    dim(diag.sums) <- c(length(diag.sums)/(2^d), rep(2, d))
    if (boolean) {
        answer <- all(diag.sums)
    }
    else {
        answer <- minmax(diag.sums)
    }
    if (give.answers) {
        return(list(answer = answer, diag.sums = drop(diag.sums)))
    }
    else {
        return(answer)
    }
}
"is.latin" <-
function (m, give.answers = FALSE) 
{
    is.latinhypercube(a = m, give.answers = give.answers)
}
"is.latinhypercube" <-
function (a, give.answers = FALSE) 
{
    f <- function(x) {
        minmax(c(1, diff(sort(x))))
    }
    is.consecutive <- is.semimagichypercube(a, FUN = f, give.answers = TRUE)$rook.sums
    answer <- all(is.consecutive)
    if (give.answers) {
        return(list(answer = answer, is.consecutive))
    }
    else {
        return(answer)
    }
}
"is.magic" <-
function (m, give.answers = FALSE, FUN = sum,  boolean = FALSE) 
{
    sums <- allsums(m, FUN = FUN)
    jj <- c(sums$rowsums, sums$colsums, sums$majors[1], sums$minors[1])
    if (boolean) {
        answer <- all(jj)
    }
    else {
        answer <- minmax(jj)
    }
    if (give.answers) {
        return(c(answer = answer, sums))
    }
    else {
        return(answer)
    }
}
"is.magichypercube" <-
function (a, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    diag.sums <- is.diagonally.correct(a, give.answers = TRUE, 
        FUN = FUN, boolean = boolean)$diag.sums
    jj.semi <- is.semimagichypercube(a, give.answers = TRUE, 
        FUN = FUN, boolean = boolean)
    answer <- minmax(diag.sums) & jj.semi$answer
    if (give.answers) {
        return(list(answer = answer, rook.sums = jj.semi$rook.sums, 
            diag.sums = diag.sums))
    }
    else {
        return(answer)
    }
}
"is.mostperfect" <-
function (m, give.answers = FALSE) 
{
    if (give.answers) {
        ibc <- is.bree.correct(m, give = TRUE)
        i2c <- is.2x2.correct(m, give = TRUE)
        ipd <- is.panmagic(m, give = TRUE)
        return(list(answer = ibc$answer & i2c$answer, rowsums = ipd$rowsums, 
            colsums = ipd$colsums, majors = ipd$majors, minors = ipd$minors, 
            diag.sums = ibc$diag.sums, tbt.sums = i2c$tbt.sums))
    }
    else {
        return(is.bree.correct(m) & is.2x2.correct(m))
    }
}
"is.normal" <-
function (m) 
{
    minmax(c(1, diff(sort(m))))
}
"is.ok" <-
function (vec, n = length(vec), d = 2) 
{
    return(sum(vec) == magic.constant(n, d = d))
}
"is.panmagic" <-
function (m, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    sums <- allsums(m, FUN = FUN)
    jj <- c(sums$rowsums, sums$colsums, sums$majors, sums$minors)
    if (boolean) {
        answer <- all(jj)
    }
    else {
        answer <- minmax(jj)
    }
    if (give.answers) {
        return(c(answer = answer, sums))
    }
    else {
        return(answer)
    }
}
"is.perfect" <-
function (a, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    n <- dim(a)[1]
    d <- length(dim(a))
    putative.magic.constant <- FUN(do.call("[", c(list(a), alist(a = )$a, 
        rep(1, d - 1))))
    jj.is.ok <- function(jj, jj.give) {
        if (length(dim(jj)) == 1) {
            if (boolean) {
                return(FUN(jj))
            }
            else {
                if (jj.give) {
                  return(FUN(jj))
                }
                else {
                  return(FUN(jj) == putative.magic.constant)
                }
            }
        }
        else {
            return(is.semimagichypercube(jj, FUN = FUN, boolean = boolean, 
                give.answers = jj.give))
        }
    }
    semi.stuff <- is.semimagichypercube(a, give.answers = TRUE, 
        FUN = FUN, boolean = boolean)
    diag.stuff <- unlist(lapply(allsubhypercubes(a), jj.is.ok, 
        jj.give = FALSE))
    answer <- semi.stuff$answer & all(diag.stuff)
    if (give.answers) {
        diag.sums <- lapply(allsubhypercubes(a), jj.is.ok, jj.give = TRUE)
        return(list(answer = answer, rook.sums = semi.stuff$rook.sums, 
            diag.sums = unlist(diag.sums, recursive = FALSE)))
    }
    else {
        return(answer)
    }
}
"is.persymmetric" <-
function (m)
{
    jj <- m[,nrow(m):1]
    all(jj==t(jj))
}
"is.semimagic" <-
function (m, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    sums <- allsums(m, FUN = FUN)
    jj <- c(sums$rowsums, sums$colsums)
    if (boolean) {
        answer <- all(jj)
    }
    else {
        answer <- minmax(jj)
    }
    if (give.answers) {
        return(c(answer = answer, sums))
    }
    else {
        return(answer)
    }
}
"is.semimagichypercube" <-
function (a, give.answers = FALSE, FUN = sum, boolean = FALSE) 
{
    d <- length(dim(a))
    f <- function(i) {
        apply(a, (1:d)[-i], FUN)
    }
    jj <- sapply(1:d, f)
    if (minmax(dim(a))) {
        dim(jj) <- c(dim(a)[-1], d)
        if (boolean) {
            answer <- all(jj)
        }
        else {
            answer <- minmax(jj)
        }
    }
    else {
        if (boolean) {
            answer <- all(unlist(jj))
        }
        else {
            answer <- minmax(unlist(jj))
        }
    }
    if (give.answers) {
        return(list(answer = answer, rook.sums = jj))
    }
    else {
        return(answer)
    }
}
"is.standard" <-
function (a) 
{
    minmax(a - as.standard(a))
}
"le" <-
function (m1, m2) 
{
    return(eq(m1, m2) || lt(m1, m2))
}
"%le%" <-
function (m1, m2) 
{
    return(le(m1, m2))
}
"lozenge" <-
function (m) 
{
    n <- 2 * m + 1
    out <- matrix(NA, n, n)
    jj <- cbind(m:-m, 0:(2 * m)) + 1
    odd.a <- jj[1:(1 + m), ]
    odd.b <- odd.a
    odd.b[, 2] <- odd.b[, 2] + 1
    odd.b <- odd.b[-(m + 1), ]
    odd.coords <- rbind(odd.a, odd.b)
    even.a <- jj[(m + 2):(2 * m + 1), ]
    even.b <- jj[(m + 1):(2 * m + 1), ]
    even.b[, 2] <- even.b[, 2] + 1
    even.coords <- rbind(even.a, even.b)
    f <- function(a, x) {
        x + a
    }
    all.odd.coords <- do.call("rbind", sapply(0:m, f, x = odd.coords, 
        simplify = FALSE))
    all.even.coords <- do.call("rbind", sapply(0:m, f, x = even.coords, 
        simplify = FALSE))
    all.even.coords <- process(all.even.coords, n)
    diam.odd <- 1:(1 + 2 * m * (1 + m))
    out[all.odd.coords[diam.odd, ]] <- 2 * diam.odd - 1
    diam.even <- 1:(2 * m * (1 + m))
    out[all.even.coords[diam.even, ]] <- 2 * diam.even
    return(force.integer(out))
}
"lt" <-
function (m1, m2) 
{
    jj <- t(as.standard(m1) - as.standard(m2))
    return(ne(m1, m2) && jj[min(which(jj != 0))] < 0)
}
"%lt%" <-
function (m1, m2) 
{
    return(lt(m1, m2))
}
"magic" <-
function (n) 
{
    n <- round(n)
    if (n == 2) {
        stop("Normal magic squares of order 2 do not exist")
    }
    if (n%%2 == 1) {
        return(as.standard(magic.2np1(floor(n/2))))
    }
    if (n%%4 == 0) {
        return(as.standard(magic.4n(round(n/4))))
    }
    if (n%%4 == 2) {
        return(as.standard(magic.4np2(round((n - 2)/4))))
    }
    stop("This cannot happen")
}
"magic.2np1" <-
function (m, ord.vec = c(-1, 1), break.vec = c(1, 0), start.point = NULL) 
{
    n <- 2 * m + 1
    if (is.null(start.point)) {
        start.row <- 0
        start.col <- n + 1
    }
    else {
        start.row <- start.point[1] - 1
        start.col <- m + start.point[2] + 1
    }
    f <- function(n) {
        ord.row <- seq(from = start.row, by = ord.vec[1], length = n)
        ord.col <- seq(from = start.col, by = ord.vec[2], length = n)
        out <- cbind(rep(ord.row, n) - (n - 1), rep(ord.col, 
            n) + m)
        break.row <- ord.vec[1] - break.vec[1]
        break.col <- ord.vec[2] - break.vec[2]
        adjust <- cbind(rep(seq(from = 0, by = break.row, len = n), 
            each = n), rep(seq(from = 0, by = break.col, len = n), 
            each = n))
        return(process(out - adjust, n))
    }
    a <- matrix(NA, n, n)
    a[f(n)] <- 1:(n * n)
    return(a)
}
"magic.4n" <-
function (m) 
{
    n <- 4 * m
    a <- matrix(1:(n^2), n, n)
    jj.1 <- kronecker(diag(2), matrix(1, 2, 2))
    jj <- as.logical(kronecker(matrix(1, m + 1, m + 1), jj.1)[2:(n + 
        1), 2:(n + 1)])
    a[jj] <- rev(a[jj])
    return(force.integer(a))
}
"magic.4np2" <-
function (m) 
{
    n <- 4 * m + 2
    f <- function(n) {
        if (n == 1) {
            return(matrix(c(4, 2, 1, 3), 2, 2))
        }
        if (n == 2) {
            return(matrix(c(1, 2, 4, 3), 2, 2))
        }
        if (n == 3) {
            return(matrix(c(1, 3, 4, 2), 2, 2))
        }
        return(NULL)
    }
    lux.n <- function(m) {
        lux <- matrix(1, 2 * m + 1, 2 * m + 1)
        lux[(m + 2), ] <- 2
        if (m > 1) {
            lux[(m + 3):(2 * m + 1), ] <- 3
        }
        lux[m + 1, m + 1] <- 2
        lux[m + 2, m + 1] <- 1
        return(lux)
    }
    i <- function(a, r) {
        jj <- which(a == r, arr.ind = TRUE)
        indices <- (cbind(jj[, 1] + (jj[, 3] - 1) * 2, jj[, 2] + 
            (jj[, 4] - 1) * 2))
        o <- order(indices[, 1] * nrow(jj) + indices[, 2])
        return(indices[o, ])
    }
    a <- apply(lux.n(m), 1:2, FUN = f)
    dim(a) <- c(2, 2, 2 * m + 1, 2 * m + 1)
    out <- matrix(NA, n, n)
    sequ <- as.vector(t(magic.2np1(m))) * 4 - 4
    out[i(a, 1)] <- sequ + 1
    out[i(a, 2)] <- sequ + 2
    out[i(a, 3)] <- sequ + 3
    out[i(a, 4)] <- sequ + 4
    return(force.integer(out))
}
"magic.8" <-
function (...) 
{
    f <- function(...) {
        0:1
    }
    j <- array(t(do.call("expand.grid", lapply(1:16, FUN = f))), 
        c(4, 4, 65536))
    all.rowsums.eq.2 <- apply(apply(j, c(1, 3), sum) == 2, 2, 
        all)
    all.colsums.eq.2 <- apply(apply(j, c(2, 3), sum) == 2, 2, 
        all)
    both.sums.eq.2 <- all.rowsums.eq.2 & all.colsums.eq.2
    j <- j[c(1:4, 4:1), c(1:4, 4:1), both.sums.eq.2] > 0
    n <- dim(j)[3]
    magics <- array(1:64, c(8, 8, n))
    ref <- function(magics, j) {
        magics[j] <- rev(magics[j])
        return(magics)
    }
    f <- function(i) {
        ref(magics[, , i], j[, , i])
    }
    return(array(sapply(1:n, f), c(8, 8, n)))
}
"magic.constant" <-
function (n, d = 2, start = 1) 
{
    if (is.array(n)) {
        return(Recall(n = dim(n)[1], d = length(dim(n))))
    }
    n * (start + (n^d - 1)/2)
}
"magiccube.2np1" <-
function (m) 
{
    n <- 2 * m + 1
    jj <- array(1:n, rep(n, 3))
    x <- dimension(jj, 1)
    y <- dimension(jj, 2)
    z <- dimension(jj, 3)
    return(force.integer(((x - y + z - 1) - n * floor((x - y + z - 1)/n)) * 
                         n * n + ((x - y - z) - n * floor((x - y - z)/n)) * n + 
                         ((x + y + z - 2) - n * floor((x + y + z - 2)/n)) + 1))
}
"magichypercube.4n" <-
function (m, d = 3) 
{
    n <- 4 * m
    a <- array(0, rep(2, d))
    jj.f <- function(i) {
        dimension(a, i)
    }
    x <- apply(sapply(1:d, jj.f, simplify = TRUE), 1, sum)
    dim(x) <- rep(2, d)
    a[x%%2 == 1] <- 1
    i <- kronecker(array(1, rep(m + 1, d)), kronecker(a, array(1, 
        rep(2, d)))) == 1
    i <- do.call("[", c(list(i), lapply(1:d, function(jj.i) {
        2:(n + 1)
    })))
    j <- array(1:(n^d), rep(n, d))
    j[i] <- rev(j[i])
    return(j)
}
"magicplot" <-
function (m, number = TRUE, do.circuit = FALSE, ...) 
{
    par(pch = 16)
    n <- nrow(m)
    jj <- sort(t(m[n:1, ]), index.return = TRUE)$ix
    x <- process(jj, n)
    y <- (jj - 1)%/%n
    par(pty = "s", xaxt = "n", yaxt = "n")
    plot(x, y, type = "n", asp = 1, xlab = "", ylab = "", frame = FALSE)
    if (number == TRUE) {
        text(x, y, as.character(1:(n * n)))
        if (missing(...)) {
            points(x, y, type = "l")
        }
        else {
            points(x, y, cex = 0, ...)
        }
    }
    else {
        if (missing(...)) {
            points(x, y, type = "o")
        }
        else {
            points(x, y, ...)
        }
    }
    if (do.circuit == TRUE) {
        lines(c(x[1], x[n * n]), c(y[1], y[n * n]), ...)
    }
}
"magic.prime" <-
function (n, i = 2, j = 3) 
{
    a <- matrix(0, n, n)
    return(force.integer(n * (col(a) - i * row(a) + i - 1)%%n + 
        (col(a) - j * row(a) + j - 1)%%n + 1))
}
"magic.product" <-
function (a, b, mat = NULL) 
{
    if (length(a) == 1) {
        a <- magic(a)
    }
    if (length(b) == 1) {
        b <- magic(b)
    }
    if (is.null(mat)) {
        mat <- a * 0
    }
    if (any(dim(mat) != dim(a))) {
        stop("third argument must be same size as a")
    }
    ra <- nrow(a)
    ca <- ncol(a)
    rb <- nrow(b)
    cb <- ncol(b)
    aa <- a
    aa[aa] <- 1:length(a)
    out <- sapply(mat[aa], transf, a = b)
    out <- sweep(out, 2, length(b) * (0:(length(a) - 1)), FUN = "+")
    out <- out[, a]
    dim(out) <- c(rb, cb, ra, ca)
    out <- aperm(out, c(1, 3, 2, 4))
    dim(out) <- c(ra * rb, ca * cb)
    return(force.integer(out))
}
"magic.product.fast" <-
function (a, b) 
{
    if ((length(a) == 1) & (length(b) == 1)) {
        return(Recall(magic(a), magic(b)))
    }
    a.l <- nrow(a)
    b.l <- nrow(b)
    return(force.integer(b.l * b.l * (kronecker(a, matrix(1, 
        b.l, b.l)) - 1) + kronecker(matrix(1, a.l, a.l), b)))
}
"minmax" <-
function (x, tol=1e-6) 
{   
    if(is.integer(x)){
      return(identical(max(x), min(x)))
    }
    if(is.real(x)){
      return(abs(max(x)-min(x))/max(abs(x)) < tol)
    } else {
      return(
             abs(max(Re(x))-min(Re(x)))/max(abs(x)) < tol &
             abs(max(Im(x))-min(Im(x)))/max(abs(x)) < tol)
    }
    
}
"ne" <-
function (m1, m2) 
{
    any(as.standard(m1) != as.standard(m2))
}
"%ne%" <-
function (m1, m2) 
{
    return(ne(m1, m2))
}
"notmagic.2n" <-
function (m) 
{
    options(warn = -1)
    n <- 2 * m
    a <- matrix(NA, n, n)
    s <- seq(from = 2, by = 2, to = m)
    jj.down <- kronecker(rep(1, m), rbind(1:n, n:1))[, 1:m]
    jj.down[, s] <- jj.down[n:1, s]
    jj.down <- cbind(c(1:n, n:1), as.vector(jj.down))
    jj.up <- jj.down
    jj.up[, 2] <- (m + jj.up[, 2])%%n
    jj.up[jj.up == 0] <- n
    jj.both <- rbind(jj.down, jj.up)
    a[jj.both] <- 1:(n^2)
    return(a)
}
"panmagic.4" <-
function (vals = 2^(0:3)) 
{
    a <- rep(1, 2)
    S <- kronecker(a, kronecker(diag(a), t(a)))
    A <- diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 2)]
    N <- t(S)
    C <- t(A)
    jj <- array(c(S, A, N, C), c(4, 4, 4))
    return(force.integer(1 + apply(sweep(jj, 3, vals, "*"), 1:2, 
        sum)))
}
"panmagic.8" <-
function (chosen = 1:6, vals = 2^(0:5)) 
{
    a <- rep(1, 2)
    a.01 <- kronecker(matrix(1, 2, 2), kronecker(diag(a), t(a)))[c(1:4, 
        4:1), ]
    a.03 <- kronecker(t(a), diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 
        2)])[c(1:4, 4:1), ]
    a.05 <- kronecker(a, kronecker(kronecker(a, kronecker(diag(a), 
        t(a))), t(a)))
    a.07 <- kronecker(kronecker(a, diag(a)[c(1, 2, 1, 2), c(2, 
        1, 1, 2)]), t(a))
    a.09 <- kronecker(a, kronecker(kronecker(diag(a), t(c(a, 
        a))), a))
    a.11 <- kronecker(diag(a)[c(1, 2, 1, 2), c(2, 1, 1, 2)], 
        matrix(1, 2, 2))
    a.02 <- t(a.01)
    a.04 <- t(a.03)
    a.06 <- t(a.05)
    a.08 <- t(a.07)
    a.10 <- t(a.09)
    a.12 <- t(a.11)
    jj <- array(c(a.01, a.02, a.03, a.04, a.05, a.06, a.07, a.08, 
        a.09, a.10, a.11, a.12), c(8, 8, 12))
    jj <- jj[, , chosen, drop = FALSE]
    return(force.integer(1 + apply(sweep(jj, 3, vals, "*"), 1:2, 
        sum)))
}
"process" <-
function (x, n) 
{
    x <- x%%n
    x[x == 0] <- n
    return(x)
}
"recurse" <-
function (start = 1:n, perm, i) 
{
    invert <- function(perm) {
        perm[perm] <- 1:length(perm)
        return(perm)
    }
    n <- length(perm)
    i <- as.integer(i)
    if (i < 0) {
        return(Recall(start = start, invert(perm), -i))
    }
    perm.final <- 1:n
    while (i != 0) {
        perm.final <- perm[perm.final]
        i <- i - 1
    }
    return(start[perm.final])
}
"shift" <-
function (x, i) 
{
    n <- length(x)
    i <- i%%n
    if (i == 0) {
        return(x)
    }
    return(x[c((n - i + 1):n, 1:(n - i))])
}
"strachey" <-
function (m, square = magic.2np1(m)) 
{
    m <- round(m)
    n <- 4 * m + 2
    r <- 2 * m + 1
    out <- kronecker(matrix(c(0, 3, 2, 1), 2, 2), matrix(1, r, 
        r)) * r^2 + kronecker(matrix(1, 2, 2), square)
    coords.top <- as.matrix(expand.grid(1:r, 1:m))
    coords.top[m + 1, 2] <- m + 1
    if (m > 1) {
        coords.top <- rbind(coords.top, as.matrix(expand.grid(1:r, 
            n:(n - m + 2))))
    }
    coords.low <- sweep(coords.top, 2, c(r, 0), "+")
    jj <- out[coords.top]
    out[coords.top] <- out[coords.low]
    out[coords.low] <- jj
    return(force.integer(out))
}
"subsums" <-
function (a, p, FUN = "sum", wrap = TRUE, pad.value = 0) 
{
    if (wrap == FALSE) {
        jj <- a.block.diag(array(pad.value, p - 1), a)
        return(Recall(jj, p, FUN = FUN, pad.value = pad.value, 
            wrap = TRUE))
    }
    if (is.vector(p)) {
        sub.coords <- 1 - as.matrix(expand.grid(sapply(p, function(i) {
            1:i
        }, simplify = FALSE)))
    }
    else {
        sub.coords <- 1 - p
    }
    out <- apply(sub.coords, 1, function(v) {
        ashift(a, v)
    })
    dim(out) <- c(dim(a), nrow(sub.coords))
    if (nchar(FUN) == 0) {
        return(out)
    }
    else {
        return(apply(out, 1:length(dim(a)), FUN))
    }
    return(out)
}
"transf" <-
function (a, i) 
{
    i <- as.integer(i%%8)
    if (i%%2) {
        a <- t(a)
    }
    if ((i%/%2)%%2) {
        a <- a[nrow(a):1, ]
    }
    if ((i%/%4)%%2) {
        a <- a[, ncol(a):1]
    }
    return(a)
}
