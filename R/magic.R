"allsubhypercubes" <-
function (a) 
{
    require(abind)
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
        print(string)
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
    return(lapply(1:number.of.subhypercubes, dummy))
}
"allsums" <-
function (m) 
{
    n <- nrow(m)
    rowsums <- apply(m, 1, sum)
    colsums <- apply(m, 2, sum)
    f1 <- function(i) {
        sum(diag.off(m, i, nw.se = TRUE))
    }
    f2 <- function(i) {
        sum(diag.off(m, i, nw.se = FALSE))
    }
    majors <- sapply(0:(n - 1), f1)
    minors <- sapply(0:(n - 1), f2)
    return(list(rowsums = rowsums, colsums = colsums, majors = majors, 
        minors = minors))
}
"as.standard" <-
function (a) 
{
    n <- dim(a)[1]
    d <- length(dim(a))
    j <- c(1, n)
    repeated.string <- paste(rep(",j", d - 1), collapse = " ")
    string <- paste("corner.coords <- expand.grid(j", repeated.string, 
        ")")
    eval(parse(text = string))
    corner.coords <- as.matrix(corner.coords)
    corners <- a[corner.coords]
    pos.of.min <- rbind(corner.coords[which.min(a[corner.coords]), 
        ])
    f <- 1:n
    b <- n:1
    jj <- rep(NA, d)
    jj[pos.of.min == 1] <- "f"
    jj[pos.of.min > 1] <- "b"
    string <- paste("a <- a[", paste(jj, collapse = ","), "]")
    eval(parse(text = string))
    return(aperm(a, order(-a[1 + diag(nrow = d)])))
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
    return(aperm(process(a, d), permute))
}
"eq" <-
function (m1, m2) 
{
    all(as.standard(m1) == as.standard(m2))
}
"%eq%" <-
function (m1, m2) 
{
    all(as.standard(m1) == as.standard(m2))
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
    return((m1 %eq% m2) || (m1 %gt% m2))
}
"%ge%" <-
function (m1, m2) 
{
    return((m1 %eq% m2) || (m1 %gt% m2))
}
"gt" <-
function (m1, m2) 
{
    jj <- t(as.standard(m1) - as.standard(m2))
    return((m1 %ne% m2) && jj[min(which(jj != 0))] > 0)
}
"%gt%" <-
function (m1, m2) 
{
    jj <- t(as.standard(m1) - as.standard(m2))
    return((m1 %ne% m2) && jj[min(which(jj != 0))] > 0)
}
"is.associative" <-
function (m) 
{
    is.magic(m) & minmax(c(m + rev(m), nrow(m)^2 + 1))
}
"is.magic" <-
function (m, give.answers = FALSE) 
{
    jj <- allsums(m)
    answer <- minmax(c(jj$rowsums, jj$colsums, jj$majors[1], 
        jj$minors[1]))
    if (give.answers) {
        return(c(answer = answer, jj))
    }
    else {
        return(answer)
    }
}
"is.magichypercube" <-
function (a) 
{
    n <- dim(a)[1]
    d <- length(dim(a))
    b <- c(-1, 1)
    f <- function(dir, v) {
        if (dir > 0) {
            return(v)
        }
        else {
            return(rev(v))
        }
    }
    g <- function(jj) {
        sum(a[sapply(jj, f, v = 1:n)])
    }
    repeated.string <- paste(rep(",b", n - 1), collapse = " ")
    string <- paste("ans <- apply(expand.grid(b", repeated.string, 
        "),1,g)")
    eval(parse(text = string))
    return(minmax(ans) & is.semimagichypercube(a))
}
"is.normal" <-
function (m) 
{
    minmax(c(1, diff(sort(m))))
}
"is.panmagic" <-
function (m, give.answers = FALSE) 
{
    jj <- allsums(m)
    answer <- minmax(c(jj$rowsums, jj$colsums, jj$majors, jj$minors))
    if (give.answers == FALSE) {
        return(answer)
    }
    else {
        return(c(answer = answer, jj))
    }
}
"is.perfect" <-
function (a) 
{
    require(abind)
    n <- dim(a)[1]
    d <- length(dim(a))
    magic.constant <- apply(a, 1:(d - 1), sum)[1]
    is.ok <- function(jj) {
        if (length(dim(jj)) == 1) {
            return(sum(jj) == magic.constant)
        }
        else {
            return(is.semimagichypercube(jj) & apply(jj, 1:(length(dim(jj)) - 
                1), sum)[1] == magic.constant)
        }
    }
    return(all(unlist(lapply(allsubhypercubes(a), is.ok))))
}
"is.semimagic" <-
function (m, give.answers = FALSE) 
{
    jj <- allsums(m)
    answer <- minmax(c(jj$rowsums, jj$colsums))
    if (give.answers) {
        return(c(answer = answer, jj))
    }
    else {
        return(answer)
    }
}
"is.semimagichypercube" <-
function (a) 
{
    n <- dim(a)[1]
    d <- length(dim(a))
    f <- function(i) {
        apply(a, (1:d)[-i], sum)
    }
    return(minmax(sapply(1:d, f)))
}
"is.standard" <-
function (a) 
{
    minmax(a - as.standard(a))
}
"le" <-
function (m1, m2) 
{
    return((m1 %eq% m2) || (m1 %lt% m2))
}
"%le%" <-
function (m1, m2) 
{
    return((m1 %eq% m2) || (m1 %lt% m2))
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
    all.even.coords <- process(all.even.coords,n)
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
    return((m1 %ne% m2) && jj[min(which(jj != 0))] < 0)
}
"%lt%" <-
function (m1, m2) 
{
    jj <- t(as.standard(m1) - as.standard(m2))
    return((m1 %ne% m2) && jj[min(which(jj != 0))] < 0)
}
"magic" <-
function (n) 
{
    n <- round(n)
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
function (m, method = "a") 
{
    n <- 4 * m
    a <- matrix(1:(n^2), n, n)
    switch(method, a = {
        jj <- as.logical(kronecker(matrix(1, m + 1, m + 1), kronecker(diag(2), 
            matrix(1, 2, 2)))[2:(n + 1), 2:(n + 1)])
    }, b = {
        majors <- kronecker(1:n, cbind(0, 4 * 1:m), "+")
        majors <- process(majors, n)
        minors <- majors
        minors[, 2] <- n - minors[, 2] + 1
        minors <- process(minors, n)
        jj <- matrix(FALSE, n, n)
        jj[majors] <- TRUE
        jj[minors] <- TRUE
    })
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
    b <- 0:1
    j <- array((t(expand.grid(b, b, b, b, b, b, b, b, b, b, b, 
        b, b, b, b, b))), c(4, 4, 65536))
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
"magic.composite" <-
function (a, b) 
{
    if ((length(a) == 1) & (length(b) == 1)) {
        return(Recall(magic(a), magic(b)))
    }
    else {
        a.l <- nrow(a)
        b.l <- nrow(b)
        return(force.integer(b.l * b.l * (kronecker(a, matrix(1, 
            b.l, b.l)) - 1) + kronecker(matrix(1, a.l, a.l), 
            b)))
    }
}
"magiccube.2np1" <-
function (m) 
{
    n <- 2 * m + 1
    jj <- array(1:n, rep(n, 3))
    x <- dimension(jj, 1)
    y <- dimension(jj, 2)
    z <- dimension(jj, 3)
    return(((x - y + z - 1) - n * floor((x - y + z - 1)/n)) * 
        n * n + ((x - y - z) - n * floor((x - y - z)/n)) * n + 
        ((x + y + z - 2) - n * floor((x + y + z - 2)/n)) + 1)
}
"magichypercube.4n" <-
function (m, d = 3) 
{
    n <- 4 * m
    a <- array(0, rep(2, d))
    f <- function(i) {
        dimension(a, i)
    }
    x <- apply(sapply(1:d, f, simplify = TRUE), 1, sum)
    dim(x) <- rep(2, d)
    a[x%%2 == 1] <- 1
    i <- kronecker(array(1, rep(m + 1, d)), kronecker(a, array(1, 
        rep(2, d)))) == 1
    subs <- " , 2:(n+1)"
    string <- paste("i <- i[2:(n+1)", paste(rep(subs, d - 1), 
        collapse = ""), "]")
    eval(parse(text = string))
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
    plot(x, y, type = "n")
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
"minmax" <-
function (x) 
{
    max(x) == min(x)
}
"ne" <-
function (m1, m2) 
{
    any(as.standard(m1) != as.standard(m2))
}
"%ne%" <-
function (m1, m2) 
{
    any(as.standard(m1) != as.standard(m2))
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
    jj <- jj[, , chosen]
    if (length(chosen) == 1) {
        dim(jj) <- c(8, 8, 1)
    }
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
