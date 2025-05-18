#####################################################################
## Compute tidy versions of ks function
#####################################################################

## default auxiliary function to compute ks object
## y = data frame
compute.ks <- function(y, args_ks, fun_ks)
{
    x <- as.matrix(dplyr::ungroup(y))
    if (length(args_ks)==0) args_ks <- list(x=x)
    else args_ks <- c(list(x=x), args_ks)
    fhat <- do.call(fun_ks, args_ks)
    if (is.matrix(fhat[["x"]])) 
        if (ncol(fhat[["x"]])==1) fhat[["x"]] <- as.vector(fhat[["x"]])
    fhat[["group"]] <- dplyr::group_vars(y)

    return(fhat)
}

## default auxiliary function to convert ks object to tidy format
## x=untidy ks object or tidy ks object
compute.tidy <- function(x, y, d, tidy)
{
    if (missing(tidy)) tidy <- inherits(x, "tbl")
    if (tidy) x <- getElement(dplyr::pull(x, "ks"),1)
    d <- as.integer(d)
    if (d==1)
        gg <- data.frame(x=x$eval.points, estimate=x$estimate)
    else if (d==2)
        gg <- data.frame(expand.grid(x=x$eval.points[[1]], y=x$eval.points[[2]]), estimate=as.vector(x$estimate))
    else if (d>2)
        stop("Data with d>2 not yet implemented.")

    ## keep copy of ks object in first row of ks column
    ## otherwise the data dimension
    gg$ks <- list(d)
    if (missing(y)) gg$ks[[1]] <- x else gg$ks[[1]] <- y
    gg <- dplyr::as_tibble(gg)

    return(gg)
}

## default function to compute tidy versions of ks::k* functions
## data = data frame
tidy_ks <- function(data, fun_ks, rename=TRUE, label, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))

    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(.data=data, .f=~dplyr::tibble(ks=list(compute.ks(.x, args_ks=args_ks, fun_ks=fun_ks))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy(.x, d=d))

    ## add label for type of kernel estimate
    if (missing(label)) label <- fun_ks
    gg <- dplyr::mutate(gg, tks=fun_ks, label=label, .after="ks")

    ## rename gg$x, gg$y to data variable names
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)

    ## move grouping variable to last column
    gg <- move_group_vars(gg, y=data)

    ## create "tidy_ks" class
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kde
## data = data frame
tidy_kde <- function(data, ...)
{
    ## compute KDE
    gg <- tidy_ks(data=data, fun_ks="kde", label="Density", ...)

    return(gg)
}

## tidy version of ks::kdcde
## data = data frame
tidy_kdcde <- function(data, ...)
{
    ## compute deconvolved KDE
    gg <- tidy_ks(data=data, fun_ks="kdcde", label="Deconv density", ...)

    return(gg)
}

## tidy version of ks::kcde
## data = data frame
tidy_kcde <- function(data,  ...)
{
    ## compute KCDE
    gg <- tidy_ks(data=data, fun_ks="kcde", label="Distribution", ...)

    return(gg)
}

## tidy version of ks::kcopula
## data = data frame
tidy_kcopula <- function(data, ...)
{
    ## compute KCDE
    gg <- tidy_ks(data=data, fun_ks="kcopula", label="Copula", ...)

    return(gg)
}

## tidy version of ks::kde.boundary
## data = data frame
tidy_kde_boundary <- function(data,  ...)
{
    ## compute boundary corrected KDE
    gg <- tidy_ks(data=data, fun_ks="kde.boundary", label="Bound density", ...)

    return(gg)
}

## tidy version of ks::kde.balloon
## data = data frame
tidy_kde_balloon <- function(data, ...)
{
    ## compute balloon variable KDE
    gg <- tidy_ks(data=data, fun_ks="kde.balloon", label="Balloon density", ...)

    return(gg)
}

## tidy version of ks::kde.sp
## data = data frame
tidy_kde_sp <- function(data, ...)
{
    ## compute sample point variable KDE
    gg <- tidy_ks(data=data, fun_ks="kde.sp", label="SP density", ...)

    return(gg)
}

## tidy version of ks::kde.truncate
## data = data frame
tidy_kde_truncate <- function(data, boundary, ...) { .tidy_kde_truncate(data, boundary, ...) }

.tidy_kde_truncate <- function(data, boundary, rename=TRUE, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    ## compute truncated KDE
    compute.ks.local <- function(y, args_ks, boundary)
    {
        x <- as.matrix(dplyr::select(dplyr::ungroup(y), dplyr::all_of(vars)))
        boundary <- as.matrix(boundary)
        if (length(args_ks)==0) args_ks <- list(x=x)
        else args_ks <- c(list(x=x), args_ks)
        fhat <- do.call("kde", args_ks)
        fhat <- do.call("kde.truncate", list(fhat=fhat, boundary=boundary))
        if (is.matrix(fhat[["x"]])) 
            if (ncol(fhat[["x"]])==1) fhat[["x"]] <- as.vector(fhat[["x"]])
        fhat[["group"]] <- dplyr::group_vars(y)
        return(fhat)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(.data=data, .f=~dplyr::tibble(ks=list(compute.ks.local(.x, args_ks=args_ks, boundary=boundary))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy(.x, d=d))
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    gg <- dplyr::mutate(gg, tks="kde.truncate", label="Trunc density", .before=dplyr::last_col())
    gg <- move_group_vars(gg, y=data)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kda
## data = data frame
tidy_kda <- function(data, ...) { .tidy_kda(data, ...) }

.tidy_kda <- function(data, rename=TRUE, ...)
{
    if (is.vector(data)) data <- data.frame(data)
    if (length(dplyr::group_vars(data))==0)
        stop("Requires input data to be grouped")
    else if (!all(is.factor(dplyr::pull(data, dplyr::group_vars(data)))))
        stop("Requires grouping variable to be a factor")

    ## compute KDA
    compute.ks.local <- function(y, args_ks)
    {
        x <- as.matrix(dplyr::select(dplyr::ungroup(y), dplyr::all_of(vars)))
        x.group <- dplyr::pull(y, dplyr::group_vars(y))
        if (length(args_ks)==0) args_ks <- list(x=x, x.group=x.group)
        else args_ks <- c(list(x=x, x.group=x.group), args_ks)

        fhat <- do.call("kda", args_ks)
        return(fhat)
    }

    ## convert to tidy
    compute.tidy.local <- function(y, d, rename)
    {
        j <- match(dplyr::pull(y, dplyr::group_vars(gg1)), levels(dplyr::pull(y, dplyr::group_vars(gg1))))
        x <- dplyr::pull(y, "ks")[[1]]
        xj <- x

        ## replace KDE with weighted KDE
        pp <- xj$prior.prob[j]
        xj$estimate <- x$estimate[[j]]*pp
        xj$H <- x$H[[j]]
        fhat.tidy <- dplyr::mutate(compute.tidy(x=xj, y=x, d=d, tidy=FALSE), prior.prob=pp, rn=dplyr::row_number())
        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg1 <- dplyr::summarise(data, .groups="keep")
    gg1 <- dplyr::mutate(gg1, ks=list(compute.ks.local(y=data, args_ks=args_ks)))
    gg <- dplyr::group_modify(gg1, .f=~compute.tidy.local(.x,d=d), .keep=TRUE)
    gg <- dplyr::rename(gg, prior_prob=dplyr::all_of("prior.prob"))

    ## create tidy version of partition (x_group)
    eval.points.group <- dplyr::summarise(dplyr::group_by(gg, .data$rn), label=levels(dplyr::pull(gg1,dplyr::group_vars(gg1)))[which.max(.data$estimate)], .groups="keep")
    gg <- dplyr::left_join(gg, eval.points.group, by="rn")
    gg <- dplyr::select(gg, -dplyr::one_of("rn"))
    ggv <- dplyr::group_vars(gg)[1]
    gg <- dplyr::mutate(gg, label=ifelse(.data$label==!!sym(ggv), .data$label, NA))

    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    gg <- dplyr::relocate(gg, "prior_prob", .after="estimate")
    gg <- dplyr::mutate(gg, tks="kda", .before="label")
    gv <- dplyr::group_vars(gg1)[1]
    gg$label <- factor(gg$label, levels=levels(getElement(gg1,gv)))
    gg <- dplyr::relocate(gg, dplyr::group_vars(gg), .after=dplyr::last_col())
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kdde
## data = data frame
tidy_kdde <- function(data, deriv_order=1, ...) { .tidy_kdde(data, deriv_order=deriv_order, ...) }

.tidy_kdde <- function(data, deriv_order=1, rename=TRUE, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    ## compute KDDE
    compute.ks.local <- function(y, args_ks)
    {
        x <- as.matrix(dplyr::ungroup(y))
        if (length(args_ks)==0) args_ks <- list(x=x, deriv.order=deriv_order)
        else args_ks <- c(list(x=x, deriv.order=deriv_order), args_ks)
        fhat <- do.call("kdde", args_ks)
        if (is.matrix(fhat[["x"]])) 
            if (ncol(fhat[["x"]])==1) fhat[["x"]] <- as.vector(fhat[["x"]])
        return(fhat)
    }

    ## convert to tidy
    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        deriv.order <- as.integer(x$deriv.order)
        deriv.ind <- x$deriv.ind
        if (d==1)
        {
            x$deriv.ind.scalar <- 1
            fhat.tidy <- dplyr::mutate(compute.tidy(x=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=1L, deriv_group=deriv.order)
        }
        else if (d==2)
        {
            for (j in 1:nrow(deriv.ind))
            {
                xj <- x
                x$deriv.ind.scalar <- j
                deriv.indj <- paste0("(",paste0(deriv.ind[j,], collapse=","),")")
                xj$estimate <- x$estimate[[j]]
                if (j==1) fhat.tidy <- dplyr::mutate(compute.tidy(xj, y=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=j, deriv_group=deriv.indj)
                else fhat.tidy <- rbind(fhat.tidy, dplyr::mutate(compute.tidy(xj, y=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=j, deriv_group=deriv.indj))
            }
        }

        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(data, ~dplyr::tibble(ks=list(compute.ks.local(y=.x, args_ks=args_ks))))
    gg <- dplyr::group_modify(gg, ~compute.tidy.local(y=.x, d=d))
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    gg$deriv_group <- factor(paste("deriv", gg$deriv_group))
    gg <- dplyr::mutate(gg, tks="kdde", label="Density deriv", .before=dplyr::last_col())
    gg <- move_group_vars(gg, data)
    gg <- dplyr::relocate(gg, "deriv_group", .after=dplyr::last_col())
    gg <- dplyr::relocate(gg, "deriv_order", "deriv_ind", .after="estimate")
    gg <- dplyr::group_by(gg, .data$deriv_group, .add=TRUE)
    gg$deriv_group <- factor(gg$deriv_group , levels=rev(levels(gg$deriv_group)))
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kde_local_test
## data1, data2 = data frames
tidy_kde_local_test <- function(data1, data2, labels, ...) 
{ 
    tklt <- .tidy_kde_local_test(data1=data1, data2=data2, rename=TRUE, ...) 
    if (missing(labels)) 
    {
        mc <- as.character(match.call())[2:3]
        labels <- gsub("f1", mc[1], levels(tklt$label))
        labels <- gsub("f2", mc[2], labels)
    }
    else labels <- c(paste0(labels[1],"<",labels[2]), NA, paste0(labels[1],">",labels[2]))
    levels(tklt$label) <- labels

    count <- dplyr::summarise(tklt, count=dplyr::n())
    count <- head(c(1, cumsum(count$count)+1), n=-1)
    for (i in count)
    {
        fhat <- tklt$ks[[i]]
        fhat$label <- levels(tklt$label)
        tklt$ks[[i]] <- fhat
    }

    return(tklt)
}

.tidy_kde_local_test <- function(data1, data2, rename=TRUE, ...)
{
    if (is.vector(data1)) data1 <- data.frame(data1)
    if (is.vector(data2)) data2 <- data.frame(data2)

    ## compute KDDE
    compute.ks.local <- function(y, args_ks)
    {
        y1 <- dplyr::filter(y, .data$.group=="data1")
        y2 <- dplyr::filter(y, .data$.group=="data2")
        y1 <- dplyr::select(y1, -dplyr::all_of(".group"))
        y2 <- dplyr::select(y2, -dplyr::all_of(".group"))
        x1 <- as.matrix(dplyr::select(dplyr::ungroup(y1), dplyr::all_of(vars)))
        x2 <- as.matrix(dplyr::select(dplyr::ungroup(y2), dplyr::all_of(vars)))
        if (ncol(x1)==1) x1 <- as.vector(x1)
        if (ncol(x2)==1) x2 <- as.vector(x2)
        if (length(args_ks)==0) args_ks <- list(x1=x1, x2=x2)
        else args_ks <- c(list(x1=x1, x2=x2), args_ks)
        fhat <- do.call("kde.local.test", args_ks)
        fhat$names <- names(y1)
        return(fhat)
    }

    ## convert to tidy
    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat.diff <- x$fhat1
        fhat.diff$estimate <- x$fhat.diff
        fhat.diff.lab <- x$fhat1
        fhat.diff.lab$estimate <- x$fhat.diff.pos$estimate - x$fhat.diff.neg$estimate
        fhat.tidy <- compute.tidy(fhat.diff, d=d, tidy=FALSE)
        diff.lab <- compute.tidy(fhat.diff.lab, d=d, tidy=FALSE)$estimate
        diff.lab <- factor(diff.lab, levels=c(-1,0,1), labels=c("f1<f2", NA, "f1>f2"))
        fhat.tidy <- dplyr::mutate(fhat.tidy, label=diff.lab)
        
        fhat.tidy[1,]$ks <- list(x)
        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data1), dplyr::groups(data1))
    d <- ncol(dplyr::select(dplyr::ungroup(data1), dplyr::all_of(vars)))
    ddata <- dplyr::bind_rows(dplyr::mutate(data1, .group="data1"), dplyr::mutate(data2, .group="data2"))
    gg <- dplyr::group_modify(ddata, ~dplyr::tibble(ks=list(compute.ks.local(y=.x, args_ks=args_ks))))
    gg <- dplyr::group_modify(gg, ~compute.tidy.local(.x, d=d))
    gg <- dplyr::mutate(gg, tks="kde.loctest", .after="ks")
    gg <- dplyr::relocate(gg, "ks", .after="estimate")
    if (rename) gg <- rename_ks(data=data1, gg=gg, d=d)
    gg <- move_group_vars(gg, data1)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kroc
## data1, data2 = data frames
tidy_kroc <- function(data1, data2, ...)
{
    if (is.vector(data1)) data1 <- data.frame(data1)
    if (is.vector(data2)) data2 <- data.frame(data2)

    ## compute KDDE
    compute.ks.local <- function(y1, y2, args_ks)
    {
        x1 <- as.matrix(dplyr::select(dplyr::ungroup(y1), dplyr::all_of(vars)))
        x2 <- as.matrix(dplyr::select(dplyr::ungroup(y2), dplyr::all_of(vars)))
        if (length(args_ks)==0) args_ks <- list(x1=x1, x2=x2)
        else args_ks <- c(list(x1=x1, x2=x2), args_ks)
        fhat <- do.call("kroc", args_ks)
        return(fhat)
    }

    ## convert to tidy
    compute.tidy.local <- function(y, d=1)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat.tidy <- compute.tidy(x,d=d)
        fhat.tidy[1,]$ks <- list(x)

        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data1), dplyr::groups(data1))
    gg <- dplyr::group_modify(data1, ~dplyr::tibble(ks=list(compute.ks.local(y1=.x, y2=data2, args_ks=args_ks))))
    gg <- dplyr::group_modify(gg, ~compute.tidy(.x, d=1L))
    gg <- dplyr::rename(gg, fpr=dplyr::all_of("x"))
    gg <- dplyr::mutate(gg, tks="kroc", label="ROC", .before=dplyr::last_col())
    gg <- move_group_vars(gg, data1)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kfs
## data = data frame
tidy_kfs <- function(data, ...) 
{ 
    labels <- as.character(match.call())[2]
    .tidy_kfs(data, labels=labels, ...) 
}

.tidy_kfs <- function(data, rename=TRUE, labels, ...)
{
    ## convert to tidy
    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat2 <- x
        if (d==1)
            fhat2$estimate <- ks::kde(x=x$x, eval.points=x$eval.points)$estimate
        else
            fhat2$estimate <- ks::kde(x=x$x, eval.points=expand.grid(x$eval.points))$estimate
        class(fhat2) <- "kde"
        fhat.tidy <- compute.tidy(fhat2, d=d, tidy=FALSE)
        signif.group <- factor(x$estimate, levels=c(0,1), labels=c(NA, labels)) ##factor(x$estimate, levels=c(0,1), labels=c(NA,"Signif curv"))
        fhat.tidy <- dplyr::mutate(fhat.tidy, label=signif.group)
        x$cont <- 0.5; names(x$cont) <- "50%"
        x$label <- labels #"Signif curv"
        fhat.tidy[1,]$ks <- list(x)
        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(data, ~dplyr::tibble(ks=list(compute.ks(.x, fun_ks="kfs", args_ks=args_ks))))
    gg <- dplyr::group_modify(gg, ~compute.tidy.local(.x, d=d))
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    gg <- dplyr::relocate(gg, "ks", .after="estimate")
    gg <- dplyr::mutate(gg, tks="kfs", .after="ks")
    gg <- move_group_vars(gg, data)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kdr
## data = data frame
tidy_kdr <- function(data, dTolerance, ...) 
{ 
    labels <- as.character(match.call())[2]
    .tidy_kdr(data, dTolerance=dTolerance, labels=labels, ...) 
}

.tidy_kdr <- function(data, rename=TRUE, dTolerance, labels, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    compute.tidy.local <- function(y, d)
    { 
        x <- dplyr::pull(y, "ks")[[1]]
        end.points <- data.frame(x$end.points)
        names(end.points)[1:d] <- c("x","y","z")[1:d]
        fhat.tidy <- dplyr::as_tibble(end.points)
        fhat.tidy$estimate <- 0
        fhat.tidy$ks <- list(d)
        fhat.tidy$ks[[1]] <- x
        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
    d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(.data=data, .f=~dplyr::tibble(ks=list(compute.ks(.x, args_ks=args_ks, fun_ks="kdr"))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy.local(.x, d=d))
    gv <- dplyr::group_vars(data) 
    gg <- dplyr::group_by(gg, dplyr::across("segment"), .add=TRUE) 

    ## reduce density of points on KDR
    if (missing(dTolerance))
    {
        ## median distance of consecutive point in segmentsgg
        gg3 <- dplyr::filter(gg, dplyr::n()>4)
        gg3 <- dplyr::group_modify(gg3, .f=~sf::st_as_sf(.x, coords=c("x","y")))
        gg3 <- sf::st_sf(gg3)
        gg3 <- dplyr::group_modify(gg3, .f=~dplyr::tibble(geom1=head(.x$geometry, n=-1), geom2=tail(.x$geometry, n=-1)))
        gg3 <- dplyr::summarise(gg3, dist=median(as.numeric(sf::st_distance(.data$geom1, .data$geom2, by_element=TRUE))))
        dTolerance <- median(gg3$dist)
        dTolerance <- min(10^(round(log10(mean(dTolerance))-1,0)), 10)
    }

    gg2 <- dplyr::group_modify(gg, ~data.frame(.x[1,], sf::st_sf(geometry=sf::st_sfc(sf::st_linestring(cbind(.x$x,.x$y))))))
    gg2 <- dplyr::select(gg2, -c("x","y"))
    gg2 <- sf::st_sf(gg2)
    gg2 <- gg2[sf::st_length(gg2)>0,] ## remove singleton segments
    gg2 <- sf::st_simplify(gg2, dTolerance=dTolerance)
    gg2 <- st_add_coordinates(gg2, rename=FALSE)
    gg2 <- dplyr::relocate(gg2, "X", "Y", .before=1)
    gg2 <- dplyr::rename(gg2, x="X", y="Y")
    gg <- dplyr::mutate(gg2, ks=ifelse(dplyr::row_number()>1, d, .data$ks))
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    #gg <- dplyr::mutate(gg, estimate=factor(.data$estimate), tks="kdr", label="Density ridge", .before=dplyr::last_col())
    gg <- dplyr::mutate(gg, estimate=factor(.data$estimate), tks="kdr", label=factor(labels), .before=dplyr::last_col())
    gg <- dplyr::group_by(gg, dplyr::across("segment"), .add=TRUE)
    gg <- move_group_vars(gg, data)
    gg <- dplyr::relocate(gg, "segment", .after="label") 
    sfname <- attr(gg2, "sf_column")
    gg <- dplyr::relocate(gg, dplyr::all_of(sfname), .after=dplyr::last_col())
    gg <- as_tidy_ks(gg)
    
    return(gg)
}

## tidy version of ks::kdr.segment
## data = tidy KDR object
tidy_kdr_segment <- function(data, dTolerance, coords=1:2, ...)
{
    gv <- dplyr::group_vars(data) 
    gv1 <- setdiff(gv, "segment")
    
    ## compute KDR segments
    compute.ks.local <- function(y, args_ks)
    {
        x <- untidy_ks(y)
        if (length(args_ks)==0) args_ks <- list(x=x)
        else args_ks <- c(list(x=x), args_ks)
        fhat <- do.call("kdr.segment", args_ks)

        return(fhat)
    }
    
    compute.tidy.local <- function(y, d=1)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat.tidy <- dplyr::tibble(x$end.points)
        fhat.tidy$estimate <- NA
        fhat.tidy$ks <- list(d)
        fhat.tidy$ks[[1]] <- x

        return(fhat.tidy)
    }

    if (length(gv1)>0) gg <- dplyr::group_by(data, dplyr::across(dplyr::all_of(gv1))) else gg <- dplyr::ungroup(data)
    gg <- dplyr::slice_head(gg, n=1)
    args_ks <- list(...)
    gg <- dplyr::group_modify(.data=gg, .f=~dplyr::tibble(ks=list(compute.ks.local(.x, args_ks=args_ks))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy.local(.x, d=dim_ks(data)))
    gg <- dplyr::group_by(gg, dplyr::across(dplyr::all_of(gv)))

    if (any(names(data) %in% "dTolerance")) dTolerance <- head(data$dTolerance,n=1)
    if (missing(dTolerance))
    {
        ## median distance of consecutive point in segmentsgg
        gg3 <- dplyr::filter(gg, dplyr::n()>4)
        gg3 <- dplyr::group_modify(gg3, .f=~sf::st_as_sf(.x, coords=c("x","y")))
        gg3 <- sf::st_sf(gg3)
        gg3 <- dplyr::group_modify(gg3, .f=~dplyr::tibble(geom1=head(.x$geometry, n=-1), geom2=tail(.x$geometry, n=-1)))
        gg3 <- dplyr::summarise(gg3, dist=median(geos::geos_distance(.data$geom1, .data$geom2)))
        dTolerance <- median(gg3$dist)
        dTolerance <- min(10^(round(log10(mean(dTolerance))-1,0)), 10)
    }
    
    gg <- dplyr::relocate(gg, "segment", .after="estimate")
    gg$segment <- factor(gg$segment)   
    gg <- rename_ks(data=data, gg=gg, d=dim_ks(data))
    gg <- dplyr::mutate(gg, tks="kdr", label="Density ridge", .before=dplyr::last_col())
    gg <- move_group_vars(gg, data)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks::kms
## data = data frame
tidy_kms <- function(data, ...) { .tidy_kms(data, ...) }

.tidy_kms <- function(data, rename=TRUE, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat.tidy <- data.frame(cbind(x$x, x$label))
        names(fhat.tidy)[1:d] <- c("x","y","z")[1:d]
        names(fhat.tidy)[d+1] <- "estimate"

        fhat.tidy <- dplyr::as_tibble(fhat.tidy)
        fhat.tidy <- dplyr::mutate(fhat.tidy, estimate=factor(.data$estimate), label=.data$estimate)
        fhat.tidy$ks <- list(d)
        fhat.tidy$ks[[1]] <- x
        return(fhat.tidy)
    }

    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    vars <- setdiff(names(data), dplyr::groups(data))
   	d <- ncol(dplyr::select(dplyr::ungroup(data), dplyr::all_of(vars)))
    gg <- dplyr::group_modify(.data=data, .f=~dplyr::tibble(ks=list(compute.ks(.x, args_ks=args_ks, fun_ks="kms"))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy.local(.x, d=d))
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    gg <- dplyr::mutate(gg, tks="kms", label=.data$estimate, .after="ks")
    gg <- dplyr::relocate(gg, "label", .after=dplyr::last_col())
    gg <- move_group_vars(gg, data)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of kernel support estimate
## data = tidy kde object
tidy_ksupp <- function(data, cont=95, convex_hull=TRUE, ...) { .tidy_ksupp(data, cont=cont, convex_hull=convex_hull, ...) }

.tidy_ksupp <- function(data, cont=95, convex_hull=TRUE, rename=TRUE, ...)
{
    if (is.vector(data)) data <- data.frame(data)

    compute.ks.local <- function(y, args_ks, cont, convex_hull)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        if (length(args_ks)==0) args_ks <- list(fhat=x, cont=cont, convex.hull=convex_hull)
        else args_ks <- c(list(fhat=x, cont=cont, convex.hull=convex_hull), args_ks)
        fhat <- do.call("ksupp", args_ks)
        names(fhat) <- names(data2)[1:2]
        return(fhat)
    }

    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        fhat.tidy <- dplyr::as_tibble(data.frame(x=x[,1], y=x[,2], estimate=0))
        fhat.tidy$ks <- list(d)
        fhat.tidy$ks[[1]] <- x
        return(fhat.tidy)
    }
    
    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
    d <- ncol(get_data_ks(data[1,]))
    data2 <- dplyr::select_if(dplyr::ungroup(data), ~!is.factor(.))

    gg <- dplyr::group_modify(data, ~dplyr::tibble(ks=list(compute.ks.local(y=.x, cont=cont, convex_hull=convex_hull, args_ks=args_ks))))
    gg <- dplyr::group_modify(gg, ~compute.tidy.local(.x, d=d))

    gg <- dplyr::mutate(gg, estimate=factor(.data$estimate), tks="ksupp", label=paste0(cont,"%"), .before=dplyr::last_col())
    if (rename) gg <- rename_ks(data=data2, gg=gg, d=d)
    gg <- move_group_vars(gg, data)
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of thinned quiver plot from first density derivative estimate
## data = tidy kdde object
tidy_kquiver <- function(data, thin=5, transf=1/4, neg.grad=FALSE) { .tidy_kquiver(data, thin=thin, transf=transf, neg.grad=neg.grad, rename=TRUE) }

.tidy_kquiver <- function(data, thin=5, transf=1/4, neg.grad=FALSE, rename=TRUE)
{
    ## compute quiver
    compute.tidy.local <- function(y)
    {
        fhat <- untidy_ks(y)
        while (!(class(fhat) %in% "kdde")) { fhat <- fhat[[1]] }
        if (fhat$deriv.order!=1) stop("quiver plots are defined only for 2D data with deriv_order=1.")

        ## create tidy version suitable for geom_quiver
        ev <- fhat$eval.points
        est <- fhat$estimate
        if (transf != 0) {
            est[[1]] <- sign(est[[1]]) * abs(est[[1]])^transf
            est[[2]] <- sign(est[[2]]) * abs(est[[2]])^transf
        }
        thin1.ind <- seq(1, length(ev[[1]]), by=thin)
        thin2.ind <- seq(1, length(ev[[2]]), by=thin)
        fx <- est[[1]][thin1.ind, thin2.ind]
        fy <- est[[2]][thin1.ind, thin2.ind]
        if (neg.grad) {
            fx <- -fx
            fy <- -fy
        }
        xy <- expand.grid(x=ev[[1]][thin1.ind], y=ev[[2]][thin2.ind])
        xyuv <- data.frame(xy, u=as.vector(fx), v=as.vector(fy))
        xyuv$ks <- list(2L)
        xyuv$ks[[1]] <- fhat
        return(xyuv)
    }

    ## compute grouped data frame of quiver estimate
    ## coalesce groups for each derivative into 1 group
    if (dplyr::is_grouped_df(data)) data1 <- data
    else data1 <- dplyr::group_by(data, .data$deriv_group)
    data1$.group <- dplyr::group_indices(data1)
    gg1 <- dplyr::filter(data1, .data$.group%%2 == 1)
    gg2 <- dplyr::slice_tail(gg1)
    gg1 <- dplyr::group_modify(gg1, ~compute.tidy.local(.x), .keep=TRUE)
    jv <- "deriv_group"
    if (dplyr::is_grouped_df(data)) jv <- dplyr::group_vars(data)
    data2 <- dplyr::select_if(dplyr::ungroup(data), ~!is.factor(.))
    if (rename) gg1 <- rename_ks(data=data2, gg=gg1, d=2)
    if (rename)
        gg <- dplyr::left_join(gg1, dplyr::select(gg2, -dplyr::one_of(names(data2)[1:2], "ks")), by=jv)
    else
        gg <- dplyr::left_join(gg1, dplyr::select(gg2, -dplyr::one_of("x","y", "ks")), by="group")
    if (length(jv)>1)
    {   
        jv <- setdiff(jv, "deriv_group")

        ##gg <- dplyr::group_by(gg, !!sym(jv))
        gg <- dplyr::group_by(gg, dplyr::across(dplyr::all_of(jv)))
    }
    else gg <- dplyr::ungroup(gg)
    gg <- dplyr::relocate(gg, dplyr::all_of(jv), .before=dplyr::last_col())
    gg <- dplyr::relocate(gg, "estimate", .before="u")
    gg <- dplyr::relocate(gg, "ks", .before="tks")
    gg <- dplyr::mutate(gg, tks="kquiver", label="Density quiver", .after="ks")
    gg <- dplyr::select(gg, -dplyr::one_of(c("deriv_ind","deriv_order","deriv_group",".group")))
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks:::kcurv
## data = tidy kdde object
tidy_kcurv <- function(data, ...) { .tidy_kcurv(data, ...) }

.tidy_kcurv <- function(data, rename=TRUE, ...)
{
    ## compute curvature
    compute.ks.local <- function(y)
    {
        fhat <- untidy_ks(y)
        while (!(class(fhat) %in% "kdde")) { fhat <- fhat[[1]] }
        if (fhat$deriv.order!=2) stop("Curvature plots are defined only for 2D data with deriv_order=2")
        fhat <- ks::kcurv(fhat=fhat)
        return(fhat)
    }

    ## compute grouped data frame of kernel estimate
    data <- dplyr::filter(data, .data$deriv_ind %in% 1)
    d <- 2
    gg <- dplyr::group_modify(.data=data, .f=~dplyr::tibble(ks=list(compute.ks.local(.x))))
    gg <- dplyr::group_modify(.data=gg, .f=~compute.tidy(.x, d=d))

    if (dplyr::is_grouped_df(data)) jv <- dplyr::group_vars(data)
    if (length(jv)>1)
    {   
        jv2 <- setdiff(jv, "deriv_group")
        gg <- dplyr::group_by(gg, dplyr::across(dplyr::all_of(jv2)))
    }
    else gg <- dplyr::ungroup(gg)
    gv <- head(names(dplyr::select_if(data, ~is.factor(.))),1)
    data2 <- dplyr::select_if(dplyr::ungroup(data), ~!is.factor(.))
    if (rename) gg <- rename_ks(data=data2, gg=gg, d=d)
    gg <- dplyr::mutate(gg, tks="kcurv", label="Density curv")
    if (length(jv)>1) gg <- dplyr::relocate(gg, dplyr::all_of(jv), .after=dplyr::last_col())
    gg <- dplyr::select(gg, -dplyr::all_of("deriv_group"))
    gg <- as_tidy_ks(gg)

    return(gg)
}

## tidy version of ks:::as.kde
## data = estimation grid + list of evaluation points
tidy_as_kde <- function(data, density, ...) 
{  
    is_tidy <- !any(class(data) %in% "list")
   .tidy_as_kde(data=data, density=density, is_tidy=is_tidy, ...) 
}

.tidy_as_kde <- function(data, density, rename=TRUE, is_tidy, ...)
{
    ## compute grouped data frame of kernel estimate
    args_ks <- list(...)
   
    compute.ks.local <- function(y, args_ks)
    {
        if (length(args_ks)==0) args_ks <- list(x=y, density=density)
        else args_ks <- c(list(x=y, density=density), args_ks)
        fhat <- do.call("as.kde", args_ks)
        if (is.matrix(fhat[["x"]])) 
            if (ncol(fhat[["x"]])==1) fhat[["x"]] <- as.vector(fhat[["x"]])
        
        return(fhat)
    }

    ## convert to tidy
    compute.tidy.local <- function(y, d)
    {
        x <- dplyr::pull(y, "ks")[[1]]
        deriv.order <- as.integer(x$deriv.order)
        deriv.ind <- x$deriv.ind
        if (d==1)
        {
            x$deriv.ind.scalar <- 1
            fhat.tidy <- dplyr::mutate(compute.tidy(x=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=1L, deriv_group=deriv.order)
        }
        else if (d==2)
        {
            for (j in 1:nrow(deriv.ind))
            {
                xj <- x
                x$deriv.ind.scalar <- j
                deriv.indj <- paste0("(",paste0(deriv.ind[j,], collapse=","),")")
                xj$estimate <- x$estimate[[j]]
                if (j==1) fhat.tidy <- dplyr::mutate(compute.tidy(xj, y=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=j, deriv_group=deriv.indj)
                else fhat.tidy <- rbind(fhat.tidy, dplyr::mutate(compute.tidy(xj, y=x, d=d, tidy=FALSE), deriv_order=deriv.order, deriv_ind=j, deriv_group=deriv.indj))
            }
        }

        return(fhat.tidy)
    }

    if (is_tidy)
    {
        d <- ncol(data)-1
        data <- as.data.frame(data)
    }
    else 
    {
        if (!is.list(data$eval.points)) d <- 1 else d <- length(data$eval.points)
        eval.points <- expand.grid(data$eval.points)
        if (all(is.null(names(data$eval.points)))) colnames(eval.points) <- c("x","y","z")[1:d]
        data <- data.frame(eval.points, estimate=as.vector(data$estimate))
    }
    ## reduce #s.f. to ensure exactly replicated values in data[,1:d] 
    data[,1:d] <- signif(data[,1:d],10)
    if (missing(density)) 
    {
        if (all(data$estimate>=0)) density <- TRUE
        else if (all(data$estimate<0)) density <- FALSE
        else density <- (abs(mean(data$estimate[data$estimate<0]))/mean(data$estimate[data$estimate>=0]) < 1e-6) 
    }
    
    gg <- dplyr::tibble(ks=list(compute.ks.local(data, args_ks=args_ks)))
    
    ## add label for type of kernel estimate
    if (density) 
    {
        gg <- compute.tidy(gg, d=d)
        gg <- dplyr::mutate(gg, tks="kde", label="Density", .after="ks")
    }
    else
    {
        gg <- compute.tidy.local(gg, d=d)
        gg <- dplyr::mutate(gg, tks="kdde", label="Density", .after="ks")
        gg <- dplyr::group_by(gg, .data$deriv_group, .add=TRUE)
    }
    
    ## rename gg$x, gg$y to data variable names
    if (rename) gg <- rename_ks(data=data, gg=gg, d=d)
    
    ## move grouping variable to last column
    gg <- move_group_vars(gg, y=data)

    ## create "tidy_ks" class
    gg <- as_tidy_ks(gg)

    return(gg)
}

## data = complete sf polygon grid 
## or incomplete point tidy matrix
## convert data to tidy matrix of eval points + 
## density (to be suitable input for ks::as.kde or eks::tidy_as_kde)
tidy_intergrid <- function(data, attrib, cellsize, verbose=FALSE)
{        
    ## spatial data
    if (inherits(data, "sf"))
    {    
        if (missing(attrib)) attrib <- 1
        if (is.numeric(attrib)) attrib <- names(data)[attrib]
       
        ## add grid cell indices if missing
        if (!all(names(data) %in% c("cell_id", "cell_id1", "cell_id2")))
            data <- st_add_index(data)

        ## replace polygon by centroids 
        xgrid <- suppressWarnings(st_add_coordinates(sf::st_centroid(data)))
        xgrid <- dplyr::arrange(xgrid, .data$cell_id2, .data$cell_id1)
        rownames(xgrid) <- NULL
        xgrid <- dplyr::distinct(xgrid, .data$geometry, .keep_all=TRUE)
       
        ## create tidy data frame with lon, lat + attrib
        eval.points <- sf::st_drop_geometry(xgrid[,c("lon", "lat")])
        eval.points <- signif(eval.points, 10)
        estimate <- xgrid[[attrib]]
        estimate[is.na(estimate)] <- 0
        estimate <- data.frame(eval.points, estimate) 
    }
    else ## non-spatial data
    {
        d <- ncol(data)-1
        coords <- 1:d 
        coord.names <- colnames(data)[coords]
        attrib <- d+1
        if (is.numeric(attrib)) attrib <- names(data)[attrib]

        datasf <- sf::st_as_sf(data, coords=coords) 
        data2 <- sf::st_drop_geometry(st_add_coordinates(datasf)) 
        data2.lon <- sort(unique(signif(data2$lon,10)))
        data2.lat <- sort(unique(signif(data2$lat,10)))
        if (missing(cellsize))
        {
            dlon <- signif(diff(data2.lon),4)
            dlat <- signif(diff(data2.lat),4)
            dlon <- dlon[dlon>0]
            dlat <- dlat[dlat>0]
            cellsize <- c(as.numeric(names(which.max(table(dlon)))),  as.numeric(names(which.max(table(dlat)))))
        } 
        data2.lon <- seq(min(data2.lon), max(data2.lon), by=cellsize[1]) 
        data2.lat <- seq(min(data2.lat), max(data2.lat), by=cellsize[2])
        
        xgrid <- expand.grid(lon=data2.lon, lat=data2.lat)
        xgrid <- dplyr::arrange(xgrid, .data$lat, .data$lon)
        xgrid <- sf::st_as_sf(xgrid, coords=1:2)
        rownames(xgrid) <- NULL
        xgrid <- dplyr::distinct(xgrid, .data$geometry, .keep_all=TRUE)
        xgrid <- st_add_coordinates(xgrid)

        ## join attributes from data onto xgrid
        xgrid <- dplyr::mutate(xgrid, point_id=1:nrow(xgrid), .before=1)
        xgrid <- sf::st_join(xgrid, sf::st_buffer(datasf, dist=min(cellsize)/4))
        rownames(xgrid) <- NULL
        
        ## create tidy data frame with lon, lat + attrib
        eval.points <- sf::st_drop_geometry(xgrid[,c( "lon", "lat")])
        eval.points <- signif(eval.points, 10)
        estimate <- xgrid[[attrib]]
        estimate[is.na(estimate)] <- 0
        estimate <- data.frame(eval.points, estimate) 
        colnames(estimate)[colnames(estimate) %in% c("lon","lat")] <- coord.names
    }

    if (verbose) cat(paste0("Original grid     sum(", attrib, ") = ", sum(data[[attrib]], na.rm=TRUE)), paste0("\nInterpolated grid sum(", attrib, ") = ", sum(estimate$estimate)), "\n") 
    return(estimate)
}
