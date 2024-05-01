#####################################################################
## Compute sf versions of ks functions
#####################################################################

## default auxiliary function to convert ks object to sf object
## x = sf point object or tidy_ks object

st_ks <- function(x, fun_ks, crs, ...)
{
	## x = tidy_ks object
	if (inherits(x, "tidy_ks"))
    {
        fhat <- x
    }
    ## x = sf object
    else
    {
        ## extract coordinates
        if (missing(crs)) crs <- sf::st_crs(x)
        xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
        
        ## compute tidy kernel estimate via tidy_k* function
        fhat <- do.call(fun_ks, args=list(data=xcoord, rename=FALSE, ...))
    }
    fhat.tidy <- fhat
    gv <- dplyr::group_vars(x)

    ## compute slice_head version of kernel estimate for geometry conversion
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    fhat <- as_tidy_ks(fhat)

    ## convert ks grid to required geometry
    fhat.geom <- st_evalpoints(x=fhat)
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.geom <- st_contourline(x=fhat, cont=1:99)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
    cont.geom <- sf::st_sf(cont.geom, crs=crs)

    ## if original data is not tibble
    if (!inherits(x, "tbl_df")) 
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }

	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of ks::kde
## x = sf point geometry

st_kde <- function(x, ...)
{
    ## compute KDE
    fhat.sf <- st_ks(x=x, fun_ks="tidy_kde", ...)

    return(fhat.sf)
}

## sf version of ks::kdcde
## x = sf point geometry

st_kdcde <- function(x, ...)
{
    ## compute deconvolved KDE
    fhat.sf <- st_ks(x=x, fun_ks="tidy_kdcde", ...)

    return(fhat.sf)
}

## sf version of ks::kcde
## data = sf point geometry

st_kcde <- function(x, ...)
{
	## compute kernel CDE
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- tidy_kcde(xcoord, rename=FALSE, ...)
	fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))

	## convert ks grid to required geometry
    gv <- dplyr::group_vars(x)
	fhat.geom <- st_evalpoints(x=fhat)
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.geom <- st_contourline(x=fhat, cont=1:99)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
    cont.geom <- sf::st_sf(cont.geom, crs=crs)

    ## if original data is not tibble
    if (!inherits(x, "tbl_df")) 
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }
    fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
    class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of ks::kde.boundary
## x = sf point geometry

st_kde_boundary <- function(x, ...)
{
    ## compute boundary corrected KDE
    fhat.sf <- st_ks(x=x, fun_ks="tidy_kde_boundary", ...)

    return(fhat.sf)
}

## sf version of ks::kde.balloon
## x = sf point geometry

st_kde_balloon <- function(x, ...)
{
    ## compute balloon variable KDE
    fhat.sf <- st_ks(x=x, fun_ks="tidy_kde_balloon", ...)

    return(fhat.sf)
}

## sf version of ks::kde.sp
## x = sf point geometry

st_kde_sp <- function(x, ...)
{
    ## compute sample point variable KDE
    fhat.sf <- st_ks(x=x, fun_ks="tidy_kde_sp", ...)

    return(fhat.sf)
}

## sf version of ks::kde.truncate
## x = sf point geometry
## boundary = sf point/polygon geometry

st_kde_truncate <- function(x, boundary, ...)
{
	## compute boundary kernel density estimate
    xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
    fhat <- .tidy_kde_truncate(data=xcoord, boundary=sf::st_coordinates(boundary), rename=FALSE, ...)
    #fhat.tidy <- fhat
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))

    gv <- dplyr::group_vars(fhat)
    fhat.geom <- st_evalpoints(x=fhat)
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.geom <- st_contourline(x=fhat, cont=1:99)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
    cont.geom <- sf::st_sf(cont.geom, crs=crs)

    ## if original data is not tibble
    if (!inherits(x, "tbl_df"))
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }
    
	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of ks::kda
## x = sf point geometry

st_kda <- function(x, ...)
{
    ## compute kernel discriminant analysis
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- .tidy_kda(xcoord, rename=FALSE, ...)
	fhat.tidy <- fhat
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","prior_prob","estimate")))
    fhat <- dplyr::mutate(fhat, label="Classif")
    gv <- dplyr::group_vars(fhat)

    ## convert ks partition grid to required geometry
    fhat.geom <- st_evalpoints(x=fhat[1,])
    fhat.geom <- dplyr::mutate(fhat.geom, label=.data$estimate, .after="estimate")
    levels(fhat.geom$estimate) <- 1:length(levels(fhat.geom$estimate))
    fhat.geom$estimate <- as.integer(levels(fhat.geom$estimate))[fhat.geom$estimate]
    gv <- dplyr::group_vars(fhat.geom)
    fhat.geom <- dplyr::select(dplyr::ungroup(fhat.geom), -dplyr::all_of(gv))
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)

    ## convert contours to multipolygon geometry
	fhat <- dplyr::bind_cols(fhat, rn=1:nrow(fhat))
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kda(x=.x))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
	fhat <- dplyr::select(fhat, -dplyr::one_of("rn"))
	if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry") 
	cont.geom <- sf::st_sf(cont.geom, crs=crs)

	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of ks::kdde
## x = sf point geometry

st_kdde <- function(x, deriv_order=1, ...)
{
	## compute kernel density derivatives estimate
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- .tidy_kdde(xcoord, deriv_order=deriv_order, rename=FALSE, ...)
    fhat <- dplyr::slice_head(fhat)
    gv <- dplyr::group_vars(fhat)
    gv2 <- c("deriv_ind", gv)
    fhat <- dplyr::arrange(fhat, .data$deriv_ind)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    fhat <- as_tidy_ks(fhat)

	## convert ks grid to required geometry
	crs <- sf::st_crs(x) 
    fhatg <- dplyr::group_by(fhat, dplyr::across(dplyr::all_of(gv2)))
    
    fhat.geom <- st_evalpoints(x=fhatg)   
    fhat.geom <- dplyr::mutate(fhat.geom, estimate=get(paste0("estimate.", dplyr::first(.data$deriv_ind))), .before="geometry")
    fhat.geom <- dplyr::select(fhat.geom, -dplyr::starts_with("estimate."))
    fhat.geom <- dplyr::group_by(fhat.geom, dplyr::across(dplyr::all_of(gv2)))
    fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv2), .before="geometry")
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contours to multipolygon geometry
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, cont=1:99))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
	cont.geom <- dplyr::left_join(cont.geom, dplyr::select(fhat, dplyr::all_of(gv2)), by=gv)
    cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv2), .before="geometry") 
    cont.geom <- dplyr::relocate(cont.geom, "deriv_ind", .after="estimate") 
    cont.geom <- dplyr::arrange(cont.geom, .data$deriv_ind)
    cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x))

    fhat <- dplyr::select(fhat, -dplyr::one_of(c("deriv_order")))
    fhat <- as_tidy_ks(fhat)
	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

	return(fhat.sf)
}

## sf version of ks::kde_local_test
## x = sf point geometry

st_kde_local_test <- function(x1, x2, labels, ...)
{
    ## compute kernel local test
	xcoord1 <- dplyr::group_modify(.data=x1, .f=~as.data.frame(sf::st_coordinates(.x)))
	xcoord2 <- dplyr::group_modify(.data=x2, .f=~as.data.frame(sf::st_coordinates(.x)))

    fhat <- .tidy_kde_local_test(data1=xcoord1, data2=xcoord2, rename=FALSE, ...)
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    fhat <- dplyr::mutate(fhat, label="Density diff")
    gv <- dplyr::group_vars(fhat)

    ## convert ks grid to required geometry
    fhat.geom <- st_evalpoints(x=fhat)
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    crs <- sf::st_crs(x1)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
	## convert contour polygons to multipolygon geometry
	abs.cont <- c(-0.5, 0.5); names(abs.cont) <- c("50%", "50%")
    for (i in 1:nrow(fhat))
    {
        fhat$ks[[i]]$eval.points <- fhat$ks[[i]]$fhat1$eval.points
        fhat$ks[[i]]$estimate <- fhat$ks[[i]]$fhat.diff.pos$estimate - fhat$ks[[i]]$fhat.diff.neg$estimate
    }
    cont.geom <- dplyr::bind_rows(st_contourline(x=fhat, abs.cont=abs.cont[1]), st_contourline(x=fhat, abs.cont=abs.cont[2]))
    cont.geom <- dplyr::mutate(cont.geom, contlabel=factor(ifelse(.data$estimate<0, paste0("-", 100-as.numeric(as.character(.data$contlabel))), as.character(.data$contlabel))))
    cont.geom <- dplyr::mutate(cont.geom, label=.data$contlabel, .after="estimate")
    levels(cont.geom$label) <- c("f1<f2", "f1>f2")

    if (missing(labels)) mc <- as.character(match.call())[-1]
    else 
    {
        if (tolower(labels[1])=="default") mc <- c("f1", "f2")
        else mc <- labels
    }
    
    lev <- levels(cont.geom$label)
    lev <- gsub("f1", mc[1], lev) 
    lev <- gsub("f2", mc[2], lev)
    tk <- rep(NA, len=length(cont.geom$label))
    tk[cont.geom$label=="f1<f2"] <- lev[1]
    tk[cont.geom$label=="f1>f2"] <- lev[2]
    cont.geom$label <- factor(tk, levels=lev) 

    cont.geom <- dplyr::arrange(cont.geom, sign(.data$estimate), abs(.data$estimate))
    if (length(gv)>0) 
    {   
        cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
        cont.geom <- dplyr::arrange(cont.geom, dplyr::across(dplyr::all_of(gv)))
    }
    cont.geom <- sf::st_sf(cont.geom, crs=crs)

    ## if original data is not tibble
    if (!inherits(x1, "tbl_df"))
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }

	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

	return(fhat.sf)
}

## sf version of ks::kfs
## x = sf point geometry

st_kfs <- function(x, ...)
{
    ## compute kernel feature significance
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- .tidy_kfs(xcoord, rename=FALSE, ...)
	fhat.tidy <- fhat
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    fhat <- dplyr::mutate(fhat, label="Signif curv", .after="tks")

    ## convert ks grid to required geometry
    gv <- dplyr::group_vars(x)
    fhat.geom <- st_evalpoints(x=dplyr::select(fhat, -dplyr::one_of("tks", "label")))
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
   
	## convert contour polygons to multipolygon geometry
	abs.cont <- 0.5; names(abs.cont) <- "50%"
    cont.geom <- st_contourline(x=fhat, abs.cont=abs.cont)
    cont.geom <- dplyr::relocate(cont.geom, "contlabel", "estimate", .before=1)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
    cont.geom <- dplyr::mutate(cont.geom, label=fhat$label[1], .after="estimate")
    cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x))

    ## if original data is not tibble
    if (!inherits(x, "tbl_df")) 
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }
    
	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

	return(fhat.sf)
}

## sf version of ks::kdr
## x = sf point geometry

st_kdr <- function(x, dTolerance=100, ...)
{
    ## compute kernel density ridge
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- .tidy_kdr(xcoord, rename=FALSE, dTolerance=dTolerance, ...)
    
    ## convert density ridge to linestring 
    fhat.sf <- sf::st_as_sf(fhat, coords=c("x","y"), crs=sf::st_crs(x))
    fhat.sf <- dplyr::summarise(fhat.sf, do_union=FALSE, tks=head(.data$tks,n=1), label=head(.data$label,n=1), .groups="keep")
    fhat.sf <- sf::st_cast(fhat.sf, to="LINESTRING")
    fhat.sf <- dplyr::mutate(fhat.sf, contlabel=50L, estimate=0, .before="label")
    
    gv <- dplyr::group_vars(x)
    if (length(gv)>0) fhat <- dplyr::group_by(fhat, dplyr::across(dplyr::all_of(gv))) else fhat <- dplyr::ungroup(fhat)
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","segment", "estimate")))
    
    fhat.sf <- list(tidy_ks=fhat, grid=NULL, sf=fhat.sf)
    class(fhat.sf) <- "sf_ks"

    return(fhat.sf) 
}

## sf version of ks::kms
## x = sf point geometry

st_kms <- function(x, ...)
{
    ## compute kernel mean shift
	xcoord <- dplyr::group_modify(.data=x, .f=~as.data.frame(sf::st_coordinates(.x)))
	fhat <- .tidy_kms(xcoord, rename=FALSE, ...)
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    gv <- dplyr::group_vars(x)

    ## compute cluster labels
    fhat1 <- dplyr::group_modify(.data=fhat, .f=~dplyr::as_tibble(untidy_ks(.x, "x")))
    fhat2 <- dplyr::group_modify(.data=fhat, .f=~dplyr::rename(dplyr::as_tibble(untidy_ks(.x, "label")), label=dplyr::all_of("value")))
    
    fhat1 <- dplyr::mutate(fhat1, tks="kms")
    fhat1$label <- as.integer(fhat2$label)

    fhat.sf <- dplyr::mutate(fhat1, contlabel=50, contlabel=factor(.data$contlabel), estimate=.data$label, .before="label")
    fhat.sf <- dplyr::relocate(fhat.sf, "contlabel", "estimate", .before="tks")
    fhat.sf <- dplyr::mutate(fhat.sf, estimate=factor(.data$estimate), label=factor(.data$label))
    fhat.sf <- sf::st_as_sf(fhat.sf, coords=c("X","Y"), crs=sf::st_crs(x))
    if (length(gv)>0) 
    {
        fhat.sf <- dplyr::relocate(fhat.sf, dplyr::all_of(gv), .before="geometry")
        fhat.sf <- sf::st_sf(fhat.sf)
    }

    ## if original data is not tibble
    if (!inherits(x, "tbl_df")) fhat.sf <- sf::st_sf(as.data.frame(fhat.sf))

    fhat.sf <- list(tidy_ks=fhat, grid=NULL, sf=fhat.sf)
    class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of kernel support estimate
## x= sf_ks object (output from st_kde)

st_ksupp <- function(x, cont=95, convex_hull=TRUE, ...)
{
    ## compute kernel support estimate
    fhat <- x$tidy_ks 
    fhat.tidy <- .tidy_ksupp(fhat, cont=cont, rename=FALSE, convex_hull=convex_hull, ...)
    fhat.tidy <- dplyr::slice_head(fhat.tidy)
    fhat.tidy <- dplyr::select(fhat.tidy, -dplyr::one_of(c("x","y","estimate")))

    ## convert kernel support to polygon
    gv <- dplyr::group_vars(x$tidy_ks)
    ##gv <- names(x$tidy_ks)[seq_along(which(names(x$tidy_ks)=="x")-1)]
    fhat.sf <- dplyr::group_modify(.data=fhat.tidy, .f=~dplyr::tibble(geometry=list(sf::st_linestring(as.matrix(untidy_ks(.x))))))
    fhat.sf <- sf::st_sf(fhat.sf)
    fhat.sf <- sf::st_cast(fhat.sf, to="POLYGON")
    fhat.sf <- dplyr::mutate(fhat.sf, estimate=cont, contlabel=factor(cont), .before="geometry")
    fhat.sf <- dplyr::relocate(fhat.sf, "contlabel", "estimate", .before=1)
    fhat.sf <- dplyr::mutate(fhat.sf, tks="ksupp", label="Support", .before="geometry")
    if (length(gv)>0) fhat.sf <- dplyr::relocate(fhat.sf, dplyr::all_of(gv), .before="geometry")
    
    ## if original data is not tibble
    if (!dplyr::is_grouped_df(x$sf)) fhat.sf <- data.frame(fhat.sf)
    fhat.sf <- sf::st_sf(fhat.sf, crs=sf::st_crs(x$sf))

    fhat.sf <- list(tidy_ks=fhat.tidy, grid=NULL, sf=fhat.sf)
	class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## sf version of ks::kcurv
## x = sf_ks object (output from st_kdde)

st_kcurv <- function(x, ...)
{
    ## compute kernel curvature estimate
    fhat <- dplyr::slice_head(x$tidy_ks)
    fhat <- tidy_kcurv(fhat, rename=FALSE)
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    gv <- setdiff(dplyr::group_vars(x$tidy_ks), "deriv_group")

    ## convert ks grid to required geometry
    fhat.geom <- st_evalpoints(x=fhat)
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before="geometry")
    crs <- sf::st_crs(x$sf)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)

	## convert contours to multipolygon geometry
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(.st_contourline(x=untidy_ks(.x)))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
    cont.geom <- dplyr::bind_cols(label=fhat$label[1], cont.geom)
	cont.geom <- dplyr::relocate(cont.geom, "contlabel", "estimate", .before=1)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before="geometry")
    cont.geom <- sf::st_sf(cont.geom, crs=crs)

    ## if original data is not tibble
    if (length(gv)==0) 
    {
        fhat.geom <- sf::st_sf(as.data.frame(fhat.geom))
        cont.geom <- sf::st_sf(as.data.frame(cont.geom))
    }
	fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
	class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## spatial version of ks::kroc
## x1, x2 = sf point geometry

st_kroc <- function(x1, x2, ...)
{
    ## compute kernel ROC 
    xcoord1 <- dplyr::group_modify(.data=x1, .f=~as.data.frame(sf::st_coordinates(.x)))
    xcoord2 <- dplyr::group_modify(.data=x2, .f=~as.data.frame(sf::st_coordinates(.x)))
    fhat <- tidy_kroc(data1=xcoord1, data2=xcoord2, ...)

    return(fhat)
}

## spatial version of ks::kquiver
## x = output from st_kdde(, deriv_order=1)

st_kquiver <- function(x, thin=5, transf=1/4, neg.grad=FALSE, scale=1)
{ 
    ## compute kernel quiver estimate
    jv <- dplyr::group_vars(x$tidy_ks)
    t1 <- dplyr::group_modify(x$tidy_ks, .f=~data.frame(expand.grid(untidy_ks(.x)$eval.points), estimate1=as.vector(untidy_ks(.x)$estimate[[1]]), estimate2=as.vector(untidy_ks(.x)$estimate[[2]])))
    t1 <- dplyr::left_join(t1, x$tidy_ks, by=jv)
    t1 <-dplyr::mutate(t1, rn=dplyr::row_number(), ks=ifelse(.data$rn==1, .data$ks, list(2L)))
    t1 <- dplyr:: mutate(t1, deriv_order=1L, deriv_ind=ifelse(.data$deriv_group=="deriv (1,0)", 1, 2), estimate=ifelse(.data$deriv_ind==1, .data$estimate1, .data$estimate2))
    t1 <- dplyr::rename(t1, x=dplyr::all_of("Var1"), y=dplyr::all_of("Var2"))
    t1 <- dplyr::select(t1, dplyr::all_of(c("x", "y", "estimate", "deriv_order", "deriv_ind", "ks", "tks", "label", jv)))
    t2 <- tidy_kquiver(t1, thin=thin, transf=transf, neg.grad=neg.grad)

    ## rescale quiver arrows 
    dx <- max(diff(t2$x))/2; dy <- max(diff(t2$y))/2
    t3 <- dplyr::filter(t2, abs(.data$u)>0 & abs(.data$v)>0)   
    t3 <- dplyr::mutate(t3, u=(.data$u-min(.data$u))/(max(.data$u)-min(.data$u))*2*dx - dx, v=(.data$v-min(.data$v))/(max(.data$v)-min(.data$v))*2*dy - dy, xend=.data$x+.data$u, yend=.data$y+.data$v)   
    t4 <- dplyr::select(dplyr::mutate(sf::st_as_sf(t3, coords=c("x","y")), .group=dplyr::row_number()), -dplyr::one_of(c("xend", "yend")))          
    t5 <- dplyr::select(dplyr::mutate(sf::st_as_sf(t3, coords=c("xend","yend")), .group=dplyr::row_number()), -dplyr::one_of(c("x", "y"))) 
    t6 <- dplyr::group_by(dplyr::bind_rows(t4,t5), .data$.group)
    t6 <- dplyr::summarise(t6, do_union=FALSE)
    t6 <- sf::st_cast(t6, to="LINESTRING")        
    t6 <- dplyr::left_join(t6, sf::st_drop_geometry(t4), by=".group") 
    t6 <- sf::st_set_crs(t6, sf::st_crs(x$sf))
    t6 <- dplyr::bind_cols(t6, x=t3$x, y=t3$y, xend=t3$xend, yend=t3$yend)
    t6$len <- unit(scale*sqrt(t3$u^2+t3$v^2)/10^ceiling(log10(max(sqrt(t3$u^2+t3$v^2)))), "npc")
    t6 <- dplyr::rename(t6, lon=dplyr::all_of("x"), lat=dplyr::all_of("y"), lon_end=dplyr::all_of("xend"), lat_end=dplyr::all_of("yend"))
    if (length(jv)>=1) jv <- setdiff(jv, "deriv_group")
    if (length(jv)>=1) t6 <- dplyr::select(t6, dplyr::one_of(c("lon","lat","lon_end","lat_end","len","u","v",jv, "geometry"))) 
    else t6 <- dplyr::select(t6, dplyr::one_of(c("lon","lat","lon_end","lat_end","len","u","v","geometry"))) 
    
    fhat <- dplyr::slice_head(t2, n=1)
    if (length(jv)==0) fhat <- dplyr::ungroup(fhat) 
    else { t6 <- dplyr::group_by(t6, dplyr::across(dplyr::all_of(jv))); t6 <- dplyr::arrange(t6, dplyr::across(dplyr::all_of(jv)), .data$lat, .data$lon) } 
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","u","v","estimate")))
    fhat.sf <- list(tidy_ks=fhat, grid=NULL, sf=t6)
    class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}
