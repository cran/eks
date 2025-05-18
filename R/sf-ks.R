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
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.seq <- seq(5,95,by=5)
    cont.geom <- st_contourline(x=fhat, cont=cont.seq)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
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
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.seq <- seq(5, 95, by=5)
    cont.geom <- st_contourline(x=fhat, cont=cont.seq)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
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
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    
    ## convert contour polygons to multipolygon geometry
    cont.seq <- seq(5, 95, by=5)
    cont.geom <- st_contourline(x=fhat, cont=cont.seq)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
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
    #if (length(gv)==1)names(fhat.geom)[names(fhat.geom) %in% "label"] <- gv
    levels(fhat.geom$estimate) <- 1:length(levels(fhat.geom$estimate))
    fhat.geom$estimate <- as.integer(levels(fhat.geom$estimate))[fhat.geom$estimate]
    gv <- dplyr::group_vars(fhat.geom)
    fhat.geom <- dplyr::select(dplyr::ungroup(fhat.geom), -dplyr::all_of(gv))
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
    sfname <- attr(fhat.geom, "sf_column")

    ## convert contours to multipolygon geometry
	fhat <- dplyr::bind_cols(fhat, rn=1:nrow(fhat))
    cont.seq <- seq(5, 95, by=5)
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kda(x=.x, cont=cont.seq))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
	fhat <- dplyr::select(fhat, -dplyr::one_of("rn"))
	if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
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
    fhatg <- dplyr::group_by(fhat, dplyr::across({{gv2}}))
    
    fhat.geom <- st_evalpoints(x=fhatg)
    sfname <- attr(fhat.geom, "sf_column")   
    fhat.geom <- dplyr::mutate(fhat.geom, estimate=get(paste0("estimate.", dplyr::first(.data$deriv_ind))), .before=!!sfname)
    fhat.geom <- dplyr::select(fhat.geom, -dplyr::starts_with("estimate."))
    fhat.geom <- dplyr::group_by(fhat.geom, dplyr::across(dplyr::all_of(gv)))
    fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    fhat.geom <- dplyr::relocate(fhat.geom, "deriv_ind", .after="estimate")
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
  
    ## convert contours to multipolygon geometry
    cont.seq <- seq(5, 95, by=5)
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, cont=cont.seq))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
	cont.geom <- dplyr::left_join(cont.geom, dplyr::select(fhat, dplyr::all_of(gv2)), by=gv)
    cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv2), .before=!!sfname) 
    cont.geom <- dplyr::relocate(cont.geom, "deriv_ind", .after="estimate")
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1) 
    cont.geom <- dplyr::arrange(cont.geom, .data$deriv_ind)
    cont.geom <- st_create_label(cont.geom, is_kdde=TRUE)
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
    if (missing(labels)) 
    {
        mc <- as.character(match.call())[2:3]
        labels <- gsub("f1", mc[1], levels(fhat$label))
        labels <- gsub("f2", mc[2], labels)
    }
    else labels <- c(paste0(labels[1],"<",labels[2]), NA, paste0(labels[1],">",labels[2]))
    fhat <- dplyr::mutate(fhat, label="Density diff")
    gv <- dplyr::group_vars(fhat)

    ## convert ks grid to required geometry
    fhat.geom <- st_evalpoints(x=fhat)
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
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
    levels(cont.geom$label) <- as.vector(na.omit(labels))    
    cont.geom <- dplyr::arrange(cont.geom, sign(.data$estimate), abs(.data$estimate))
    if (length(gv)>0) 
    {   
        cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
        cont.geom <- dplyr::arrange(cont.geom, dplyr::across(dplyr::all_of(gv)))
    }
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom, is_kdde=FALSE) 
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
    labels <- as.character(match.call())[2]
	fhat <- .tidy_kfs(xcoord, rename=FALSE, labels=labels, ...)
	fhat.tidy <- fhat
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    fhat <- dplyr::mutate(fhat, label="Signif curv", .after="tks")

    ## convert ks grid to required geometry
    gv <- dplyr::group_vars(x)
    fhat.geom <- st_evalpoints(x=dplyr::select(fhat, -dplyr::one_of("tks", "label")))
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)
   
	## convert contour polygons to multipolygon geometry
	abs.cont <- 0.5; names(abs.cont) <- "50%"
    cont.geom <- st_contourline(x=fhat, abs.cont=abs.cont)

    cont.geom <- dplyr::relocate(cont.geom, "estimate", "contlabel", .before=1)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::mutate(cont.geom, label=!!labels, .after="estimate")
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
    cont.geom$contregion <- factor(cont.geom$label)
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
    labels <- as.character(match.call())[2]
	fhat <- .tidy_kdr(xcoord, rename=FALSE, dTolerance=dTolerance, labels=labels, ...)
    if (inherits(fhat, "sf")) fhat <- sf::st_drop_geometry(fhat)
    
    ## convert density ridge to linestring 
    fhat.sf <- sf::st_as_sf(fhat, coords=c("x","y"), crs=sf::st_crs(x))
    fhat.sf <- dplyr::summarise(fhat.sf, do_union=FALSE, label=head(.data$label,n=1), .groups="keep")
    fhat.sf <- sf::st_cast(fhat.sf, to="LINESTRING")
    fhat.sf <- dplyr::mutate(fhat.sf, contlabel=50L, estimate=0, .before="label")
    fhat.sf <- st_create_label(fhat.sf)
    fhat.sf$contregion <- factor(fhat.sf$label)
    fhat.sf <- dplyr::relocate(fhat.sf, "segment", .after="label")
    fhat.sf <- sf::st_as_sf(fhat.sf)

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
    fhat1$label <- as.integer(fhat2$label)

    fhat.sf <- dplyr::mutate(fhat1, contlabel=50, contlabel=factor(.data$contlabel), estimate=.data$label, .before="label")
    fhat.sf <- dplyr::mutate(fhat.sf, estimate=factor(.data$estimate), label=factor(.data$label))
    fhat.sf <- sf::st_as_sf(fhat.sf, coords=c("X","Y"), crs=sf::st_crs(x))
    fhat.sf <- dplyr::relocate(fhat.sf, "estimate", .before=1)
    fhat.sf <- st_create_label(fhat.sf)
    sfname <- attr(fhat.sf, "sf_column")
    if (length(gv)>0) 
    {
        fhat.sf <- dplyr::relocate(fhat.sf, dplyr::all_of(gv), .before=!!sfname)
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
    
    ## convert ks grid to required geometry
    gv <- dplyr::group_vars(x$tidy_ks)
    fhat.geom <- st_evalpoints(x=fhat)
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    crs <- sf::st_crs(x)
    fhat.geom <- sf::st_sf(fhat.geom, crs=sf::st_crs(x$sf))

    ## convert kernel support to polygon
    gv <- dplyr::group_vars(x$tidy_ks)
    fhat.sf <- dplyr::group_modify(.data=fhat.tidy, .f=~dplyr::tibble(geometry=list(sf::st_linestring(as.matrix(untidy_ks(.x))))))
    fhat.sf <- sf::st_sf(fhat.sf)
    fhat.sf <- sf::st_cast(fhat.sf, to="POLYGON")

    sfname <- attr(fhat.sf, "sf_column")
    contlev <- min(contourLevels(fhat,cont=cont)$breaks)
    fhat.sf <- dplyr::mutate(fhat.sf, estimate=contlev, contlabel=factor(cont), .before=!!sfname)
    fhat.sf <- dplyr::relocate(fhat.sf, "contlabel", "estimate", .before=1)
    fhat.sf <- dplyr::mutate(fhat.sf, label="Support", .before=!!sfname)
    fhat.sf <- dplyr::relocate(fhat.sf, "estimate", .before=1)
    if (length(gv)>0) fhat.sf <- dplyr::relocate(fhat.sf, dplyr::all_of(gv), .before=!!sfname)
    fhat.sf <- st_create_label(fhat.sf)
    fhat.sf$contregion <- factor(round_signif(fhat.sf$estimate))
    levels(fhat.sf$contregion) <- paste("\u2265", levels(fhat.sf$contregion))
    fhat.sf <- sf::st_sf(fhat.sf)

    ## if original data is not tibble
    if (!dplyr::is_grouped_df(x$sf)) fhat.sf <- data.frame(fhat.sf)
    fhat.sf <- sf::st_sf(fhat.sf, crs=sf::st_crs(x$sf))

    fhat.sf <- list(tidy_ks=fhat.tidy, grid=fhat.geom, sf=fhat.sf)
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
    sfname <- attr(fhat.geom, "sf_column")
    if (length(gv)>0) fhat.geom <- dplyr::relocate(fhat.geom, dplyr::all_of(gv), .before=!!sfname)
    crs <- sf::st_crs(x$sf)
    fhat.geom <- sf::st_sf(fhat.geom, crs=crs)

	## convert contours to multipolygon geometry
	cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(.st_contourline(x=untidy_ks(.x)))))
	cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
    cont.geom <- dplyr::filter(cont.geom, .data$contlabel %in% seq(5,95,by=5))
    #cont.geom <- dplyr::bind_cols(label=fhat$label[1], cont.geom)
    #cont.geom <- dplyr::relocate(cont.geom, "contlabel", "estimate", .before=1)
    if (length(gv)>0) cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv), .before=!!sfname)
    cont.geom <- dplyr::relocate(cont.geom, "estimate", .before=1)
    cont.geom <- st_create_label(cont.geom)
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
   
    if (any(names(t1) %in% "X")) t1 <- dplyr::rename(t1, x=dplyr::all_of("X"))
    if (any(names(t1) %in% "Y")) t1 <- dplyr::rename(t1, y=dplyr::all_of("Y"))   
    if (any(names(t1) %in% "Var1")) t1 <- dplyr::rename(t1, x=dplyr::all_of("Var1"))
    if (any(names(t1) %in% "Var2")) t1 <- dplyr::rename(t1, y=dplyr::all_of("Var2"))   
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
    sfname <- attr(t6, "sf_column")
    if (length(jv)>=1) jv <- setdiff(jv, "deriv_group")
    if (length(jv)>=1) t6 <- dplyr::select(t6, dplyr::one_of(c("lon","lat","lon_end","lat_end","len","u","v",jv, sfname))) 
    else t6 <- dplyr::select(t6, dplyr::one_of(c("lon","lat","lon_end","lat_end","len","u","v",sfname))) 
    
    fhat <- dplyr::slice_head(t2, n=1)
    if (length(jv)==0) fhat <- dplyr::ungroup(fhat) 
    else 
    { 
        t6 <- dplyr::group_by(t6, dplyr::across(dplyr::all_of(jv)))
        t6 <- dplyr::arrange(t6, dplyr::across(dplyr::all_of(jv)), .data$lat, .data$lon) 
    } 
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","u","v","estimate")))
    fhat.sf <- list(tidy_ks=fhat, grid=NULL, sf=t6)
    class(fhat.sf) <- "sf_ks"

    return(fhat.sf)
}

## spatial version of ks::as.kde
## x = complete sf polygon grid
st_as_kde <- function(x, attrib=1, density, ...)
{
    if (is.numeric(attrib)) attrib <- names(x)[attrib]
    x2 <- tidy_intergrid(x, attrib=attrib) 
    crs <- sf::st_crs(x)

    if (missing(density)) 
    {
        if (all(x2$estimate>=0)) density <- TRUE
        else if (all(x2$estimate<0)) density <- FALSE
        else density <- (abs(mean(x2$estimate[x2$estimate<0]))/mean(x2$estimate[x2$estimate>=0]) < 1e-6) 
    }

    ## convert to tidy KDE
    ## fhat is based on ks::kde so is evaluated at vertices of M1 x M2 estimation grid
    fhat <- .tidy_as_kde(data=x2, is_tidy=TRUE, density=density, rename=FALSE, ...)
    fhat <- dplyr::slice_head(fhat)
    fhat <- dplyr::select(fhat, -dplyr::one_of(c("x","y","estimate")))
    if (!density)
    {
        fhat <- dplyr::select(fhat, -dplyr::one_of(c("deriv_order")))
        fhat <- as_tidy_ks(fhat)
    }
    
    ## convert ks grid to required geometry
    ## fhat.geom is evaluated at (M1-1) x (M2-1) centroids of M1 x M2 estimation grid
    ## i.e. fhat.geom is (M1-1) x (M2-1) rectangle polygon grid
    fhat.geom <- dplyr::transmute(x, estimate=.data[[attrib]])
    if (!density)
    {
        fhat.geom <- dplyr::mutate(fhat.geom, deriv_order=1L, deriv_ind=1L, deriv_group="(1)", .after="estimate")
        fhat.geom <- dplyr::group_by(fhat.geom, .data$deriv_group, .add=TRUE)
    }
    
    ## convert contour polygons to multipolygon geometry
    cont.seq <- seq(5, 95, by=5)
    if (density) 
    {
        cont.geom <- st_contourline(x=fhat, cont=cont.seq)
        cont.geom <- sf::st_sf(cont.geom, crs=crs)
        cont.geom <- suppressWarnings(sf::st_crop(cont.geom, sf::st_bbox(x)))
        cont.geom <- st_create_label(cont.geom)
    }
    else 
    {
        fhat <- dplyr::group_by(fhat, .data$deriv_group, .add=TRUE)
        gv <- dplyr::group_vars(fhat)
        gv2 <- c("deriv_ind", gv)
        crs <- sf::st_crs(x) 
        
        ## convert contours to multipolygon geometry
        cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, cont=cont.seq))))
        cont.geom <- dplyr::group_modify(.data=cont.geom, .f=~data.frame(.x$geometry))
        cont.geom <- dplyr::left_join(cont.geom, dplyr::select(fhat, dplyr::all_of(gv2)), by=gv)
        cont.geom <- dplyr::arrange(cont.geom, .data$deriv_ind)
        cont.geom <- sf::st_sf(cont.geom, crs=crs)
        sfname <- attr(cont.geom, "sf_column")
        cont.geom <- dplyr::relocate(cont.geom, "deriv_ind", .after="estimate") 
        cont.geom <- dplyr::relocate(cont.geom, dplyr::all_of(gv2), .before=!!sfname) 
        cont.geom <- sf::st_sf(cont.geom, crs=crs)
        cont.geom <- suppressWarnings(sf::st_crop(cont.geom, sf::st_bbox(x)))
        cont.geom <- st_create_label(cont.geom)

        #fhat.sf <- list(tidy_ks=fhat, grid=fhat.geom, sf=cont.geom)
        #class(fhat.sf) <- "sf_ks"
    }

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


## extract or compute new contours
## x = ks_sf object
st_get_contour <- function(x, cont=c(25,50,75), breaks, disjoint=TRUE, digits)
{
    ## some defaults
    crop <- TRUE      ## crop to bounding box of x$grid
    cont100 <- TRUE   ## compute 100% contour = bounding box \ union {all other contour regions}
    oct <- object_class(x, type="tks")
    if (missing(digits)) digits <- ifelse(oct %in% c("kcde", "kcopula", "kde.loctest", "kfs", "kdr", "kms"), 2, 4)
    fhat <- dplyr::slice_head(x$tidy_ks)
    missing_breaks <- missing(breaks)

    if (oct %in% c("kdde", "kde.loctest")) cont <- unique(c(cont, -cont))
    if (oct %in% c("kms", "kfs"))
    {
        xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
    }
    else if (oct %in% c("kdr", "kms", "kquiver", "ksupp"))
    {
        xc <- x$sf
    }
    else if (oct %in% "kdde")
    { 
        ## absolute contour levels 
        if (!missing(breaks))
        {
            xc <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, abs.cont=breaks, digits=digits)), deriv_ind=.x$deriv_ind))
            xc <- dplyr::group_modify(.data=xc, .f=~data.frame(deriv_ind=.x$deriv_ind, .x$geometry))
            xc <- dplyr::arrange(xc, .data$deriv_group, .data$contlabel)
            xc <- sf::st_as_sf(xc)
        }
        ## add non-integer contour levels
        else
        {
            xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
            xc <- dplyr::arrange(xc, .data$deriv_group, .data$estimate)
            cont2 <- cont[!(cont %in% x$sf$contlabel)]
            cont2 <- cont2[cont2>=0 & cont2<=100]
            cont2 <- cont2[round(cont2,0)!=cont2]
            sfname <- attr(xc, "sf_column")
            if (length(cont2)>0)
            {   
                cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, cont=cont2, digits=digits)), deriv_ind=.x$deriv_ind))
                cont.geom <- lapply(cont.geom$deriv_ind, function(.) cbind(dplyr::select(cont.geom[.,], -!!sfname), cont.geom[[sfname]][[.]]))
                cont.geom <- do.call(rbind, cont.geom)
                cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x$sf))
                #xc <- rbind(xc, cont.geom)
                xc <- rbind(xc, st_create_label(cont.geom, is_kdde=TRUE))
            }
            cont.pos <- cont[cont>=0]
            xc$contlabel <- factor(xc$contlabel, levels=c(-sort(cont.pos), sort(cont.pos, decreasing=TRUE)))
        }
         
        sfname <- attr(xc, "sf_column")
        if (is.factor(xc$estimate)) xc$estimate <- unfactor(xc$estimate)

        xc <- dplyr::mutate(xc, pos=sign(.data$estimate), .before=!!sfname)
        xc <- dplyr::arrange(xc, .data$deriv_ind, .data$contlabel)
        xc <- dplyr::group_by(xc, .data$deriv_group, .data$pos, .add=TRUE)
        xc <- dplyr::group_modify(xc, .f=~dplyr::arrange(.x, .data$estimate))
        gv <- dplyr::group_vars(xc)
        gv <- setdiff(gv, "pos")
        xc <- dplyr::group_by(xc, dplyr::across({{gv}}))
        xc <- dplyr::select(xc, -"pos")
    } 
    else 
    {
        ## absolute contour levels 
        if (!missing(breaks))
        {
            if (!inherits(breaks, c("matrix","data.frame","tbl_df"))) breaks <- dplyr::as_tibble(data.frame(breaks=breaks))
            if (oct %in% "kda") stop("Explicit contour breaks not supported for kda objects")
            xc <- st_contourline(x=fhat, abs.cont=breaks$breaks, digits=digits)
            xc <- sf::st_sf(xc, crs=sf::st_crs(x$sf))
            xc <- dplyr::arrange(xc, .data$contlabel)
        }
        ## add non-integer contour levels
        else
        {
            xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
            cont2 <- cont[!(cont %in% x$sf$contlabel)]
            cont2 <- cont2[cont2>=0 & cont2<=100]
            cont2 <- cont2[round(cont2,0)!=cont2]
            if (length(cont2)>0)
            {   
                cont.geom <- st_contourline(x=fhat, cont=cont2, digits=digits)
                cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x$sf))
                xc <- rbind(xc, st_create_label(cont.geom))
                xc <- dplyr::arrange(xc, .data$contlabel)
                xc <- dplyr::arrange(xc, sign(.data$estimate), abs(.data$estimate))
            }
        }
    }
    xc <- sf::st_as_sf(xc)
    sf::st_crs(xc) <- sf::st_crs(x$grid)
    gv <- dplyr::group_vars(xc)
    
    ## crop contours to bounding box
    if (crop & !is.null(x$grid)) 
    {
        xbbox <- sf::st_as_sfc(sf::st_bbox(x$grid))
        xc <- suppressWarnings(dplyr::group_modify(xc, .f=~sf::st_crop(.x, xbbox))) 
        xc <- sf::st_as_sf(xc)
    }

    ## add 100% contour for KDDE 
    if (oct %in% "kdde" & cont100) xc <- st_add_contour100(xc, xbbox=xbbox, missing_breaks=missing_breaks) 
    
    ## create disjoint contours
    xc <- dplyr::group_by(xc, dplyr::across({{gv}}))
    if (!(oct %in% c("kdr", "kms")) & disjoint) xc <- st_contour_disjoint(xc) 

    ## add more informative labels
    xc <- st_create_label(xc, digits=digits, missing_breaks=missing_breaks, is_kdde=oct %in% "kdde")
    if (oct %in% "ksupp") xc$contregion <- x$sf$contregion
    
    ## factor for $estimate has fixed number of signif digits displayed 
    ## slightly different from $contlabel
    xc$estimate <- round_signif(xc$estimate, digits=digits)
    xc$estimate <- factor(xc$estimate, levels=unique(xc$estimate))

    ## if original data is not tibble
    if (!inherits(x$sf, "tbl_df"))  xc <- sf::st_sf(data.frame(xc))
    xc <- sf::st_sf(xc, crs=sf::st_crs(x$sf))

    return(xc)
}

## add 100% contour
st_add_contour100 <- function(x, xbbox, missing_breaks)
{
    if (missing(xbbox)) xbbox <- sf::st_as_sfc(sf::st_bbox(x))
    xc <- dplyr::group_modify(x, .f=~.st_add_contour100(.x, xbbox=xbbox, missing_breaks=missing_breaks))
    xc <- dplyr::select(xc, names(x))
    xc <- sf::st_sf(xc)

    return(xc)
}

.st_add_contour100 <- function(x, xbbox, missing_breaks)
{
    if (any(x$estimate==0)) xc <- x
    else
    {
        xc100.geom <- sf::st_union(sf::st_geometry(x))
        xc100.geom <- sf::st_difference(xbbox, xc100.geom)
        xc100 <- x[1,]
        sf::st_geometry(xc100) <-  sf::st_geometry(xc100.geom) 
        if (missing_breaks) xc100$contlabel <- "100" else xc100$contlabel <- "0"
        xc100$estimate <- 0;
        xc <- rbind(x[x$estimate<0,], xc100, x[x$estimate>0,])
        xc$contlabel <- factor(xc$contlabel)
    }    
    return(xc)
}

## x = contour regions as multipolygons (e.g. $sf from sf_ks objects)
st_create_label <- function(x, digits=4, is_kdde=FALSE, missing_breaks=TRUE)
{      
    ## create labels for contour regions    
    ## contlabel is deprecated but kept backwards compatibility < 1.1.0
    gv <- dplyr::group_vars(x)
    if (is.factor(x$estimate)) x$estimate <- unfactor(x$estimate)
     
    ## create contregion label
    if (any(names(x) %in% "label")) { xc <- x; xc$contregion <- factor(xc$label) }
    else 
    {
        if (is_kdde)
            xc <- dplyr::mutate(x, contregion=paste0(ifelse(.data$estimate>0, ">", "<"), round_signif(.data$estimate, digits=digits)), .after="estimate")
        else
            xc <- dplyr::mutate(x, contregion=paste0(">", round_signif(.data$estimate, digits=digits)), .after="estimate")
        xc$contregion[xc$contlabel %in% c("0","100")] <- gsub("<", ">", xc$contregion[which(xc$contlabel %in% c("0","100"))-1])
        xc <- dplyr::arrange(xc, .data$estimate, .by.group=TRUE)
        xc$contregion <- gsub("-", "\u2013", gsub(">", "\u2265 ", gsub("<", "\u2264 ", xc$contregion)))
    }
    if (is.factor(xc$contlabel)) 
        xc$contlabel <- factor(droplevels(xc$contlabel), levels=unique(xc$contlabel))
    else xc$contlabel <- factor(xc$contlabel, levels=unique(xc$contlabel))
    xc$contregion <- factor(xc$contregion, levels=unique(xc$contregion))

    ## create contline label
    if (is_kdde)
    {
        xc <- dplyr::mutate(xc, contline=ifelse(.data$estimate==0, NA, 1L), .after="contregion")
        xc$contline <- factor(xc$contline)
    }

    ## create cont percentage label  
    if (missing_breaks) 
    {   
        xc$contperc <- xc$contlabel
        levels(xc$contperc) <- paste0(levels(xc$contperc), "%")
    }
    else 
    {
        xc$contperc <- NA  
        xc$contperc <- factor(xc$contperc)
    }
    
    sfname <- attr(xc, "sf_column")
    xc <- dplyr::relocate(xc, "contlabel", "contperc", .after="estimate")
    xc <- dplyr::relocate(xc, dplyr::starts_with("cont"), .after="estimate")
    xc <- dplyr::relocate(xc, dplyr::all_of(gv), .before=!!sfname)
    if (any(names(xc) %in% "deriv_ind") & any(names(xc) %in% "deriv_group"))  xc <- dplyr::relocate(xc, "deriv_ind", .before="deriv_group")
    xc <- dplyr::group_by(xc, dplyr::across({{gv}}))
    xc <- dplyr::arrange(xc, dplyr::across({{gv}}), .data$contlabel)
    xc <- sf::st_as_sf(xc)
     
    return(xc)
}

## interpolate sf polygon grid and fill in missing values on x 
st_intergrid <- function(x, attrib, cellsize, verbose=FALSE)
{
    sf::st_geometry(x) <- "geometry"
    sfname <- attr(x, "sf_column")
    if (missing(attrib)) attrib <- 1
    if (is.numeric(attrib)) attrib <- names(x)[attrib] 
    xcrs <- sf::st_crs(x)
    if (sf::st_is_longlat(x)) sf::st_crs(x) <- NA
    
    ## default grid cell size - based on grid cell with mean area 
    if (missing(cellsize))
    {
        xarea <- geos::geos_area(x)
        xarea.mn <- mean(xarea, na.rm=TRUE)
        xseg1 <- st_rectangle_segments(x[which.min(abs(xarea.mn-xarea)),])
        cellsize.default <- geos::geos_length(xseg1)[1:2]
        cellsize <- signif(cellsize.default,3)
    }
    else if (!missing(cellsize)) { if (length(cellsize)<2) cellsize <- rep(cellsize, times=2)[1:2] }
  
    ## create bounding box which fits exactly ncell gridcells of dimension cellsize
    xbbox <- sf::st_as_sfc(sf::st_bbox(x))
    xbbox.seg <- st_bbox_segments(x)
    ncell <- round(geos::geos_length(xbbox.seg)[1:2]/cellsize)
    xgrid.bbox <- sf::st_bbox(x)
    xgrid.bbox[3:4] <- xgrid.bbox[1:2] + ncell*cellsize
    xgrid.bbox <- sf::st_as_sfc(xgrid.bbox)
    sf::st_crs(xgrid.bbox) <- sf::st_crs(x)
    
    ## for grid cells with varying areas, subdivide grid cells before interpolation     
    xarea <- geos::geos_area(x)    
    nsubdiv <- ifelse(all(abs(xarea - mean(xarea)) <= 0.1*mean(xarea)), 1, 2) 
    xgrid.subdiv <- sf::st_sf(geometry=sf::st_make_grid(xgrid.bbox, cellsize=cellsize/nsubdiv, crs=sf::st_crs(x), square=TRUE))
    sf::st_geometry(xgrid.subdiv) <- sfname
    xgrid.subdiv <- st_add_index(xgrid.subdiv)

    ## assign those sub-divided cells which cover >=2 different original 
    ## grid cells in xgrid to single value that is closest to weighted mean
    ## search for more efficient operation than st_intersection?
    xgrid.subdiv.int <- suppressWarnings(sf::st_intersection(x, xgrid.subdiv))
    xgrid.subdiv.int <- xgrid.subdiv.int[sf::st_is(xgrid.subdiv.int, "POLYGON"),]
    xgrid.subdiv.int <- dplyr::mutate(xgrid.subdiv.int, area=geos::geos_area(xgrid.subdiv.int), .before=!!sfname)
    ## remove small intersections <= 0.01 * grid cell area
    ## as these may unduly affect the attrib weighted mean
    ## if these intersections have attrib value=0
    xgrid.subdiv.int <- dplyr::filter(xgrid.subdiv.int, .data$area>prod(!!cellsize)*1e-2)
    xgrid.subdiv.int <- dplyr::arrange(xgrid.subdiv.int, dplyr::desc(.data$area*.data[[attrib]]))
    xgrid.subdiv.int <- dplyr::group_by(xgrid.subdiv.int, .data$cell_id)
    xgrid.subdiv.int <- dplyr::summarise(sf::st_drop_geometry(xgrid.subdiv.int), attribm=weighted.mean(.data[[attrib]], .data$area, na.rm=TRUE), .attrib=.data[[attrib]][which.min(abs(.data$attribm-.data[[attrib]]))[1]])
    #xgrid.subdiv.int <- dplyr::summarise(sf::st_drop_geometry(xgrid.subdiv.int), .attrib=weighted.mean(.data[[attrib]], .data$area, na.rm=TRUE))
    xgrid.subdiv.int <- dplyr::rename_with(xgrid.subdiv.int, function(.) attrib, ".attrib")
    xgrid.subdiv.int$attribm <- NULL
   
    ## fill in other grid cells that are not intersected 
    xgrid.subdiv <- dplyr::right_join(xgrid.subdiv.int, xgrid.subdiv, by=c("cell_id"))
    xgrid.subdiv <- sf::st_sf(xgrid.subdiv)
    xgrid.subdiv <- dplyr::relocate(xgrid.subdiv, dplyr::all_of(attrib), .before=1)
    xgrid.subdiv <- dplyr::arrange(xgrid.subdiv, .data$cell_id, .data[[attrib]])
    xgrid.subdiv <- dplyr::distinct(xgrid.subdiv, .data$cell_id, .keep_all=TRUE)

    ## combine subdivided grid cells into original grid cells
    xgrid <- st_union_grid(xgrid.subdiv, attrib=attrib, n=nsubdiv, .f=mean)
    xgrid <- st_add_index(xgrid)

    ## keep track of NA attribute values
    # xgrid[["attrib_na"]] <- xgrid[[attrib]]
    xgrid[[attrib]][is.na(xgrid[[attrib]])] <- 0
    # xgrid <- dplyr::relocate(xgrid, "attrib_na", .after={{attrib}})
    # xgrid <- dplyr::rename_with(xgrid, function(.) paste0(attrib,"_na"), "attrib_na")
    if (is.na(sf::st_crs(xgrid))) sf::st_crs(xgrid) <- xcrs

    if (verbose) cat(paste0("Original grid     sum(", attrib, ") = ", sum(x[[attrib]], na.rm=TRUE)), paste0("\nInterpolated grid sum(", attrib, ") = ", sum(xgrid[[attrib]])), "\n") 

    return(xgrid)
}

## convert raster to sf polygon grid
st_raster_as_grid <- function(x)
{
    dims <- dim(x)-1
    xsf <- sf::st_as_sf(x)
    xsf <- cbind(expand.grid(cell_id1=1:dims[1], cell_id2=dims[2]:1), xsf)
    xsf <- dplyr::as_tibble(xsf)
    xsf <- dplyr::arrange(xsf, .data$cell_id2, .data$cell_id1)
    xsf <- sf::st_sf(xsf) 
    sfname <- attr(xsf, "sf_column")
    xsf <- dplyr::relocate(xsf, "cell_id1", "cell_id2", .before=!!sfname)
    
    return(xsf)
}

## prepare kde grid for raster conversion
## x = ks_sf object (output from st_kde/st_as_kde) 
st_add_contour_label <- function(x, cont=c(25,50,75))
{
    if (inherits(x, "sf_ks")) xgrid <- sf::st_drop_geometry(st_add_index(x$grid))
    else xgrid <- sf::st_drop_geometry(st_add_index(x))
    raster.dim <- as.vector(apply(xgrid[,c("cell_id1", "cell_id2")], 2, max))
    breaks <- contourLevels(x, cont=cont)
    contlev <- sort(c(min(x$grid$estimate)-0.1*abs(max(x$grid$estimate)), breaks$breaks))
    breaks$cont <- sort(c(0,head(breaks$cont,n=-1)))
    xcontlabel <- cut(x$grid$estimate, contlev, right=TRUE, labels=breaks$cont)
    xcontperc <- xcontlabel 
    levels(xcontperc) <- paste0(levels(xcontperc), "%")
    xg <- dplyr::mutate(x$grid, contlabel=!!xcontlabel, contperc=!!xcontperc,  nx=raster.dim[1], ny=raster.dim[2], .after="estimate")
    
    return(xg)
}

## add indices, starting from SW corner, to rectangular grid
st_add_index <- function(x)
{
    xcrs <- sf::st_crs(x)
    if (!is.na(sf::st_crs(x)))
    {
        if (sf::st_is_longlat(x)) sf::st_crs(x) <- NA
    } 
    ## centroids of first horizontal row and first vertical col
    xbbox <- st_bbox_segments(x)
    h1 <- xbbox[1]; v1 <- xbbox[2]
    xh.cent <- sf::st_filter(x, h1) 
    xv.cent <- sf::st_filter(x, v1)    
    xh.cent <- suppressWarnings(sf::st_centroid(xh.cent))
    xv.cent <- suppressWarnings(sf::st_centroid(xv.cent))

    sfname <- attr(x, "sf_column")
    xcent <- suppressWarnings(sf::st_centroid(x))

    ## compute complete aligned grid with row,col indices
    xgrid <- dplyr::mutate(x, cell_id1=sf::st_nearest_feature(xcent, xh.cent), cell_id2=sf::st_nearest_feature(xcent, xv.cent), .before=!!sfname)
    xgrid$cell_id1 <- factor(xgrid$cell_id1) 
    levels(xgrid$cell_id1) <- 1:(length(levels(xgrid$cell_id1)))
    xgrid$cell_id1 <- unfactor(xgrid$cell_id1)
    xgrid$cell_id2 <- factor(xgrid$cell_id2) 
    levels(xgrid$cell_id2) <- 1:(length(levels(xgrid$cell_id2)))
    xgrid$cell_id2 <- unfactor(xgrid$cell_id2)
    xgrid <- dplyr::arrange(xgrid, .data$cell_id2, .data$cell_id1)
    xgrid <- dplyr::mutate(xgrid, cell_id=1:nrow( xgrid), .before="cell_id1")
    
    if (is.na(sf::st_crs(xgrid))) sf::st_crs(xgrid) <- xcrs

    return(xgrid)
}

median_closest <- function(x) { md <- median(x, na.rm=TRUE); ifelse(is.na(md), NA, x[which.min(abs(md-x))]) }
mean_closest <- function(x) { mn <- mean(x, na.rm=TRUE); ifelse(is.na(mn), NA, x[which.min(abs(mn-x))]) }

## aggregate grid by combining n x n grid cells  
st_union_grid <- function(x, attrib, n=1, .f)
{
    if (missing(.f)) .f <- function(x) (mean(x, na.rm=TRUE))
    xu <- dplyr::mutate(x, celln_id1=(.data$cell_id1-1)%/%n+1, celln_id2=(.data$cell_id2-1)%/%n+1, .after="cell_id2")
    xu <- dplyr::group_by(xu, .data$celln_id1, .data$celln_id2)
    xu2 <- dplyr::filter(xu, dplyr::n()>1)
    xu2.attrib <- dplyr::summarise(xu2, .attrib=.f(.data[[attrib]]), do_union=FALSE, .groups="drop")
    xu2.attrib <- sf::st_convex_hull(xu2.attrib)
    xu2.attrib <- dplyr::rename(xu2.attrib, cell_id1="celln_id1", cell_id2="celln_id2")
    xu2.attrib <- dplyr::rename_with(xu2.attrib, function(.) attrib, ".attrib")
    xu2.attrib <- dplyr::relocate(xu2.attrib, {{attrib}}, .before=1)
    
    xu.attrib <- rbind(xu2.attrib, dplyr::select(dplyr::filter(xu, dplyr::n()<=1), names(xu2.attrib)))
    xu.attrib <- dplyr::arrange(xu.attrib, .data$cell_id2, .data$cell_id1)
    xu.attrib <- dplyr::mutate(xu.attrib, cell_id=1:nrow(xu.attrib), .before="cell_id1")
    
    return(xu.attrib)
}
