#####################################################################
## Auxiliary functions for st_* functions
#####################################################################

## plot function for sf_ks object

ggplot.sf_ks <- function(data=NULL, mapping=aes(), ..., which_geometry="sf")
{
    if (!is.null(data))
    {
        g <- match.arg(which_geometry, c("sf","grid"))
        data1 <- data[[g]]
    }
    gg <- ggplot2::ggplot(data=data1, mapping=mapping, ...) + labs_ks(data) + guides_ks(data) 
    return(gg)
}

plot.sf_ks <- function(x, ...) { plot_sf_ks(x=x, ...) }

plot_sf_ks <- function(x, which_geometry="sf", cont=c(25,50,75), abs_cont=breaks, breaks, which_deriv_ind=1, main="", pal, col, pos="bottomleft", key.pos=NULL, alpha, legend=TRUE, ...) 
{
    g <- match.arg(which_geometry, c("sf","grid"))
    y <- x[[g]]
    oc <- head(x$tidy_ks$tks,1)

    if (missing(alpha)) { if (g=="sf") alpha <-1 else alpha <- 0.1 }
    if (missing(col)) 
    {
        if (any(oc %in% "kdr")) col <- ggplot2::alpha(6, alpha=alpha)
        else if (any(oc %in% "kfs")) col <- ggplot2::alpha(7, alpha=alpha)
        else if (any(oc %in% "ksupp")) col <- NA
        else if (any(oc %in% "kquiver")) col <- ggplot2::alpha(1, alpha=alpha)
 
        pal.missing <- FALSE
        if (missing(pal)) 
        {
            pal.missing <- TRUE
            if (any(oc %in% "kdde")) pal <- function(.) { colorspace::diverging_hcl(n=., palette="Blue-Red", alpha=alpha) }
            else if (any(oc %in% "kcurv")) pal <- function(.) { colorspace::sequential_hcl(n=., h1=30, c1=360, c2=60, alpha=alpha, rev=TRUE) }
            else if (any(oc %in% "kde.loctest")) pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Dark2", alpha=alpha, rev=TRUE) }
            else if (any(oc %in% "kcde")) pal <- function(.) { colorspace::sequential_hcl(n=., rev=FALSE, palette="Viridis", alpha=alpha) }
            else if (any(oc %in% c("kroc","kms"))) pal <- function(.) { colorspace::qualitative_hcl(n=., rev=TRUE, palette="Set2", alpha=alpha) }
            else if (any(oc %in% "kda")) 
            {
                if (g=="sf") 
                {
                    pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Dark2", alpha=alpha) } 
                    ng <- length(table(dplyr::group_indices(y)))
                    cols <- pal(ng)
                    if (missing(breaks)) nc <- length(cont) else nc <- length(breaks)
                    col <- rep(cols, each=nc)
                } 
                else pal <- function(.) { colorspace::qualitative_hcl(n=., palette="Dark2", alpha=alpha) }
            }
            else pal <- function(.) { colorspace::sequential_hcl(n=., rev=TRUE, palette="Heat2", alpha=alpha) }
        }
    }
    else 
    {
        col <- ggplot2::alpha(col, alpha=alpha)
        pal.missing <- TRUE
    }

    if (g=="sf")
	{   
        if (!(oc %in% c("kms", "kdr"))) yd <- dplyr::mutate(y, contlabel=as.character(.data$contlabel))
           
        if (any(oc %in% c("kdr", "kms")))
        {
            yd <- dplyr::select(y, dplyr::all_of("label")) 
        }
        else if (any(oc %in% "ksupp"))
        {
            yd <- dplyr::select(y, dplyr::all_of("contlabel")) 
        }
        else if (any(oc %in% "kdde"))
        {
            y <- st_get_contour(x, cont=cont, breaks=abs_cont, which_deriv_ind=which_deriv_ind) 
            y <- y[y$deriv_ind == which_deriv_ind,]
            y$estimate <- droplevels(y$estimate)
            yd <- dplyr::select(y, dplyr::all_of("estimate")) 
        }
        else
        {
            y <- st_get_contour(x, cont=cont, breaks=abs_cont)
            yd <- dplyr::select(y, dplyr::all_of("estimate"))
        }
        
        if (missing(col)) plot(yd, main=main, pal=pal, ...)
        else 
        { 
            if (any(oc %in% "kda")) plot(yd, main=main, col=NA, ...)
            else plot(yd, main=main, col=col, ...)
        }
	}
	else if (g=="grid") 
	{
        if (any(oc %in% "kdde")) y <- y[y$deriv_ind == which_deriv_ind,]
        yd <- dplyr::select(y, dplyr::all_of("estimate"))
        
        if (missing(col)) plot(yd, main=main, pal=pal, key.pos=key.pos, ...)
        else plot(yd, main=main, col=col, key.pos=key.pos, ...)
    }

    ## add legend
    ## mapsf legends don't allow for line types in legend boxes so st_ksupp plots
    ## aren't supported
    legend <- legend & !(oc %in% c("ksupp"))
    if (legend)
    {
        if (!requireNamespace("mapsf", quietly=TRUE)) stop("Install the mapsf package as it is required.", call.=FALSE)
    
        forms <- list(...)
        forms <- forms[names(forms) %in% names(formals(mapsf::mf_legend_t))]
        forms <- forms[!(names(forms) %in% "border")]
        gu <- guides_ks(dplyr::add_row(dplyr::ungroup(x$tidy_ks), ks=list(2L)))
        gu.title <- gu$guides$fill$params$title

        if (!(oc %in% "kda"))
        {
            contlabel <- sort(unique(y$contlabel), decreasing=TRUE)
            if (missing(breaks)) contlabel <- paste0(contlabel,"%")
            nc <- length(contlabel)
        }

        if (g=="grid")
        {
            if (oc %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv, pal=pal(length(gv)), pos=pos, title=gu.title), forms))
            }
            else
            {
                y2 <- st_get_contour(x, cont=cont, breaks=abs_cont)
                yd2 <- dplyr::select(y2, dplyr::all_of("estimate"))
                gv <- yd2$estimate
                do.call(mapsf::mf_legend, args=c(list(type="choro", val=gv, pal=pal(length(gv)), pos=pos, title=gu.title), forms))
            }
        }
        else if (!is.null(forms$border)) 
        {
            if (oc %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), dplyr::all_of(dplyr::group_vars(y))))
                ng <- length(gv)
                do.call(mapsf::mf_legend, args=c(list(type="symb", val=gv, pal=unique(forms$border), pos=pos, title=gu.title, cex=rep(3,ng), pch=rep("-", ng)), forms))
            }
            else 
            {
                if (!all(is.na(forms$border)) & pal.missing) do.call(mapsf::mf_legend, args=c(list(type="symb", val=contlabel, pal=rev(forms$border), pos=pos, title=gu.title, cex=3, pch="-"), forms))
                else do.call(mapsf::mf_legend,, args=c(list(type="typo", val=contlabel, pal=rev(pal(nc)), pos=pos, title=gu.title), forms)) 
            }
        }
        else if (!missing(col))
        {
            if (oc %in% "kdr") do.call(mapsf::mf_legend, args=c(list(type="symb", val="Density ridge", pal=col,  pos=pos, title=gu.title, cex=3, pch="-"), forms))
            else if (oc %in% "kfs") do.call(mapsf::mf_legend, args=c(list(type="typo", val="Signif curv", pos=pos, pal=col, title=gu.title), forms))
            else if (oc %in% "ksupp") do.call(mapsf::mf_legend, args=c(list(type="symb", val="Support\nconvex hull", pos=pos, pal=col, title=gu.title), forms))
            else if (oc %in% "kda")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), dplyr::all_of(dplyr::group_vars(y))))
                ng <- length(gv)
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv, pal=pal(ng), pos=pos, title=gu.title), forms))
            }
            else 
            {
                if (!missing(breaks)) { col <- rev(col) }
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pal=col, pos=pos, title=gu.title), forms))
            }
        }
        else if (!missing(pal)) 
        {        
            if (oc %in% "kms") 
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                nc <- length(gv)
                formsg <- list(...)
                if (is.null(formsg$pch)) pch <- 1 else pch <- unique(formsg$pch) 
                do.call(mapsf::mf_legend, args=c(list(type="symb", val=gv, pal=pal(nc), pos=pos, title=gu.title, cex=rep(1,nc), pch=pch), forms))
            }
            else if (oc %in% "kde.loctest")
            {
                gv <- levels(dplyr::pull(sf::st_drop_geometry(y), .data$label))
                ind <- order(gv, decreasing=TRUE)
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=gv[ind], pal=pal(length(gv))[ind], pos=pos, title=gu.title), forms))
            }
            else 
                do.call(mapsf::mf_legend, args=c(list(type="typo", val=contlabel, pal=rev(pal(nc)), pos=pos, title=gu.title), forms))
        }
    }
}

contourLevels.sf_ks <- function(x, cont=c(25,50,75), group=FALSE, ...) { ks::contourLevels(x=x$tidy_ks, cont=cont, group=group, ...) }

## extract or compute new contours
## x = ks_sf object

st_get_contour <- function(x, cont=c(25,50,75), breaks, which_deriv_ind, disjoint=TRUE, as_point=FALSE)
{
    oc <- head(x$tidy_ks$tks,1)
    fhat <- dplyr::slice_head(x$tidy_ks)
    if (oc %in% c("kdde", "kde.loctest")) cont <- unique(c(cont, -cont))
    
    if (oc %in% c("kms", "kfs"))
    {
        xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
    }
    else if (oc %in% "kdr")
    {
        xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
        disjoint <- FALSE
    }
    else if (oc %in% "ksupp")
    {
        xc <- x$sf
    }
    else if (oc %in% "kdde")
    {
        ## absolute contour levels 
        if (!missing(breaks))
        {
            xc <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, abs.cont=breaks)), deriv_ind=.x$deriv_ind))
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
            cont2 <- cont2[cont2>0 & cont2<100]
            cont2 <- cont2[round(cont2,0)!=cont2]
            if (length(cont2)>0)
            {   
                cont.geom <- dplyr::group_modify(.data=fhat, .f=~dplyr::tibble(geometry=list(st_contourline_kdde(x=untidy_ks(.x), which_deriv_ind=.x$deriv_ind, cont=cont2)), deriv_ind=.x$deriv_ind))
                cont.geom <- lapply(cont.geom$deriv_ind, function(.) cbind(dplyr::select(cont.geom[.,], -"geometry"), cont.geom$geometry[[.]]))
                cont.geom <- do.call(rbind, cont.geom)
                cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x$sf))
                xc <- rbind(xc, cont.geom)
            }
        }
        if (!missing(which_deriv_ind)) xc <- dplyr::filter(xc, .data$deriv_ind %in% which_deriv_ind)
        xc <- dplyr::mutate(xc, pos=sign(.data$estimate), .before="geometry")
        xc <- dplyr::arrange(xc, .data$deriv_ind, .data$contlabel)
        xc <- dplyr::group_by(xc, .data$deriv_group, .data$pos)
        xc <- dplyr::group_modify(xc, .f=~dplyr::arrange(.x, .data$estimate))
        xc <- sf::st_as_sf(xc)
    } 
    else 
    {
        ## absolute contour levels 
        if (!missing(breaks))
        {
            if (!inherits(breaks,"tbl_df")) breaks <- dplyr::as_tibble(data.frame(breaks=breaks))
            if (oc %in% "kda") stop("Explicit contour breaks not supported for kda objects")
            
            xc <- st_contourline(x=fhat, abs.cont=breaks$breaks)
            xc <- sf::st_sf(xc, crs=sf::st_crs(x$sf))
            xc <- dplyr::arrange(xc, .data$contlabel)
        }
        ## add non-integer contour levels
        else
        {
            xc <- dplyr::filter(x$sf, .data$contlabel %in% cont)
            cont2 <- cont[!(cont %in% x$sf$contlabel)]
            cont2 <- cont2[cont2>0 & cont2<100]
            cont2 <- cont2[round(cont2,0)!=cont2]
            if (length(cont2)>0)
            {   
                cont.geom <- st_contourline(x=fhat, cont=cont2)
                cont.geom <- sf::st_sf(cont.geom, crs=sf::st_crs(x$sf))
                xc <- rbind(xc, cont.geom)
                xc <- dplyr::arrange(xc, .data$contlabel)
                xc <- dplyr::arrange(xc, sign(.data$estimate), abs(.data$estimate))
            }
        }
    }

    if (disjoint & !as_point) xc <- st_contour_disjoint(xc) 
    gv <- dplyr::group_vars(xc)
    xc <- dplyr::relocate(xc, dplyr::all_of(gv), .before="geometry")
    xc <- sf::st_sf(xc, crs=sf::st_crs(x$sf))
    
    if (!(oc %in% "kms"))
    {
        xc$contlabel <- factor(droplevels(xc$contlabel), levels=unique(xc$contlabel))
        xc$estimate <- signif(xc$estimate, 3)
        xc$estimate <- factor(xc$estimate) 
    }
     
    ## if original data is not tibble
    if (!inherits(x$sf, "tbl_df"))  xc <- sf::st_sf(data.frame(xc))

    ## convert polygons to point coordinates in tidy format 
    if (as_point)
    {
        gv <- dplyr::group_vars(xc)
        if (length(gv)>0) xc <- dplyr::group_by(xc, dplyr::across(dplyr::all_of(c(gv,"contlabel")))) 
        else xc <- dplyr::group_by(xc, .data$contlabel)
        xc <- dplyr::group_modify(xc, ~dplyr::as_tibble(sf::st_coordinates(.x)))
        xc$contlabel_group <- factor(apply(dplyr::select(dplyr::ungroup(xc), dplyr::num_range("L", range=1:1000)), 1, paste, collapse="."))
        xc <- dplyr::select(xc, -dplyr::num_range("L", range=1:1000))
        xc$contlabel_group <- factor(dplyr::group_indices(dplyr::group_by(xc, .data$contlabel, .data$contlabel_group)))
    }

    if (oc %in% "kdde")
    {
        xc <- dplyr::group_by(xc, .data$deriv_group)
        xc <- dplyr::select(xc, -"pos")
    }
     
    return(xc)
}

## create set of disjoint contour polygons
## x = multipolygon geometry of contour regions

st_contour_disjoint <- function(x)
{
    xc <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contour_disjoint(.x))) 
    xc <- sf::st_sf(xc)
    xc <- sf::st_cast(xc, to="MULTIPOLYGON")
    
    return(xc)
}

.st_contour_disjoint <- function(x)
{
    if (any(as.numeric(as.character(x$contlabel))<0)) oc <- "kdde" else oc <- "kde"
    if (any(oc %in% "kdde"))
    {
        cp <- x
        cp <- dplyr::mutate(cp, estimate=as.numeric(as.character(.data$estimate)), contlabel=as.numeric(as.character(.data$contlabel)))
        if (is.null(cp$deriv_ind)) cp <- dplyr::arrange(cp, .data$estimate) 
        else cp <- dplyr::arrange(cp, .data$deriv_ind, .data$estimate)

        cp1 <- dplyr::filter(cp, .data$contlabel<0)
        cp1 <- sf::st_geometry(cp1)
        cp1 <- st_difference_sequence(cp1, headtail=FALSE)
        cp1 <- sf::st_sfc(cp1)
        cp1 <- c(head(sf::st_geometry(dplyr::filter(cp, .data$contlabel<0)),n=1),cp1)
        cp1 <- sf::st_set_crs(cp1, sf::st_crs(x))
        
        cp2 <- dplyr::filter(cp, .data$contlabel>0)
        cp2 <- sf::st_geometry(cp2)
        cp2 <- st_difference_sequence(cp2, headtail=TRUE)
        cp2 <- sf::st_sfc(cp2)
        cp2 <- c(cp2, tail(sf::st_geometry(dplyr::filter(cp, .data$contlabel>0)),n=1))
        cp2 <- sf::st_set_crs(cp2, sf::st_crs(x))
    
        x <- sf::st_set_geometry(x, value=c(cp1,cp2))
    }
    else
    { 
        cp <- sf::st_geometry(x)
        cp <- st_difference_sequence(cp, headtail=TRUE)
        cp <- sf::st_sfc(cp)
        cp <- c(cp, tail(sf::st_geometry(x),n=1))
        cp <- sf::st_set_crs(cp, sf::st_crs(x))
        x <- sf::st_set_geometry(x, value=cp)
    }
    
    return(x)
}

## compute st_difference on pairs of adjacent rows in sf object
## x = sf object

st_difference_sequence <- function(x, headtail=TRUE) 
{ 
    for (i in 1:max(1,length(x)-1)) 
    {
        xh <- head(x, n=-1)
        xt <- tail(x, n=-1)
        if (headtail) 
        {
            if (i==1) xp <- sf::st_difference(xh[i], xt[i])
            else xp <- c(xp, sf::st_difference(xh[i], xt[i]))
        }
        else
        {
            if (i==1) xp <- sf::st_difference(xt[i], xh[i])
            else xp <- c(xp, sf::st_difference(xt[i], xh[i]))
        }
    }
    
    return(xp)
}

## convert contour polygons in ks object to multipolygon geometry
## x = ks object   

.st_contourline <- function(x, cont, abs.cont, edge_zero=TRUE)
{   
    ## flag to use values of abs.cont in contlabel rather than names(abs.cont)
    if (missing(abs.cont)) abs.cont.flag <- FALSE 
    else if (is.null(names(abs.cont))) abs.cont.flag <- TRUE
    else abs.cont.flag <- !all(grepl("[0-9]%", names(abs.cont)))
    
    if (missing(abs.cont) & missing(cont)) abs.cont <- x$cont
    if (!missing(cont)) 
    {
        if (class(x) %in% c("kfs", "kde.loctest")) abs.cont <- abs.cont
        else abs.cont <- contourLevels(x, cont=cont)
    }
    abs.cont <- sort(abs.cont)
   
    ## convert contourLines output to linestring
    cline2ls <- function(y) 
    { 
        sf::st_linestring(as.matrix(data.frame(x=y["x"], y=y["y"]))) 
    }
    
    ## add zero estimate values at edges of estimation grid 
    if (edge_zero)
    {
        x.zero <- x
        delta <- as.vector(head(sapply(x$eval.points, diff), n=1))
        a <- 1
        x.zero$eval.points[[1]] <- c(min(x.zero$eval.points[[1]]) - delta[1]*(a:1), x.zero$eval.points[[1]], max(x.zero$eval.points[[1]]) + delta[1]*(1:a))
        x.zero$eval.points[[2]] <- c(min(x.zero$eval.points[[2]]) - delta[2]*(a:1), x.zero$eval.points[[2]], max(x.zero$eval.points[[2]]) + delta[2]*(1:a))
        x.zero$estimate <- matrix(0, ncol=ncol(x.zero$estimate)+2*a, nrow=nrow(x.zero$estimate)+2*a)
        x.zero$estimate[(1:nrow(x$estimate))+a, (1:ncol(x$estimate))+a] <- x$estimate
        x <- x.zero
    }
    j <- 0   
   	for (i in 1:length(abs.cont))
    {
        if (is.list(x$estimate)) estimate <- x$estimate[[1]] else estimate <- x$estimate
    	cline <- contourLines(x=x$eval.points[[1]], y=x$eval.points[[2]], z=estimate, levels=abs.cont[i])
        cline.len <- lengths(lapply(cline, getElement, "x"))
        cline <- cline[cline.len>=4]
    	cont.polygon <- lapply(cline, cline2ls)
    	if (length(cont.polygon)==0) cont.polygon <- sf::st_polygon()
        cont.polygon <- sf::st_sfc(cont.polygon)
    	cont.polygon <- sf::st_cast(cont.polygon, to="POLYGON")
        cont.polygon <- cont.polygon[sf::st_is_valid(cont.polygon)]
    	cont.polygon <- dplyr::summarise(sf::st_sf(data.frame(x=1,cont.polygon)))
        
        if (abs.cont.flag) 
            cont.polygon <- cbind(contlabel=signif(abs.cont[i],3), estimate=as.numeric(abs.cont[i]), cont.polygon) 
    	else 
            cont.polygon <- cbind(contlabel=as.numeric(gsub("%", "", names(abs.cont[i]))), estimate=as.numeric(abs.cont[i]), cont.polygon) 
    	
        j <- j + !sf::st_is_empty(cont.polygon)
    	if (j<=1) { cont.polygon.all <- cont.polygon ; j <- j+1 }
        else if (j>1) cont.polygon.all <- rbind(cont.polygon.all, cont.polygon)
    }
     
    cont.polygon.all <- cont.polygon.all[!sf::st_is_empty(cont.polygon.all),]
    if (nrow(cont.polygon.all)>0)
    {
        cont.polygon.all <- sf::st_cast(cont.polygon.all, to="MULTIPOLYGON")
        cont.polygon.all <- sf::st_make_valid(cont.polygon.all)
        cont.polygon.all <- dplyr::distinct(cont.polygon.all)
         
        if (abs.cont.flag) 
        {
            cont.polygon.all$contlabel <- factor(cont.polygon.all$contlabel)
        }
        else
        {
            cont.polygon.all$contlabel <- factor(cont.polygon.all$contlabel)
            levels(cont.polygon.all$contlabel) <- as.character(100-as.numeric(levels(cont.polygon.all$contlabel)))
        }
    }
   
	return(cont.polygon.all)
}

## x = tibble with ks column (output from first part of st_ks)

st_contourline <- function(x, cont, abs.cont, edge_zero=TRUE, ...)
{
    ## convert contour polygons to multipolygon geometry
    if (missing(abs.cont)) cont.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contourline(untidy_ks(.x), cont=cont, edge_zero=edge_zero))) 
    else if (!missing(abs.cont)) cont.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_contourline(untidy_ks(.x), abs.cont=abs.cont, edge_zero=edge_zero))) 
    
    return(cont.geom)
    
}

## x = tibble with ks column (output from first part of st_kdde)

st_contourline_kdde <- function(x, cont=1:99, abs.cont, which_deriv_ind=1)
{
    x.temp <- x
    if (missing(abs.cont)) 
    { 
        abs.cont <- contourLevels(x, cont=cont, which.deriv.ind=which_deriv_ind)
        abs.cont.flag <- FALSE
    }
    else 
    {
        if (!is.vector(abs.cont) & !is.null(abs.cont$breaks)) abs.cont <- abs.cont$breaks
        abs.cont <- suppressWarnings(rbind(abs.cont[abs.cont<0], abs.cont[abs.cont>0]))
        abs.cont.flag <- TRUE
    }
    abs.cont <- sort(c(abs.cont[1,], abs.cont[2,])) 
    
    x.temp$estimate <- x$estimate[[which_deriv_ind]]
    cont.geom <- .st_contourline(x=x.temp, abs.cont=abs.cont) 
    if (!abs.cont.flag) 
    {
        cont.geom <- dplyr::mutate(cont.geom, contlabel=ifelse(.data$estimate<0, paste0("-", as.character(.data$contlabel)), as.character(.data$contlabel)))
        cont.geom$contlabel <- factor(cont.geom$contlabel, levels=c(-sort(cont), sort(cont, decreasing=TRUE)))
    }
    else cont.geom$contlabel <- factor(cont.geom$contlabel)
    return(cont.geom)
}

## x = tibble with ks column (output from first part of st_kda)

st_contourline_kda <- function(x, cont=1:99)
{
    xj <- untidy_ks(x)
    ## replace KDE with weighted KDE
    j <- x$rn
    pp <- xj$prior.prob[j]
    abs.cont <- contourLevels(xj, cont=cont)[[j]]*pp
    xj$estimate <- xj$estimate[[j]]*pp
    xj$H <- xj$H[[j]]
    xj$cont <- xj$cont[[j]]
    class(xj) <- "kde"
    
    cont.geom <- .st_contourline(x=xj, abs.cont=abs.cont) 
    cont.geom <- dplyr::mutate(cont.geom, contlabel=factor(ifelse(.data$estimate<0, paste0("-", 100-as.numeric(as.character(.data$contlabel))), as.character(.data$contlabel))))
    cont.geom <- dplyr::arrange(cont.geom, sign(.data$estimate), abs(.data$estimate))
    levels(cont.geom) <- as.character(cont.geom$contlabel)
       
    return(cont.geom)
}

## convert eval.points in ks object to rectangle polygon/point geometry
## x = tibble with ks column (output from first part of st_ks)

st_evalpoints <- function(x, gridtype="POLYGON")
{
    ## convert ks grid to rectangle polygon geometry
    if (inherits(untidy_ks(x[1,]), "kde.loctest"))
    {
        fhat.geom.neg <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x)$fhat.diff.neg)))
        fhat.geom.pos <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x)$fhat.diff.pos)))
        fhat.geom.neg <- dplyr::mutate(fhat.geom.neg, estimate=as.integer(.data$estimate>=0.5))
        fhat.geom.pos <- dplyr::mutate(fhat.geom.pos, estimate=as.integer(.data$estimate>=0.5))
        fhat.geom <- fhat.geom.pos
        fhat.geom$estimate <- fhat.geom.pos$estimate - fhat.geom.neg$estimate
    }
    else 
        fhat.geom <- dplyr::group_modify(.data=x, .f=~dplyr::tibble(.st_evalpoints(untidy_ks(.x))))
    
    fhat.geom <- sf::st_as_sf(fhat.geom, crs=sf::st_crs(x))
    
    ## convert if required to point  
    if (gridtype %in% "POINT")
    {
        ep <- untidy_ks(x)$eval.points
        if (is.null(ep)) ep <- untidy_ks(x)[[1]]$eval.points
        nx <- lengths(ep)[1]
        ny <- lengths(ep)[2]
    
        fhat.geom <- dplyr::group_modify(.data=fhat.geom, .f=~dplyr::tibble(st_cast_ks(.x, to=gridtype, nx=nx-1, ny=ny-1)))
    }
       
    return(fhat.geom)
}

## convert eval.points in ks object to rectangle polygon geometry
## x = ks object

.st_evalpoints <- function(x)
{
    ## convert KDE to rectangle polygon geometry
    nx <- lengths(x$eval.points)[1]
    ny <- lengths(x$eval.points)[2]
    ep1 <- x$eval.points[[1]]
    ep2 <- x$eval.points[[2]]
    fhat.point <- sf::st_as_sf(x=expand.grid(x=ep1[c(1,nx)], y=ep2[c(1,ny)]), coords=1:2)
    fhat.geom <- sf::st_make_grid(fhat.point, n=c(nx-1,ny-1))
    fhat.geom.centroid <- sf::st_centroid(sf::st_geometry(fhat.geom))
    fhat.geom <- data.frame(estimate=predict(x, x=sf::st_coordinates(fhat.geom.centroid)), fhat.geom)
  
    return(fhat.geom)
}

## convert polygon geometry to point geometry/raster
## x = polygon geometry

st_cast_ks <- function(x, to, nx, ny)
{
    to <- match.arg(toupper(to), c("POINT","RASTER"))
	if (to %in% "POINT") 
 	{
 		## convert to point geometry
 		sf::st_geometry(x) <- sf::st_centroid(sf::st_geometry(x))
 	}
 	else if (to %in% "RASTER")
    { 
        ## convert to raster 
        ## x <- stars::st_rasterize(x, nx=nx, ny=ny)	
        stop("Not yet implemented")
    } 
    
	return(x)
}

## remove segments in linestring where distance between two consecutive points > len 
## (metres) useful to neaten appearance of polygons which are split to separate edges of 
## world map after projection  
## x = linestring or multipoint geometry

st_remove_long_segment <- function(x, len=500e3)
{
    x.geom <- dplyr::group_modify(x, .f=~dplyr::as_tibble(sf::st_coordinates(sf::st_geometry(.))))
    gvar <- setdiff(names(x.geom), c("X","Y")) 
    x.geom <- dplyr::group_by(x.geom, dplyr::across(dplyr::all_of(gvar)))
    x.geom <- dplyr::mutate(x.geom, L1=cumsum(c(0, sqrt(diff(.data$X)^2 + diff(.data$Y)^2) > len)))
    x.geom <- dplyr::group_by(x.geom, dplyr::across(dplyr::all_of(gvar)))
    valid.LS <- names(table(x.geom$L1))[table(x.geom$L1)>2]
    
    if (length(valid.LS)>0) 
    {
        x.geom <- dplyr::filter(x.geom, .data$L1 %in% valid.LS)
        x.geom <- dplyr::summarise(sf::st_as_sf(x.geom, coords=c("X","Y"), crs=sf::st_crs(x)), do_union=FALSE, .groups="drop")
        x.geom <- sf::st_cast(x.geom, to="LINESTRING")
        x.geom <- dplyr::filter(x.geom, as.numeric(sf::st_length(.data$geometry))>0)
        x.geom <- dplyr::group_by(x.geom, dplyr::across(dplyr::all_of(gvar)))
        x.geom <- dplyr::summarise(sf::st_sf(x.geom), do_union=FALSE, .groups="drop")
    
        x2 <- dplyr::right_join(sf::st_drop_geometry(x), x.geom, by=dplyr::group_vars(x))
        x2 <- dplyr::relocate(x2, "geometry", .after=dplyr::last_col())
        x <- sf::st_sf(x2)
    }
    
    return(x)
}

## concatenate sf objects

c_sf <- function(..., labels)
{
    args <- list(...)
    nargs <- length(args)
    if (nargs>=1) 
    {
        tt <- dplyr::mutate(args[[1]], group=1L)
        if (nargs>=2)
        {
            for (i in 2:nargs) 
                tt <- dplyr::bind_rows(tt, dplyr::mutate(args[[i]], group=i))
        }
        if (!missing(labels)) tt$group <- factor(tt$group, labels=labels)
        tt <- dplyr::relocate(tt, "geometry", .after=dplyr::last_col())
        tt <- sf::st_sf(tt)
    }
    else tt <- NULL

    return(tt)
}

## concatenate sf_ks objects

c.sf_ks <- function(..., labels) { c_sf_ks(..., labels=labels) }

c_sf_ks <- function(..., labels)
{
    args <- list(...)
    nargs <- length(args)
    if (nargs>=1) 
    {
        csf <- list()
        csf$tidy_ks <- do.call("c_tidy_ks", lapply(args, getElement, "tidy_ks"))
        csf$grid <- do.call("c_sf", lapply(args, getElement, "grid"))
        csf$sf <- do.call("c_sf", lapply(args, getElement, "sf"))

         if (!missing(labels)) 
        {
            csf$tidy_ks$group <- factor(csf$tidy_ks$group, labels=labels)
            csf$grid$group <- factor(csf$grid$group, labels=labels)
            csf$sf$group <- factor(csf$sf$group, labels=labels)
        }

        csf$tidy_ks <- dplyr::group_by(csf$tidy_ks, .data$group)
        csf$grid <- dplyr::group_by(csf$grid, .data$group)
        csf$sf <- dplyr::group_by(csf$sf, .data$group)
        csf$tidy_ks <- as_tidy_ks(csf$tidy_ks)
    }
    else csf <- NULL
    class(csf) <- "sf_ks"

    return(csf)
}

## add coordinates as attributes

st_add_coordinates <- function(x, as_sf=FALSE, as_tibble=FALSE, rename=TRUE)
{
    xc <- dplyr::mutate(x, .coord=sf::st_coordinates(.data$geometry), .before="geometry")
    xc$X <- xc$.coord[,1]
    xc$Y <- xc$.coord[,2]
    xc <- dplyr::select(xc, -dplyr::all_of(".coord"))

    if (!as_sf)  xc <- sf::st_drop_geometry(xc)
    else  xc <- dplyr::relocate(xc, dplyr::all_of(c("X", "Y")), .before="geometry")
    if (as_tibble)  xc <- dplyr::as_tibble(xc)
    if (rename) xc <- dplyr::rename(xc, lon=dplyr::all_of("X"), lat=dplyr::all_of("Y"))

    return(xc)
}

## transform point geometries to simplified line geometry

st_simplify_point <- function(x, dTolerance)
{
    gv <- dplyr::group_vars(x)
    xs <- dplyr::summarise(x, n=dplyr::n(), do_union=FALSE, .groups="keep")
    xs <- sf::st_cast(xs, to="LINESTRING")
    xs <- dplyr::filter(xs, .data$n>1)  ## remove segments with singleton points
    xs <- sf::st_simplify(xs, dTolerance=dTolerance)
    xs <- dplyr::select(xs, -dplyr::all_of("n"))
    xs <- dplyr::arrange(xs, .data$segment)
    if (length(gv)>0) xs <- dplyr::group_by(xs, dplyr::across(dplyr::all_of(gv)))

    return(xs)
}

## add "%" suffix  

label_percent <- function(y)
{
    factor(y, labels=paste0(levels(y),"%"))
}
