## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::opts_chunk$set(global.par=TRUE, collapse=TRUE, comment="#>", fig.width=5, fig.height=5, fig.align="center", dpi=96)
options(tibble.print_min=3L, tibble.print_max=3L)

## ----message=FALSE------------------------------------------------------------
library(eks)
library(colorspace)
library(ggplot2)
library(dplyr)

## crabs data 
data(crabs, package="MASS")
crabs2 <- select(crabs, FL, CW)
xlab <- "Frontal lobe size (mm)"
ylab <- "Carapace width (mm)"

## -----------------------------------------------------------------------------
## KDE contour plot + scatter plot
tkde2 <- tidy_kde(crabs2)
gkde2 <- ggplot(tkde2, aes(x=FL, y=CW)) + labs(x=xlab, y=ylab)
gkde2 + geom_point_ks(colour=8) + geom_contour_ks(colour=1)

## -----------------------------------------------------------------------------
## geom_density_2d KDE contour plot + scatter plot
mkde2 <- ggplot(crabs2, aes(x=FL, y=CW))
mkde2 + geom_point(colour=8) + geom_density_2d(colour=1, bins=4)

## -----------------------------------------------------------------------------
## KDE filled contour plot
gkde2 + geom_contour_filled_ks(colour=1) + scale_fill_discrete_sequential(h1=275)

## -----------------------------------------------------------------------------
## KDE continuous colour plot
ggplot(tkde2, aes(x=FL,y=CW)) + geom_raster(aes(fill=estimate), interpolate=TRUE) + 
   labs(x=xlab, y=ylab) + scale_fill_continuous_sequential(h1=275)

## -----------------------------------------------------------------------------
crabs2g <- select(crabs, FL, CW, sp)
crabs2g <- group_by(crabs2g, sp)
tkde2g <- tidy_kde(crabs2g)
gkde2g <- ggplot(tkde2g, aes(x=FL, y=CW, group=sp)) + labs(x=xlab, y=ylab, 
   colour="Species") + scale_colour_manual(values=c(4, 7)) 

## superposed KDE contour plots + scatter plots
gkde2g + geom_point_ks(colour=8) + 
   geom_contour_ks(aes(colour=sp)) + guides(colour=guide_legend(title="Species"))

## ----fig.width=7--------------------------------------------------------------
## facetted KDE contour plots + scatter plots
gkde2g + geom_point_ks(colour=8) + 
   guides(colour=guide_legend(reverse=FALSE)) + 
   geom_contour_ks(aes(colour=sp)) + facet_wrap(~sp) 

## ----fig.width=7--------------------------------------------------------------
## facetted KDE filled contour plots with fixed contour levels for all facets
bkde2g <- contour_breaks(tkde2g)
gkde2g + geom_contour_filled_ks(breaks=bkde2g, colour=1) + 
   scale_fill_discrete_sequential(h1=275) + facet_wrap(~sp)

## ----message=FALSE------------------------------------------------------------
library(sf)

## Grevillea data
data(grevilleasf, package="eks")
grevilleasf <- mutate(grevilleasf, species=factor(species))
paradoxa <- filter(grevilleasf, name %in% "Grevillea paradoxa")
eryngioides <- filter(grevilleasf, name %in% "Grevillea eryngioides")
grevillea_ep <- filter(grevilleasf, name %in% c("Grevillea eryngioides", 
   "Grevillea paradoxa"))
grevillea_ep <- group_by(grevillea_ep, name)
xlim <- c(165479.3, 1096516.3); ylim <- c(6101931, 7255991)

## WA polygon
data(wa, package="eks")
gwa <- geom_sf(data=wa, fill=NA, colour=1)

## -----------------------------------------------------------------------------
## base R scatter plot
plot(st_geometry(wa), xlim=xlim, ylim=ylim)
plot(st_geometry(eryngioides), add=TRUE, col=3, pch=16, cex=0.5)
plot(st_geometry(paradoxa), add=TRUE, col=6, pch=17, cex=0.5)
mapsf::mf_legend(type="symb", val=c("Grevillea eryngioides", "Grevillea paradoxa"), 
   pal=c(3,6), pch=16:17, cex=c(1,1), title="Species", pos="topright")
mapsf::mf_scale(size=200, lwd=4)
## geom_sf scatter plot 
gsc <- ggspatial::annotation_scale(data=data.frame(name="Grevillea paradoxa"), 
   location="br", width_hint=0.2, bar_cols=1)
theme_set(ggthemes::theme_map())
theme_update(legend.position=c(0.99,0.99), legend.justification=c(1,1)) 
ggplot() + gwa + gsc + 
   geom_sf(data=grevillea_ep, aes(colour=name, shape=name), size=0.5) + 
   coord_sf(xlim=xlim, ylim=ylim) + scale_colour_manual(values=c(3, 6)) + 
   guides(colour=guide_legend(title="Species"), shape=guide_legend(title="Species"))

## -----------------------------------------------------------------------------
skde1 <- st_kde(paradoxa)

## base R contour plot
plot(st_geometry(wa), xlim=xlim, ylim=ylim)
plot(st_geometry(paradoxa), add=TRUE, pch=16, col=8, cex=0.5)
plot(skde1, add=TRUE, col=NA, border=1, legend=FALSE)
mapsf::mf_scale(size=200, lwd=4)

## -----------------------------------------------------------------------------
## geom_sf contour plot
gs <- ggplot(skde1) + gwa + gsc
gs + geom_sf(data=paradoxa, col=8, size=0.5) + 
   geom_sf(data=st_get_contour(skde1), colour=1, fill=NA, show.legend=FALSE) + 
   coord_sf(xlim=xlim, ylim=ylim) + ggthemes::theme_map()

## -----------------------------------------------------------------------------
## R base filled contour plot
plot(st_geometry(wa), xlim=xlim, ylim=ylim)
plot(skde1, add=TRUE, pal=function(.) { sequential_hcl(n=., h1=275, rev=TRUE) })
mapsf::mf_scale(size=200, lwd=4)
## geom_sf filled contour
gs + geom_sf(data=st_get_contour(skde1), aes(fill=label_percent(contlabel))) + 
  scale_fill_discrete_sequential(h1=275)  + 
  coord_sf(xlim=xlim, ylim=ylim) + ggthemes::theme_map()

## ----fig.width=7--------------------------------------------------------------
skde1g <- st_kde(grevillea_ep)

## facetted geom_sf filled contour
gsg <- ggplot(skde1g) + gwa + gsc
gsg + geom_sf(data=st_get_contour(skde1g), aes(fill=label_percent(contlabel))) + 
   scale_fill_discrete_sequential(h1=275) + 
   facet_wrap(~name) + coord_sf(xlim=xlim, ylim=ylim) + ggthemes::theme_map()
## facetted geom_sf filled contour with fixed contour levels for all facets
bkde1g <- contour_breaks(skde1g)
gsg + geom_sf(data=st_get_contour(skde1g, breaks=bkde1g), aes(fill=contlabel)) +
   scale_fill_discrete_sequential(h1=275) + facet_wrap(~name) + 
   coord_sf(xlim=xlim, ylim=ylim) + ggthemes::theme_map()

