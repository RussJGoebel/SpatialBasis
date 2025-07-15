#' Title
#'
#' @param map
#' @param position
#' @param pal
#' @param values
#' @param na.label
#' @param bins
#' @param colors
#' @param opacity
#' @param labels
#' @param labFormat
#' @param title
#' @param className
#' @param layerId
#' @param group
#' @param data
#' @param decreasing
#'
#' @return
#' @export
#'
#' @examples
addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft",
                                                    "topleft"), pal, values, na.label = "NA", bins = 7, colors,
                                  opacity = 0.5, labels = NULL, labFormat = leaflet::labelFormat(),
                                  title = NULL, className = "info legend", layerId = NULL,
                                  group = NULL, data = leaflet::getMapData(map), decreasing = FALSE) {
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors))
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula"))
      title <- deparse(values[[2]])
    values <- leaflet::evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] ==
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins))
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1)
        pretty(values, bins)
      else bins

      if (length(bins) > 2)
        if (!all(abs(diff(bins, differences = 2)) <=
                 sqrt(.Machine$double.eps)))
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")

    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }

    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2,
                       na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values)))
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels))
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)),
                 na_color = na.color, na_label = na.label, opacity = opacity,
                 position = position, type = type, title = title, extra = extra,
                 layerId = layerId, className = className, group = group)
  leaflet::invokeMethod(map, data, "addLegend", legend)
}


#' make_sif_color_palette
#'
#' A function for creating a leaflet numeric color palette
#'
#' @param domain (see `leaflet::colorNumeric`): this can be a simple numeric range (e.g. c(0,100)).
#' @param na.color (see `leaflet::colorNumeric`): The color to return for NA values. Note that na.color = NA is valid.
#'
#' @return
#' @export
#'
#' @examples
make_sif_color_palette <- function(domain,
                                   na.color = "transparent"){

  sif_color_palette <- leaflet::colorNumeric(c( '#FEFBE9',
                                   '#FCF7D5',
                                   '#F5F3C1',
                                   '#EAF0B5',
                                   '#DDECBF',
                                   '#D0E7CA',
                                   '#C2E3D2',
                                   '#B5DDD8',
                                   '#A8D8DC',
                                   '#9BD2E1',
                                   '#8DCBE4',
                                   '#81C4E7',
                                   '#7BBCE7',
                                   '#7EB2E4',
                                   '#88A5DD',
                                   '#9398D2',
                                   '#9B8AC4',
                                   '#9D7DB2',
                                   '#9A709E',
                                   '#906388',
                                   '#805770',
                                   '#684957',
                                   '#46353A'),
                                domain,
                                na.color="transparent")

  return(sif_color_palette)

}




#' Title
#'
#' @param sf_object An sf object with at least one numeric column.
#' @param column Either the name of a column in the dataset or a number representing which column stores the values.
#' @param pal A function (as returned by leaflet::colorNumeric) that works as a color palette
#' @param group_name A character string denoting the name of the layer for the sf_object
#' @param generate_plot If TRUE, render the final plot.
#'
#' @return
#' @export
#'
#' @examples
#'
#' leaflet_visualization(soundings,"SIF_757nm",group_name = "Soundings")
#' leaflet_visualization(target_grid,"mean_value", title = "Mean Albedo")
#' leaflet_visualization(target_grid,"dominant_class")
#'
leaflet_visualization <- function(sf_object,
                                  column,
                                  pal = make_sif_color_palette(range(sf_object[[column]])),
                                  group_name = "Target",
                                  title = "SIF<sub>757</sub> [mW m<sup>-2</sup>nm<sup>-1</sup>sr<sup>-1</sup>]",
                                  generate_plot = TRUE
                                  ){


  m <- leaflet::leaflet() %>%
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron, group = 'Positron') %>%
    leaflet::addProviderTiles(leaflet::providers$OpenStreetMap, group = 'OpenStreetMap') %>%
    leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = 'Esri')

  m <- m %>%
    leaflet::addPolygons(data=sf_object,
                         weight = .04,
                         color = "grey",
                         fillOpacity = 0.6,
                         fillColor = pal(sf_object[[column]]),
                         group = group_name) %>%
    #addRectangles(boston_bbox[1],boston_bbox[2],boston_bbox[3],boston_bbox[4], fill = F, color = "grey") %>%
    addLegend_decreasing("topright",
                         pal = pal,
                         values = sf_object[[column]],
                         opacity =  0.8,
                         decreasing = T,
                         title = title,
                         group="All") %>%
    leaflet::addScaleBar(position="bottomright")%>%
    leaflet::addLayersControl(baseGroups = c('Positron', 'OpenStreetMap','Esri'), overlayGroups = c(group_name))
  # addControl(title, position = "bottomleft", className="map-title")

  if(generate_plot){
    print(m)
  }


  m


}


#interactive_sif_plot_with_landcover <- function()




