#' Create a SpatialExperiment object
#'
#' @description This function creates a SpatialExperiment object from a SingleCellExperiment object and a spatial barcode file.
#' @param sce The SingleCellExperiment object obtained from running the \code{\link{sc_long_pipeline}} function.
#' @param spatial_barcode_file The path to the spatial barcode file, e.g. \code{"spaceranger-2.1.1/lib/python/cellranger/barcodes/visium-v2_coordinates.txt"}.
#' @param mannual_align_json The path to the mannual alignment json file.
#' @param image 'DataFrame' containing the image data. See \code{?SpatialExperiment::readImgData} and \code{?SpatialExperiment::SpatialExperiment}.
#' @param tissue_positions_file The path to Visium positions file, e.g. \code{"spaceranger-2.1.1/lib/python/cellranger/barcodes/visium-v2_tissue_positions_list.csv"}.
#' @return A SpatialExperiment object.
#' @importFrom readr read_table
#' @importFrom jsonlite fromJSON
#' @importFrom tidyr as_tibble
#' @importFrom dplyr mutate left_join select
#' @importFrom SummarizedExperiment assays colData rowData rowRanges rowRanges<- colData<-
#' @importFrom SpatialExperiment SpatialExperiment readImgData imgData imgData<-
#' @importFrom methods is
#' @export
create_spe <- function(sce, spatial_barcode_file, mannual_align_json, image, tissue_positions_file) {
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object")
  }

  # Read the full list file
  full_list <- readr::read_table(
    spatial_barcode_file,
    col_names = c("barcode", "col", "row"),
    col_types = "cii"
  )

  # use mannual alignment of image and spots if provided
  if (!missing(mannual_align_json)) {
    align_df <- jsonlite::fromJSON(mannual_align_json)$oligo |>
      tidyr::as_tibble() |>
      dplyr::mutate(row = row + 1, col = col + 1)
    full_list <- dplyr::left_join(align_df, full_list, by = c("row", "col"))
  }

  # add spatial info to colData
  coldata <- full_list |>
    dplyr::select(-barcode) |>
    as.data.frame()
  rownames(coldata) <- full_list$barcode
  coldata <- coldata[colnames(sce), ]
  coldata <- cbind(SummarizedExperiment::colData(sce), coldata)

  # Create a SpatialExperiment object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = SummarizedExperiment::assays(sce),
    colData = coldata,
    rowData = SummarizedExperiment::rowData(sce),
    spatialCoordsNames = c("imageX", "imageY")
  )
  SummarizedExperiment::rowRanges(spe) <- SummarizedExperiment::rowRanges(sce)

  if (!missing(tissue_positions_file)) {
    SummarizedExperiment::colData(spe)$in_tissue <- NULL
    message(sprintf("Reading tissue positions from %s", tissue_positions_file))
    tissue_positions <- readr::read_csv(
      tissue_positions_file,
      col_names = c(
        "barcode", "in_tissue", "array_row",
        "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"
      )
    ) |>
      dplyr::mutate(barcode = stringr::str_remove(barcode, "-1$"))
    SummarizedExperiment::colData(spe)$in_tissue <- SummarizedExperiment::colData(spe) |>
      as.data.frame() |>
      tibble::as_tibble(rownames = "barcode") |>
      dplyr::left_join(tissue_positions, by = "barcode") |>
      dplyr::mutate(in_tissue = in_tissue == 1) |>
      dplyr::pull(in_tissue)
  } else if ("tissue" %in% names(SummarizedExperiment::colData(spe))) {
    SummarizedExperiment::colData(spe)$in_tissue <- SummarizedExperiment::colData(spe)$tissue
    message(sprintf("in_tissue column is set to tissue column from %s", mannual_align_json))
    cat("in_tissue counts:")
    print(table(SummarizedExperiment::colData(spe)$in_tissue, useNA = "always"))
  } else {
    SummarizedExperiment::colData(spe)$in_tissue <- TRUE
    warning("No tissue column in full list file, all spots are considered to be in tissue")
  }

  # add image file if provided
  if (!missing(image)) {
    if (is.character(image)) {
      img_path <- file.path(image, "tissue_hires_image.png")
      if (!file.exists(img_path)) {
        stop(sprintf("Image file not found: %s", img_path))
      }
      image <- SpatialExperiment::readImgData(
        image,
        sample_id = SummarizedExperiment::colData(spe)$sample_id[1],
        imageSources = img_path
      )
    }
    SpatialExperiment::imgData(spe) <- image
  }

  return(spe)
}

#' Plot spatial pie chart
#'
#' @description This function plots a spatial pie chart for given features.
#' @param spe The SpatialExperiment object.
#' @param features The features to plot.
#' @param assay_type The assay that contains the given features.
#' @param opacity The opacity of the background tissue image.
#' @param grayscale Whether to convert the background image to grayscale.
#' @param pie_scale The size of the pie charts.
#' @param color_palette Named vector of colors for each feature.
#' @importFrom ggplot2 ggplot annotation_raster coord_fixed theme_void aes
#' @importFrom scatterpie geom_scatterpie
#' @importFrom magick image_read image_colorize
#' @importFrom grDevices as.raster
#' @importFrom dplyr mutate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment assay
#' @keywords internal
#' @return A ggplot object.
plot_spatial_pie <- function(
    spe, features, assay_type = "counts", color_palette,
    opacity = 50, grayscale = TRUE, pie_scale = 0.8) {
  if (missing(color_palette)) {
    color_palette <- RColorBrewer::brewer.pal(8, "Set2") |>
      head(length(features)) |>
      setNames(features)
  }

  feature <- SummarizedExperiment::assay(spe, assay_type)[features, , drop = FALSE] |>
    as.matrix() |>
    t() |>
    as.data.frame()
  plot_spatial(spe,
    opacity = opacity, grayscale = grayscale,
    feature = feature,
    gglayerFunc = scatterpie::geom_scatterpie,
    aes = ggplot2::aes(x = imageX, y = imageY),
    cols = features, pie_scale = pie_scale, color = NA
  ) +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
}

#' Plot spatial pie chart of isoforms
#'
#' @description This function plots a spatial pie chart for given features.
#' @param spe The SpatialExperiment object.
#' @param isoforms The isoforms to plot.
#' @param assay_type The assay that contains the given features. E.g. 'counts', 'logcounts'.
#' @param color_palette Named vector of colors for each isoform.
#' @param ... Additional arguments to pass to \code{\link{plot_spatial_pie}}, including \code{opacity}, \code{grayscale}, \code{pie_scale}.
#' @return A ggplot object.
#' @importFrom cowplot plot_grid
#' @export
plot_spatial_isoform <- function(spe, isoforms, assay_type = "counts", color_palette, ...) {
  if (missing(color_palette)) {
    color_palette <- RColorBrewer::brewer.pal(8, "Set2") |>
      head(length(isoforms))
  }
  isoform_plot <- plot_isoforms(spe, transcript_ids = isoforms, colors = color_palette)
  pie_plot <- plot_spatial_pie(spe, isoforms,
    assay_type = assay_type,
    color_palette = color_palette, ...
  )
  cowplot::plot_grid(pie_plot, isoform_plot, ncol = 1, rel_heights = c(4, 1))
}

#' Plot feature on spatial image
#'
#' @description This function plots a spatial point plot for given feature
#' @param spe The SpatialExperiment object.
#' @param feature The feature to plot. Could be either a feature name or index 
#' present in the assay or a numeric vector of length nrow(spe).
#' @param assay_type The assay that contains the given features. E.g. 'counts', 'logcounts'.
#' @param opacity The opacity of the background tissue image.
#' @param grayscale Whether to convert the background image to grayscale.
#' @param size The size of the points.
#' @param color The maximum color for the feature. Minimum color is transparent.
#' @param ... Additional arguments to pass to \code{\link[ggplot2]{geom_point}}.
#' @return A ggplot object.
#' @importFrom cowplot plot_grid
#' @export
plot_spatial_feature <- function(
    spe, feature, opacity = 50, grayscale = TRUE, size = 1,
    assay_type = "counts", color = "red", ...) {
  stopifnot("feature must be either length 1 or ncol(spe)" = length(feature) == 1 || length(feature) == ncol(spe))
  if (length(feature) == 1) {
    if (is.character(feature) || is.numeric(feature)) {
      feature <- SummarizedExperiment::assay(spe, assay_type)[feature, , drop = TRUE]
    } else {
      stop(sprintf("Invalid feature type: %s", class(feature)))
    }
  }
  # othwerwise, use the feature as is

  plot_spatial(spe = spe, opacity = opacity, grayscale = grayscale,
    feature = feature, gglayerFunc = ggplot2::geom_point,
    aes = ggplot2::aes(x = imageX, y = imageY, alpha = feature), col = color, size = size
  ) +
    ggplot2::geom_point(
      data = tibble(
        x = SpatialExperiment::spatialCoords(spe)[1, "imageX"],
        y = SpatialExperiment::spatialCoords(spe)[1, "imageY"],
        feature = feature
      ),
      ggplot2::aes(x = x, y = y, col = feature), alpha = 0
    ) +
    guides(alpha = "none") +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_colour_gradient(low = "white", high = color)
}

plot_spatial <- function(
    spe, opacity, grayscale = TRUE, feature, gglayerFunc = ggplot2::geom_point,
    aes = ggplot2::aes(x = imageX, y = imageY, color = feature), ...) {
  stopifnot("No image data found in the SpatialExperiment object" = nrow(imgData(spe)) > 0)
  background_img <- SpatialExperiment::imgData(spe)$data[[1]] |>
    SpatialExperiment::imgRaster() |>
    magick::image_read()
  if (grayscale) {
    background_img <- background_img |>
      magick::image_quantize(colorspace = "gray")
  }
  if (!missing(opacity)) {
    background_img <- background_img |>
      magick::image_colorize(opacity = opacity, color = "white")
  }
  background_img <- grDevices::as.raster(background_img)

  maxX <- dim(background_img)[1]
  maxY <- dim(background_img)[2]
  p1 <- ggplot2::ggplot(mapping = ggplot2::aes(1:maxX, 1:maxY)) +
    ggplot2::annotation_raster(background_img,
      xmin = 1, xmax = maxX, ymin = 1, ymax = maxY
    )

  plot_d <- SpatialExperiment::spatialCoords(spe) |>
    as.data.frame() |>
    dplyr::mutate(imageY = maxY - imageY)
  if (is.null(dim(feature))) {
    plot_d <- cbind(plot_d, feature = feature)
  } else {
    plot_d <- cbind(plot_d, feature)
  }

  p1 +
    gglayerFunc(
      data = plot_d, aes, ...
    )
}
