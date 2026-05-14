base_dir <- getwd()
png_files <- list.files(base_dir, pattern = "\\.png$", full.names = TRUE)
if (requireNamespace("magick", quietly = TRUE)) {
  for (f in png_files) {
    out <- sub("\\.png$", ".jpg", f, ignore.case = TRUE)
    img <- magick::image_read(f)
    img <- magick::image_background(img, "white", flatten = TRUE)
    magick::image_write(img, path = out, format = "jpeg", quality = 95)
  }
} else if (requireNamespace("png", quietly = TRUE)) {
  for (f in png_files) {
    out <- sub("\\.png$", ".jpg", f, ignore.case = TRUE)
    img <- png::readPNG(f)
    if (length(dim(img)) == 3 && dim(img)[3] == 4) {
      alpha <- img[, , 4]
      rgb <- img[, , 1:3, drop = FALSE]
      img <- rgb * alpha + (1 - alpha)
    }
    jpeg(out, width = dim(img)[2], height = dim(img)[1], quality = 95)
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    plot.new()
    rasterImage(img, 0, 0, 1, 1)
    dev.off()
  }
}