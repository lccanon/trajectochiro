# Ubuntu dependencies:
# sudo apt install libcgal-dev libglu1-mesa-dev libglu1-mesa-dev
# (maybe not required for WebGL)

library(rgl)
require(manipulateWidget)

plot_trajectory <- function(positions, recorders, speeds,
                            min_distance_m, min_speed_kmh,
                            filename = NULL) {
  # Sort according to time for representation
  positions <- positions %>% arrange(id)

  # Hide animation when two successive points seem unrelated
  positions.alpha <- rep(1, nrow(positions))
  if (!is.null(speeds)) {
    positions.slow <- speeds %>%
      filter(speed <= min_speed_kmh | (destination - source) *
             speed / 3.6 >= min_distance_m)
    positions.alpha[positions$id %in% positions.slow$source |
                    positions$id %in% positions.slow$destination] <- 0
  }

  # Recorder and call positions
  recid <- plot3d(recorders %>% select(-id), size = 10)["data"]
  lineid <- plot3d(positions %>% select(-id), type = "l", size = 5,
                   col = rainbow(nrow(positions)),
                   alpha = positions.alpha, add = TRUE)["data"]
  pointid <- plot3d(positions %>% select(-id), type = "p", size = 5,
                    col = rainbow(nrow(positions)), add = TRUE)["data"]
  sphereid <- spheres3d(positions %>% select(-id) %>% slice(1),
                        radius = 0.2)

  # Real-time rendering
  w <- rglwidget() %>%
    playwidget(list(ageControl(births = 0, ages = positions$id,
                               vertices = positions %>% select(-id),
                               objids = sphereid,
                               colors = rainbow(nrow(positions)),
                               alpha = positions.alpha)),
               start = floor(min(positions$id)),
               stop = ceiling(max(positions$id)), rate = 1,
               components = c("Reverse", "Play", "Slower", "Faster",
                              "Reset", "Slider", "Label"),
               loop = TRUE)

  if (!is.null(filename)) {
    htmlwidgets::saveWidget(w, filename, selfcontained = FALSE,
                            libdir = "lib")
    close3d()
    browseURL(filename)
  }
}
