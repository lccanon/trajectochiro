library(dplyr)
library(purrr)
library(tidyr)

# distance must contain the following columns: source, destination and
# length.
# Returns a vector of ids (sorted by order of appearance).
get_ids <- function(distance) {
  require(dplyr)
  require(purrr)
  require(tidyr)
  distance %>%
    pivot_longer(cols = c("source", "destination"), values_to = "id") %>%
    select(id) %>%
    distinct() %>%
    as_vector()
}

# Transforms tibble in simple vector for associated ids.
get_distance <- function(distance, id1, id2) {
  distance %>%
    filter(source == id1 & destination == id2 |
           source == id2 & destination == id1) %>%
    slice(1) %>%
    select(distance) %>%
    as.numeric()
}

# Transforms tibble in simple vector for associated ids.
get_position <- function(position, iden) {
  require(purrr)
  position %>%
    filter(id == iden) %>%
    slice(1) %>%
    select(x, y, z) %>%
    as_vector()
}

# Remove half of the most extreme values.
get_most_common <- function(val, length.out = length(val) / 2) {
  if (length(val) == 1)
    return(val)
  val <- sort(val)
  CS <- cumsum(val)
  CS2 <- cumsum(val ^ 2)
  i_start <- head(seq_along(val), -length.out + 1)
  i_end <- tail(seq_along(val), -length.out + 1)
  sd <- sqrt((CS2[i_end] - c(0, CS2)[i_start]) / length.out -
             ((CS[i_end] - c(0, CS)[i_start]) / length.out) ^ 2 + 1e-10)
  i <- which.min(sd)
  val[i_start[i]:i_end[i]]
}

# position1 must contain the following columns: id, x, y and z.
# Either compute the distances between each location pairwise or the
# distances between the locations from position1 and the ones from
# position2 (if specified).
# Returns a tibble with columns: source, destination and distance.
position2distance <- function(position1, position2 = NULL) {
  require(dplyr)
  distance <- list()
  pos1 <- position1 %>% select(-id)
  internal <- is.null(position2)
  if (internal)
    position2 <- position1
  pos2 <- position2 %>% select(-id)
  for (i in seq_len(nrow(position1)))
    for (j in seq_len(nrow(position2)))
      if (!internal || internal && i < j) {
        s <- position1[i,"id"][[1]]
        d <- position2[j,"id"][[1]]
        diff <- pos1[i,] - pos2[j,]
        dist <- sqrt(sum(diff ^ 2))
        distance[[length(distance) + 1]] <- list(source = s,
                                                 destination = d,
                                                 distance = dist)
      }
  bind_rows(distance)
}

# Inspired from https://math.stackexchange.com/a/2969614/137036
trilateration <- function(p1, p2, p3, d1, d2, d3) {
  e1 <- p2 - p1
  h <- sqrt(sum(e1 ^ 2))
  e1 <- e1 / h

  i <- sum(e1 * (p3 - p1))
  e2 <- (p3 - p1) - i * e1
  t <- sqrt(sum(e2 ^ 2))
  e2 <- e2 / t

  e3 <- c(x = e1[2] * e2[3] - e1[3] * e2[2],
          y = e1[3] * e2[1] - e1[1] * e2[3],
          z = e1[1] * e2[2] - e1[2] * e2[1])

  j <- sum(e2 * (p3 - p1))
  u <- (d1 ^ 2 - d2 ^ 2 + h ^ 2) / (2 * h)
  v <- (d1 ^ 2 - d3 ^ 2 + i * (i - 2 * u) + j ^ 2) / (2 * j)
  if (d1 ^ 2 - u ^ 2 - v ^ 2 < 0)
    return(NaN)
  w <- sqrt(d1 ^ 2 - u ^ 2 - v ^ 2)

  r1 <- p1 + u * e1 + v * e2 + w * e3
  r2 <- p1 + u * e1 + v * e2 - w * e3
  list(r1 = r1, r2 = r2)
}

# distance must contain the following columns: source, destination and
# distance.
# The first id appearing in the distances correspond to the origin,
# the second to a point on the x-axis, the third to the plane x-y and
# the fourth to a point with positive z-value. Each distance must be
# known between each of the first 3 points. Each additional point must
# contain a distance to at least 4 previous points.
distance2position <- function(distance) {
  require(dplyr)
  ids <- get_ids(distance)

  # Get distances to compute first four points
  d12 <- get_distance(distance, ids[1], ids[2])
  d13 <- get_distance(distance, ids[1], ids[3])
  d23 <- get_distance(distance, ids[2], ids[3])
  d14 <- get_distance(distance, ids[1], ids[4])
  d24 <- get_distance(distance, ids[2], ids[4])
  d34 <- get_distance(distance, ids[3], ids[4])
  # Position for first point
  p1 <- c(x = 0, y = 0, z = 0)
  # Position for second point
  p2 <- c(x = d12, y = 0, z = 0)
  # Position for third point
  p3 <- c(x = (d12 ^ 2 + d13 ^ 2 - d23 ^ 2) / (2 * d12),
          y = sqrt(- d23 ^ 4 + d23 ^ 2 * (2 * d12 ^ 2 + 2 * d13 ^ 2) -
               d13 ^ 4 + 2 * d12 ^ 2 * d13 ^ 2 - d12 ^ 4) / (2 * d12), z = 0)
  # Position for fourth point (by convention, it has a positive z)
  p4 <- trilateration(p1, p2, p3, d14, d24, d34)[[1]]

  # Compute all remaining points
  position <- tibble(id = ids[1:4], bind_rows(p1, p2, p3, p4))
  for (i in 5:length(ids)) {
    ids_curr <- head(ids, i)
    dist <- distance %>% filter(source == ids[i] & destination %in% head(ids, i) |
                                destination == ids[i] & source %in% head(ids, i))
    stopifnot(nrow(dist) >= 3)
    comb <- combn(nrow(dist), 3)
    candidate <- list()
    for (j in seq_len(ncol(comb))) {
      co <- comb[,j]
      co.ids <- dist[co,c("source", "destination")] %>% t() %>% as.vector()
      co.ids <- co.ids[co.ids != ids[i]]
      p1 <- get_position(position, co.ids[1])
      p2 <- get_position(position, co.ids[2])
      p3 <- get_position(position, co.ids[3])
      d1 <- get_distance(distance, dist[co[1],"source"][[1]],
                         dist[co[1],"destination"][[1]])
      d2 <- get_distance(distance, dist[co[2],"source"][[1]],
                         dist[co[2],"destination"][[1]])
      d3 <- get_distance(distance, dist[co[3],"source"][[1]],
                         dist[co[3],"destination"][[1]])
      ps <- trilateration(p1, p2, p3, d1, d2, d3)
      if (!any(is.na(ps)))
        candidate <- c(candidate, ps)
      else
        print(paste0("Warning: inconsistent distances for ", ids[i]))
    }
    candidate <- bind_rows(candidate)
    # If there are only 3 distance for a value, we keep the position
    # with positive z-value (assuming this is above ground)
    if (nrow(candidate) == 2)
      candidate <- candidate[which.max(candidate$z),]
    # Otherwise, we combine all the resulting positions
    x <- mean(get_most_common(candidate$x))
    y <- mean(get_most_common(candidate$y))
    z <- mean(get_most_common(candidate$z))
    position[nrow(position) + 1,] <- list(id = ids[i], x = x, y = y, z = z)
  }
  position
}

# Failed attempt: algorithm from wikipedia seems incorrect
## Inspired from https://en.wikipedia.org/wiki/Multilateration#Cartesian_solution_with_limited_computational_resources
#multilaterate <- function(p0, p1, p2, p3, p4, T1, T2, T3, T4,
#                          sound_speed_ms) {
#  p1 <- p1 - p0
#  p2 <- p2 - p0
#  p3 <- p3 - p0
#  p4 <- p4 - p0
#  cT1 <- sound_speed_ms * T1
#  cT2 <- sound_speed_ms * T2
#  cT3 <- sound_speed_ms * T3
#  cT4 <- sound_speed_ms * T4
#  ABC2 <- 2 * p2 / cT2 - 2 * p1 / cT1
#  ABC3 <- 2 * p3 / cT3 - 2 * p1 / cT1
#  ABC4 <- 2 * p4 / cT4 - 2 * p1 / cT1
#  D2 <- cT2 - cT1 - sum(p2 ^ 2) / cT2 + sum(p1 ^ 2) / cT1
#  D3 <- cT3 - cT1 - sum(p3 ^ 2) / cT3 + sum(p1 ^ 2) / cT1
#  D4 <- cT4 - cT1 - sum(p4 ^ 2) / cT4 + sum(p1 ^ 2) / cT1
#  solve(rbind(ABC2, ABC3, ABC4), -c(D2, D3, D4)) + p0
#}
#
## TDOA must contain the following columns: id, recorder and
## TDOA.
## recorder must contain the following columns: id, x, y and z.
## The TDOA must be of size at least 5 for each measure.
#TDOA2position <- function(TDOA, recorder, sound_speed_ms = 343) {
#  require(dplyr)
#  require(purrr)
#  ids <- TDOA %>%
#    select(id) %>%
#    distinct() %>%
#    as_vector()
#  pos <- list()
#  for (iden in ids) {
#    diff_iden <- TDOA %>% filter(id == iden)
#    stopifnot(nrow(diff_iden) >= 5)
#    comb <- combn(nrow(diff_iden), 5)
#    candidate <- list()
#    for (j in seq_len(ncol(comb))) {
#      diff <- diff_iden %>%
#        filter(recorder %in% diff_iden$recorder[comb[,j]]) %>%
#        mutate(TDOA = TDOA - min(TDOA)) %>%
#        arrange(TDOA)
#      p0 <- get_position(recorder, diff$recorder[1])
#      p1 <- get_position(recorder, diff$recorder[2])
#      p2 <- get_position(recorder, diff$recorder[3])
#      p3 <- get_position(recorder, diff$recorder[4])
#      p4 <- get_position(recorder, diff$recorder[5])
#      T1 <- diff$TDOA[2]
#      T2 <- diff$TDOA[3]
#      T3 <- diff$TDOA[4]
#      T4 <- diff$TDOA[5]
#      ps <- multilaterate(p0, p1, p2, p3, p4, T1, T2, T3, T4, sound_speed_ms)
#      candidate[[length(candidate) + 1]] <- ps
#    }
#    candidate <- bind_rows(candidate)
#    x = median(candidate$x)
#    y = median(candidate$y)
#    z = median(candidate$z)
#    pos[[length(pos) + 1]] <- list(id = iden, x = x, y = y, z = z)
#  }
#  bind_rows(pos)
#}

# TDOA must contain the following columns: id, recorder and
# TDOA.
# recorder must contain the following columns: id, x, y and z.
# The TDOA must be of size at least 5 for each measure.
#TODO return confidence value or variance for each position
TDOA2position <- function(TDOA, recorders, sound_speed_ms = 343) {
  require(dplyr)
  require(purrr)
  ids <- TDOA %>%
    select(id) %>%
    distinct() %>%
    as_vector()
  call <- list()
  pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
  setTxtProgressBar(pb, 0)
  for (iden in ids) {
    diff <- TDOA %>%
      filter(id == iden)
    closest <- (diff %>% filter(TDOA == min(TDOA)))[1,"recorder"][[1]]
    pos <- recorders %>%
      filter(id == closest) %>%
      select(-id) %>%
      as_vector()
    # TODO rely on CPP to optimize this optimization
    error_call <- function(pos, iden, TDOA, recorders, sound_speed_ms = 343) {
      require(dplyr)
      position <- tibble(id = iden, x = pos[1], y = pos[2], z = pos[3])
      position2distance(position, recorders) %>%
        mutate(id = source) %>%
        mutate(recorder = destination) %>%
        right_join(TDOA, by = c("id", "recorder")) %>%
        mutate(diff = distance / sound_speed_ms) %>%
        group_by(source) %>%
        mutate(diff = diff - min(diff)) %>%
        mutate(TDOA = TDOA - min(TDOA)) %>%
        ungroup() %>%
        mutate(err = (TDOA - diff) ^ 2) %>%
        select(err) %>%
        sum() %>%
        sqrt()
    }
    res <- optim(par = pos, fn = error_call, gr = NULL,
                 iden, diff, recorders, sound_speed_ms)
    call[[length(call) + 1]] <- list(id = iden, x = res$par[1],
                                     y = res$par[2], z = res$par[3])
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  bind_rows(call)
}

# Compute speeds between successions of positions
#TODO avoid computing all distances pairwise
speed_locations <- function(positions) {
  if (nrow(positions) <= 1)
    return(NULL)
  require(dplyr)
  positions %>%
    arrange(id) %>%
    position2distance() %>%
    group_by(source) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(speed = distance / (destination - source) * 3.6) %>%
    select(-distance)
}

# Remove positions that are too far away or too fast between two
# consecutive locations.
#TODO also remove when variance is high
clean_locations <- function(positions, recorders, speeds,
                            max_distance_m, max_speed_kmh) {
  require(dplyr)
  positions.close <- position2distance(positions, recorders) %>%
    filter(distance <= max_distance_m)
  if (!is.null(speeds))
    positions.slow <- speeds %>%
      filter(speed <= max_speed_kmh)
  else
    positions.slow <- positions.close
  positions %>%
    filter(id %in% positions.close$source) %>%
    filter(id %in% positions.slow$source |
           id %in% positions.slow$destination)
}
