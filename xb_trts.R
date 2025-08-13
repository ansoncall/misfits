# identifying multi-polygon xb treatments

# load packages ####
library(mapview)
library(sf)
library(tidyverse)
library(nngeo)
library(tictoc)

# load data ####
# treatments
treatments_raw <- st_read("raw/LarimerCoTreatments.gdb",
                          fid_column_name = "OBJECTID") %>%
  st_make_valid
# conservation areas
comap_manager_detail <- st_read("raw/misfits.gdb",
                                layer = "COMaP_LarimerCo_Manager")

# make lines for COMaP "Manager Detail" attribute
comap_manager_lines <- st_cast(comap_manager_detail, "MULTIPOLYGON") %>%
  st_cast("MULTILINESTRING")

# wrangle ####
# id single polygon xb treatments

treatments <- treatments_raw
comap_lines <- comap_manager_lines
buffer_distance <- -15
# step 1. buffer treatments by a negative buffer_distance
treatments_buffered <- st_buffer(treatments_raw, dist = buffer_distance)
# step 2. identify which *buffered* treatments intersect comap boundaries
treatments_comap <- st_intersection(treatments_buffered,
                                    comap_lines)
# id single poly xb treatments
treatments_comap_polys <- treatments_buffered %>%
  filter(OBJECTID %in% treatments_comap$OBJECTID)

treatments <- treatments %>%
  mutate(single_poly_xb = case_when(
    OBJECTID %in% treatments_comap$OBJECTID ~ 1,
    .default = 0
  ))


# assign treatments to "home" comap polygon
st_centroid_within_poly <- function (poly) {

  # check if centroid is in polygon
  ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
  in_poly <- diag(st_within(ctrd, poly, sparse = F))

  # replace geometries that are not within polygon with st_point_on_surface()
  st_geometry(ctrd[!in_poly,]) <- st_geometry(st_point_on_surface(poly[!in_poly,]))

  ctrd
}

trts_centroids <- map(seq_len(nrow(treatments)), \(x) {
  y <- treatments[x, ]
  st_centroid_within_poly(y)
  }) %>%
  do.call(rbind, .)

comap <- comap_manager_detail %>%
  mutate(id =  row_number()) %>%
  st_cast("MULTIPOLYGON")


trts_centroids_join <- st_join(trts_centroids, comap, st_intersects)

treatments_join <- left_join(treatments,
                             trts_centroids_join %>% select(OBJECTID, id) %>% st_drop_geometry,
                             by = "OBJECTID")

# TODO fix edge case: one multipart polygon not assigned to comap region because
# part of it falls outside the county
# mapview(comap, col.regions = "red") +
#   mapview(treatments_join[1607, ]) +
#   mapview(trts_centroids_join[1607, ])
treatments_join$id[1607] <- 57

# # find NEAREST neighboring treatment that's not in the same comap boundary
# treatments_join$closest_OBJECTID <- NULL
# treatments_join$closest_dist <- NULL
# treatments_join[1607, ]
# for (i in 1:nrow(treatments_join)) {
#   # i <- 1
#   home_comap_poly <- treatments_join$id[i]
#   target_trt <- treatments_join[i, ]
#   filtered <- treatments_join %>% filter(id != home_comap_poly)
#   closest <- st_nearest_feature(target_trt, filtered)
#   c_dist <- st_distance(target_trt, filtered[closest, ])
#   treatments_join$closest_OBJECTID[i] <- filtered[closest, "OBJECTID"][[1]]
#   treatments_join$closest_dist[i] <- c_dist
#   print(i)
# }
#
# mapview(comap, col.regions = "red") + mapview(treatments_join)
#
# treatments_join %>% filter(is.na(id))

treatments_join$closest_same_yr_OBJECTID_1 <- NULL
treatments_join$closest_same_yr_dist_1 <- NULL
treatments_join$closest_same_yr_OBJECTID_2 <- NULL
treatments_join$closest_same_yr_dist_2 <- NULL
treatments_join$closest_same_yr_OBJECTID_3 <- NULL
treatments_join$closest_same_yr_dist_3 <- NULL
treatments_join$closest_plus1_OBJECTID_1 <- NULL
treatments_join$closest_plus1_dist_1 <- NULL
treatments_join$closest_plus1_OBJECTID_2 <- NULL
treatments_join$closest_plus1_dist_2 <- NULL
treatments_join$closest_plus1_OBJECTID_3 <- NULL
treatments_join$closest_plus1_dist_3 <- NULL
treatments_join$closest_plus2_OBJECTID_1 <- NULL
treatments_join$closest_plus2_dist_1 <- NULL
treatments_join$closest_plus2_OBJECTID_2 <- NULL
treatments_join$closest_plus2_dist_2 <- NULL
treatments_join$closest_plus2_OBJECTID_3 <- NULL
treatments_join$closest_plus2_dist_3 <- NULL


tic()
for (i in 1:nrow(treatments_join)) {
  # i <- 1
  home_comap_poly <- treatments_join$id[i]
  target_trt <- treatments_join[i, ]
  target_yr <- target_trt$YEAR_COMP

  filtered_same_yr <- treatments_join %>% filter(id != home_comap_poly,
                                                 YEAR_COMP == target_yr)
  filtered_plus1_yr <- treatments_join %>% filter(id != home_comap_poly,
                                                  YEAR_COMP == target_yr + 1 |
                                                    YEAR_COMP == target_yr - 1)
  filtered_plus2_yr <- treatments_join %>% filter(id != home_comap_poly,
                                                  YEAR_COMP == target_yr + 2 |
                                                    YEAR_COMP == target_yr - 2)

  matches_yr <- nrow(filtered_same_yr)
  matches_plus1 <- nrow(filtered_plus1_yr)
  matches_plus2 <- nrow(filtered_plus2_yr)

  if (matches_yr == 0) {
    treatments_join$closest_same_yr_OBJECTID_1[i] <- NA
    treatments_join$closest_same_yr_dist_1[i] <- NA
    treatments_join$closest_same_yr_OBJECTID_2[i] <- NA
    treatments_join$closest_same_yr_dist_2[i] <- NA
    treatments_join$closest_same_yr_OBJECTID_3[i] <- NA
    treatments_join$closest_same_yr_dist_3[i] <- NA

  } else if (matches_yr > 0 && matches_yr < 3) {
    closest_n_same_yr <- st_nn(target_trt,
                               filtered_same_yr,
                               k = matches_yr,
                               returnDist = TRUE)
    # hacky way to return OIDs instead of row idx
    # TODO clean up
    closest_n_same_yr$nn[[1]] <- st_drop_geometry(filtered_same_yr[closest_n_same_yr$nn[[1]], "OBJECTID"])[[1]]
  } else {
    closest_n_same_yr <- st_nn(target_trt,
                               filtered_same_yr,
                               k = 3,
                               returnDist = TRUE)
    closest_n_same_yr$nn[[1]] <- st_drop_geometry(filtered_same_yr[closest_n_same_yr$nn[[1]], "OBJECTID"])[[1]]
  }

  if (matches_plus1 == 0) {
    treatments_join$closest_plus1_OBJECTID_1[i] <- NA
    treatments_join$closest_plus1_dist_1[i] <- NA
    treatments_join$closest_plus1_OBJECTID_2[i] <- NA
    treatments_join$closest_plus1_dist_2[i] <- NA
    treatments_join$closest_plus1_OBJECTID_3[i] <- NA
    treatments_join$closest_plus1_dist_3[i] <- NA
  } else if (matches_plus1 > 0 && matches_plus1 < 3) {
    closest_n_plus1 <- st_nn(target_trt,
                             filtered_plus1_yr,
                             k = matches_plus1,
                             returnDist = TRUE)
    closest_n_plus1$nn[[1]] <- st_drop_geometry(filtered_plus1_yr[closest_n_plus1$nn[[1]], "OBJECTID"])[[1]]
  } else {
    closest_n_plus1 <- st_nn(target_trt,
                             filtered_plus1_yr,
                             k = 3,
                             returnDist = TRUE)
    closest_n_plus1$nn[[1]] <- st_drop_geometry(filtered_plus1_yr[closest_n_plus1$nn[[1]], "OBJECTID"])[[1]]
  }

  if (matches_plus2 == 0) {
    treatments_join$closest_plus2_OBJECTID_1[i] <- NA
    treatments_join$closest_plus2_dist_1[i] <- NA
    treatments_join$closest_plus2_OBJECTID_2[i] <- NA
    treatments_join$closest_plus2_dist_2[i] <- NA
    treatments_join$closest_plus2_OBJECTID_3[i] <- NA
    treatments_join$closest_plus2_dist_3[i] <- NA

  } else if (matches_plus2 > 0 && matches_plus2 < 3){
    closest_n_plus2 <- st_nn(target_trt,
                             filtered_plus2_yr,
                             k = matches_plus2,
                             returnDist = TRUE)
    closest_n_plus2$nn[[1]] <- st_drop_geometry(filtered_plus2_yr[closest_n_plus2$nn[[1]], "OBJECTID"])[[1]]
  } else {
    closest_n_plus2 <- st_nn(target_trt,
                             filtered_plus2_yr,
                             k = 3,
                             returnDist = TRUE)
    closest_n_plus2$nn[[1]] <- st_drop_geometry(filtered_plus2_yr[closest_n_plus2$nn[[1]], "OBJECTID"])[[1]]
  }

  treatments_join$closest_same_yr_OBJECTID_1[i] <- closest_n_same_yr$nn[[1]][[1]]
  treatments_join$closest_same_yr_dist_1[i] <- closest_n_same_yr$dist[[1]][[1]]
  treatments_join$closest_same_yr_OBJECTID_2[i] <- closest_n_same_yr$nn[[1]][[2]]
  treatments_join$closest_same_yr_dist_2[i] <- closest_n_same_yr$dist[[1]][[2]]
  treatments_join$closest_same_yr_OBJECTID_3[i] <- closest_n_same_yr$nn[[1]][[3]]
  treatments_join$closest_same_yr_dist_3[i] <- closest_n_same_yr$dist[[1]][[3]]

  treatments_join$closest_plus1_OBJECTID_1[i] <- closest_n_plus1$nn[[1]][[1]]
  treatments_join$closest_plus1_dist_1[i] <- closest_n_plus1$dist[[1]][[1]]
  treatments_join$closest_plus1_OBJECTID_2[i] <- closest_n_plus1$nn[[1]][[2]]
  treatments_join$closest_plus1_dist_2[i] <- closest_n_plus1$dist[[1]][[2]]
  treatments_join$closest_plus1_OBJECTID_3[i] <- closest_n_plus1$nn[[1]][[3]]
  treatments_join$closest_plus1_dist_3[i] <- closest_n_plus1$dist[[1]][[3]]

  # TODO fix subscript out of bounds when matches > 0 && matches < 3
  treatments_join$closest_plus2_OBJECTID_1[i] <- closest_n_plus2$nn[[1]][[1]]
  treatments_join$closest_plus2_dist_1[i] <- closest_n_plus2$dist[[1]][[1]]
  treatments_join$closest_plus2_OBJECTID_2[i] <- closest_n_plus2$nn[[1]][[2]]
  treatments_join$closest_plus2_dist_2[i] <- closest_n_plus2$dist[[1]][[2]]
  treatments_join$closest_plus2_OBJECTID_3[i] <- closest_n_plus2$nn[[1]][[3]]
  treatments_join$closest_plus2_dist_3[i] <- closest_n_plus2$dist[[1]][[3]]

  print(i)
}
toc()
test <- treatments_join %>% filter(OBJECTID %in% c(1, 1244, 1245, 1246))
mapview(test)
# TODO graph

# graph: all data
treatments_join %>%
  pivot_longer(contains("_dist_"), names_to = "type", values_to = "distance") %>%
  mutate(
    year_slice = case_when(
      str_detect(type, "_same_yr_") ~ "same",
      str_detect(type, "_plus1_") ~ "plus1",
      str_detect(type, "_plus2_") ~ "plus2",
      .default = "ERROR"
    ),
    neighbor = case_when(
      str_detect(type, "_dist_1") ~ "first",
      str_detect(type, "_dist_2") ~ "second",
      str_detect(type, "_dist_3") ~ "third",
      .default = "ERROR_neighbor"
    )
  ) %>%
  ggplot(aes(x = distance, fill = neighbor)) +
  geom_histogram(alpha = 0.6, position = "dodge") +
  facet_wrap(~year_slice, nrow = 3) +
  scale_fill_manual(values = c("#162521", "#4F7CAC", "#C0E0DE")) +
  theme_minimal() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 14))

# graph: less than 2.5 kms
treatments_join %>%
  pivot_longer(contains("_dist_"), names_to = "type", values_to = "distance") %>%
  mutate(
    year_slice = case_when(
      str_detect(type, "_same_yr_") ~ "same",
      str_detect(type, "_plus1_") ~ "plus1",
      str_detect(type, "_plus2_") ~ "plus2",
      .default = "ERROR"
    ),
    neighbor = case_when(
      str_detect(type, "_dist_1") ~ "first",
      str_detect(type, "_dist_2") ~ "second",
      str_detect(type, "_dist_3") ~ "third",
      .default = "ERROR_neighbor"
    )
  ) %>%
  filter(distance < 2500) %>%
  ggplot(aes(x = distance, fill = neighbor)) +
  geom_histogram(alpha = 0.6, position = "dodge") +
  facet_wrap(~year_slice, nrow = 3) +
  scale_fill_manual(values = c("#162521", "#4F7CAC", "#C0E0DE")) +
  theme_minimal() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 14))

# map out some treatments where the neighbor is really close
potential_xb <- treatments_join %>%
  filter(closest_same_yr_dist_1 < 20 |
           closest_plus1_dist_1 < 20 |
           closest_plus2_dist_1 < 20)

mapview(comap_manager_lines, color = "red") +
  mapview(potential_xb,
          zcol = "OBJECTID",
          col.regions =
            c(
            rep(hcl.colors(palette = "Dynamic", n = 5), 32),
            hcl.colors(palette = "Dynamic", n = 3)
            )
          )


# TODO see how many trts were xb by PRIVATE/PRIVATE CONSERVATION specifically
# (using comap LEGEND)

# TODO see how we can leverage the NAME attribute to further confirm xb nature
# of treatment pairs

# TODO refactor/tidy put into rmd
