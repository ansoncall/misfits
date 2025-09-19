# identifying multi-polygon xb treatments

# load packages ####
library(mapview)
library(sf)
library(tidyverse)
library(nngeo)
library(tictoc)

# load data ####
# treatments
treatments <- st_read("raw/LarimerCoTreatments.gdb",
                          fid_column_name = "OBJECTID") %>%
  st_make_valid
# conservation areas
comap_manager_detail <- st_read("raw/misfits.gdb",
                                layer = "COMaP_LarimerCo_Manager")

# make lines for COMaP "Manager Detail" attribute
comap_lines <- st_cast(comap_manager_detail, "MULTIPOLYGON") %>%
  st_cast("MULTILINESTRING")

# wrangle ####
# id single polygon xb treatments
# step 1. buffer treatments by a negative buffer_distance
buffer_distance <- -15
treatments_buffered <- st_buffer(treatments, dist = buffer_distance)
# step 2. identify which *buffered* treatments intersect comap boundaries
treatments_comap <- st_intersection(treatments_buffered,
                                    comap_lines)
# id single poly xb treatments
treatments_comap_polys <- treatments_buffered %>%
  filter(OBJECTID %in% treatments_comap$OBJECTID)

# add column to flag single-polygon xb treatments
treatments <- treatments %>%
  mutate(
    single_poly_xb = if_else(OBJECTID %in% treatments_comap$OBJECTID, 1, 0)
  )

# function to identify polygon centroid or, if centroid is not in the polygon, a
# point that is definitely within the polygon.
st_centroid_within_poly <- \(poly) {
  # calculate polygon centroid
  ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
  # check if centroid is in polygon
  in_poly <- diag(st_within(ctrd, poly, sparse = F))
  # replace centroids that are not within polygon with st_point_on_surface().
  st_geometry(ctrd[!in_poly,]) <-
    st_geometry(st_point_on_surface(poly[!in_poly,]))
  # return the point
  ctrd
}

# map centroid-finding function over each treatment, and rbind the result into
# an sf dataframe.
trts_centroids <- map(seq_len(nrow(treatments)), \(x) {
  y <- treatments[x, ]
  st_centroid_within_poly(y)
  }) %>%
  do.call(rbind, .)

# add id col and cast as multipolygon to facilitate the spatial join that
# follows.
comap <- comap_manager_detail %>%
  mutate(id =  row_number()) %>%
  st_cast("MULTIPOLYGON")

# spatial join the treatment centroids to the comap regions. This is the key
# step in assigning each treatment to a "home" region.
trts_centroids_join <- trts_centroids %>%
  # join
  st_join(comap, st_intersects) %>%
  # keep only OBJECTID (of treatment) and id (of comap region). This is all the
  # info we need.
  select(OBJECTID, id) %>%
  # drop the spatial data.
  st_drop_geometry

# join the "home" comap region id to the treatments data
treatments_join <- left_join(treatments, trts_centroids_join, by = "OBJECTID")

# TODO fix edge case: one multipart polygon not assigned to comap region because
# part of it falls outside the county. See:
# mapview(comap, col.regions = "red") +
#   mapview(treatments_join[1607, ]) +
#   mapview(trts_centroids_join[1607, ])
treatments_join$id[1607] <- 57
# Since there is only one instance, a one-off hard-coded fix is okay. Moving to
# full-state analysis may require a programmatic fix.

# rather than init empty columns and fill with for loop, write a function that
# does what you want and map it over the df with mutate?

# function that takes a "target" treatment, an sf dataframe of "candidate"
# treatments (treatments within the right year slice), and a k_max parameter
# (the max number of neighbors to return). returns the OBJECTIDs of the nearest
# neighbors and the distances to them as tibbles.
find_neighbors <- \(target, candidates, k_max = 3) {
  # if there are no other treatments (of any distance) with the right completion
  # year, return tibble with NAs
  if (nrow(candidates) == 0) {
    return(tibble(OBJECTID = rep(NA, k_max),
                  dist = rep(NA, k_max)))
  }
  # set the number of neighbors to return. this handles cases when there are
  # less than k_max other treatments with the right completion year.
  k <- min(nrow(candidates), k_max)
  # first, hide pesky "lines or polygons" output with suppressMessages().
  # st_nn() prints its "mode" every time it is called, but we don't care to see
  # that printed 1000x.
  suppressMessages(
  # identify the neighbors and the distances to them.
    nn <- st_nn(target, candidates, k = k, returnDist = TRUE, progress = FALSE)
  )
  # return output as tibble.
  tibble(
    # st_nn has returned the index of the nearest neighbor(s). here, use the
    # index to subset the candidates dataframe in order to grab the correct
    # OBJECTID(s) in a list-col.
    OBJECTID = st_drop_geometry(candidates[nn$nn[[1]], "OBJECTID"])[[1]] %>%
      c(rep(NA, k_max - k)),  # pad with NA when there are <k_max neighbors
    # return the distance to those nearest neighbors as the "dist" list-col
    dist     = nn$dist[[1]] %>% c(rep(NA, k_max - k))
  )
}

# helper that applies find_neighbors to one treatment (row index i)
neighbors_for_i <- function(i, df, k_max = 3) {
  target_trt <- df[i, ]                       # keep geometry intact
  this_id <- target_trt$id[[1]]               # scalar id
  this_yr <- target_trt$YEAR_COMP[[1]]        # scalar year

  candidates_same  <- df %>% filter(id != this_id, YEAR_COMP == this_yr)
  candidates_plus1 <- df %>% filter(id != this_id,
                                    YEAR_COMP %in% c(this_yr - 1, this_yr + 1))
  candidates_plus2 <- df %>% filter(id != this_id,
                                    YEAR_COMP %in% c(this_yr - 2, this_yr + 2))

  list(
    same  = find_neighbors(target_trt, candidates_same,  k_max = k_max),
    plus1 = find_neighbors(target_trt, candidates_plus1, k_max = k_max),
    plus2 = find_neighbors(target_trt, candidates_plus2, k_max = k_max)
  )
}


# find neighbors for all treatments
res_list <- map(seq_len(nrow(treatments_join)),
                neighbors_for_i,
                df = treatments_join,
                k_max = 3,
                .progress = TRUE)

# save/set aside geometry
tj_geo <- st_geometry(treatments_join)
# attach results as list-columns
treatments_join_out <- treatments_join %>%
  mutate(
    same_yr = map(res_list, "same"),
    plus1   = map(res_list, "plus1"),
    plus2   = map(res_list, "plus2")
  ) %>%
  # tidy output: ungroup and expand list-cols to separate columns
  ungroup %>%
  # expand top-level list to "_OBJECTID" and "_dist" list_cols
  # note: unnest_wider is dropping the geometry and turning the sf object into a
  # tibble. be sure to reattach!
  unnest_wider(starts_with(c("same_yr", "plus")), names_sep = "_") %>%
  # expand "_OBJECTID" and "_dist" list cols to single values for each neighbor
  unnest_wider(starts_with(c("same_yr", "plus")), names_sep = "_n") %>%
  # reset geometry
  st_set_geometry(tj_geo)

# graph: all data
multi_poly_xb_all_dists <- treatments_join_out %>%
  pivot_longer(contains("_dist"), names_to = "type", values_to = "distance") %>%
  mutate(
    year_slice = case_when(
      str_detect(type, "same_yr_") ~ "same",
      str_detect(type, "plus1_") ~ "plus1",
      str_detect(type, "plus2_") ~ "plus2",
      .default = "ERROR"
    ),
    neighbor = case_when(
      str_detect(type, "_dist_n1") ~ "first",
      str_detect(type, "_dist_n2") ~ "second",
      str_detect(type, "_dist_n3") ~ "third",
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
multi_poly_xb_close <- treatments_join_out %>%
  pivot_longer(contains("_dist"), names_to = "type", values_to = "distance") %>%
  mutate(
    year_slice = case_when(
      str_detect(type, "same_yr_") ~ "same",
      str_detect(type, "plus1_") ~ "plus1",
      str_detect(type, "plus2_") ~ "plus2",
      .default = "ERROR"
    ),
    neighbor = case_when(
      str_detect(type, "_dist_n1") ~ "first",
      str_detect(type, "_dist_n2") ~ "second",
      str_detect(type, "_dist_n3") ~ "third",
      .default = "ERROR_neighbor"
    )
  ) %>%
  filter(distance < 1200) %>%
  ggplot(aes(x = distance, fill = neighbor)) +
  geom_histogram(alpha = 0.6, position = "dodge") +
  facet_wrap(~year_slice, nrow = 3) +
  scale_fill_manual(values = c("#162521", "#4F7CAC", "#C0E0DE")) +
  theme_minimal() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 14))

# TODO map out some treatments where the neighbor is really close
potential_xb <- treatments_join_out %>%
  filter(if_any(contains("_dist"), \(x)  x < 20))

# nolint start
# mapview(comap_lines, color = "red") +
#   mapview(potential_xb,
#           zcol = "OBJECTID",
#           col.regions =
#             c(
#             rep(hcl.colors(palette = "Dynamic", n = 5), 32),
#             hcl.colors(palette = "Dynamic", n = 3)
#             )
#           )
# nolint end

treatments_xb <- treatments_join_out %>%
  mutate(multi_poly_xb = case_when(
    OBJECTID %in% potential_xb$OBJECTID ~ 1,
    .default = 0
  ),
  any_xb = case_when(
    single_poly_xb == 1 | multi_poly_xb == 1 ~ 1,
    .default = 0
  ))

# calculate proportion xb area in Larimer Co region
xb_area <- st_area(treatments_xb %>% filter(any_xb == 1)) %>% sum
non_xb_area <- st_area(treatments_xb %>% filter(any_xb != 1)) %>% sum

xb_proportion <- xb_area / (xb_area + non_xb_area)

# using the name attribute ####

extract_neighbors <- \ (data, target_OID, pattern) {
  data %>%
    filter(OBJECTID == target_OID) %>%
    select(contains(pattern)) %>%
    st_drop_geometry() %>%
    pivot_longer(everything(), values_to = "value") %>%
    pull(value)
}

neighbor_name_dist <- \ (i) {
  # i = 125
  target_name <- potential_xb[i, ] %>% pull(PRJ_NAME)
  target_OID <- potential_xb[i, ] %>% pull(OBJECTID)

  neighbor_ids <- extract_neighbors(treatments_xb, target_OID, "_OBJECTID_n")
  neighbor_dists <- extract_neighbors(treatments_xb, target_OID, "_dist_n")

  names_df = treatments_xb %>%
    filter(OBJECTID %in% neighbor_ids) %>%
    st_drop_geometry %>%
    select(OBJECTID, PRJ_NAME)

  nears <- cbind(neighbor_ids, neighbor_dists) %>%
    as_tibble %>%
    set_names(c("OBJECTID", "dists")) %>%
    filter(as.numeric(dists) < 20) %>%
    left_join(names_df, by = "OBJECTID") %>%
    mutate(lev_dist = adist(PRJ_NAME,
                            target_name,
                            costs = list("insertions" = 1,
                                         "deletions" = 1)))

  list(target_OID = target_OID,
       mean_lev_dist = mean(nears$lev_dist),
       name = target_OID,
       details = nears)

}

names_analysis <- map(seq_len(nrow(potential_xb)),
                      neighbor_name_dist,
                      .progress = TRUE)

# map over analysis results and extract mean lev_dist for each potential_xb
# treatment
mean_lev_xbs <- map(names_analysis, \(trt) {trt$mean_lev_dist}) %>%
  unlist %>%
  tibble %>%
  set_names("lev_dist") %>%
  mutate(source = "near xb neighbors")

# compare distribution of lev distances between near neighbors and between any
# two random treatments.

random_names_1 <- treatments_xb$PRJ_NAME %>% sample(1000)
random_names_2 <- treatments_xb$PRJ_NAME %>% sample(1000)
levs <- map2(random_names_1, random_names_2,
                    \(x,y) adist(x, y, costs = list("insertions" = 1,
                                                    "deletions" = 1))) %>%
  unlist %>%
  tibble %>%
  set_names("lev_dist") %>%
  mutate(source = "random") %>%
  rbind(mean_lev_xbs)

lev_density <- ggplot(levs, aes(x = lev_dist, color = source)) +
  geom_density() +
  theme_minimal(base_size = 16) +
  labs(x = "Levenshtein distance", y = "Density") +
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.8, 0.8))

# summary: if names were a good indicator of cross-boundary planning, we might
# expect nearby xb treatments to have similar names. however, the similarity of
# names is not really any different between near neighbors and random pairs.
# Also, there's a big density spike around ~38 distance, which reflects the
# common occurrence of the name "None" in Colorado State Forest Service records.
# I think that these "None"-named projects are often XB in nature, but they just
# don't seem to use that attribute in their system, at least not regularly.

# yet, there are some instances where lev distance is low among neighbors, and
# this might in fact confirm that they are xb. just know that this low lev
# distance is not exceptionally low:

# identify which in names_analysis have common names in neighbors (lev < 2)
low_lev <- names_analysis %>%
  # count number of neighbors with lev_dist < 2
  map(\(x) {
    df <- x$details %>% filter(lev_dist < 2)
    x$details <- df
    x$num_low_lev <- nrow(df)
    x
  }) %>%
  # filter out target treatments whose xb neighbors do not have matching names
  keep(\(x) {
    x$num_low_lev > 0
  })

# check the names of the remaining treatments
low_lev_names <-  treatments_xb %>%
  # OBJECTID should match one of the $name values from low_lev
  filter(OBJECTID %in% map(low_lev,
                           \(x) x$name[[1]]) %>%
           unlist) %>%
  select(OBJECTID, PRJ_NAME) %>%
  st_drop_geometry

# we can see that, of the 18 xb trts with name matches, there are really only
# two groups here. The rest of the matches are due to missing data.

# bottom line: leveraging the name attribute seems like a good idea, but the
# data just isn't there.


# TODO refactor below
# # private/private conservation ####
# # conservation areas
# comap_legend_detail <- st_read("raw/misfits.gdb",
#                                layer = "COMaP_LarimerCo_Legend")
#
# # make lines for COMaP "Legend Detail" attribute
# comap_legend_lines <- st_cast(comap_legend_detail, "MULTIPOLYGON") %>%
#   st_cast("MULTILINESTRING")
#
# # assign id col to comap polys
# comap_l <- comap_legend_detail %>%
#   mutate(id =  row_number()) %>%
#   st_cast("MULTIPOLYGON")
#
# comap_p_pc <- rbind(pc_lands,  p_lands %>% select(legend:id))
#
# # starting over with treatments, id single-poly private-private conservation
# # xb treatments
# treatments <- treatments
#
# # id centroids
# trts_centroids <- map(seq_len(nrow(treatments)), \(x) {
#   y <- treatments[x, ]
#   st_centroid_within_poly(y)
# }) %>%
#   do.call(rbind, .)
#
# trts_centroids_join <- st_join(trts_centroids, comap_p_pc, st_intersects)
#
# # build private set
# treatments_join_private <- left_join(treatments,
#                                      trts_centroids_join %>% select(OBJECTID, id) %>% st_drop_geometry,
#                                      by = "OBJECTID") %>%
#   filter(id == 5)
#
# treatments_join_private_b <- st_buffer(treatments_join_private, -15)
#
# # build private conservation set
# treatments_join_private_cons <- left_join(treatments,
#                                      trts_centroids_join %>% select(OBJECTID, id) %>% st_drop_geometry,
#                                      by = "OBJECTID") %>%
#   filter(id == 6)
#
#
# treatments_join_private_cons_b <- st_buffer(treatments_join_private_cons, -15)
#
# mapview(treatments_join_private, col.regions = "blue") +
# mapview(treatments_join_private_b, col.regions = "red")
#
# # id private lands treatments that overlap private conservation
# pc_lands <- comap_l %>% filter(id == 6)
# pc_landlines <- pc_lands %>% st_cast("MULTILINESTRING")
# private_to_private_cons <- st_intersection(treatments_join_private_b, pc_lands)
#
# # and vice versa
# p_lands <- comap_l %>% filter(id == 5) %>% st_make_valid %>% st_difference(pc_lands)
# p_landlines <- p_lands %>% st_cast("MULTILINESTRING")
# private_cons_to_private <- st_intersection(treatments_join_private_cons_b, p_lands)
#
#
# output <- treatments %>% filter(OBJECTID %in% c(private_cons_to_private$OBJECTID, private_to_private_cons$OBJECTID))
#
# test <- private_cons_to_private %>% filter(OBJECTID == 629) %>% st_buffer(-15)
#
# mapview(p_lands, col.regions = "yellow") +
#   mapview(pc_lands, col.regions = "blue") +
#   mapview(output, col.regions = "red") +
#   mapview(test, col.regions = "hotpink")
#
#
#
#
#
# buffer_distance <- -15
# # step 1. buffer treatments by a negative buffer_distance
# treatments_buffered <- st_buffer(treatments, dist = buffer_distance)
# # step 2. identify which *buffered* treatments intersect comap boundaries
# treatments_comap <- st_join(treatments_buffered, comap_lines, st_intersection)
# # id single poly xb treatments
# treatments_comap_polys <- treatments_buffered %>%
#   filter(OBJECTID %in% treatments_comap$OBJECTID)
#
# mapview(treatments_comap_polys) + mapview(comap_legend_lines)
#
# treatments <- treatments %>%
#   mutate(single_poly_xb = case_when(
#     OBJECTID %in% treatments_comap$OBJECTID ~ 1,
#     .default = 0
#   ))
#
#
# trts_centroids_legend_join <- st_join(trts_centroids, comap_l, st_intersects)
#
# treatments_legend_join <- left_join(treatments,
#                              trts_centroids_legend_join %>% select(OBJECTID, id) %>% st_drop_geometry,
#                              by = "OBJECTID")
#
# # filter to private/private conservation only
# treatments_private <- treatments_legend_join %>%
#   filter(id %in% c(5, 6))
#
# # TODO fix edge case: one multipart polygon not assigned to comap region because
# # part of it falls outside the county
# mapview(comap_l, col.regions = "red") +
#   mapview(treatments_legend_join[1607, ]) +
#   mapview(trts_centroids_legend_join[1607, ])
# treatments_legend_join$id[1607] <- 8
#
#
# # id nearest neighbors
#
# treatments_private$closest_same_yr_OBJECTID_1 <- NULL
# treatments_private$closest_same_yr_dist_1 <- NULL
# treatments_private$closest_same_yr_OBJECTID_2 <- NULL
# treatments_private$closest_same_yr_dist_2 <- NULL
# treatments_private$closest_same_yr_OBJECTID_3 <- NULL
# treatments_private$closest_same_yr_dist_3 <- NULL
# treatments_private$closest_plus1_OBJECTID_1 <- NULL
# treatments_private$closest_plus1_dist_1 <- NULL
# treatments_private$closest_plus1_OBJECTID_2 <- NULL
# treatments_private$closest_plus1_dist_2 <- NULL
# treatments_private$closest_plus1_OBJECTID_3 <- NULL
# treatments_private$closest_plus1_dist_3 <- NULL
# treatments_private$closest_plus2_OBJECTID_1 <- NULL
# treatments_private$closest_plus2_dist_1 <- NULL
# treatments_private$closest_plus2_OBJECTID_2 <- NULL
# treatments_private$closest_plus2_dist_2 <- NULL
# treatments_private$closest_plus2_OBJECTID_3 <- NULL
# treatments_private$closest_plus2_dist_3 <- NULL
#
#
# tic()
# for (i in 1:nrow(treatments_private)) {
#   # i <- 1
#   home_comap_poly <- treatments_private$id[i]
#   target_trt <- treatments_private[i, ]
#   target_yr <- target_trt$YEAR_COMP
#
#   filtered_same_yr <- treatments_private %>% filter(id != home_comap_poly,
#                                                  YEAR_COMP == target_yr)
#   filtered_plus1_yr <- treatments_private %>% filter(id != home_comap_poly,
#                                                   YEAR_COMP == target_yr + 1 |
#                                                     YEAR_COMP == target_yr - 1)
#   filtered_plus2_yr <- treatments_private %>% filter(id != home_comap_poly,
#                                                   YEAR_COMP == target_yr + 2 |
#                                                     YEAR_COMP == target_yr - 2)
#
#   matches_yr <- nrow(filtered_same_yr)
#   matches_plus1 <- nrow(filtered_plus1_yr)
#   matches_plus2 <- nrow(filtered_plus2_yr)
#
#   if (matches_yr == 0) {
#     treatments_private$closest_same_yr_OBJECTID_1[i] <- NA
#     treatments_private$closest_same_yr_dist_1[i] <- NA
#     treatments_private$closest_same_yr_OBJECTID_2[i] <- NA
#     treatments_private$closest_same_yr_dist_2[i] <- NA
#     treatments_private$closest_same_yr_OBJECTID_3[i] <- NA
#     treatments_private$closest_same_yr_dist_3[i] <- NA
#
#   } else if (matches_yr > 0 && matches_yr < 3) {
#     closest_n_same_yr <- st_nn(target_trt,
#                                filtered_same_yr,
#                                k = matches_yr,
#                                returnDist = TRUE)
#     # hacky way to return OIDs instead of row idx
#     # TODO clean up
#     closest_n_same_yr$nn[[1]] <- st_drop_geometry(filtered_same_yr[closest_n_same_yr$nn[[1]], "OBJECTID"])[[1]]
#   } else {
#     closest_n_same_yr <- st_nn(target_trt,
#                                filtered_same_yr,
#                                k = 3,
#                                returnDist = TRUE)
#     closest_n_same_yr$nn[[1]] <- st_drop_geometry(filtered_same_yr[closest_n_same_yr$nn[[1]], "OBJECTID"])[[1]]
#   }
#
#   if (matches_plus1 == 0) {
#     treatments_private$closest_plus1_OBJECTID_1[i] <- NA
#     treatments_private$closest_plus1_dist_1[i] <- NA
#     treatments_private$closest_plus1_OBJECTID_2[i] <- NA
#     treatments_private$closest_plus1_dist_2[i] <- NA
#     treatments_private$closest_plus1_OBJECTID_3[i] <- NA
#     treatments_private$closest_plus1_dist_3[i] <- NA
#   } else if (matches_plus1 > 0 && matches_plus1 < 3) {
#     closest_n_plus1 <- st_nn(target_trt,
#                              filtered_plus1_yr,
#                              k = matches_plus1,
#                              returnDist = TRUE)
#     closest_n_plus1$nn[[1]] <- st_drop_geometry(filtered_plus1_yr[closest_n_plus1$nn[[1]], "OBJECTID"])[[1]]
#   } else {
#     closest_n_plus1 <- st_nn(target_trt,
#                              filtered_plus1_yr,
#                              k = 3,
#                              returnDist = TRUE)
#     closest_n_plus1$nn[[1]] <- st_drop_geometry(filtered_plus1_yr[closest_n_plus1$nn[[1]], "OBJECTID"])[[1]]
#   }
#
#   if (matches_plus2 == 0) {
#     treatments_private$closest_plus2_OBJECTID_1[i] <- NA
#     treatments_private$closest_plus2_dist_1[i] <- NA
#     treatments_private$closest_plus2_OBJECTID_2[i] <- NA
#     treatments_private$closest_plus2_dist_2[i] <- NA
#     treatments_private$closest_plus2_OBJECTID_3[i] <- NA
#     treatments_private$closest_plus2_dist_3[i] <- NA
#
#   } else if (matches_plus2 > 0 && matches_plus2 < 3){
#     closest_n_plus2 <- st_nn(target_trt,
#                              filtered_plus2_yr,
#                              k = matches_plus2,
#                              returnDist = TRUE)
#     closest_n_plus2$nn[[1]] <- st_drop_geometry(filtered_plus2_yr[closest_n_plus2$nn[[1]], "OBJECTID"])[[1]]
#   } else {
#     closest_n_plus2 <- st_nn(target_trt,
#                              filtered_plus2_yr,
#                              k = 3,
#                              returnDist = TRUE)
#     closest_n_plus2$nn[[1]] <- st_drop_geometry(filtered_plus2_yr[closest_n_plus2$nn[[1]], "OBJECTID"])[[1]]
#   }
#
#   treatments_private$closest_same_yr_OBJECTID_1[i] <- NA
#   try(treatments_private$closest_same_yr_OBJECTID_1[i] <- closest_n_same_yr$nn[[1]][[1]])
#   treatments_private$closest_same_yr_dist_1[i] <- NA
#   try(treatments_private$closest_same_yr_dist_1[i] <- closest_n_same_yr$dist[[1]][[1]])
#   treatments_private$closest_same_yr_OBJECTID_2[i] <- NA
#   try(treatments_private$closest_same_yr_OBJECTID_2[i] <- closest_n_same_yr$nn[[1]][[2]])
#   treatments_private$closest_same_yr_dist_2[i] <- NA
#   try(treatments_private$closest_same_yr_dist_2[i] <- closest_n_same_yr$dist[[1]][[2]])
#   treatments_private$closest_same_yr_OBJECTID_3[i] <- NA
#   try(treatments_private$closest_same_yr_OBJECTID_3[i] <- closest_n_same_yr$nn[[1]][[3]])
#   treatments_private$closest_same_yr_dist_3[i] <- NA
#   try(treatments_private$closest_same_yr_dist_3[i] <- closest_n_same_yr$dist[[1]][[3]])
#
#   treatments_private$closest_plus1_OBJECTID_1[i] <- NA
#   try(treatments_private$closest_plus1_OBJECTID_1[i] <- closest_n_plus1$nn[[1]][[1]])
#   treatments_private$closest_plus1_dist_1[i] <- NA
#   try(treatments_private$closest_plus1_dist_1[i] <- closest_n_plus1$dist[[1]][[1]])
#   treatments_private$closest_plus1_OBJECTID_2[i] <- NA
#   try(treatments_private$closest_plus1_OBJECTID_2[i] <- closest_n_plus1$nn[[1]][[2]])
#   treatments_private$closest_plus1_dist_2[i] <- NA
#   try(treatments_private$closest_plus1_dist_2[i] <- closest_n_plus1$dist[[1]][[2]])
#   treatments_private$closest_plus1_OBJECTID_3[i] <- NA
#   try(treatments_private$closest_plus1_OBJECTID_3[i] <- closest_n_plus1$nn[[1]][[3]])
#   treatments_private$closest_plus1_dist_3[i] <- NA
#   try(treatments_private$closest_plus1_dist_3[i] <- closest_n_plus1$dist[[1]][[3]])
#
#   # TODO fix subscript out of bounds when matches > 0 && matches < 3
#   treatments_private$closest_plus2_OBJECTID_1[i] <- closest_n_plus2$nn[[1]][[1]]
#   treatments_private$closest_plus2_dist_1[i] <- closest_n_plus2$dist[[1]][[1]]
#
#   treatments_private$closest_plus2_OBJECTID_2[i] <- NA
#   try(treatments_private$closest_plus2_OBJECTID_2[i] <- closest_n_plus2$nn[[1]][[2]])
#   treatments_private$closest_plus2_OBJECTID_2[i] <- NA
#   try(treatments_private$closest_plus2_OBJECTID_2[i] <- closest_n_plus2$nn[[1]][[2]])
#   treatments_private$closest_plus2_dist_2[i] <- NA
#   try(treatments_private$closest_plus2_dist_2[i] <- closest_n_plus2$dist[[1]][[2]])
#   treatments_private$closest_plus2_OBJECTID_3[i] <- NA
#   try(treatments_private$closest_plus2_OBJECTID_3[i] <- closest_n_plus2$nn[[1]][[3]])
#   treatments_private$closest_plus2_dist_3[i] <- NA
#   try(treatments_private$closest_plus2_dist_3[i] <- closest_n_plus2$dist[[1]][[3]])
#
#   print(i)
# }
# toc()
#
# # plots
#
# # graph: all data
# treatments_private %>%
#   pivot_longer(contains("_dist_"), names_to = "type", values_to = "distance") %>%
#   mutate(
#     year_slice = case_when(
#       str_detect(type, "_same_yr_") ~ "same",
#       str_detect(type, "_plus1_") ~ "plus1",
#       str_detect(type, "_plus2_") ~ "plus2",
#       .default = "ERROR"
#     ),
#     neighbor = case_when(
#       str_detect(type, "_dist_1") ~ "first",
#       str_detect(type, "_dist_2") ~ "second",
#       str_detect(type, "_dist_3") ~ "third",
#       .default = "ERROR_neighbor"
#     )
#   ) %>%
#   ggplot(aes(x = distance, fill = neighbor)) +
#   geom_histogram(alpha = 0.6, position = "dodge") +
#   facet_wrap(~year_slice, nrow = 3) +
#   scale_fill_manual(values = c("#162521", "#4F7CAC", "#C0E0DE")) +
#   theme_minimal() +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 14))
#
# # graph: less than 2.5 kms
# treatments_private %>%
#   pivot_longer(contains("_dist_"), names_to = "type", values_to = "distance") %>%
#   mutate(
#     year_slice = case_when(
#       str_detect(type, "_same_yr_") ~ "same",
#       str_detect(type, "_plus1_") ~ "plus1",
#       str_detect(type, "_plus2_") ~ "plus2",
#       .default = "ERROR"
#     ),
#     neighbor = case_when(
#       str_detect(type, "_dist_1") ~ "first",
#       str_detect(type, "_dist_2") ~ "second",
#       str_detect(type, "_dist_3") ~ "third",
#       .default = "ERROR_neighbor"
#     )
#   ) %>%
#   filter(distance < 1200) %>%
#   ggplot(aes(x = distance, fill = neighbor)) +
#   geom_histogram(alpha = 0.6, position = "dodge") +
#   facet_wrap(~year_slice, nrow = 3) +
#   scale_fill_manual(values = c("#162521", "#4F7CAC", "#C0E0DE")) +
#   theme_minimal() +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 14))
#
# # map out some treatments where the neighbor is really close
# potential_private_xb <- treatments_private %>%
#   filter(closest_same_yr_dist_1 < 20 |
#            closest_plus1_dist_1 < 20 |
#            closest_plus2_dist_1 < 20)
#
# mapview(comap_legend_lines, color = "red") +
#   mapview(potential_private_xb,
#           zcol = "OBJECTID",
#           col.regions =
#             c(
#               rep(hcl.colors(palette = "Dynamic", n = 5), 32),
#               hcl.colors(palette = "Dynamic", n = 3)
#             )
#   )
#
#
