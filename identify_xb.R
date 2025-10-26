# load packages ####
library(mapview)
library(sf)
library(tidyverse)
library(nngeo)
library(tictoc)

# load data ####
# treatments
treatments <- st_read("raw/misfits.gdb",
                      layer = "Forest_Tracker_v1_2024_ALL",
                      fid_column_name = "OBJECTID") %>%
  st_cast("MULTIPOLYGON") %>%
  st_make_valid
# conservation areas
comap_manager_detail <- st_read("raw/misfits.gdb",
                                layer = "COMaP_Manager_Dissolve") %>%
  mutate(id =  row_number()) %>%
  st_cast("MULTIPOLYGON") %>%
  st_make_valid
# make lines for COMaP "Manager Detail" attribute
comap_lines <- st_cast(comap_manager_detail, "MULTIPOLYGON") %>%
  st_cast("MULTILINESTRING")

# single poly xb ####
# buffer treatments by a negative buffer_distance
treatments_buffered <- st_buffer(treatments, dist = -15)
# identify which *buffered* treatments intersect comap boundaries
treatments_comap <- st_intersection(treatments_buffered,
                                    comap_lines)
# identify the same records among the original, *unbuffered* treatments
treatments <- treatments %>%
  mutate(
    single_poly_xb = if_else(OBJECTID %in% treatments_comap$OBJECTID, 1, 0)
  )

# multi-poly xb ####
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
# map centroid-finding function over each treatment. returns long list, each
# list element is a 1-feature sf object.
trts_centroids_list <- map(seq_len(nrow(treatments)),
                           \(x) st_centroid_within_poly(treatments[x, ]),
                           .progress = TRUE)

# map over 1000-element subsets of the output to rbind. Turns 1000 list elements
# (each element a 1-feature sf object) into 1 list element (a 1000-feature sf
# object). note: chunking is required because rbinding into a single +25,000 sf
# element all at once is not memory efficient.
d <- c(1:length(trts_centroids_list))
chunks <- split(d, ceiling(seq_along(d)/1000))
centroids_chunks <- map(chunks,
                        \(x) do.call(rbind, trts_centroids_list[x]),
                        .progress = TRUE)
# now, rbind the 26 ~1000-feature chunks to create one +25,000-feature sf
# object.
trt_centroids <- do.call(rbind, centroids_chunks)

# spatial join the treatment centroids to the comap regions. This is the key
# step in assigning each treatment to a "home" region.
trt_centroids_join <- trt_centroids %>%
  # join
  st_join(comap_manager_detail, st_intersects) %>%
  # keep only OBJECTID (of treatment) and id (of comap region). This is all the
  # info we need.
  select(OBJECTID, id) %>%
  # drop the spatial data.
  st_drop_geometry

# CoMAP id == 292 is "private". This mega-polygon overlaps some other comap
# polygons, creating a problem where treatment centroids are assigned to
# multiple comap types. Filter out these instances of treatments having two
# "home" polygons by removing the erroneous assignment to "private".
get_single_home <- \(oid) {
  homes <- trt_centroids_join %>% filter(OBJECTID == oid)
  if (nrow(homes) > 1) {
    homes %>% filter(id != 262)
  } else {
    homes
  }
}

# get list of unique oids to map function over. should have ONE home polygon for
# each unique oid.
unique_oids <- trt_centroids_join %>% pull(OBJECTID) %>% unique

# map function over each unique id and bind the results into a dataframe.
centroids_home <- map(unique_oids,
                      \(x) get_single_home(x),
                      .progress = TRUE) %>%
  do.call(rbind, .)

# join the "home" comap region id to the treatments data
treatments_join <- left_join(treatments,
                             centroids_home, by = "OBJECTID")

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

potential_xb <- treatments_join_out %>%
  filter(if_any(contains("_dist"), \(x)  x < 20))

# to identify WHICH year slices the target treatment has matches in:
same_year_xb <- treatments_join_out %>%
  filter(if_any(contains("same_yr_dist"), \(x)  x < 20))
plus1_year_xb <- treatments_join_out %>%
  filter(if_any(contains("plus1_dist"), \(x)  x < 20))
plus2_year_xb <- treatments_join_out %>%
  filter(if_any(contains("plus2_dist"), \(x)  x < 20))

treatments_xb <- treatments_join_out %>%
  mutate(
    # is the target xb by adjacency to another treatment in ANY year slice?
    multi_poly_xb = case_when(OBJECTID %in% potential_xb$OBJECTID ~ 1,
                              .default = 0),
    # is the target xb by adjacency to another treatment in the same year?
    xb_multi_0 = case_when(OBJECTID %in% same_year_xb$OBJECTID ~ 1,
                           .default = 0),
    # is the target xb by adjacency to another treatment in the +/-1 year slice?
    xb_multi_1 = case_when(OBJECTID %in% plus1_year_xb$OBJECTID ~ 1,
                           .default = 0),
    # is the target xb by adjacency to another treatment in the +/-2 year slice?
    xb_multi_2 = case_when(OBJECTID %in% plus2_year_xb$OBJECTID ~ 1,
                           .default = 0),
    # is the target xb on its own OR by adjacency with another treatment (any
    # year)?
    any_xb = case_when(single_poly_xb == 1 | multi_poly_xb == 1 ~ 1,
                       .default = 0),
    # if the target is xb by adjacency, what is the minimum year slice in which
    # adjacency occurs? (not xb by adjacency in any year == NA). values assigned
    # in order. If any target == TRUE for more than one condition, the value for
    # the first TRUE condition is assigned.
    xb_multi_any_min = case_when(OBJECTID %in% same_year_xb$OBJECTID ~ 0,
                                 OBJECTID %in% plus1_year_xb$OBJECTID ~ 1,
                                 OBJECTID %in% plus2_year_xb$OBJECTID ~ 2,
                                 .default = NA),
    # xb_multi_*_min: dummy variable for xb_multi_any_min. == 1 if that year
    # slice is the MINIMUM year slice in which adjacency occurs.
    xb_multi_0_min = if_else(xb_multi_any_min == 0, 1, 0, 0),
    xb_multi_1_min = if_else(xb_multi_any_min == 1, 1, 0, 0),
    xb_multi_2_min = if_else(xb_multi_any_min == 2, 1, 0, 0)
  )

# rename columns for clarity and brevity (.shp has 10-character limit for
# attribute names). just explicitly name everything for clarity.
treatments_xb %>% glimpse

out <- treatments_xb %>%
  rename(
    OID_OLD = "OBJECTID",
    PRJ_NAME = "PRJ_NAME",
    AGENCY = "AGENCY",
    AGENCY_C = "AGENCY_C",
    FUNDING = "FUNDING",
    LANDOWNER = "LANDOWNER",
    MGT_TYPE = "MGT_TYPE",
    RX_MGT = "RXFIRE_MGT",
    CAN_MGT = "CANOPY_MGT",
    SURF_MGT = "SURF_MGT",
    REFOREST = "REFOREST",
    TREE_COUNT = "TREE_COUNT",
    SPECIES = "SPECIES",
    PRJ_OBJECT = "PRJ_OBJECT",
    YEAR_COMP = "YEAR_COMP",
    ACRES_MGT = "ACRES_MGT",
    ACRES_GIS = "ACRES_GIS",
    NOTES = "NOTES",
    ORGFILE = "ORGFILE",
    UPDATED = "UPDATED",
    MODIFY_BY = "MODIFY_BY",
    FOR_TYPE = "FOR_TYPE",
    home_id = "id",
    sn1_oid = "same_yr_OBJECTID_n1",
    sn2_oid = "same_yr_OBJECTID_n2",
    sn3_oid = "same_yr_OBJECTID_n3",
    sn1_dist = "same_yr_dist_n1",
    sn2_dist = "same_yr_dist_n2",
    sn3_dist = "same_yr_dist_n3",
    p1n1_oid = "plus1_OBJECTID_n1",
    p1n2_oid = "plus1_OBJECTID_n2",
    p1n3_oid = "plus1_OBJECTID_n3",
    p1n1_dist = "plus1_dist_n1",
    p1n2_dist = "plus1_dist_n2",
    p1n3_dist = "plus1_dist_n3",
    p2n1_oid = "plus2_OBJECTID_n1",
    p2n2_oid = "plus2_OBJECTID_n2",
    p2n3_oid = "plus2_OBJECTID_n3",
    p2n1_dist = "plus2_dist_n1",
    p2n2_dist = "plus2_dist_n2",
    p2n3_dist = "plus2_dist_n3",
    xb_single = "single_poly_xb",
    xb_multi = "multi_poly_xb",
    xb_any = "any_xb",
    xb_m_0 = "xb_multi_0",
    xb_m_1 = "xb_multi_1",
    xb_m_2 = "xb_multi_2",
    xb_m_min = "xb_multi_any_min",
    xb_m_min_0 = "xb_multi_0_min",
    xb_m_min_1 = "xb_multi_1_min",
    xb_m_min_2 = "xb_multi_2_min",
    Shape = "Shape",
    Shape_Len = "Shape_Length",
    Shape_Area = "Shape_Area") %>%
  # note: no OID created here, because that is created in st_write below
  # need to coerce aly POLYGON to MULTIPOLYGON for write out
  st_cast

out
st_geometry_type(out)

# write out ####
st_write(obj = out,
         dsn = "raw/misfits.gdb",
         layer = "treatments_xb",
         append = FALSE)

write_csv(out %>% st_drop_geometry, "treatments_xb.csv")
