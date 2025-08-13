# wrangle data for analysis

# load packages ####
library(tidyverse)
library(mapview)
library(sf)
library(magrittr)
library(spatialEco)
# load data ####
# fuel treatments
treatments_raw <- st_read("raw/LarimerCoTreatments.gdb", fid_column_name = "OBJECTID")

treatments_raw %<>% st_make_valid
# conservation areas
comap_raw <- st_read("raw/misfits.gdb", layer = "COMaP_Raw")
comap_legend <- st_read("raw/misfits.gdb", layer = "COMaP_LarimerCo_Legend")
comap_manager_detail <- st_read("raw/misfits.gdb", layer = "COMaP_LarimerCo_Manager")

# defining xb stuff ####
# TODO move this to the defining_xb.Rmd file
# TODO run for entire state of CO
# do -10m buffer on treatments
treatments_neg10 <- st_buffer(treatments_raw, dist = -10)

# convert comap to lines
comap_legend_lines <- st_cast(comap_legend, "MULTIPOLYGON") %>%
  st_cast("MULTILINESTRING")

# buffered treatments that intersect a comap boundary
treatments_comap <- st_intersection(treatments_neg10,
                                    comap_legend_lines)
treatments_comap_polys <- treatments_neg10 %>%
  filter(OBJECTID %in% treatments_comap$OBJECTID)

names(treatments_raw)

mapview(treatments_comap_polys, col.regions = "red") +
  mapview(comap_legend_lines)

count_xb <- function (treatments, comap_lines, buffer_distance) {
  # step 1. buffer treatments by buffer_distance
  treatments_buffered <- st_buffer(treatments_raw, dist = buffer_distance)
  # step 2. identify which buffered treatments intersect comap boundaries
  treatments_comap <- st_intersection(treatments_buffered,
                                      comap_lines)
  # step 3. identify which polygons in the original treatments intersected comap
  # boundaries
  treatments_comap_polys <- treatments_buffered %>%
    filter(OBJECTID %in% treatments_comap$OBJECTID)
  # step 4. count the number of treatments that intersect comap boundaries
  n_xb <- nrow(treatments_comap_polys)
  return(n_xb)
}

buffer_distances <- seq(from = -100, 0, by = 10)

xb_counts <- map_dbl(buffer_distances,
                      \(x) count_xb(treatments_raw, comap_legend_lines, x),
                     .progress = TRUE)


dist_df <- tibble(buffer_distance = buffer_distances,
                     xb_count = xb_counts)
# plot the results
ggplot(dist_df, aes(x = buffer_distance, y = xb_count)) +
  geom_line() +
  geom_point() +
  labs(title = "Cross-boundary treatments by buffer distance",
       x = "Buffer distance (m)",
       y = "Number of cross-boundary treatments, Larimer Co.") +
  theme_minimal() +
  geom_vline(xintercept = -15, linetype = "dashed", color = "red")

# examine effect of COMaP divisions ####
# make lines for COMaP "Manager Detail"
comap_manager_lines <- st_cast(comap_manager_detail, "MULTIPOLYGON") %>%
  st_cast("MULTILINESTRING")

legend_manager_comp <- tibble(buffer_distance = seq(from = -25,
                                                    to = -10,
                                                    by = 5))

legend_manager_comp$xb_count_legend <- map_dbl(
  legend_manager_comp$buffer_distance,
  \(x) count_xb(treatments_raw, comap_legend_lines, x),
  .progress = TRUE
  )

legend_manager_comp$xb_count_manager <- map_dbl(
  legend_manager_comp$buffer_distance,
  \(x) count_xb(treatments_raw, comap_manager_lines, x),
  .progress = TRUE
)

legend_manager_comp
# plot the results
ggplot(legend_manager_comp, aes(x = buffer_distance)) +
  geom_line(aes(y = xb_count_legend, color = "Legend")) +
  geom_line(aes(y = xb_count_manager, color = "Manager Detail")) +
  geom_point(aes(y = xb_count_legend, color = "Legend")) +
  geom_point(aes(y = xb_count_manager, color = "Manager Detail")) +
  labs(title = "Cross-boundary treatments by buffer distance",
       x = "Buffer distance (m)",
       y = "Number of cross-boundary treatments, Larimer Co.",
       color = "COMaP Layer") +
  theme_minimal() +
  geom_vline(xintercept = -15, linetype = "dashed", color = "red")



# check for comap polygons where the manager detail spans more than one legend
# category
comap_raw %>%
  as_tibble() %>%
  group_by(MANAGER_DETAIL) %>%
  summarize(
    unique_legend = n_distinct(legend)
  ) %>%
  arrange(desc(unique_legend))




# TODO resolve why I'm not seeing more PRIVATE land parcels on this map. This
# result is not trustworthy until we figure this out!
comap_private <- comap_raw %>%
  filter(MANAGER_DETAIL == "Private")
plot(comap_raw)
mapview(comap_raw)
mapview(comap_legend_lines, col.regions = "blue") +
  mapview(comap_manager_lines, col.regions = "green") +
  mapview(treatments_comap_polys, col.regions = "red") +
  mapview(comap_private, col.regions = "purple")


xb_tab <- function (treatments, comap_lines, buffer_distance) {
  # step 1. buffer treatments by buffer_distance
  treatments_buffered <- st_buffer(treatments_raw, dist = buffer_distance)
  # step 2. identify which buffered treatments intersect comap boundaries
  treatments_comap <- st_intersection(treatments_buffered,
                                      comap_lines)
  # step 3. identify which polygons in the original treatments intersected comap
  # boundaries
  treatments_comap_polys <- treatments_buffered %>%
    filter(OBJECTID %in% treatments_comap$OBJECTID)
  # step 4. return polys
  return(treatments_comap_polys)
}

legend_buffer_out <- xb_tab(treatments_raw, comap_legend_lines, -15)
nrow(legend_buffer_out)

manager_buffer_out <- xb_tab(treatments_raw, comap_manager_lines, -15)
nrow(manager_buffer_out)

diff <- legend_buffer_out %>%
  filter(!OBJECTID %in% manager_buffer_out$OBJECTID)

diff2 <- manager_buffer_out %>%
  filter(!OBJECTID %in% legend_buffer_out$OBJECTID)
diff2



nrow(legend_buffer_out) - nrow(manager_buffer_out)

mapview(diff) + mapview(diff2, col.regions = "purple") +
  mapview(comap_legend_lines, col.regions = "blue") +
  mapview(comap_manager_lines, col.regions = "green")

# one trt, two polys ####
# TODO:
# next steps:

# treatments can be recorded as separate polygons on either side of a comap
# boundary, but functionally "one" cross-boundary treatment.

# check to see if the number of xb treatments changes when such treatments are
# "consolidated".

# In other words: Two treatments are on separately owned lands AND are adjacent
# (say within 10 meters of each other) AND have completion dates within N years
# (maybe 2 years?) of each other.

# step 1. filter out treatments that aren't already cross-boundary
treatments_notxb <- treatments_raw %>%
  # set appropriate buffer!
  filter(st_intersects(treatments_neg10, comap_legend_lines, sparse = FALSE))


# step 2. spatial join to add COMaP detail to treatments
treatments_notxb_detail <- st_join(treatments_notxb,
                                   comap_legend,
                                   join = st_intersects)

# step 3. spatial join treatments with their near neighbors
treatments_near <- st_join(treatments_raw,
                           treatments_raw,
                           join = st_is_within_distance,
                           dist = 10)

# step 4. filter to only those joins that are cross-boundary
treatments_near_xb <- treatments_near %>%
  filter(OBJECTID.x != OBJECTID.y) %>%
  filter(st_intersects(treatments_neg10[OBJECTID.x, ],
                       comap_legend_lines, sparse = FALSE)) %>%
  select(OBJECTID.x, OBJECTID.y)
