# treatment area by fireshed
# load packages ####
library(mapview)
library(sf)
library(tidyverse)
library(nngeo)
library(tictoc)

# load data ####
# treatments
all_treatments <- st_read("raw/misfits.gdb",
                          layer = "treatments_xb",
                          fid_column_name = "OBJECTID")
# firesheds
firesheds <- st_read("raw/firesheds/Firesheds_Clip_USFS.shp") %>%
  # mutate(id =  row_number()) %>%
  st_cast("MULTIPOLYGON") %>%
  st_make_valid %>%
  st_transform(st_crs(all_treatments))

# slim down dataframes for efficiency
t_slim <- all_treatments %>%
  select(OBJECTID, starts_with("xb_")) %>%
  # make xb_no column (not xb in any way)
  mutate(xb_no = if_else(xb_any == 1, 0, 1))
f_slim <- firesheds %>%
  select(Fireshed_I, Fireshed_N)

# calculate treatment area by fireshed
# 1) total treatment area
treatments_joined <- st_intersection(
  t_slim,
  f_slim  # assumes first column uniquely identifies firesheds
) %>%
  # area in hectares
  mutate(area = st_area(.)/1e4) %>%
  st_drop_geometry()

# summarize total treated area and areas for each type of xb treatment (as a
# ratio *:total area)
treatment_summary <- treatments_joined %>%
  group_by(Fireshed_I, Fireshed_N) %>%
  summarise(total_area = sum(as.numeric(area), na.rm = TRUE),
            .groups = "drop") %>%
  left_join(
    treatments_joined %>%
      pivot_longer(cols = starts_with("xb_"),
                   names_to = "treatment_type",
                   values_to = "is_treated") %>%
      filter(is_treated == 1) %>%
      group_by(Fireshed_I, Fireshed_N, treatment_type) %>%
      summarise(treat_area = sum(as.numeric(area), na.rm = TRUE),
                .groups = "drop"),
    by = c("Fireshed_I", "Fireshed_N")
  ) %>%
  mutate(ratio = treat_area / total_area)

# wide output
treatment_summary_wide <- treatment_summary %>%
  select(Fireshed_I, Fireshed_N, treatment_type, ratio, total_area) %>%
  pivot_wider(names_from = treatment_type,
              values_from = ratio,
              values_fill = 0) %>%
  arrange(total_area)

# show
glimpse(treatment_summary_wide)
treatment_summary_wide

# check
treatments_joined %>% filter(Fireshed_I == 2152) %>% summarize(area= sum(area))
treatments_joined %>% filter(Fireshed_I == 2152)
treatments_joined %>% filter(Fireshed_I == 2152 & xb_no == 1) %>%
  summarize(area= sum(area))

treatment_summary %>% filter(Fireshed_I == 2257)
treatments_joined %>% filter(Fireshed_I == 2257) %>% summarize(area= sum(area))
treatments_joined %>% filter(Fireshed_I == 2257)
treatments_joined %>% filter(Fireshed_I == 2257 & xb_no == 1) %>%
  summarize(area= sum(area))
treatments_joined %>% filter(Fireshed_I == 2257 & xb_single == 1) %>%
  summarize(area= sum(area))

# looks good, but need to drop one meaningless col and rearrange
treatment_summary_wide <- treatment_summary_wide %>%
  select(Fireshed_I, Fireshed_N, total_area, xb_no, xb_any, xb_single, xb_multi,
         xb_m_0,
         xb_m_1,
         xb_m_2,
         xb_m_min_0,
         xb_m_min_1,
         xb_m_min_2,
         )
write_csv(treatment_summary_wide, "treatments_by_fireshed.csv")

