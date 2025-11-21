# additional fireshed statistics

library(tidyverse)

df <- read_csv("raw/t_f_intersect.csv")
names(df)

# treatment*fireshed stats ####
# number of treatments/treatment area in each fireshed, by treatment type
t_f_summary <- df %>%
  group_by(Fireshed_I, xb_any) %>%
  summarize(n = n(),
            area = sum(SHAPE_Area)) # in map units (sq. meters)
# this number is slightly inflated because some treatmens were cut in half by
# the intersect:
sum(t_f_summary$n)

# number of treatments/treatment area in each fireshed, all treatment types
f_summary <- df %>%
  group_by(Fireshed_I) %>%
  summarize(n = n(),
            area = sum(SHAPE_Area)) # in map units (sq. meters)

# manager*fireshed stats ####
df_c <- read_csv("raw/comap_firesheds_intersect.csv")
names(df_c)

f_summary2 <- df_c %>%
  group_by(Fireshed_I) %>%
  summarize(n_distinct_managers = n_distinct(MANAGER_DETAIL))

c_f_summary <- df_c %>%
  group_by(Fireshed_I, MANAGER_DETAIL) %>%
  summarize(n_distinct_managers = n_distinct(MANAGER_DETAIL),
            area_managed = sum(Shape_Area)) # in map units
