# Nuclei counts from optimised protocol

setwd("/Users/sd21/Documents/WORK/02_PROJECTS/single_cell-nuclei_transcriptomics/")

library(tidyverse)
options(scipen = 999)

counts <- tribble(
  ~species,         ~stage,         ~nuclei,
  "H. contortus",   "Adult female", c(30200, 31820, 33370),
  "H. contortus",   "Adult male",   c(60450, 56250, 97500),
  "H. contortus",   "L3",           c(442500, 642000, 566000),
  "T. muris",       "Adult female", c(622650, 740000, 638000),
  "T. muris",       "Adult male",   c(731120, 935935, 506570)
) |>
  unnest(nuclei) |>
  mutate(stage = fct_inorder(stage), species = fct_inorder(species))

p <- ggplot(counts, aes(stage, nuclei)) +
  stat_summary(fun = mean, geom = "col",
               fill = "skyblue", alpha = 0.7, colour = "darkgrey",
               linewidth = 0.3, width = 0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", colour = "orange",
               width = 0.25, linewidth = 0.4) +
  geom_jitter(width = 0.15, colour = "grey40", size = 0.8) +
  facet_grid(~ species, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k"),
                     limits = c(0, 1e6)) +
  labs(x = NULL, y = "Nuclei counts") +
  theme_light(base_size = 7) +
  theme(strip.text = element_text(face = "italic"))

ggsave("fig3a.pdf", plot = p, width = 90, height = 70,
       units = "mm", device = cairo_pdf)