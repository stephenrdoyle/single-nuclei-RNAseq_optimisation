# Nuclei counts from optimised protocol

setwd("/Users/sd21/Documents/WORK/02_PROJECTS/single_cell-nuclei_transcriptomics/")


# Figure 1C

library(tidyverse)
options(scipen = 999)

protocol <- tribble(
  ~method,        ~nuclei,
  "Celegans",     c(8800, 5940, 8360),
  "HCA",          c(15300, 20280, 22100),
  "Frankenstein", c(14620, 14175, 12240)
) |>
  unnest(nuclei) |>
  mutate(method = fct_inorder(method))

## Welch's t-tests 
brackets <- tribble(
  ~a,         ~b,             ~x, ~xend, ~y,
  "Celegans", "HCA",           1,  2,    24500,
  "Celegans", "Frankenstein",  1,  3,    26500
) |>
  rowwise() |>
  mutate(
    p     = t.test(protocol$nuclei[protocol$method == a],
                   protocol$nuclei[protocol$method == b])$p.value,
    label = if_else(p < 0.001, "p < 0.001",
                    paste0("p = ", format(round(p, 3), nsmall = 3)))
  ) |>
  ungroup()

p <- ggplot(protocol, aes(method, nuclei)) +
  stat_summary(fun = mean, geom = "col",
               fill = "skyblue", alpha = 0.7, colour = "darkgrey",
               linewidth = 0.3, width = 0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", colour = "orange",
               width = 0.25, linewidth = 0.4) +
  geom_jitter(width = 0.12, colour = "grey40", size = 0.8) +
  geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.25) +
  geom_text(data = brackets,
            aes(x = (x + xend) / 2, y = y + 400, label = label),
            inherit.aes = FALSE, size = 2, vjust = 0) +
  scale_y_continuous(limits = c(0, 28500), breaks = seq(0, 20000, 5000),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Nuclei counts") +
  theme_light(base_size = 7)
p
ggsave("fig1c.pdf", plot = p, width = 90, height = 70,
       units = "mm", device = cairo_pdf)


# statistical comparison between protocols
library(tidyverse)

protocol <- tribble(
  ~method,        ~nuclei,
  "Celegans",     c(8800, 5940, 8360),
  "HCA",          c(15300, 20280, 22100),
  "Frankenstein", c(14620, 14175, 12240)
) |>
  unnest(nuclei)

combn(unique(protocol$method), 2, simplify = FALSE) |>
  map_dfr(function(pair) {
    a <- filter(protocol, method == pair[1])$nuclei
    b <- filter(protocol, method == pair[2])$nuclei
    tt <- t.test(a, b)                    # Welch by default
    tibble(
      comparison = paste(pair[1], "vs", pair[2]),
      mean_a     = mean(a),
      mean_b     = mean(b),
      fold       = mean(b) / mean(a),
      t          = unname(tt$statistic),
      df         = unname(tt$parameter),
      p          = tt$p.value,
      ci_low     = tt$conf.int[1],
      ci_high    = tt$conf.int[2]
    )
  }) |>
  mutate(across(c(t, df, fold), \(x) round(x, 3)),
         p = signif(p, 3)) |>
  print(n = Inf)




# Figure 3A
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