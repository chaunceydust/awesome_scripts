setwd("C:/Users/93740/Desktop")

library(tidyverse)
library(ggpubr)
library(patchwork)

c <- read_csv("jibing_abandance.csv")
c <- c[2:4]
colnames(c) <- c("Country", "Groups", "abundance")
head(c)
c$Groups[which(c$Groups == "C")] <- "Healthy"
c$Groups[which(c$Groups == "D")] <- "Disease"
c$Groups <- factor(
  c$Groups,
  
  levels = c("Healthy", "Disease")
)

c <- c |> 
  separate(
    col = "Country",
    sep = "-",
    into = c("Country", "Disease")
  )
c$abundance[which(is.na(c$abundance))] <- 0

data.prevalence <- c
data.prevalence$abundance <- as.logical(data.prevalence$abundance)
data.prevalence2 <- data.prevalence |> 
  group_by(
    Country, Disease, Groups
  ) |> 
  summarise(
    prevalence = mean(abundance) * 100
  ) |> 
  ungroup()

p2 <- ggplot(data.prevalence2,
       aes(Country, prevalence, fill = Groups)) +
  geom_bar(
    stat = "identity",
    width = .6,
    # alpha = .85,
    position = position_dodge()
  ) +
  facet_grid(
    .~Disease,
    # nrow = 1,
    scales = "free_x",
    space = 'free'
  ) +
  labs(y = "prevalence (%)") +
  scale_fill_manual(
    values = c("#68C5DB", "#EA6653")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      # size = 1
    ),
    strip.text = element_text(
      face = "bold",
      size = rel(.7)
    )
  )


# 2 -----------------------------------------------------------------------

p1 <- ggplot(
  filter(c, abundance != 0),
  aes(Country, abundance, fill = Groups)
) +
  geom_boxplot(outlier.alpha = 1) +
  scale_y_log10(
    limits = c(0.0001, 1)
  ) +
  facet_grid(
    .~Disease,
    # nrow = 1,
    scales = "free_x",
    space = 'free'
  ) +
  labs(y = "relative abundace") +
  scale_fill_manual(
    values = c("#68C5DB", "#EA6653"),
  ) +
  stat_compare_means(
    aes(group = Groups),
    method = "t.test", 
    label = "p.signif",
    label.y = -0.35,
    size = 6.5,
    hide.ns = TRUE
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      # size = 1
    ),
    strip.text = element_text(
      face = "bold",
      size = rel(.7)
    )
    
  )

p1 / p2 +
  plot_annotation(tag_levels = "A") 

ggsave("diff.pdf",
       width = 8, height = 6)
ggsave("diff.tiff",
       width = 8, height = 6)
