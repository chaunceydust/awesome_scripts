
library(ggplot2)
# 3.1 Making a Basic Bar Graph --------------------------------------------

library(gcookbook)  # Load gcookbook for the pg_mean data set
ggplot(pg_mean, aes(x = group, y = weight)) +
  geom_col()

# There's no entry for Time == 6
BOD
#>   Time demand
#> 1    1    8.3
#> 2    2   10.3
#> 3    3   19.0
#> 4    4   16.0
#> 5    5   15.6
#> 6    7   19.8

# Time is numeric (continuous)
str(BOD)
#> 'data.frame':    6 obs. of  2 variables:
#>  $ Time  : num  1 2 3 4 5 7
#>  $ demand: num  8.3 10.3 19 16 15.6 19.8
#>  - attr(*, "reference")= chr "A1.4, p. 270"

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_col()

# Convert Time to a discrete (categorical) variable with factor()
ggplot(BOD, aes(x = factor(Time), y = demand)) +
  geom_col()

ggplot(pg_mean, aes(x = group, y = weight)) +
  geom_col(fill = "lightblue", colour = "black")

# 3.2 Grouping Bars Together ----------------------------------------------

library(gcookbook)  # Load gcookbook for the cabbage_exp data set
cabbage_exp
#>   Cultivar Date Weight        sd  n         se
#> 1      c39  d16   3.18 0.9566144 10 0.30250803
#> 2      c39  d20   2.80 0.2788867 10 0.08819171
#> 3      c39  d21   2.74 0.9834181 10 0.31098410
#> 4      c52  d16   2.26 0.4452215 10 0.14079141
#> 5      c52  d20   3.11 0.7908505 10 0.25008887
#> 6      c52  d21   1.47 0.2110819 10 0.06674995

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(position = "dodge")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_brewer(palette = "Pastel1")

ce <- cabbage_exp[1:5, ]
ce
#>   Cultivar Date Weight        sd  n         se
#> 1      c39  d16   3.18 0.9566144 10 0.30250803
#> 2      c39  d20   2.80 0.2788867 10 0.08819171
#> 3      c39  d21   2.74 0.9834181 10 0.31098410
#> 4      c52  d16   2.26 0.4452215 10 0.14079141
#> 5      c52  d20   3.11 0.7908505 10 0.25008887

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_brewer(palette = "Pastel1")

# 3.3 Making a Bar Graph of Counts ----------------------------------------

# Equivalent to using geom_bar(stat = "bin")
ggplot(diamonds, aes(x = cut)) +
  geom_bar()

# 3.4 Using Colors in a Bar Graph -----------------------------------------

library(gcookbook) # Load gcookbook for the uspopchange data set
library(dplyr)

upc <- uspopchange %>%
  arrange(desc(Change)) %>%
  slice(1:10)

upc
#>             State Abb Region Change
#> 1          Nevada  NV   West   35.1
#> 2         Arizona  AZ   West   24.6
#> 3            Utah  UT   West   23.8
#>  ...<4 more rows>...
#> 8         Florida  FL  South   17.6
#> 9        Colorado  CO   West   16.9
#> 10 South Carolina  SC  South   15.3

ggplot(upc, aes(x = Abb, y = Change, fill = Region)) +
  geom_col()

ggplot(upc, aes(x = reorder(Abb, - Change), y = Change, fill = Region)) +
  geom_col(colour = "black") +
  scale_fill_manual(values = c("#669933", "#FFCC66")) +
  xlab("State")

# 3.5 Coloring Negative and Positive Bars Differently ---------------------

library(gcookbook) # Load gcookbook for the climate data set
library(dplyr)

climate_sub <- climate %>%
  filter(Source == "Berkeley" & Year >= 1900) %>%
  mutate(pos = Anomaly10y >= 0)

climate_sub
#>       Source Year Anomaly1y Anomaly5y Anomaly10y Unc10y   pos
#> 1   Berkeley 1900        NA        NA     -0.171  0.108 FALSE
#> 2   Berkeley 1901        NA        NA     -0.162  0.109 FALSE
#> 3   Berkeley 1902        NA        NA     -0.177  0.108 FALSE
#>  ...<99 more rows>...
#> 103 Berkeley 2002        NA        NA      0.856  0.028  TRUE
#> 104 Berkeley 2003        NA        NA      0.869  0.028  TRUE
#> 105 Berkeley 2004        NA        NA      0.884  0.029  TRUE

# Notice that we use position=“identity” with the bars. This will prevent a warning message about stacking not being well defined for negative numbers:
ggplot(climate_sub, aes(x = Year, y = Anomaly10y, fill = pos)) +
  geom_col(position = "identity")

ggplot(climate_sub, aes(x = Year, y = Anomaly10y, fill = pos)) +
  geom_col(position = "identity", colour = "black", size = 0.25) +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = FALSE)

# 3.6 Adjusting Bar Width and Spacing -------------------------------------

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(width = 0.5, position = position_dodge(0.7))

# 3.7 Making a Stacked Bar Graph ------------------------------------------

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  guides(fill = guide_legend(reverse = TRUE))

# 3.8 Making a Proportional Stacked Bar Graph -----------------------------

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1")

library(gcookbook)
library(dplyr)

cabbage_exp
#>   Cultivar Date Weight        sd  n         se
#> 1      c39  d16   3.18 0.9566144 10 0.30250803
#> 2      c39  d20   2.80 0.2788867 10 0.08819171
#> 3      c39  d21   2.74 0.9834181 10 0.31098410
#> 4      c52  d16   2.26 0.4452215 10 0.14079141
#> 5      c52  d20   3.11 0.7908505 10 0.25008887
#> 6      c52  d21   1.47 0.2110819 10 0.06674995

# Do a group-wise transform(), splitting on "Date"
ce <- cabbage_exp %>%
  group_by(Date) %>%
  mutate(percent_weight = Weight / sum(Weight) * 100)

ce
#> # A tibble: 6 x 7
#> # Groups:   Date [3]
#>   Cultivar Date  Weight    sd     n     se percent_weight
#>   <fct>    <fct>  <dbl> <dbl> <int>  <dbl>          <dbl>
#> 1 c39      d16     3.18 0.957    10 0.303            58.5
#> 2 c39      d20     2.8  0.279    10 0.0882           47.4
#> 3 c39      d21     2.74 0.983    10 0.311            65.1
#> 4 c52      d16     2.26 0.445    10 0.141            41.5
#> 5 c52      d20     3.11 0.791    10 0.250            52.6
#> 6 c52      d21     1.47 0.211    10 0.0667           34.9

# 3.9 Adding Labels to a Bar Graph ----------------------------------------

library(gcookbook) # Load gcookbook for the cabbage_exp data set

# Below the top
ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
  geom_col() +
  geom_text(aes(label = Weight), vjust = 1.5, colour = "white")

# Above the top
ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
  geom_col() +
  geom_text(aes(label = Weight), vjust = -0.2)


ggplot(mtcars, aes(x = factor(cyl))) +
  geom_bar() +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white")


# Adjust y limits to be a little higher
ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
  geom_col() +
  geom_text(aes(label = Weight), vjust = -0.2) +
  ylim(0, max(cabbage_exp$Weight) * 1.05)

# Map y positions slightly above bar top - y range of plot will auto-adjust
ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
  geom_col() +
  geom_text(aes(y = Weight + 0.1, label = Weight))


ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = Weight),
    colour = "white", size = 5,
    vjust = 1.5, position = position_dodge(.9)
  )


library(dplyr)

# Sort by the Date and Cultivar columns
ce <- cabbage_exp %>%
  # arrange(Date, Cultivar)
  arrange(Date, rev(Cultivar))


# Get the cumulative sum
ce <- ce %>%
  group_by(Date) %>%
  mutate(label_y = cumsum(Weight))

ce
#> # A tibble: 6 x 7
#> # Groups:   Date [3]
#>   Cultivar Date  Weight    sd     n     se label_y
#>   <fct>    <fct>  <dbl> <dbl> <int>  <dbl>   <dbl>
#> 1 c52      d16     2.26 0.445    10 0.141     2.26
#> 2 c39      d16     3.18 0.957    10 0.303     5.44
#> 3 c52      d20     3.11 0.791    10 0.250     3.11
#> 4 c39      d20     2.8  0.279    10 0.0882    5.91
#> 5 c52      d21     1.47 0.211    10 0.0667    1.47
#> 6 c39      d21     2.74 0.983    10 0.311     4.21

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col() +
  geom_text(aes(y = label_y, label = Weight), vjust = 1.5, colour = "white")



ce <- cabbage_exp %>%
  arrange(Date, rev(Cultivar))

# Calculate y position, placing it in the middle
ce <- ce %>%
  group_by(Date) %>%
  mutate(label_y = cumsum(Weight) - 0.5 * Weight)

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col() +
  geom_text(aes(y = label_y, label = Weight), colour = "white")

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
  geom_col(colour = "black") +
  geom_text(aes(y = label_y, label = paste(format(Weight, nsmall = 2), "kg")), size = 4) +
  scale_fill_brewer(palette = "Pastel1")




# 3.10 Making a Cleveland Dot Plot ----------------------------------------

library(gcookbook) # Load gcookbook for the tophitters2001 data set
tophit <- tophitters2001[1:25, ] # Take the top 25 from the tophitters data set

ggplot(tophit, aes(x = avg, y = name)) +
  geom_point()

tophit[, c("name", "lg", "avg")]
#>             name lg    avg
#> 1   Larry Walker NL 0.3501
#> 2  Ichiro Suzuki AL 0.3497
#> 3   Jason Giambi AL 0.3423
#>  ...<19 more rows>...
#> 23  Jeff Cirillo NL 0.3125
#> 24   Jeff Conine AL 0.3111
#> 25   Derek Jeter AL 0.3111

ggplot(tophit, aes(x = avg, y = reorder(name, avg))) +
  geom_point(size = 3) +  # Use a larger dot
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
  )


ggplot(tophit, aes(x = reorder(name, avg), y = avg)) +
  geom_point(size = 3) +  # Use a larger dot
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey60", linetype = "dashed"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


# Get the names, sorted first by lg, then by avg
nameorder <- tophit$name[order(tophit$lg, tophit$avg)]

# Turn name into a factor, with levels in the order of nameorder
tophit$name <- factor(tophit$name, levels = nameorder)


ggplot(tophit, aes(x = avg, y = name)) +
  geom_segment(aes(yend = name), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = lg)) +
  scale_colour_brewer(palette = "Set1", limits = c("NL", "AL")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5)
  )

ggplot(tophit, aes(x = avg, y = name)) +
  geom_segment(aes(yend = name), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = lg)) +
  scale_colour_brewer(palette = "Set1", limits = c("NL", "AL"), guide = FALSE) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  facet_grid(lg ~ ., scales = "free_y", space = "free_y")
