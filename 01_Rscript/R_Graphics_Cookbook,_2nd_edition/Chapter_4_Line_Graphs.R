
# 4.1 Making a Basic Line Graph -------------------------------------------

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line()


BOD1 <- BOD  # Make a copy of the data
BOD1$Time <- factor(BOD1$Time)

ggplot(BOD1, aes(x = Time, y = demand, group = 1)) +  #  When the x variable is a factor, you must also use aes(group=1)
  geom_line()

# These have the same result
ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line() +
  ylim(0, max(BOD$demand))

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line() +
  expand_limits(y = 0)

# 4.2 Adding Points to a Line Graph ---------------------------------------

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line() +
  geom_point()

library(gcookbook) # Load gcookbook for the worldpop data set

ggplot(worldpop, aes(x = Year, y = Population)) +
  geom_line() +
  geom_point()

# Same with a log y-axis
ggplot(worldpop, aes(x = Year, y = Population)) +
  geom_line() +
  geom_point() +
  scale_y_log10()

# 4.3 Making a Line Graph with Multiple Lines -----------------------------

library(gcookbook) # Load gcookbook for the tg data set

# Map supp to colour
ggplot(tg, aes(x = dose, y = length, colour = supp)) +
  geom_line()

# Map supp to linetype
ggplot(tg, aes(x = dose, y = length, linetype = supp)) +
  geom_line()

ggplot(tg, aes(x = factor(dose), y = length, colour = supp, group = supp)) +
  geom_line()

ggplot(tg, aes(x = factor(dose), y = length, colour = supp)) + geom_line()
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?

ggplot(tg, aes(x = dose, y = length)) +
  geom_line()

ggplot(tg, aes(x = dose, y = length, shape = supp)) +
  geom_line() +
  geom_point(size = 4)  # Make the points a little larger

ggplot(tg, aes(x = dose, y = length, fill = supp)) +
  geom_line() +
  geom_point(size = 4, shape = 21)  # Also use a point with a color fill

ggplot(tg, aes(x = dose, y = length, shape = supp)) +
  geom_line(position = position_dodge(0.2)) +           # Dodge lines by 0.2
  geom_point(position = position_dodge(0.2), size = 4)  # Dodge points by 0.2

# 4.4 Changing the Appearance of Lines ------------------------------------

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line(linetype = "dashed", size = 1, colour = "blue")

library(gcookbook)  # Load gcookbook for the tg data set

ggplot(tg, aes(x = dose, y = length, colour = supp)) +
  geom_line() +
  scale_colour_brewer(palette = "Set1")

# If both lines have the same properties, you need to specify a variable to
# use for grouping
ggplot(tg, aes(x = dose, y = length, group = supp)) +
  geom_line(colour = "darkgreen", size = 1.5)

# Since supp is mapped to colour, it will automatically be used for grouping
ggplot(tg, aes(x = dose, y = length, colour = supp)) +
  geom_line(linetype = "dashed") +
  geom_point(shape = 22, size = 3, fill = "white")

# 4.5 Changing the Appearance of Points -----------------------------------

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line() +
  geom_point(size = 4, shape = 22, colour = "darkred", fill = "pink")

ggplot(BOD, aes(x = Time, y = demand)) +
  geom_line() +
  geom_point(size = 4, shape = 21, fill = "white")

library(gcookbook)  # Load gcookbook for the tg data set

# Save the position_dodge specification because we'll use it multiple times
pd <- position_dodge(0.2)

ggplot(tg, aes(x = dose, y = length, fill = supp)) +
  geom_line(position = pd) +
  geom_point(shape = 21, size = 3, position = pd) +
  scale_fill_manual(values = c("black","white"))

# 4.6 Making a Graph with a Shaded Area -----------------------------------

# Convert the sunspot.year data set into a data frame for this example
sunspotyear <- data.frame(
  Year     = as.numeric(time(sunspot.year)),
  Sunspots = as.numeric(sunspot.year)
)

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) +
  geom_area()

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) +
  geom_area(colour = "black", fill = "blue", alpha = .2)

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) +
  geom_area(fill = "blue", alpha = .2) +
  geom_line()

# 4.7 Making a Stacked Area Graph -----------------------------------------

library(gcookbook) # Load gcookbook for the uspopage data set

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) +
  geom_area()

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) +
  geom_area(colour = "black", size = .2, alpha = .4) +
  scale_fill_brewer(palette = "Blues")

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = dplyr::desc(AgeGroup))) +
  geom_area(colour = NA, alpha = .4) +
  scale_fill_brewer(palette = "Blues") +
  geom_line(position = "stack", size = .2)

# 4.8 Making a Proportional Stacked Area Graph ----------------------------

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = .4) +
  scale_fill_brewer(palette = "Blues")

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = .4) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = scales::percent)

# 4.9 Adding a Confidence Region ------------------------------------------

library(gcookbook) # Load gcookbook for the climate data set
library(dplyr)

# Grab a subset of the climate data
climate_mod <- climate %>%
  filter(Source == "Berkeley") %>%
  select(Year, Anomaly10y, Unc10y)

climate_mod
#>     Year Anomaly10y Unc10y
#> 1   1800     -0.435  0.505
#> 2   1801     -0.453  0.493
#> 3   1802     -0.460  0.486
#>  ...<199 more rows>...
#> 203 2002      0.856  0.028
#> 204 2003      0.869  0.028
#> 205 2004      0.884  0.029

# Shaded region
ggplot(climate_mod, aes(x = Year, y = Anomaly10y)) +
  geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = 0.2) +
  geom_line()

# With a dotted line for upper and lower bounds
ggplot(climate_mod, aes(x = Year, y = Anomaly10y)) +
  geom_line(aes(y = Anomaly10y - Unc10y), colour = "grey50", linetype = "dotted") +
  geom_line(aes(y = Anomaly10y + Unc10y), colour = "grey50", linetype = "dotted") +
  geom_line()

