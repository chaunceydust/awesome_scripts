
# 5.1 Making a Basic Scatter Plot -----------------------------------------

library(gcookbook) # Load gcookbook for the heightweight data set
library(dplyr)

# Show the head of the two columns we'll use in the plot
heightweight %>%
  select(ageYear, heightIn)

ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point(shape = 21, size = 1.5)


# 5.2 Grouping Points Together using Shapes or Colors ---------------------

library(gcookbook) # Load gcookbook for the heightweight data set

# Show the head of the three columns we'll use
heightweight %>%
  select(sex, ageYear, heightIn)

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) +
  geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) +
  geom_point() +
  scale_shape_manual(values = c(1,2)) +
  scale_colour_brewer(palette = "Set1")


# 5.3 Using Different Point Shapes ----------------------------------------

library(gcookbook) # Load gcookbook for the heightweight data set

ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point(shape = 3)

# Use slightly larger points and use custom values for the shape scale
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1, 4))

# Using the heightweight data set, create a new column that indicates if the
# child weighs < 100 or >= 100 pounds. We'll save this modified dataset as 'hw'.
hw <- heightweight %>%
  mutate(weightgroup = ifelse(weightLb < 100, "< 100", ">= 100"))

# Specify shapes with fill and color, and specify fill colors that includes an empty (NA) color
ggplot(hw, aes(x = ageYear, y = heightIn, shape = sex, fill = weightgroup)) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(
    values = c(NA, "black"),
    guide = guide_legend(override.aes = list(shape = 21))
  )

# 5.4 Mapping a Continuous Variable to Color or Size ----------------------

library(gcookbook) # Load gcookbook for the heightweight data set

# Show the head of the four columns we'll use
heightweight %>%
  select(sex, ageYear, heightIn, weightLb)

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = weightLb)) +
  geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb)) +
  geom_point()
