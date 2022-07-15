
library(ggplot2)
# 2.1 Creating a Scatter Plot -------------------------------------------------

plot(mtcars$wt, mtcars$mpg)



ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point()

ggplot(data = NULL, aes(x = mtcars$wt, y = mtcars$mpg)) +
  geom_point()

# 2.2 Creating a Line Graph -----------------------------------------------

plot(pressure$temperature, pressure$pressure, type = "l")

plot(pressure$temperature, pressure$pressure, type = "l")
points(pressure$temperature, pressure$pressure)

lines(pressure$temperature, pressure$pressure/2, col = "red")
points(pressure$temperature, pressure$pressure/2, col = "red")


ggplot(pressure, aes(x = temperature, y = pressure)) +
  geom_line()

ggplot(pressure, aes(x = temperature, y = pressure)) +
  geom_line() +
  geom_point()

# 2.3 Creating a Bar Graph ------------------------------------------------

barplot(BOD$demand, names.arg = BOD$Time)


####  geom_bar(stat = "identity")  == geom_col()

# Bar graph of values. This uses the BOD data frame, with the
# "Time" column for x values and the "demand" column for y values.
ggplot(BOD, aes(x = Time, y = demand)) +
  geom_col()

# Convert the x variable to a factor, so that it is treated as discrete
ggplot(BOD, aes(x = factor(Time), y = demand)) +
  geom_col()


# Bar graph of counts This uses the mtcars data frame, with the "cyl" column for
# x position. The y position is calculated by counting the number of rows for
# each value of cyl.
ggplot(mtcars, aes(x = cyl)) +
  geom_bar()

# Bar graph of counts
ggplot(mtcars, aes(x = factor(cyl))) +
  geom_bar()

# 2.4 Creating a Histogram ------------------------------------------------

hist(mtcars$mpg)

# Specify approximate number of bins with breaks
hist(mtcars$mpg, breaks = 10)


ggplot(mtcars, aes(x = mpg)) +
  geom_histogram()
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

# With wider bins
ggplot(mtcars, aes(x = mpg)) +
  geom_histogram(binwidth = 4)

# 2.5 Creating a Box Plot -------------------------------------------------

plot(ToothGrowth$supp, ToothGrowth$len)

# Formula syntax
boxplot(len ~ supp, data = ToothGrowth)

# Put interaction of two variables on x-axis
boxplot(len ~ supp + dose, data = ToothGrowth)


# library(ggplot2)
ggplot(ToothGrowth, aes(x = supp, y = len)) +
  geom_boxplot()

ggplot(ToothGrowth, aes(x = interaction(supp, dose), y = len)) +
  geom_boxplot()

# 2.6 Plotting a Function Curve -------------------------------------------

curve(x^3 - 5*x, from = -4, to = 4)

# Plot a user-defined function
myfun <- function(xvar) {
  1 / (1 + exp(-xvar + 10))
}
curve(myfun(x), from = 0, to = 20)
# Add a line:
curve(1 - myfun(x), add = TRUE, col = "red")

library(ggplot2)
# This sets the x range from 0 to 20
ggplot(data.frame(x = c(0, 20)), aes(x = x)) +
  stat_function(fun = myfun, geom = "line")
