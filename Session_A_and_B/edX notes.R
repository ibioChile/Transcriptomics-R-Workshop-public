
# 1.1

# installing the dslabs package
install.packages("dslabs")
# loading the dslabs package into the R session
library(dslabs)



install.packages("dslabs")  # to install a single package
install.packages(c("tidyverse", "dslabs")) # to install two packages at the same time
installed.packages() # to see the list of all installed packages



library(tidyverse)
library(dslabs)
data(murders)

murders %>%
  ggplot(aes(population, total, label=abb, color=region)) +
  geom_label()


natey <- function(a,b=10) {
  a * b
}


mean(c(1,2,4))





# Playing with course data ------------------------------------------------



stop("these data are not available fo you all, this example will only work on my computer")
df <- read.delim("data/students.txt")
names(df)

df$Dirección.de.correo.electrónico <- NULL
df$X.Por.qué.quieres.tomar.este.workshop. <- NULL
df$Nombre <- NULL

tab <- table(df$Nivel.de.conocimiento.de.R)
barplot(tab, col=c('green','orange','grey'),
        main="Experience with R",
        ylab="Students")

tab <- table(df$Laboratorio, df$Tipo.de.investigador)
rownames(tab)[rownames(tab) == ''] <- "No Lab Listed"
par(mar=c(10,5,2,2))
barplot(t(tab), las=2, ylab='Students', main="Laboratory", beside=T, col=rainbow(3), legend.text=colnames(tab))

head(df)
df$date <- as.Date(df$Marca.temporal, "%d/%m/%Y %H:%M")
df <- df[order(df$date),]
df$date_rank <- 1:nrow(df)


graphics.off()
plot(df$date, df$date_rank, type='l',
     xlab='Date', ylab='Students',
     main="Sign up dates")


tab <- table(df$X.Con.qué.computador.vas.a.trabajar.)
barplot(tab)





nate <- function(a, b) {
  return(a  + b)
}

nate(1,2)

df



head(df)

df$Laboratorio
df[["Laboratorio"]]


get_lab <- function() {
  # getting the laboratory field of the data frame (trivial)
  
  log.vec <- names(df) == "Laboratorio"
  print(log.vec)
  
  idx     <- which(log.vec)
  print(idx)
  
  return(df[,idx])
}
get_lab()





df[[idx]]
df[, "Laboratorio"]

# which ~ order


idx <- which(df$Laboratorio == "R. Gutiérrez" )
df$Laboratorio[idx]





## Types of logical evaluation
1 == 3
1 == 1
1 != 3
1 != 1

1 >  0
1 >= 1
1 >= 2
1 <= 2

## Vector evaluaation

c(1,2,3,4,5) == 5
c(1,2,3) == c(1,2,3)

## Comparing entire vectors
identical(1, 1)
identical(c(1,2,3), c(1,2,3))


## "And" "Y" and/y "Or" "O" 
# | <- o or
# & <- y and

1 == 2 | 1 == 1
1 == 2 & 1 == 1


## Simple loop
for (i in 10:1) {
  print(i)
}


## if else evaluation

if (1==1) { print("This worked!!") }

if (1==1) { 
  print('A')
} else if (1 == 2) {
  print("B")
} else {
  print("C")
}




test.df <- data.frame(a=1:10, b=10:19)





# 1.2 R basics ------------------------------------------------------------

# assigning values to variables
a <-1
b <-1
c <--1

# solving the quadratic equation
(-b + sqrt(b^2 - 4*a*c))/(2*a)
(-b - sqrt(b^2 - 4*a*c))/(2*a)


# 1.2 data types ----------------------------------------------------------

# loading the dslabs package and the murders dataset
library(dslabs)
data(murders)

# determining that the murders dataset is of the "data frame" class
class(murders)
# finding out more about the structure of the object
str(murders)
# showing the first 6 lines of the dataset
head(murders)

# using the accessor operator to obtain the population column
murders$population
# displaying the variable names in the murders dataset
names(murders)
# determining how many entries are in a vector
pop <- murders$population
length(pop)
# vectors can be of class numeric and character
class(pop)
class(murders$state)

# logical vectors are either TRUE or FALSE
z <- 3 == 2
z
class(z)

# factors are another type of class
class(murders$region)
# obtaining the levels of a factor
levels(murders$region)



# 2.1 vectors -------------------------------------------------------------


# We may create vectors of class numeric or character with the concatenate function
codes <- c(380, 124, 818)
country <- c("italy", "canada", "egypt")

# We can also name the elements of a numeric vector
# Note that the two lines of code below have the same result
codes <- c(italy = 380, canada = 124, egypt = 818)
codes <- c("italy" = 380, "canada" = 124, "egypt" = 818)

# We can also name the elements of a numeric vector using the names() function
codes <- c(380, 124, 818)
country <- c("italy","canada","egypt")
names(codes) <- country

# Using square brackets is useful for subsetting to access specific elements of a vector
codes[2]
codes[c(1,3)]
codes[1:2]

# If the entries of a vector are named, they may be accessed by referring to their name
codes["canada"]
codes[c("egypt","italy")]



## coercion
x <- c(1, "canada", 3)



# 2.2 sorting -------------------------------------------------------------

library(dslabs)
data(murders)
sort(murders$total)

x <- c(31, 4, 15, 92, 65)
x
sort(x)    # puts elements in order

index <- order(x)    # returns index that will put x in order
x[index]    # rearranging by this index puts elements in order
order(x)

murders$state[1:10]
murders$abb[1:10]

index <- order(murders$total)
murders$abb[index]    # order abbreviations by total murders

max(murders$total)    # highest number of total murders
i_max <- which.max(murders$total)    # index with highest number of murders
murders$state[i_max]    # state name with highest number of total murders

x <- c(31, 4, 15, 92, 65)
x
rank(x)    # returns ranks (smallest to largest)



# 2.3 vector arithmatic ---------------------------------------------------

# The name of the state with the maximum population is found by doing the following
murders$state[which.max(murders$population)]

# how to obtain the murder rate
murder_rate <- murders$total / murders$population * 100000

# ordering the states by murder rate, in decreasing order
murders$state[order(murder_rate, decreasing=TRUE)]










# 3.1 indexing ------------------------------------------------------------



# defining murder rate as before
murder_rate <- murders$total / murders$population * 100000
# creating a logical vector that specifies if the murder rate in that state is less than or equal to 0.71
index <- murder_rate <= 0.71
# determining which states have murder rates less than or equal to 0.71
murders$state[index]
# calculating how many states have a murder rate less than or equal to 0.71
sum(index)

# creating the two logical vectors representing our conditions
west <- murders$region == "West"
safe <- murder_rate <= 1
# defining an index and identifying states with both conditions true
index <- safe & west
murders$state[index]




x <- c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE)
which(x)    # returns indices that are TRUE

# to determine the murder rate in Massachusetts we may do the following
index <- which(murders$state == "Massachusetts")
index
murder_rate[index]

# to obtain the indices and subsequent murder rates of New York, Florida, Texas, we do:
index <- match(c("New York", "Florida", "Texas"), murders$state)
index
murders$state[index]
murder_rate[index]

x <- c("a", "b", "c", "d", "e")
y <- c("a", "d", "f")
y %in% x

# to see if Boston, Dakota, and Washington are states
c("Boston", "Dakota", "Washington") %in% murders$state

## %in% is basically a much more simple way of writing:
!is.na(match(y,x))



# 3.2 basic data wrangling ------------------------------------------------


# installing and loading the dplyr package
install.packages("dplyr")
library(dplyr)

# adding a column with mutate
library(dslabs)
data("murders")
murders <- mutate(murders, rate = total / population * 100000)

# subsetting with filter
filter(murders, rate <= 0.71)

# selecting columns with select
new_table <- select(murders, state, region, rate)

# using the pipe
murders %>% select(state, region, rate) %>% filter(rate <= 0.71) %>% filter(region == "Northeast")



# creating a data frame with stringAsFactors = FALSE
grades <- data.frame(names = c("John", "Juan", "Jean", "Yao"), 
                     exam_1 = c(95, 80, 90, 85), 
                     exam_2 = c(90, 85, 85, 90),
                     stringsAsFactors = FALSE)





# 3.3 basic plots ---------------------------------------------------------

murders$rate <- murders$total / murders$population * 100000


# a simple scatterplot of total murders versus population
x <- murders$population /10^6
y <- murders$total
plot(x, y)

# a histogram of murder rates
hist(murders$rate)
?hist

# boxplots of murder rates by region
boxplot(rate~region, data = murders)




# 4.2 Conditionals ------------------------------------

# an example showing the general structure of an if-else statement
a <- 0
if(a!=0){
  print(1/a)
} else{
  print("No reciprocal for 0.")
}

# an example that tells us which states, if any, have a murder rate less than 0.5
library(dslabs)
data(murders)
murder_rate <- murders$total / murders$population*100000
ind <- which.min(murder_rate)
if(murder_rate[ind] < 0.5){
  print(murders$state[ind]) 
} else{
  print("No state has murder rate that low")
}

# changing the condition to < 0.25 changes the result
if(murder_rate[ind] < 0.25){
  print(murders$state[ind]) 
} else{
  print("No state has a murder rate that low.")
}

# the ifelse() function works similarly to an if-else conditional
a <- 0
ifelse(a > 0, 1/a, NA)

# the ifelse() function is particularly useful on vectors
a <- c(0,1,2,-4,5)
result <- ifelse(a > 0, 1/a, NA)

# the ifelse() function is also helpful for replacing missing values
data(na_example)
no_nas <- ifelse(is.na(na_example), 0, na_example) 
sum(is.na(no_nas))

# the any() and all() functions evaluate logical vectors
z <- c(TRUE, TRUE, FALSE)
any(z)
all(z)

# !any(z) is this different than all(z) ? not sure




# 4.3 functions -----------------------------------------------------------



# example of defining a function to compute the average of a vector x
avg <- function(x){
  s <- sum(x)
  n <- length(x)
  s/n
}

# we see that the above function and the pre-built R mean() function are identical
x <- 1:100
identical(mean(x), avg(x))

# variables inside a function are not defined in the workspace
s <- 3
avg(1:10)
s

# the general form of a function
# my_function <- function(VARIABLE_NAME){
#   perform operations on VARIABLE_NAME and calculate VALUE
#   VALUE
# }

# functions can have multiple arguments as well as default values
avg <- function(x, arithmetic = TRUE){
  n <- length(x)
  ifelse(arithmetic, sum(x)/n, prod(x)^(1/n))
}








# 4.4 for loops -----------------------------------------------------------


# creating a function that computes the sum of integers 1 through n
compute_s_n <- function(n){
  x <- 1:n
  return(sum(x))
  print("will this print?")
}

compute_s_n(10)


# a very simple for-loop
for (i in 1:5){
  print(i)
}



# a for-loop for our summation
m <- 25
s_n <- vector(length = m) # create an empty vector
for(n in 1:m){
  s_n[n] <- compute_s_n(n)
}

# creating a plot for our summation function
n <- 1:m
plot(n, s_n)

# a table of values comparing our function to the summation formula
head(data.frame(s_n = s_n, formula = n*(n+1)/2))

# overlaying our function with the summation formula
plot(n, s_n)
lines(n, n*(n+1)/2)


# Dplyr help --------------------------------------------------------------

# There are multiple ways to do everything in R. I find indicies much more simple than dplyr, so i'll explain them again here. Use whatever way makes more sense to you!

# remember, dataframes are indexed by row and column:  df[row,column]

data("murders")

# Mutate defines variables in a dataframe. We learned to do this in the first session.
# The following lines are equivalent:
murders <- mutate(murders, rate = total / population * 100000)
murders$rate <- murders$total / murders$population * 100000


# Filter selects rows of a dataframe, based on a logical or indicies
# For this example, we use the logial vector produced by: murders$rate <= 0.71
# The following lines are equivalent:
filter(murders, rate <= 0.71)
murders[murders$rate <= 0.71,]

# selecting columns with select, based on their names
# The following lines are equivalent:
select(murders, state, region, rate)
murders[,c('state', 'region', 'rate')]


# using the pipe, we can string multiple together
# remember: df[rows, columns]
# The following lines are equivalent:
murders %>% select(state, region, rate) %>% filter(rate <= 0.71) %>% filter(region == "Northeast")
murders[murders$rate <= 0.71 & murders$region == "Northeast", c("state", "region", "rate")]

# we could also break these filters out, no need to do this in one line!

filter1 = murders$rate <= 0.71
filter2 = murders$region == "Northeast"
columns = c("state", "region", "rate")

murders[filter1 & filter2, columns]



# Dataframe help ----------------------------------------

data("heights")

# dataframes have rows and columns, which can be accessed by name or index using the following format:
# df[row, column]

heights[c(1,2,3),c(1,2)] # extracts the first 3 rows

# notice that empty indicies are treated as "all"

new_heights <- heights[,]
identical(heights, new_heights)


# also, if you extract only one column, something weird happens... it returns as a vector (or factor)
# the following are the equivalent:
heights[,2]
heights$height 

# we can avoid this and keep the df shape through one more option
heights[,2, drop=F]
# however, most of the time you probably just want the vector or factor...






## Functions to have a look at the shape of a dataframe

nrow(heights)
ncol(heights)
dim(heights) # returns rows and cols
str(heights)
head(heights)
tail(heights)
length(heights) # does not work! dataframes have dimensions, not a length. This should only work on vectors/lists.


# Loop help ---------------------------------------------------------------

for (i in 1:5) {
  print(i)
  if (i > 4) {
    # this locates breaks out of the first enclosing loop. Note, it doesn't mean anything to the conditional clause that calls it (i > 4) {}, it only matters to the loop (for ()...).
    break
  }
}

i
# i == 4


# Function help -----------------------------------------------------------

# return can be used to stop a function and return something. This is more explicit than using the last line of a function to return.

my_funct <- function(j) { 
  print(j)
  return(j)
  print("this won't ever print")
}

my_funct()



# functions can be used to organize your code, and don't need to have inputs

plot_murder_rates <- function() { # notice, I gave no input variables. This function is entirely self-contained.
  data(murders)
  
  murders$murders_per_million <- murders$total / murders$population * 1000000
  
  boxplot(murders$murders_per_million ~ murders$region) 
  
  # I also highlight the murder rates in the states I've lived
  michigan <- murders[murders$abb == "MI",]
  points(3, michigan$murders_per_million, pch=19, col='red')
  text(3, michigan$murders_per_million, 'Michigan', pos=4, col='red')
  
  pennsylvania <- murders[murders$abb == "PA",]
  points(3, pennsylvania$murders_per_million, pch=19, col='darkseagreen4')
  text(3, pennsylvania$murders_per_million, 'Pennsylvania', pos=4, col='darkseagreen4')
  
}

plot_murder_rates() # note, my 'murders' object outside of the function is not changed.






# Pheatmap example --------------------------------------------------------

# Just showing how a headplotting package can be installed and tested using random data produced in R.


install.packages("pheatmap")

library(pheatmap)

# the following all comes from the help documentation, easily accessed with the following command:
?pheatmap

test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)












