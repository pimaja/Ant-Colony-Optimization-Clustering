# K-Means Clustering

# Importing the dataset
dataset = read.csv('Wine.csv')
dataset = dataset[2:14]

# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
# library(caTools)
# set.seed(123)
# split = sample.split(dataset$DependentVariable, SplitRatio = 0.8)
# training_set = subset(dataset, split == TRUE)
# test_set = subset(dataset, split == FALSE)

# Feature Scaling
# training_set = scale(training_set)
# test_set = scale(test_set)

# Using the elbow method to find the optimal number of clusters

# Fitting K-Means to the dataset
set.seed(29)
kmeans = kmeans(x = dataset, centers = 3)
y_kmeans = kmeans$cluster

#Print to file
sink("output_wine.txt")
print(y_kmeans)
sink()

# Visualising data
boxplot(dataset[,2:8])
