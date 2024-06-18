plot(iris$Petal.Length, iris$Petal.Width)
result <- lm(Sepal.Width ~ Sepal.Length, data = iris)
summary(result)
#三番目のコミット