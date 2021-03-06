#mi=200,400,600,800, 1000, 1200, 1400, 1600, 1800, 2000
#r=10
x <- c(200,400,600,800, 1000, 1200, 1400, 1600)
y1 <- c(664.973, 659.386, 656.069, 658.279, 651.19, 653.6, 647.412, 644.176)

#r=30
y2 <- c(622.738, 626.64, 636.993, 640.369, 643.742, 646.574, 646.458,642.047)

#r=50
y3 <- c(366.319, 157.386, 100.638, 82.7069, 81.1712, 83.301, 77.2944, 77.2944)

plot(x, y3, xlab="Broj iteracija", ylab="ACO fitness", type = "b", col = "green", ylim=c(75, 665), lwd = 4, main="Iris")
lines(x, y2, col="blue", lwd = 4, type = "b")
lines(x, y1, col="red", lwd = 4, type = "b")
legend(1600, 500, legend=c("r=10", "r=30", "r=50"), col=c("red", "blue", "green"),lty=1,  cex=1, text.font=2, lwd=3)