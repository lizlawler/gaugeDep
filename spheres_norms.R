x <- sort(runif(500, min = 0, max = 1))
y2 <- sqrt(1 - x^2)

plot(x, y2, type = "l", xlim = c(0,2), ylim = c(0,2), lwd = 2)
abline(v = 0, h = 0, lwd = 2)
points(1.25, 1, col = "red", pch = 19)
segments(x0 = 0, y0 = 0, x1 = 1.25, y1=1, lwd = 2, col = "red")
points(1.25/sqrt(1.25^2 + 1), 1/sqrt(1.25^2 + 1), col = "purple", pch = 19)
segments(x0 = 0, y0 = 0, x1 = 1.25/sqrt(1.25^2 + 1), y1=1/sqrt(1.25^2 + 1), 
         lwd = 2, col = "purple", lty = "dashed")
segments(x0 = 1.25/sqrt(1.25^2 + 1), y0 = 0, x1 = 1.25/sqrt(1.25^2 + 1), y1=1/sqrt(1.25^2 + 1), 
         lwd = 2, col = "purple", lty = "dashed")

y1_1 <- 1 - x[x >=0]
y1_2 <- 1 + x[x <= 0]
y1_3 <- -(1 + x[x <= 0])
y1_4 <- x[x>=0] - 1
lines(x[x >=0], y1_1, lwd = 2, col = "blue")
# lines(x[x <= 0], y1_2, lwd = 2, col = "blue")
# lines(x[x <= 0], y1_3, lwd = 2, col = "blue")
# lines(x[x >=0], y1_4, lwd = 2, col = "blue")
points(1.25/2.25, 1/2.25, pch = 19, col = "orange")
segments(x0 = 0, y0 = 0, x1 = 1.25/2.25, y1=1/2.25, 
         lwd = 2, col = "orange", lty = "longdash")
segments(x0 = 1.25/2.25, y0 = 0, x1 = 1.25/2.25, y1=1/2.25, 
         lwd = 2, col = "orange")

y1_1prime <- 2.25 - x[x >=0]
y1_2prime <- 2.25 + x[x <= 0]
y1_3prime <- -(2.25 + x[x <= 0])
y1_4prime <- x[x>=0] - 2.25

lines(x[x >=0], y1_1prime, lwd = 2, col = "darkgreen")
lines(x[x <= 0], y1_2prime, lwd = 2, col = "darkgreen")
lines(x[x <= 0], y1_3prime, lwd = 2, col = "darkgreen")
lines(x[x >=0], y1_4prime, lwd = 2, col = "darkgreen")

r2 <- sqrt(1.25^2 + 1)
w2_1 <- 1.25/r2
w2_2 <- 1/r2
sin(acos(w2_1))
asin(w2_2)

r1 <- 1.25 + 1
w1_1 <- 1.25/r1
w1_2 <- 1/r1
atan(w1_2/w1_1)
