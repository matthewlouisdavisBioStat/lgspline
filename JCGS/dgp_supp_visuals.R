set.seed(1234)

cols <- hcl.colors(80, "Spectral", rev = TRUE)

## Predictor transforms.
#  The surfaces are evaluated before the final range transform; the observed
#  predictors t1, ..., tq have different ranges, skewed margins, and mild
#  mixture structure.
u_to_t2 <- function(u1, u2) {
  t1 <- 36 * u1^1.65
  t2 <- -3 + 12 * u2^0.62
  cbind(t1, t2)
}

u_to_t4 <- function(u1, u2, u3, u4) {
  t1 <- 55 * u1^1.55
  t2 <- -4 + 16 * u2^0.70
  t3 <- 7 * u3^1.25
  t4 <- -2.2 + 7.0 * u4^0.82 + 0.35 * pmax(u4 - 0.58, 0)^2
  cbind(t1, t2, t3, t4)
}

generate_u2 <- function(n) {
  u1 <- rbeta(n, 1.4, 2.4)
  keep <- runif(n) < 0.25
  u1[keep] <- rbeta(sum(keep), 4.5, 1.8)

  u2 <- rbeta(n, 0.95, 2.2)
  keep <- runif(n) < 0.22
  u2[keep] <- rbeta(sum(keep), 5.5, 2.0)

  cbind(u1, u2)
}

generate_u4 <- function(n) {
  z1 <- rnorm(n)
  z2 <- 0.55 * z1 + sqrt(1 - 0.55^2) * rnorm(n)
  z3 <- -0.35 * z1 + 0.45 * z2 + sqrt(1 - 0.35^2) * rnorm(n)
  z4 <- 0.40 * z1 - 0.25 * z3 + sqrt(1 - 0.40^2) * rnorm(n)

  u1 <- pnorm(z1)
  u2 <- pnorm(z2)
  u3 <- pnorm(z3)
  u4 <- pnorm(z4)

  keep <- runif(n) < 0.20
  u2[keep] <- rbeta(sum(keep), 1.0, 3.2)
  keep <- runif(n) < 0.22
  u4[keep] <- rbeta(sum(keep), 6.0, 2.2)

  cbind(u1, u2, u3, u4)
}

## Response surfaces.
#  The nonsmooth terms are light kinks rather than discontinuities.
f_q2 <- function(u1, u2) {
  0.75 * sin(2 * pi * u1) * cos(1.5 * pi * u2) +
    1.30 * exp(-55 * (u2 - 0.58 - 0.16 * sin(2.4 * pi * u1))^2) +
    1.35 * exp(-38 * ((u1 - 0.23)^2 + 0.50 * (u2 - 0.78)^2)) -
    1.20 * exp(-70 * ((u1 - 0.72)^2 + (u2 - 0.28)^2)) +
    0.75 * (u1 - 0.5) * (u2 - 0.5) +
    0.22 * abs(u1 - 0.62) -
    0.22 * pmax(u2 - 0.72, 0) +
    0.28 * pmax(u1 + 0.55 * u2 - 1.08, 0)
}

f_q4 <- function(u1, u2, u3, u4) {
  z1 <- u1 - 0.5
  z2 <- u2 - 0.5
  z3 <- u3 - 0.5
  z4 <- u4 - 0.5

  0.65 * sin(2 * pi * u1) +
    0.55 * cos(2 * pi * u2) +
    0.45 * sin(2 * pi * (u3 + 0.25 * u4)) +
    0.65 * z2^2 - 0.45 * z4^2 +
    1.15 * z1 * z2 -
    0.95 * z2 * z3 +
    1.70 * z1 * z3 * z4 +
    1.25 * exp(-55 * ((u1 - 0.75)^2 + (u4 - 0.25)^2)) -
    1.10 * exp(-45 * ((u2 - 0.25)^2 + (u3 - 0.70)^2)) +
    0.18 * abs(u3 - 0.55) +
    0.28 * pmax(u1 - 0.72, 0) * pmax(0.35 - u2, 0) +
    0.16 * pmax(u3 + u4 - 1.15, 0)
}

## Faceted DGP figure.
g <- seq(0, 1, length.out = 180)
g4 <- seq(0, 1, length.out = 150)

slices <- list(
  "q = 2: t1 by t2, unit scale" =
    outer(g, g, Vectorize(f_q2)),
  "q = 4: t1 by t2; t3=.25,t4=.25" =
    outer(g4, g4, Vectorize(function(u1, u2) f_q4(u1, u2, .25, .25))),
  "q = 4: t1 by t2; t3=.75,t4=.75" =
    outer(g4, g4, Vectorize(function(u1, u2) f_q4(u1, u2, .75, .75))),
  "q = 4: t2 by t3; t1=.25,t4=.75" =
    outer(g4, g4, Vectorize(function(u2, u3) f_q4(.25, u2, u3, .75))),
  "q = 4: t2 by t3; t1=.75,t4=.25" =
    outer(g4, g4, Vectorize(function(u2, u3) f_q4(.75, u2, u3, .25))),
  "q = 4: t1 by t4; t2=.25,t3=.75" =
    outer(g4, g4, Vectorize(function(u1, u4) f_q4(u1, .25, .75, u4))),
  "q = 4: t1 by t4; t2=.75,t3=.25" =
    outer(g4, g4, Vectorize(function(u1, u4) f_q4(u1, .75, .25, u4)))
)

png("dgp_surfaces.png", width = 1500, height = 2300,
    res = 150)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 0), nrow = 4, byrow = TRUE))
par(mar = c(4, 4, 3, 1))
for (i in seq_along(slices)) {
  gg <- if (i == 1) g else g4
  image(gg, gg, slices[[i]], col = cols, xlab = "", ylab = "",
        main = names(slices)[i])
  contour(gg, gg, slices[[i]], add = TRUE, drawlabels = FALSE,
          nlevels = 12, col = rgb(0, 0, 0, 0.45), lwd = 0.8)
}
dev.off()
