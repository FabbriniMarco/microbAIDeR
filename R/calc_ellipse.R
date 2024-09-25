calc_ellipse <- function(coord_beta, level = 0.99) {
  n <- nrow(coord_beta)
  mean_x <- mean(coord_beta$PC1)
  mean_y <- mean(coord_beta$PC2)
  cov_matrix <- cov(coord_beta[, c("PC1", "PC2")])
  eig_vals <- eigen(cov_matrix)$values
  eig_vecs <- eigen(cov_matrix)$vectors
  angle <- atan2(eig_vecs[2, 1], eig_vecs[1, 1])
  chisq_val <- qchisq(level, df = 2)
  scale_factor <- sqrt(chisq_val / n)
  radii <- sqrt(eig_vals) * scale_factor
  theta <- seq(0, 2 * pi, length = 100)
  ellipse_x <- mean_x + radii[1] * cos(theta) * cos(angle) - radii[2] * sin(theta) * sin(angle)
  ellipse_y <- mean_y + radii[1] * cos(theta) * sin(angle) + radii[2] * sin(theta) * cos(angle)
  data.frame(PC1 = ellipse_x, PC2 = ellipse_y, Group = unique(coord_beta$Group))
}