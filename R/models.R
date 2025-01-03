# This file sets the following global variables:
#     spl: a list governing which spline bases to apply to the dataset
#     b9, b10, etc: spline bases
#     par_init: initial parameter values (for optimization)
#     par_y: a character vector naming the parameters corresponding to the outcome model; used by f_y() and prob()
#     terms_y: data values corresponding to par_y; used by f_y()
#     terms_y2: a function that generates the "data values" (for predicted probabilities); used by prob()
#     par_x, terms_x, etc.: analogous but for the seroconversion discrete hazard model
#     par_s, terms_s, etc.: analogous but for the initial status model

# Construct spline bases
if (cfg$model_version==7) {
  # (none)
} else if (cfg$model_version %in% c(29:33,36)) {
  spl <- list(list(name="b9", var="w_1", df=4),
              list(name="b10", var="t_end", df=4))
  b9 <- construct_basis("age (13,20,30,40,60)")
  b10 <- construct_basis("year (10,13,16,19,22)")
} else if (cfg$model_version %in% c(34:35)) {
  b9 <- construct_basis("age (13,20,30,40,60)")
  b12 <- construct_basis("year (17,...,22)")
} else if (cfg$model_version==37) {
  b13 <- construct_basis("age (13,20,30,40,60)", linear=T)
  b14 <- construct_basis("year (10,13,16,19,22)", linear=T)
} else if (cfg$model_version %in% c(38:42)) {
  spl <- list(list(name="b13", var="w_1", df=4),
              list(name="b14", var="t_end", df=4),
              list(name="b15", var="w_1", df=3))
  b13 <- construct_basis("age (13,20,30,40,60)", linear=T)
  b14 <- construct_basis("year (10,13,16,19,22)", linear=T)
  b15 <- construct_basis("age (13,30,40,60)", linear=T)
}

# Set parameter initial values
if (cfg$model_version==7) {
  par_init <- c(a_x=-6.2967, g_x1=-0.1535, g_x2=0.9796, t_x=0.5343, a_s=-2.3111, g_s1=-0.5649, g_s2=0.6198, t_s=0.4245, beta_x=1.401, a_y=-5.5786, g_y1=0.3278, g_y2=4.2046, t_y=-0.7198)
} else if (cfg$model_version==29) {
  par_init <- c(a_x=-6.087, g_x1=1.283, g_x2=-0.983, g_x3=-1.032, g_x4=-0.888, t_x1=0.602, t_x2=-1.874, t_x3=-1.472, t_x4=-0.530, a_s=-3.209, g_s1=3.658, g_s2=2.351, g_s3=4.186, g_s4=0.942, beta_x1=1.981, beta_x2=-0.784, beta_x3=0.0149, beta_x4=-1.884, beta_x5=-2.042, beta_x6=-0.565, a_y=-8.153, g_y1=2.255, g_y2=1.771, g_y3=4.729, g_y4=2.989, t_y1=-0.174, t_y2=-0.103, t_y3=-0.288, t_y4=0.126)
} else if (cfg$model_version==30) {
  par_init <- c(a_x=-6.087, g_x1=1.283, g_x2=-0.983, g_x3=-1.032, g_x4=-0.888, t_x1=0.602, t_x2=-1.874, t_x3=-1.472, t_x4=-0.530, a_s=-3.209, g_s1=3.658, g_s2=2.351, g_s3=4.186, g_s4=0.942, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.153, g_y1=2.255, g_y2=1.771, g_y3=4.729, g_y4=2.989, t_y1=-0.174, t_y2=-0.103, t_y3=-0.288, t_y4=0.126)
} else if (cfg$model_version==31) {
  par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=1.632, beta_x2=-0.406, beta_x3=-0.511, beta_x4=-0.525, beta_x5=-1.337, beta_x6=-0.472, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
} else if (cfg$model_version==32) {
  par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
} else if (cfg$model_version==33) {
  par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
} else if (cfg$model_version==34) {
  par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
} else if (cfg$model_version==35) {
  par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=0, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
} else if (cfg$model_version==36) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.6, g_x1=0.69, g_x2=-1.27, g_x3=1.04, g_x4=-2.38, t_x1=1.21, t_x2=-1.64, t_x3=-1.04, t_x4=-0.51, a_s=-3.18, g_s1=3.51, g_s2=2.13, g_s3=4.35, g_s4=0.57, beta_x1=2.11, beta_x2=-0.03, beta_x3=-0.39, beta_x4=-0.15, a_y=-8.02, g_y1=2.13, g_y2=1.66, g_y3=4.37, g_y4=2.8, t_y1=-0.17, t_y2=-0.17, t_y3=-0.37, t_y4=0.06) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.61, g_x1=2.41, g_x2=-1.89, g_x3=-0.42, g_x4=-1.73, t_x1=0.48, t_x2=-3.19, t_x3=-0.52, t_x4=-0.43, a_s=-3.75, g_s1=3.84, g_s2=3.03, g_s3=4.03, g_s4=2.09, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.15, beta_x4=-0.17, a_y=-7.18, g_y1=1.91, g_y2=1.9, g_y3=3.95, g_y4=2.71, t_y1=-0.27, t_y2=-0.03, t_y3=-0.12, t_y4=0.02) }
} else if (cfg$model_version==37) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.6, g_x1=25, g_x2=-32, g_x3=-5, g_x4=7, t_x1=0.02, t_x2=0.2, t_x3=-0.8, t_x4=0.7, a_s=-3.18, g_s1=30, g_s2=-17, g_s3=-18, g_s4=0, beta_x1=2.11, beta_x2=-0.03, beta_x3=-0.39, beta_x4=-0.15, a_y=-8.02, g_y1=20, g_y2=-13, g_y3=-5, g_y4=4, t_y1=-0.07, t_y2=0.075, t_y3=-0.01, t_y4=0.04) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.61, g_x1=13, g_x2=-3, g_x3=-30, g_x4=13, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.3, a_s=-3.75, g_s1=15, g_s2=10, g_s3=-25, g_s4=-4, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.15, beta_x4=-0.17, a_y=-7.18, g_y1=15, g_y2=-7, g_y3=-4, g_y4=1, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.02) }
} else if (cfg$model_version==38) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.47, g_x1=10.01, g_x2=-27.66, g_x3=13.78, t_x1=0.01, t_x2=0.2, t_x3=-0.81, t_x4=0.7, a_s=-3.17, g_s1=30.03, g_s2=-16.89, g_s3=-17.81, g_s4=0.04, beta_x1=2.1, beta_x2=-0.03, beta_x3=-0.4, beta_x4=-0.15, a_y=-8.03, g_y1=19.99, g_y2=-12.99, g_y3=-4.94, g_y4=4.24, t_y1=-0.07, t_y2=0.08, t_y3=-0.02, t_y4=0.03) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.59, g_x1=11.58, g_x2=-30.05, g_x3=13.13, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.32, a_s=-3.75, g_s1=14.95, g_s2=9.95, g_s3=-24.96, g_s4=-3.94, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.18, beta_x4=-0.18, a_y=-7.18, g_y1=15, g_y2=-7.04, g_y3=-4.05, g_y4=1.23, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.01) }
} else if (cfg$model_version==39) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.47, g_x1=25, g_x2=-32, g_x3=-5, g_x4=7, t_x1=0.01, t_x2=0.2, t_x3=-0.81, t_x4=0.7, a_s=-3.17, g_s1=30.03, g_s2=-16.89, g_s3=-17.81, g_s4=0.04, beta_x1=2.1, beta_x2=-0.03, beta_x3=-0.4, beta_x4=-0.15, a_y=-8.03, g_y1=19.99, g_y2=-12.99, g_y3=-4.94, g_y4=4.24, t_y1=-0.07, t_y2=0.08, t_y3=-0.02, t_y4=0.03) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.59, g_x1=11.58, g_x2=-30.05, g_x3=13.13, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.32, a_s=-3.75, g_s1=14.95, g_s2=9.95, g_s3=-24.96, g_s4=-3.94, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.18, beta_x4=-0.18, a_y=-7.18, g_y1=15, g_y2=-7.04, g_y3=-4.05, g_y4=1.23, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.01) }
} else if (cfg$model_version==40) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.29, g_x1=26.2, g_x2=-32.84, g_x3=-10.42, g_x4=14.75, t_x1=-0.12, t_x2=0.27, t_x3=-0.61, t_x4=0.27, a_s=-3.21, g_s1=30.25, g_s2=-16.79, g_s3=-17.27, g_s4=-1.94, t_s1=0, beta_x1=2.17, beta_x2=-0.02, beta_x3=-0.88, beta_x4=-0.14, a_y=-7.98, g_y1=20.33, g_y2=-13.34, g_y3=-5.15, g_y4=4.74, t_y1=-0.13, t_y2=0.2, t_y3=-0.14, t_y4=0.11) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.55, g_x1=12, g_x2=-32.08, g_x3=13.73, t_x1=0.19, t_x2=-0.4, t_x3=-0.5, t_x4=1.29, a_s=-3.75, g_s1=14.87, g_s2=9.92, g_s3=-24.81, g_s4=-3.88, t_s1=0, beta_x1=2.23, beta_x2=-0.03, beta_x3=-1.21, beta_x4=-0.17, a_y=-7.16, g_y1=14.95, g_y2=-7.09, g_y3=-4.13, g_y4=1.51, t_y1=-0.05, t_y2=0.01, t_y3=0.07, t_y4=-0.02) }
} else if (cfg$model_version==41) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.29, g_x1=26.2, g_x2=-32.84, g_x3=-10.42, g_x4=14.75, t_x1=-0.12, t_x2=0.27, t_x3=-0.61, t_x4=0.27, a_s=-3.21, g_s1=30.25, g_s2=-16.79, g_s3=-17.27, g_s4=-1.94, t_s1=0, beta_x1=2.17, beta_x2=-0.02, beta_x3=0, beta_x4=-0.88, beta_x5=-0.14, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-7.98, g_y1=20.33, g_y2=-13.34, g_y3=-5.15, g_y4=4.74, t_y1=-0.13, t_y2=0.2, t_y3=-0.14, t_y4=0.11) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.55, g_x1=12, g_x2=-32.08, g_x3=13.73, t_x1=0.19, t_x2=-0.4, t_x3=-0.5, t_x4=1.29, a_s=-3.75, g_s1=14.87, g_s2=9.92, g_s3=-24.81, g_s4=-3.88, t_s1=0, beta_x1=2.23, beta_x2=-0.03, beta_x3=0, beta_x4=-1.21, beta_x5=-0.17, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-7.16, g_y1=14.95, g_y2=-7.09, g_y3=-4.13, g_y4=1.51, t_y1=-0.05, t_y2=0.01, t_y3=0.07, t_y4=-0.02) }
} else if (cfg$model_version==42) {
  if (cfg$model_sex=="Female") { par_init <- c(a_x=-5.25, g_x1=26.43, g_x2=-32.67, g_x3=-9.94, g_x4=15.47, g_x5=0, t_x1=-0.12, t_x2=0.24, t_x3=-0.61, t_x4=0.27, a_s=-3.2, g_s1=30.27, g_s2=-16.75, g_s3=-17.2, g_s4=-1.81, g_s5=0, t_s1=-0.01, beta_x1=2.18, beta_x2=-0.02, beta_x3=-0.89, beta_x4=-0.15, a_y=-7.99, g_y1=20.31, g_y2=-13.3, g_y3=-5.1, g_y4=5, g_y5=0, t_y1=-0.13, t_y2=0.19, t_y3=-0.12, t_y4=0.11) }
  if (cfg$model_sex=="Male") { par_init <- c(a_x=-6.49, g_x1=10.22, g_x2=-30.32, g_x3=17.79, g_x4=0, t_x1=0.22, t_x2=-0.33, t_x3=-0.57, t_x4=1.03, a_s=-3.69, g_s1=14.91, g_s2=9.88, g_s3=-24.08, g_s4=-5.29, g_s5=0, t_s1=-0.03, beta_x1=2.29, beta_x2=-0.03, beta_x3=-1.43, beta_x4=-0.15, a_y=-7.1, g_y1=15, g_y2=-7.17, g_y3=-5.46, g_y4=3.79, g_y5=0, t_y1=-0.04, t_y2=-0.01, t_y3=0.09, t_y4=-0.02) }
}

# Outcome model
if (cfg$model_version==7) {
  par_y_F <- c("a_y", "t_y", "g_y1", "g_y2", "beta_x")
  terms_y_F <- function(r, x) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]], x) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(1, j, w_1, w_3, x) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==29) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*max(r[["w_1"]]-0.3,0), x*max(r[["w_1"]]-0.3,0)*r[["j"]], x*max(r[["w_1"]]-0.45,0), x*max(r[["w_1"]]-0.45,0)*r[["j"]], 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*max(w_1-0.3,0), x*j*max(w_1-0.3,0), x*max(w_1-0.45,0), x*j*max(w_1-0.45,0), 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b10(j,1), b10(j,2), b10(j,3), b10(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version %in% c(30,36)) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["w_1"]], x*r[["j"]]*r[["w_1"]], 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*w_1, x*j*w_1, 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b10(j,1), b10(j,2), b10(j,3), b10(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==31) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*max(r[["w_1"]]-0.2,0), x*max(r[["w_1"]]-0.2,0)*r[["j"]], x*min(max(r[["w_1"]]-0.4,0), 0.5), x*min(max(r[["w_1"]]-0.4,0), 0.5)*r[["j"]], 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*max(w_1-0.2,0), x*j*max(w_1-0.2,0), x*min(max(w_1-0.4,0),0.5), x*j*min(max(w_1-0.4,0),0.5), 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b10(j,1), b10(j,2), b10(j,3), b10(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==32) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7", "beta_x8", "beta_x9", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*max(r[["j"]]-0.6,0), x*r[["w_1"]], x*r[["w_1"]]*r[["j"]], x*r[["w_1"]]*max(r[["j"]]-0.6,0), x*max(r[["w_1"]]-0.4,0), x*max(r[["w_1"]]-0.4,0)*r[["j"]], x*max(r[["w_1"]]-0.4,0)*max(r[["j"]]-0.6,0), 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x*1, x*j, x*max(j-0.6,0), x*w_1, x*w_1*j, x*w_1*max(j-0.6,0), x*max(w_1-0.4,0), x*max(w_1-0.4,0)*j, x*max(w_1-0.4,0)*max(j-0.6,0), 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b10(j,1), b10(j,2), b10(j,3), b10(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==33) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7", "beta_x8", "beta_x9", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["j"]]^2, x*r[["w_1"]], x*r[["w_1"]]*r[["j"]], x*r[["w_1"]]*r[["j"]]^2, x*r[["w_1"]]^2, x*r[["w_1"]]^2*r[["j"]], x*r[["w_1"]]^2*r[["j"]]^2, 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x*1, x*j, x*j^2, x*w_1, x*w_1*j, x*w_1*j^2, x*w_1^2, x*w_1^2*j, x*w_1^2*j^2, 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b10(j,1), b10(j,2), b10(j,3), b10(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version %in% c(34:35)) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["w_1"]], x*r[["j"]]*r[["w_1"]], 1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]], r[["b12_1"]], r[["b12_2"]], r[["b12_3"]], r[["b12_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*w_1, x*j*w_1, 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4), b12(j,1), b12(j,2), b12(j,3), b12(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version %in% c(37:40)) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["w_1"]], x*r[["j"]]*r[["w_1"]], 1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*w_1, x*j*w_1, 1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), b14(j,1), b14(j,2), b14(j,3), b14(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==41) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6", "beta_x7", "beta_x8", "beta_x9", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["j"]]^2, x*r[["w_1"]], x*r[["j"]]*r[["w_1"]], x*r[["j"]]^2*r[["w_1"]], x*r[["w_1"]]^2, x*r[["j"]]*r[["w_1"]]^2, x*r[["j"]]^2*r[["w_1"]]^2, 1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*j^2, x*w_1, x*j*w_1, x*j^2*w_1, x*w_1^2, x*j*w_1^2, x*j^2*w_1^2, 1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), b14(j,1), b14(j,2), b14(j,3), b14(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
} else if (cfg$model_version==42) {
  par_y_F <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3", "t_y4")
  terms_y_F <- function(r, x) { c(x, x*r[["j"]], x*r[["w_1"]], x*r[["j"]]*r[["w_1"]], 1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["w_2"]], r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]]) }
  terms_y2_F <- function(x, j, w_1, w_2) { c(x, x*j, x*w_1, x*j*w_1, 1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), w_2, b14(j,1), b14(j,2), b14(j,3), b14(j,4)) }
  par_y_M <- par_y_F; terms_y_M <- terms_y_F; terms_y2_M <- terms_y2_F;
}

# Seroconversion model
if (cfg$model_version==7) {
  par_x_F <- c("a_x", "t_x1", "g_x1", "g_x2")
  terms_x_F <- function(r) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, j, w_1, w_3) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version %in% c(29:33,36)) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["b10_1"]], r[["b10_2"]], r[["b10_3"]], r[["b10_4"]], r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b10(j,1), b10(j,2), b10(j,3), b10(j,4), b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4)) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version==34) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["b12_1"]], r[["b12_2"]], r[["b12_3"]], r[["b12_4"]], r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b12(j,1), b12(j,2), b12(j,3), b12(j,4), b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4)) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version==35) {
  par_x_F <- c("a_x", "t_x1", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["j"]], r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, j, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4)) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version==37) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4)) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version==38) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3")
  terms_x_F <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b15_1"]], r[["b15_2"]], r[["b15_3"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b15(w_1,1), b15(w_1,2), b15(w_1,3)) }
  par_x_M <- par_x_F; terms_x_M <- terms_x_F; terms_x2_M <- terms_x2_F;
} else if (cfg$model_version %in% c(39:41)) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_F <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4)) }
  par_x_M <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3")
  terms_x_M <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b15_1"]], r[["b15_2"]], r[["b15_3"]]) }
  terms_x2_M <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b15(w_1,1), b15(w_1,2), b15(w_1,3)) }
} else if (cfg$model_version==42) {
  par_x_F <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4", "g_x5")
  terms_x_F <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["w_2"]]) }
  terms_x2_F <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), w_2) }
  par_x_M <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3", "g_x4")
  terms_x_M <- function(r) { c(1, r[["b14_1"]], r[["b14_2"]], r[["b14_3"]], r[["b14_4"]], r[["b15_1"]], r[["b15_2"]], r[["b15_3"]], r[["w_2"]]) }
  terms_x2_M <- function(j, w_1, w_2) { c(1, b14(j,1), b14(j,2), b14(j,3), b14(j,4), b15(w_1,1), b15(w_1,2), b15(w_1,3), w_2) }
}

# Initial status model
if (cfg$model_version==7) {
  par_s_F <- c("a_s", "t_s1", "g_s1", "g_s2")
  terms_s_F <- function(r) { c(1, r[["j"]], r[["w_1"]], r[["w_3"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, j, w_1, w_3) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
} else if (cfg$model_version %in% c(29:36)) {
  par_s_F <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4")
  terms_s_F <- function(r) { c(1, r[["b9_1"]], r[["b9_2"]], r[["b9_3"]], r[["b9_4"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4)) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
} else if (cfg$model_version %in% c(37:39)) {
  par_s_F <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4")
  terms_s_F <- function(r) { c(1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4)) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
} else if (cfg$model_version %in% c(40:41)) {
  par_s_F <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4", "t_s1")
  terms_s_F <- function(r) { c(1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["j"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), j) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
} else if (cfg$model_version==42) {
  par_s_F <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4", "g_s5", "t_s1")
  terms_s_F <- function(r) { c(1, r[["b13_1"]], r[["b13_2"]], r[["b13_3"]], r[["b13_4"]], r[["w_2"]], r[["j"]]) }
  terms_s2_F <- function(j, w_1, w_2) { c(1, b13(w_1,1), b13(w_1,2), b13(w_1,3), b13(w_1,4), w_2, j) }
  par_s_M <- par_s_F; terms_s_M <- terms_s_F; terms_s2_M <- terms_s2_F;
}

# Helper code to construct par_init(...) statement
# Construct param init vector
if (F) {
  ests <- readRDS("objs/ests_40_full_F_20250102.rds")
  str <- "par_init <- c("
  for (par in names(ests$opt$par)) {
    str <- paste0(str, par, "=", round(ests$opt$par[[par]], 2), ", ")
  }
  str <- paste0(substr(str, 1, nchar(str)-2), ")")
  print(str)
}
