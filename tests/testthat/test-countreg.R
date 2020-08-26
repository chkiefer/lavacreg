# ---------------------------------------------------
# TEST 1 - intercept only - one group
# ---------------------------------------------------
test_that("intercept-only Poisson", {
  fit <- countreg(forml='dv ~ 1',
                  group = NULL,
                  data = example01,
                  family = 'poisson')
  par <- fit@fit$pt$par
  comp <- c(6.769661, 2.555460, 0.000000)
  expect_equal(length(par), 3)
  expect_equal(par, comp, tolerance = 1e-5)
})

test_that("intercept-only negative binomial", {
  fit <- countreg(forml='dv ~ 1',
                  group = NULL,
                  data = example01,
                  family = 'negbin')
  par <- fit@fit$pt$par
  comp <- c(6.769614, 2.555453, 8.626471)
  expect_equal(length(par), 3)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 2 - intercept-only - two groups
# ---------------------------------------------------
test_that("two-group intercept-only Poisson", {
  fit <- countreg(forml='dv ~ 1',
                  group = "treat",
                  data = example01,
                  family = 'poisson')
  pt <- fit@fit$pt
  par <- pt$par
  comp <- c(6.059399, 2.445395, 0.000000, 6.093341, 2.651439, 0.000000)
  expect_equal(length(par), 6)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group intercept-only negative binomial", {
  fit <- countreg(forml='dv ~ 1',
                  group = "treat",
                  data = example01,
                  family = 'negbin')
  pt <- fit@fit$pt
  par <- pt$par
  comp <- c(6.060023, 2.445382, 9.247623, 6.092560, 2.651660, 9.600843)
  expect_equal(length(par), 6)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 3 - one manifest covariate - two groups
# ---------------------------------------------------
test_that("two-group one manifest covariate Poisson", {
  fit <- countreg(forml='dv ~ z12',
                  group = "treat",
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  pt <- fit@fit$pt
  par <- pt$par
  comp <- c(6.0589305,  2.6549102, -0.1697253,  1.3644560,  1.5933374,  0.0000000,  
            6.0943167,  2.8358178, -0.1414390,  1.4037590,  1.4801428,  0.0000000)
  expect_equal(length(par), 12)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one manifest covariate negative binomial", {
  fit <- countreg(forml='dv ~ z12',
                  group = "treat",
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  pt <- fit@fit$pt
  par <- pt$par
  comp <- c(6.0592596,  2.6546687, -0.1693490,  1.3646616,  1.5933417, 14.8400276,  
            6.0935426,  2.8361803, -0.1416315,  1.4056919,  1.4803146, 14.7545861)
  expect_equal(length(par), 12)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 4 - three manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group three manifest covariates Poisson", {
  fit <- countreg(forml='dv ~ z12 + z11 + z21',
                  group = "treat",
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 30)
  
  # Correct parameter estimates?
  comp <- c(6.05973365,  2.34250754, -0.08741083, -0.07882434,  0.08191071, 
            1.36238200,  1.59241136,  1.59772962,  1.68492870,  3.90908989,  
            0.95103335,  1.02370271, -0.63570761, -0.50580146,  0.00000000, 
            6.09268439,  2.62190851, -0.09284062, -0.04703224,  0.05435265,  
            1.40683041,  1.47773341,  1.56318303,  1.36082238,  3.99387272,
            0.83180403,  0.86118618, -0.56709844, -0.46496553,  0.00000000)
  
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group three manifest covariates negative binomial", {
  fit <- countreg(forml='dv ~ z12 + z11 + z21',
                  group = "treat",
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 30)
  
  # Correct parameter estimates?
  comp <- c(6.05979618,  2.33778266, -0.08674098, -0.07949502,  0.08308643,  
            1.36224716,  1.57193783,  1.59789694,  1.67196204,  3.90999818,  
            0.94587503,  1.00702913, -0.62502894, -0.49742089, 15.70498430,  
            6.09337544,  2.63022367, -0.09373444, -0.04660160,  0.05242476,  
            1.40506932,  1.48904095,  1.56190589,  1.36983050,  3.99515666,  
            0.83702507,  0.87047641, -0.57452029, -0.47118257, 17.83688410)
  expect_equal(par, comp, tolerance = 1e-5)
  
})


# ---------------------------------------------------
# TEST 5 - one latent covariate - two groups
# ---------------------------------------------------
test_that("two-group one latent covariate Poisson", {
  fit <- countreg(forml='dv ~ eta1',
                  lv = list(eta1 = c("z21", "z22")),
                  group = 'treat',
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0L)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 24)
  
  # Correct parameter estimates?
  comp <- c(6.0591237, 0.4236599, 0.5009449, 0.0000000, 1.0000000, 1.4556876, 
            0.7245659, 3.9347108, 0.4246533, 0.0000000, 0.5268721, 0.3569312, 
            6.0935697, 0.8144386, 0.4543410, 0.0000000, 1.0000000, 1.4556876, 
            0.7245659, 3.9715882, 0.3250218, 0.0000000, 0.5097923, 0.3563857)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one latent covariate negative binomial", {
  fit <- countreg(forml='dv ~ eta1',
                  lv = list(eta1 = c("z21", "z22")),
                  group = 'treat',
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 24)
  
  # Correct parameter estimates?
  # CAUTION: ABWEICHUNGEN VON MPLUS, VORALLEM MESSMODELL (BEI DREI INDIKATOREN WENIGER STARK, ABER VORHANDEN)
  comp <- c(6.05912320,  1.35864961,  0.27166173,  0.00000000,  1.00000000, -0.89108964,  
            1.31800757,  3.87668969,  0.43411584, 11.15025565,  0.66704596,  0.03291918, 
            6.09356990,  1.46223253,  0.29842322,  0.00000000,  1.00000000, -0.89108964,  
            1.31800757,  3.95295807,  0.22299303, 12.90454215,  0.58219365,  0.14880052)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 6 - two latent covariates - two groups
# ---------------------------------------------------
test_that("two-group two latent covariates Poisson", {
  print("This test can take up to 6 minutes.")
  fit <- countreg(forml='dv ~ eta1 + eta2',
                  lv = list(eta1 = c("z21", "z22"),
                            eta2 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0L)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 60)
  
  # Correct parameter estimates?
  comp <- c(6.05908773,  0.62897099,  0.46384518, -0.03685026,  0.00000000,
            1.00000000,  1.32681665,  0.75733638,  0.00000000,  0.00000000,  
            0.00000000,  0.00000000,  0.00000000,  0.00000000,  1.00000000, 
           -0.09736793,  1.28238655, -0.46033865,  1.36469021,  3.93975281,  
            0.42760931,  1.58285529,  1.88298793, -0.34278999,  0.00000000,  
            0.52624708,  0.33736358,  1.53538744,  1.45101893,  1.39077548,  
            6.09360142,  0.52812449,  0.51359919,  0.03015731,  0.00000000,  
            1.00000000,  1.32681665,  0.75733638,  0.00000000,  0.00000000,  
            0.00000000,  0.00000000,  0.00000000,  0.00000000,  1.00000000, 
           -0.09736793,  1.28238655, -0.46033865,  1.36469021,  3.97667345,  
            0.30845392,  1.64881638,  2.13851931, -0.40904144,  0.00000000,  
            0.52795840,  0.35534981,  1.38432994,  1.56758856,  1.00855859)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group two latent covariates negative binomial", {
  skip("This test can take up to 10 minutes.")
  fit <- countreg(forml='dv ~ eta1 + eta2',
                  lv = list(eta1 = c("z21", "z22"),
                            eta2 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 60)
  
  # Correct parameter estimates?
  comp <- c(6.05911345,  1.47195862,  0.26775418, -0.06742940,  0.00000000,
            1.00000000,  1.42094736,  0.73367468,  0.00000000,  0.00000000,  
            0.00000000,  0.00000000,  0.00000000,  0.00000000,  1.00000000, 
           -0.09546400,  1.28058375, -0.45606900,  1.36106299,  3.93994900,  
            0.54609831,  1.57854672,  1.89179740, -0.35106555, 15.39459637,  
            0.41592476,  0.28454688,  1.53720808,  1.45580357,  1.36947597,
            6.09357275,  1.67814046,  0.24971032, -0.01986893,  0.00000000,  
            1.00000000,  1.42094736,  0.73367468,  0.00000000,  0.00000000,  
            0.00000000,  0.00000000,  0.00000000,  0.00000000,  1.00000000, 
           -0.09546400,  1.28058375, -0.45606900,  1.36106299,  3.97884302,  
            0.43404231,  1.65126977,  2.16847937, -0.42078726, 17.00167035,  
            0.40983931,  0.30044684,  1.36352092,  1.57430547,  1.04101521)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 7 - one latent, one manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group one latent, one manifest covariate Poisson", {
  fit <- countreg(forml='dv ~ eta1 + z12',
                  lv = list(eta1 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0L)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 38)
  
  # Correct parameter estimates?
  comp <- c(6.05909655,  2.75913781, -0.14161341, -0.09370232,  0.00000000,  1.00000000, 
           -0.06871014,  1.26866202, -0.40189817,  1.33204893,  1.58256133,  1.90225741,  
            1.35599410,  1.58854970,  0.50967263,  0.00000000,  1.51402004,  1.46300115,  1.47144684,  
            6.09358806,  2.86544311, -0.12885742, -0.02848887,  0.00000000,  1.00000000, 
           -0.06871014,  1.26866202, -0.40189817,  1.33204893,  1.63915321,  2.18476306,  
            1.39529322,  1.49758488,  0.69848482,  0.00000000,  1.35903882,  1.51395049,  1.10138983)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one latent, one manifest covariate negative binomial", {
  skip("This test can take up to 10 minutes.")
  fit <- countreg(forml='dv ~ eta1 + z12',
                  lv = list(eta1 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 38)
  
  # Correct parameter estimates?
  # CAUTION: ABWEICHUNGEN VON MPLUS, VORALLEM MESSMODELL
  comp <- c(6.05913042,  2.74324019, -0.14916755, -0.07454461,  0.00000000,  1.00000000, 
           -0.11744052,  1.29223974, -0.85183243,  1.59771173,  1.60579135,  1.58193485,  
            1.36335131,  1.59838439,  0.45494084, 14.70721662,  1.70005958,  1.72693384,  0.94441234,  
            6.09355591,  2.86146567, -0.13330789, -0.02223673,  0.00000000,  1.00000000, 
           -0.11744052,  1.29223974, -0.85183243,  1.59771173,  1.67656884,  1.89033247,  
            1.40280054,  1.48848257,  0.60500366, 16.32127791,  1.68752564,  2.05397280,  0.26490092)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 8 - one latent, two manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group one latent, two manifest covariates Poisson", {
  fit <- countreg(forml='dv ~ eta1 + z12 + z21',
                  lv = list(eta1 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0L)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 48)
  
  # Correct parameter estimates?
  comp <- c(6.05912400,  2.36817180, -0.10841270,  0.08546189, -0.08748268,  0.00000000,  
            1.00000000, -0.06095325,  1.26402003, -0.41094544,  1.33719711,  1.58201068,  
            1.89939682,  1.36297738,  1.59049735,  3.90870908,  0.95194527, -0.63516411,  
            0.50476204, -0.26784343,  0.00000000,  1.51100435,  1.47560853,  1.44901943,  
            6.09356933,  2.57935324, -0.10665970,  0.06166213, -0.02401931,  0.00000000,  
            1.00000000, -0.06095325,  1.26402003, -0.41094544,  1.33719711,  1.64190657,  
            2.21643094,  1.39029640,  1.50214463,  4.00295743,  0.83859741, -0.58081911,  
            0.71121117, -0.40034811,  0.00000000,  1.34648320,  1.53754843,  1.06622232)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one latent, two manifest covariate negative binomial", {
  fit <- countreg(forml='dv ~ eta1 + z12 + z21',
                  lv = list(eta1 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 48)
  
  # Correct parameter estimates?
  # CAUTION: ABWEICHUNGEN VON MPLUS, VORALLEM MESSMODELL
  comp <- c(6.05912263,  2.34575133, -0.11303908,  0.08722950, -0.07244109,  0.00000000,  
            1.00000000, -0.10679834,  1.28808926, -0.48202704,  1.37623877,  1.58238357,  
            1.82771715,  1.35804073,  1.59891498,  3.91151517,  0.95254511, -0.63766547, 
            0.49038094, -0.25929977, 15.79591870,  1.55076248,  1.46841628,  1.36032350,  
            6.09357034,  2.59017300, -0.10894758,  0.05901176, -0.02210588,  0.00000000, 
            1.00000000, -0.10679834,  1.28808926, -0.48202704,  1.37623877,  1.64709335, 
            2.11067844,  1.39376522,  1.48983733,  4.00111840,  0.83615246, -0.57428941,  
            0.68313307, -0.38176448, 17.56235530,  1.39594155,  1.55199900,  1.00597486)
  expect_equal(par, comp, tolerance = 1e-5)
})


# ---------------------------------------------------
# TEST 9 - two latent, one manifest covariates - two groups
# ---------------------------------------------------
test_that("two-group two latent, one manifest covariates Poisson", {
  skip("This test can take up to 10 minutes.")
  fit <- countreg(forml='dv ~ eta1 + eta2 + z12',
                  lv = list(eta1 = c("z21", "z22"),
                            eta2 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'poisson',
                  se = FALSE)
  
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0L)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par[pt$par_free > 0L]
  expect_equal(length(par), 50)
  
  # Correct parameter estimates?
  comp <- c(6.05911879, -0.07540952,  0.07657528,  0.61043811, -0.02530174,  
            1.54644671,  0.70142013, -0.12121390,  1.29766040, -0.50652721,  
            1.39266284,  3.93871056,  0.41808751,  1.57716324,  1.84211943, 
           -0.35097437,  1.35724708,  1.60241726, -0.63245938,  0.50768859,  
            0.54249220,  0.37320559,  1.55987944,  1.43499172,  1.37814223,  
            6.09357869, -0.22505319,  0.09227664,  0.66831004,  0.03316800,  
            1.54644671,  0.70142013, -0.12121390,  1.29766040, -0.50652721,  
            1.39266284,  3.98300003,  0.31690630,  1.63526768,  2.08584996, 
           -0.41038142,  1.39058289,  1.49239446, -0.55768012,  0.67350735,  
            0.51904462,  0.37777490,  1.40109401,  1.59260762,  0.94147950)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group two latent, one manifest covariate negative binomial", {
  skip("This test can take up to 25 minutes.")
  fit <- countreg(forml='dv ~ eta1 + eta2 + z12',
                  lv = list(eta1 = c("z21", "z22"),
                            eta2 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'negbin',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par[pt$par_free > 0L]
  expect_equal(length(par), 52)
  
  # Correct parameter estimates?
  comp <- c(6.05912263,  2.34575133, -0.11303908,  0.08722950, -0.07244109,  0.00000000,  
            1.00000000, -0.10679834,  1.28808926, -0.48202704,  1.37623877,  1.58238357,  
            1.82771715,  1.35804073,  1.59891498,  3.91151517,  0.95254511, -0.63766547, 
            0.49038094, -0.25929977, 15.79591870,  1.55076248,  1.46841628,  1.36032350,  
            6.09357034,  2.59017300, -0.10894758,  0.05901176, -0.02210588,  0.00000000, 
            1.00000000, -0.10679834,  1.28808926, -0.48202704,  1.37623877,  1.64709335, 
            2.11067844,  1.39376522,  1.48983733,  4.00111840,  0.83615246, -0.57428941,  
            0.68313307, -0.38176448, 17.56235530,  1.39594155,  1.55199900,  1.00597486)
  expect_equal(par, comp, tolerance = 1e-5)
})




