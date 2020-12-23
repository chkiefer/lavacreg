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
                  family = 'nbinom')
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
                  family = 'nbinom')
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
                  family = 'nbinom',
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
                  se = TRUE)
  
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 30)
  
  # Correct parameter estimates?
  comp <- c(6.05971523, 2.34246617, -0.08740382, -0.0788266, 0.08191891, 
            1.36239451, 1.59239172, 1.59773197, 1.68488975, 3.9090709, 
            0.95099035, 1.02372983, -0.63568355, -0.50579888, 0, 6.09270393, 
            2.62186177, -0.0928436, -0.04702289, 0.05436236, 1.4068395, 1.47774334, 
            1.56318072, 1.36078961, 3.99384264, 0.83179534, 0.86117941, -0.56710374, 
            -0.46501115, 0)
  
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group three manifest covariates negative binomial", {
  fit <- countreg(forml='dv ~ z12 + z11 + z21',
                  group = "treat",
                  data = example01,
                  family = 'nbinom',
                  se = FALSE)
  # Converged?
  conv <- fit@fit$fit$convergence
  expect_equal(conv, 0)
  
  # Right amount of parameters?
  pt <- fit@fit$pt
  par <- pt$par
  expect_equal(length(par), 30)
  
  # Correct parameter estimates?
  comp <- c(6.05820964, 2.33641826, -0.08667106, -0.07961202, 0.08353954, 1.3629161, 
            1.58195313, 1.59863642, 1.68118476, 3.90915155, 0.94687096, 1.01833197, 
            -0.62852949, -0.50192426, 15.58364617, 6.09436922, 2.63362154, -0.09385172, 
            -0.04683473, 0.05173418, 1.40531856, 1.48049197, 1.56203441, 1.36630404, 3.994657, 
            0.83546495, 0.86506982, -0.57100751, -0.46916909, 17.72343536)
  
  # Currently, the results differ slightly between 32 and 64bit systems
  # A higher tolerance would be desirable, but cannot be achieved at the moment
  expect_equal(par, comp, tolerance = 1e-2)
  
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
  comp <- c(6.05912482, 0.4235168, 0.50097813, 0, 1, 1.45579428, 0.72453059, 
            3.9347506, 0.42455089, 0, 0.52693647, 0.35698009, 6.09356923, 
            0.81433536, 0.45436361, 0, 1, 1.45579428, 0.72453059, 3.97163325, 
            0.32495448, 0, 0.50982612, 0.35643488)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one latent covariate negative binomial", {
  fit <- countreg(forml='dv ~ eta1',
                  lv = list(eta1 = c("z21", "z22")),
                  group = 'treat',
                  data = example01,
                  family = 'nbinom',
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
  comp <- c(6.05912318, 1.35949557, 0.2714377, 0, 1, -0.88976068, 
            1.31764183, 3.8767443, 0.43434811, 11.15574397, 0.66700314, 
            0.03291596, 6.09356999, 1.46169048, 0.29855387, 0, 1, -0.88976068, 
            1.31764183, 3.95304398, 0.22302417, 12.90964003, 0.5821887, 0.14893486)
  expect_equal(par, comp, tolerance = 1e-2)
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
  comp <- c(6.05911601, 0.63021608, 0.46370007, -0.03725704, 0, 1, 1.28783918, 
            0.76712164, 0, 0, 0, 0, 0, 0, 1, -0.10893721, 1.29106256, -0.47703288, 
            1.37432815, 3.94215875, 0.42635444, 1.57320128, 1.8619543, -0.34283515, 
            0, 0.52696251, 0.33559719, 1.5334956, 1.43583402, 1.39811958, 6.09358191, 
            0.48523348, 0.52380971, 0.03182268, 0, 1, 1.28783918, 0.76712164, 0, 0, 0, 
            0, 0, 0, 1, -0.10893721, 1.29106256, -0.47703288, 1.37432815, 3.97530262, 
            0.30173802, 1.65009244, 2.12180181, -0.40660376, 0, 0.53166143, 0.35623691, 
            1.38672523, 1.57462925, 1.0176448)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group two latent covariates negative binomial", {
  skip("This test can take up to 10 minutes.")
  fit <- countreg(forml='dv ~ eta1 + eta2',
                  lv = list(eta1 = c("z21", "z22"),
                            eta2 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'nbinom',
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
  comp <- c(6.05909654, 2.75913773, -0.1416134, -0.09370237, 0, 1, -0.06870994, 1.26866175, 
            -0.40189797, 1.33204864, 1.58256123, 1.90225746, 1.35599432, 1.58854965, 0.50967322, 
            0, 1.51401949, 1.46300016, 1.47144681, 6.09358806, 2.86544299, -0.12885743, -0.02848884, 
            0, 1, -0.06870994, 1.26866175, -0.40189797, 1.33204864, 1.63915325, 2.18476374, 
            1.39529371, 1.49758501, 0.69848563, 0, 1.35903901, 1.51394989, 1.10139006)
  expect_equal(par, comp, tolerance = 1e-5)
})


test_that("two-group one latent, one manifest covariate negative binomial", {
  fit <- countreg(forml='dv ~ eta1 + z12',
                  lv = list(eta1 = c("z41", "z42", "z43")),
                  group = 'treat',
                  data = example01,
                  family = 'nbinom',
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
  expect_equal(par, comp, tolerance = 1e-2)
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
  comp <- c(6.059124, 2.36817557, -0.10841153, 0.08546017, -0.08748471, 0, 1, -0.06092308, 
            1.26400824, -0.41092405, 1.3371816, 1.58195192, 1.89942736, 1.36298464, 
            1.59047871, 3.90870323, 0.95192072, -0.63512512, 0.50475972, -0.26781445, 0, 
            1.51100331, 1.47559415, 1.44904759, 6.09356947, 2.57935422, -0.10665812, 
            0.06166142, -0.02401982, 0, 1, -0.06092308, 1.26400824, -0.41092405, 1.3371816, 
            1.64191737, 2.21649531, 1.39027783, 1.50211821, 4.00297221, 0.83856677, -0.5807701, 
            0.7111785, -0.40031997, 0, 1.34647959, 1.53750437, 1.06622238)
  
  # Currently, the results differ slightly between 32 and 64bit systems
  # A higher tolerance would be desirable, but cannot be achieved at the moment
  expect_equal(par, comp, tolerance = 1e-3)
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


# ---------------------------------------------------
# TEST 10 - one latent variable - one group
# ---------------------------------------------------
test_that("one latent variable in one group - Poisson", {
  fit <- countreg(forml='dv ~ eta',
    lv = list(eta = c("z41", "z42", "z43")),
    group = NULL,
    data = example01,
    family = 'poisson')
  par <- fit@fit$pt$par
  comp <- c(6.769642, 2.72414, -0.109569, 0, 1, -0.054137, 
            1.256452, -0.343201, 1.29337, 1.61712, 2.004875, 0, 
            1.393845, 1.523669, 1.430018)
  expect_equal(length(par), 15)
  expect_equal(par, comp, tolerance = 1e-5)
})

