d_bern <- list(N = 10, y = c(1, 1, 1, 0, 1, 1, 1, 0, 1, 0))

bernoulli_model = stan_model(file = 
                               here('./code/bernoulli_testing/bernoulli_test.stan'))
fit_bern <- sampling(bernoulli_model, data = d_bern, seed = 10)
