ps_400 = read_csv(here('data/processed-data/ps400_post_data.csv'))
ps_2000 = read_csv(here('data/processed-data/ps2000_post_data.csv'))
ps_10000 = read_csv(here('data/processed-data/ps10000_post_data.csv'))

ps_400_mod = readRDS(here('output/intermediate-objects/gen_0_ps400_onerep_fit.rds'))
ps_2000_mod = readRDS(here('output/intermediate-objects/gen_0_ps2000_onerep_fit.rds'))
ps_10000_mod = readRDS(here('output/intermediate-objects/gen_0_ps10000_onerep_fit.rds'))

parms = c("theta_ll[1]", 
          #"theta_cq[1]", 
          "cstar", "NEC",
          "Lp", "Rm", "Lm", "tau_l", "tau_r")
gen_0_ps400_out = rstan::extract(ps_400_mod,
                                     permuted=TRUE,include=TRUE)
Z_df_ps400 = data.frame(gen_0_ps400_out$r_y_rep)
Z_means_ps400 = Z_df_ps400 %>% summarize(across(`X1`:`X22`, mean))
Z_means_ps400 = t(Z_means_ps400)



basinofattractionID = read.csv('C:/Users/brookson/Documents/Github/Coral-Resotration-Modeling/data/intermediate-files/basins_output/basinofattraction_a0.30_g0.00_z0.25.csv')


basinofattractionID <- data.frame(InitCond=rep(1:ntrajectory), Equilibrium = NA, 
                                  initM1 = initM[1:ntrajectory], 
                                  initC1 = initC[1:ntrajectory], initT1 = initT[1:ntrajectory])
basinofattractionID = basinofattractionID %>% 
  rename(InitCond = init_cond, Equilibrium = equilibrium, initM1 = init_M,
         initC1 = init_C, initT1 = init_T)

basinofattractionID$Equilibrium = c(rep(7, 105), rep(16, 105))

#change path to suit
png(paste0("C:/...Github/Coral-Restoration-Modeling/graphs/
           basinsplots/essentialcombinations/basinsplot_a",a,"z",z,"g"g '.png'))
a = 0.25;g = 0.3; z = 0
cols = c("#719D06", "#FF4040")
cols = cols[as.numeric(as.factor(basinofattractionID$Equilibrium))]
plot(basinofattractionID$initM1, basinofattractionID$initC1, 
     col = cols, 
     pch = 19,
     cex = 1.2,
     xlab = "Percent Cover Macroalgae",
     ylab = "Percent Cover Coral",
     main = paste("competiton = ", a, ", grazing = ", g, ", dispersal = ", z))
legend("topright", legend = levels(as.factor(
                                   as.character(basinofattractionID$Equilibrium))),
       col = c("#719D06", "#FF4040"),
       pch = c(19, 19))
dev.off()
