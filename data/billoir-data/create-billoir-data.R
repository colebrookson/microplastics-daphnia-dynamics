
# growth data ==================================================================

# growth data for control (7, 14, 21) with 59 replicates 
day7con = c(2.61,2.82,2.73,2.58,2.81,2.56,2.52,2.81,2.58,2.80,
            2.64,2.80,2.59,2.81,2.60,2.79,2.78,2.58,2.65,2.71,3.78,2.79,
            2.74,2.69,2.66,2.53,2.68,2.62,2.63,2.74,2.76,2.86,2.82,2.64,
            2.73,2.65,2.49,2.71,2.88,2.57,2.84,2.63,2.79,2.56,2.66,2.73,
            2.59,2.76,2.54,2.59,2.57,2.74,2.65,2.76,2.73,2.65,2.73,2.63,2.54)
day14con = c(3.71,3.58,3.71,3.66,3.53,3.65,3.66,3.54,3.69,3.54,
             3.61,3.65,3.61,3.53,3.38,3.55,3.50,3.59,NA,NA,3.58,3.63,3.63,
             3.71,3.66,3.71,3.78,3.66,3.65,3.81,3.82,3.69,3.70,3.81,3.67,
             3.77,NA,NA,NA,NA,3.47,3.59,3.70,3.64,3.49,3.54,3.51,3.53,3.80,
             3.65,3.70,3.68,3.64,3.55,3.73,3.60,3.64,3.60,NA)
day21con = c(3.82,4.01,3.75,4.00,3.91,3.88,3.79,3.63,3.95,4.01,
             3.96,3.90,3.93,3.79,3.82,3.72,3.79,3.96,NA,NA,4.00,3.89,3.90,
             3.89,3.66,3.85,3.92,3.81,3.87,3.96,3.90,4.05,3.85,3.94,3.85,
             3.83,3.82,3.91,3.65,NA,3.70,3.79,3.82,3.76,3.83,3.75,4.01,
             3.86,3.96,3.85,3.98,3.82,3.90,3.93,3.90,3.97,3.85,3.79,3.76)
# concat all together 
Y0 = data.frame(day7con, day14con, day21con)

# growth data for concentration 1
day7c1 = c(2.55,2.65,2.67,2.70,2.68,2.63,2.64,2.70,2.61,2.73,2.66,
           2.72,2.69,2.48,2.69,2.71,2.65,2.66,2.70,2.68,2.57,2.60,
           2.83,2.51,2.62,2.62,2.47,2.59,2.70,2.73,2.73,2.49,2.66,
           2.71,2.71,2.53,2.73,2.75,2.69,2.63,2.78,2.89,2.68,2.61,
           2.61,2.60,2.75,2.75,2.80,2.62,2.82,2.61,2.74,2.63,2.65,
           2.71,2.80,2.74,NA)
day14c1 = c(3.65,3.64,3.53,3.68,3.61,3.50,3.49,3.59,3.61,3.50,3.56,
            3.53,3.66,3.57,3.54,3.68,3.59,NA,NA,NA,3.61,3.68,3.55,
            3.57,3.42,3.60,3.58,3.52,3.57,3.66,3.50,3.15,3.50,3.52,
            3.55,3.59,3.68,3.64,3.46,3.56,3.53,3.43,3.60,3.54,3.52,
            3.57,3.47,3.66,3.61,3.75,3.53,3.52,3.59,3.57,3.58,3.64,
            3.45,3.71,NA)
day21c1 = c(3.69,3.76,3.97,3.95,3.50,3.77,3.80,3.77,3.85,3.75,3.80,
            3.97,4.03,3.77,3.85,3.64,3.76,3.60,3.81,NA,3.87,3.70,
            3.88,3.81,3.71,3.53,3.57,3.66,3.82,3.67,3.87,3.88,3.89,
            3.85,3.83,3.84,3.91,3.88,3.68,NA,3.70,3.53,3.71,3.62,3.63,
            3.53,3.57,3.65,3.60,3.70,4.00,3.62,3.69,3.43,3.57,3.53,
            3.70,3.64,NA)
# concat all together
Y1 = data.frame(day7c1, day14c1, day21c1)

# growth data for concentration 2
day7c2 = c(2.35,2.39,2.36,2.31,2.31,2.23,2.38,2.36,2.29,2.14,2.41,2.34,
           2.362.28,2.23,2.18,NA,NA,NA,NA,2.26,2.44,2.53,1.86,2.28,2.26,
           2.32,2.41,2.30,2.33,2.33,2.20,2.50,2.15,2.35,NA,NA,NA,NA,NA,
           2.02,2.20,2.19,2.25,2.36,2.24,1.87,2.42,1.92,2.18,2.28,2.37,
           2.39,2.43,2.39,2.14,1.52,2.27,NA)
day14c2 = c(2.82,2.73,2.74,2.60,2.75,2.79,2.96,2.93,2.81,2.60,2.78,2.61,
            3.01,2.74,NA,NA,NA,NA,NA,NA,2.80,2.86,2.77,2.77,2.95,2.75,2.77,
            2.45,2.83,2.69,2.51,2.74,2.74,2.90,2.87,NA,NA,NA,NA,NA,2.66,
            2.66,2.81,2.86,2.78,2.77,2.67,2.85,2.70,2.77,2.98,NA,NA,NA,
            NA,NA,NA,NA,NA)
day21c2 = c(2.88,2.98,3.09,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
            NA,NA,2.99,2.90,2.72,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
            NA,NA,NA,NA,2.90,2.94,2.78,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
            NA,NA,NA,NA,NA)
# concat all together
Y1 = data.frame(day7c2, day14c2, day21c2)

## reproduction data 

# reproduction data for controls (days 1-21 with three replicates)
repcon = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.79,6.32,7.10,
           6.63,7.58,8.75,7.84,8.68,8.75,16.51,20.26,19.10,19.56,21.95,
           20.35,22.62,21.95,20.40,28.68,32.21,32.40,31.62,34.68,33.95,
           33.56,36.42,33.95,40.45,49.11,47.40,45.62,52.32,49.80,48.73,
           53.16,50.15,54.56,65.95,60.83,63.51,71.32,67.57)

# reproduction data for concentration 1
repc1 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.63,6.58,6.63,7.63,
          7.37,7.32,7.63,7.37,7.32,19.58,18.00,18.43,19.58,18.00,19.09,
          19.58,19.89,19.93,33.37,33.11,33.70,35.68,33.11,34.82,35.68,
          34.68,34.82,43.68,49.47,50.15,47.47,49.47,51.09,47.47,51.53,
          51.09,59.26,65.42,62.70,67.21,69.26,66.04)

# reproduction data for concentration 2
repc2 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.47,4.94,3.45,6.37,
          6.61,5.35,6.37,6.61,5.85,14.32,14.67,15.00,15.58,16.61,15.05,
          16.37,17.28,16.30,26.74,27.61,27.00,29.63,31.17,28.45,30.32,
          31.94,29.70,40.43,41.94,41.00,44.87,45.82,42.80,45.70,46.82,
          44.12,51.29,55.76,51.43,57.42,59.01,58.26)

# reproduction data for concentration 3
repc3 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00,0.00,0.00,
          0.00,0.00,0.00,1.06,0.00,0.00,3.81,0.15,0.90,6.00,0.70,1.40,
          6.86,0.70,1.78,8.08,0.70,1.78,10.15,0.70,2.78,10.15,0.70,2.78,
          10.15,1.50,3.28,10.15,1.90,3.28,10.15,1.90,3.28,10.15,3.23,
          3.28,10.15,4.90,3.28)

# note ####################################################################
# since the dimensions of the reproduction data are weird (it's 21,3) so it
# goes day 7, 14, 21, 7 14, 21, etc, we need to fix this so the data are 
# organized properly 
# end note ################################################################

split_reproduction_data = function(data) {
     
            # go through some gross inefficient loop to pull out the values
            # and put them in three objects then make three objects into
            # df columns then return the df
}










