df_test_set <- NULL

for(lambda_H_real in c(0.1, 0.55, 1, 1.45, 1.9)){
  for (lambda_S_real in c(0.1, 0.55, 1, 1.45, 1.9)) {
    for (lambda_C_real in c(0.1, 0.55, 1, 1.45, 1.9)){
      for (exp_H_real in c(0.1, 0.55, 1, 1.45, 1.9)){
        
        # here we set an upper limit for the sum of all four lambdas
        
        df_test_set <- rbind(
          df_test_set,
          c(
            lambda_H_real, lambda_S_real, lambda_C_real, exp_H_real,
            lambda_H_real*1.  + lambda_S_real*1.  + lambda_C_real*(1. +1. )/2 + exp_H_real*1.  <= 5,
            lambda_H_real*0.7 + lambda_S_real*0.7 + lambda_C_real*(0.7+0.7)/2 + exp_H_real*0.7 <= 5,
            lambda_H_real*1.  + lambda_S_real*0.3 + lambda_C_real*(1. +0.3)/2 + exp_H_real*0.3 <= 5,
            lambda_H_real*0.3 + lambda_S_real*1.  + lambda_C_real*(0.3+1. )/2 + exp_H_real*1.  <= 5,
            lambda_H_real*0.3 + lambda_S_real*0.3 + lambda_C_real*(0.3+0.3)/2 + exp_H_real*0.3 <= 5
          )
          )
        
      }
    }
  }
}

# save the table
write.csv(df_test_set, "test_set.csv")
