
# Q-matrix specifications for different model structures

qmat_list <- list(
  base_model = matrix(c(1, 1, 1, 1,
                        1, 1, 1, 1,
                        0, 0, 0, 0,
                        0, 0, 0, 0), 
                      nrow = 4, ncol = 4, byrow = TRUE,
                      dimnames = list(c("M", "S", "D", "R"),
                                      c("M", "S", "D", "R"))),
  
  mod_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5, byrow = TRUE,
                 dimnames = list(c("M1", "M2", "S", "D", "R"),
                                 c("M1", "M2", "S", "D", "R"))),
  
  mod_3 = matrix(c(1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1, 
                   0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 0), 
                 nrow = 6, ncol = 6, byrow = TRUE,
                 dimnames = list(c("M1", "M2", "M3", "S", "D", "R"),
                                 c("M1", "M2", "M3", "S", "D", "R"))),
  
  sev_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5, byrow = TRUE,
                 dimnames = list(c("M", "S1", "S2", "D", "R"),
                                 c("M", "S1", "S2", "D", "R"))),
  
  reduced_trans = matrix(c(1, 1, 0, 1,
                           1, 1, 1, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0), 
                         nrow = 4, ncol = 4, byrow = TRUE,
                         dimnames = list(c("M", "S", "D", "R"),
                                         c("M", "S", "D", "R"))),
  
  hx_sev = matrix(c(1, 0, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE,
                  dimnames = list(c("M", "MS", "S", "D", "R"),
                                  c("M", "MS", "S", "D", "R")))
)
