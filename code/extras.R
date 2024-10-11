#     for(m in 1:M_count) {
#        eps_pms[s,m,c] ~ dt(0,1,1)T(0,) 
#        # create diagonal matrix
#        for(j in 1:M_count) {
#          diag.eps_pms[s,m,j,c] <- ifelse(m==j, eps_pms[s,m,c], 0)
#        }
#      }
#      Sigma.alpha_pms[s, 1:M_count,1:M_count,c] <- diag.eps_pms[s, 1:M_count,1:M_count, c]%*%Q_pms[s, 1:M_count,1:M_count, c]%*%diag.eps_pms[s, 1:M_count,1:M_count, c]
#      inv.Sigma.alpha_pms[s, 1:M_count,1:M_count,c] <- inverse(Sigma.alpha_pms[s, 1:M_count,1:M_count,c])
#   }

# for(m in 1:M_count) {
#   sd.alpha_pms[s, m, c] ~ dt(0,1,1)T(0,) # estimate sd terms
#   for(j in 1:M_count) {
#     rho.trans_alpha_pms[s,m,j,c] ~ dnorm(0, pow(sd.Z[s], -2)) # Expected corr to be 0
#     rho.alpha_pms[s,m,j,c] <- ((exp(rho.trans_alpha_pms[s,m,j,c])-1)/(exp(rho.trans_alpha_pms[s,m,j,c])+1)) # Fishers inverse-transformation of correlations
#   }
#   Sigma.alpha_pms[s,m,m,c] <- pow(sd.alpha_pms[s, m, c], 2)
#   for (j in (m+1):M_count) {
#     Sigma.alpha_pms[s,m,j,c] <- rho.alpha_pms[s,m,j,c] * sd.alpha_pms[s,m,c] * sd.alpha_pms[s,j,c]
#     Sigma.alpha_pms[s,j,m,c] <- Sigma.alpha_pms[s,m,j,c] 
#   }
# }
# Inverse for province-level params
# inv.Sigma.alpha_pms[s, 1:M_count,1:M_count,c]  <- inverse(Sigma.alpha_pms[s, 1:M_count,1:M_count,c])