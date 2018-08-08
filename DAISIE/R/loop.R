# library("snow")
# library("parallel")
# 
# cpu <- detectCores() # Number of cores requested.
# hosts <- rep("localhost",cpu)
# cl <- makeCluster(hosts, type = "SOCK")

pb <- txtProgressBar(min = 0, max = 1000, style = 3)

stt_list <- list()

stt_list[[1]] <- DAISIE_exinction_test(10,
                                       mainland_n = 1000,
                                       pars = c(0, 0, Inf, 0, 0),
                                       Apars = Apars, 
                                       Epars = Epars, 
                                       island_ontogeny = "quadratic")
plot(stt_list[[1]][1:1001,2], stt_list[[1]][1:1001,1], col = "grey", lwd = 1, type = "l",
     xlab = "Time", ylab = "Number of species", main = "Number of species through time")
for (i in 2:1000) {
  stt_list[[i]] <- DAISIE_exinction_test(10,
                                         mainland_n = 1000,
                                         pars = c(0, 0, Inf, 0, 0),
                                         Apars = Apars, 
                                         Epars = Epars, 
                                         island_ontogeny = "quadratic")
  lines(stt_list[[i]][1:1001,2], stt_list[[1]][1:1001,1], col = "grey", lwd = 1)
  setTxtProgressBar(pb, i)
}

close(pb)

legend("right", inset = 0.05, title = "Data origin", legend=c("Simulation data", "Differential equation"),
       col=c("grey", "red"), lty=1, cex=1, text.font = 2, )
