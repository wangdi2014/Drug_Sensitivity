# Mehmet Gonen (mehmet.gonen@gmail.com)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science

bayesian_multitask_multiple_kernel_learning_test <- function(Km, state) {
  T <- length(Km)
  N <- matrix(0, T, 1)
  for (o in 1:T) {
    N[o] <- dim(Km[[o]])[2]
  }
  P <- dim(Km[[1]])[3]

  G <- vector("list", T)
  for (o in 1:T) {
    G[[o]] <- list(mean = matrix(0, P, N[o]), covariance = matrix(0, P, N[o]))
    for (m in 1:P) {
      G[[o]]$mean[m,] <- crossprod(state$a[[o]]$mean, Km[[o]][,,m])
      G[[o]]$covariance[m,] <- 1 / (state$upsilon$shape[o] * state$upsilon$scale[o]) + diag(crossprod(Km[[o]][,,m], state$a[[o]]$covariance) %*% Km[[o]][,,m])
    } 
  }

  f <- vector("list", T)
  for (o in 1:T) {
    f[[o]] <- list(mean = matrix(0, N[o], 1), covariance = matrix(0, N[o], 1))
    f[[o]]$mean <- crossprod(rbind(matrix(1, 1, N[o]), G[[o]]$mean), state$be$mean[c(o, (T + 1):(T + P))])
    f[[o]]$covariance <- 1 / (state$epsilon$shape[o] * state$epsilon$scale[o]) + diag(crossprod(rbind(matrix(1, 1, N[o]), G[[o]]$mean), state$be$covariance[c(o, (T + 1):(T + P)), c(o, (T + 1):(T + P))]) %*% rbind(matrix(1, 1, N[o]), G[[o]]$mean))
  }

  prediction <- list(f = f)
}
