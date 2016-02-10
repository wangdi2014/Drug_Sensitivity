# Mehmet Gonen (mehmet.gonen@gmail.com)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science

logdet <- function(Sigma) {
    2 * sum(log(diag(chol(Sigma))))
}

repmat <- function(M, row, column) {
    kronecker(matrix(1, row, column), M)
}

bayesian_multitask_multiple_kernel_learning_train <- function(Km, y, parameters) {
  set.seed(parameters$seed)

  T <- length(Km)
  D <- matrix(0, T, 1)
  N <- matrix(0, T, 1)
  for (o in 1:T) {
    D[o] <- dim(Km[[o]])[1]
    N[o] <- dim(Km[[o]])[2]
  }
  P <- dim(Km[[1]])[3]

  log2pi <- log(2 * pi)

  lambda <- vector("list", T)
  for (o in 1:T) {
    lambda[[o]] <- list(shape = matrix(parameters$alpha_lambda + 0.5, D[o], 1), scale = matrix(parameters$beta_lambda, D[o], 1))
  }
  upsilon <- list(shape = matrix(parameters$alpha_upsilon + 0.5 * N * P, T, 1), scale = matrix(parameters$beta_upsilon, T, 1))
  a <- vector("list", T)
  for (o in 1:T) {
    a[[o]] <- list(mean = matrix(rnorm(D[o]), D[o], 1), covariance = diag(1, D[o], D[o]))
  }
  G <- vector("list", T)
  for (o in 1:T) {
    G[[o]] <- list(mean = matrix(rnorm(P * N[o]), P, N[o]), covariance = diag(1, P, P))
  }
  gamma <- list(shape = matrix(parameters$alpha_gamma + 0.5, T, 1), scale = matrix(parameters$beta_gamma, T, 1))
  omega <- list(shape = matrix(parameters$alpha_omega + 0.5, P, 1), scale = matrix(parameters$beta_omega, P, 1))
  epsilon <- list(shape = matrix(parameters$alpha_epsilon + 0.5 * N, T, 1), scale = matrix(parameters$beta_epsilon, T, 1))
  be <- list(mean = rbind(matrix(0, T, 1), matrix(1, P, 1)), covariance = diag(1, T + P, T + P))

  KmKm <- vector("list", T)
  for (o in 1:T) {
      KmKm[[o]] <- matrix(0, D[o], D[o])
      for (m in 1:P) {
          KmKm[[o]] <- KmKm[[o]] + tcrossprod(Km[[o]][,,m], Km[[o]][,,m])
      }
      Km[[o]] <- matrix(Km[[o]], D[o], N[o] * P)
  }

  if (parameters$progress == 1) {
    bounds <- matrix(0, parameters$iteration, 1)
  }

  atimesaT.mean <- vector("list", T)
  for (o in 1:T) {
    atimesaT.mean[[o]] <- tcrossprod(a[[o]]$mean, a[[o]]$mean) + a[[o]]$covariance
  }
  GtimesGT.mean <- vector("list", T)
  for (o in 1:T) {
    GtimesGT.mean[[o]] <- tcrossprod(G[[o]]$mean, G[[o]]$mean) + N[o] * G[[o]]$covariance
  }
  btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
  etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
  etimesb.mean <- matrix(0, P, T)
  for (o in 1:T) {
    etimesb.mean[,o] <- be$mean[(T + 1):(T + P)] * be$mean[o] + be$covariance[(T + 1):(T + P), o]
  }
  KmtimesGT.mean <- vector("list", T)
  for (o in 1:T) {
    KmtimesGT.mean[[o]] <- Km[[o]] %*% matrix(t(G[[o]]$mean), N[o] * P, 1)
  }
  for (iter in 1:parameters$iteration) {
    # update lambda
    for (o in 1:T) {
      lambda[[o]]$scale <- 1 / (1 / parameters$beta_lambda + 0.5 * diag(atimesaT.mean[[o]]))
    }
    # update upsilon
    for (o in 1:T) {
      upsilon$scale[o] <- 1 / (1 / parameters$beta_upsilon + 0.5 * (sum(diag(GtimesGT.mean[[o]])) - 2 * sum(matrix(crossprod(a[[o]]$mean, Km[[o]]), N[o], P) * t(G[[o]]$mean)) + sum(diag(KmKm[[o]] %*% atimesaT.mean[[o]]))))
    }
    # update a
    for (o in 1:T) {
      a[[o]]$covariance <- chol2inv(chol(diag(as.vector(lambda[[o]]$shape * lambda[[o]]$scale), D[o], D[o]) + upsilon$shape[o] * upsilon$scale[o] * KmKm[[o]]))
      a[[o]]$mean <- a[[o]]$covariance %*% (upsilon$shape[o] * upsilon$scale[o] * KmtimesGT.mean[[o]])
      atimesaT.mean[[o]] <- tcrossprod(a[[o]]$mean, a[[o]]$mean) + a[[o]]$covariance
    }
    # update G
    for (o in 1:T) {
      G[[o]]$covariance <- chol2inv(chol(diag(upsilon$shape[o] * upsilon$scale[o], P, P) + epsilon$shape[o] * epsilon$scale[o] * etimeseT.mean))
      G[[o]]$mean <- G[[o]]$covariance %*% (upsilon$shape[o] * upsilon$scale[o] * t(matrix(crossprod(a[[o]]$mean, Km[[o]]), N[o], P)) + epsilon$shape[o] * epsilon$scale[o] * (tcrossprod(be$mean[(T + 1):(T + P)], y[[o]]) - repmat(etimesb.mean[,o], 1, N[o])))
      GtimesGT.mean[[o]] <- tcrossprod(G[[o]]$mean, G[[o]]$mean) + N[o] * G[[o]]$covariance
      KmtimesGT.mean[[o]] <- Km[[o]] %*% matrix(t(G[[o]]$mean), N[o] * P, 1)
    }
    # update gamma
    gamma$scale <- 1 / (1 / parameters$beta_gamma + 0.5 * diag(btimesbT.mean))
    # update omega
    omega$scale <- 1 / (1 / parameters$beta_omega + 0.5 * diag(etimeseT.mean))
    # update epsilon
    for (o in 1:T) {
      epsilon$scale[o] <- 1 / (1 / parameters$beta_epsilon + 0.5 * as.double(crossprod(y[[o]], y[[o]]) - 2 * crossprod(y[[o]], crossprod(rbind(matrix(1, 1, N[o]), G[[o]]$mean), be$mean[c(o, (T + 1):(T + P))])) + N[o] * btimesbT.mean[o, o] + sum(diag(GtimesGT.mean[[o]] %*% etimeseT.mean)) + 2 * sum(diag(crossprod(rowSums(G[[o]]$mean), etimesb.mean[,o])))))
    }
    # update b and e
    be$covariance <- rbind(cbind(diag(as.vector(gamma$shape * gamma$scale), T, T) + diag(as.vector(N * epsilon$shape * epsilon$scale), T, T), matrix(0, T, P)), cbind(matrix(0, P, T), diag(as.vector(omega$shape * omega$scale), P, P)))
    for (o in 1:T) {
      be$covariance[(T + 1):(T + P), o] <- epsilon$shape[o] * epsilon$scale[o] * rowSums(G[[o]]$mean)
      be$covariance[o, (T + 1):(T + P)] <- epsilon$shape[o] * epsilon$scale[o] * t(rowSums(G[[o]]$mean))
      be$covariance[(T + 1):(T + P), (T + 1):(T + P)] <- be$covariance[(T + 1):(T + P), (T + 1):(T + P)] + epsilon$shape[o] * epsilon$scale[o] * GtimesGT.mean[[o]]
    }
    be$covariance <- chol2inv(chol(be$covariance))
    be$mean <- matrix(0, T + P, 1)
    for (o in 1:T) {
      be$mean[o] <- epsilon$shape[o] * epsilon$scale[o] * sum(y[[o]])
      be$mean[(T + 1):(T + P)] <- be$mean[(T + 1):(T + P)] + epsilon$shape[o] * epsilon$scale[o] * G[[o]]$mean %*% y[[o]]
    }
    be$mean <- be$covariance %*% be$mean
    btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
    etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
    for (o in 1:T) {
        etimesb.mean[,o] <- be$mean[(T + 1):(T + P)] * be$mean[o] + be$covariance[(T + 1):(T + P), o]
    }

    if (parameters$progress == 1) {
      lb <- 0

      # p(lambda)
      for (o in 1:T) {
        lb <- lb + sum((parameters$alpha_lambda - 1) * (digamma(lambda[[o]]$shape) + log(lambda[[o]]$scale)) - lambda[[o]]$shape * lambda[[o]]$scale / parameters$beta_lambda - lgamma(parameters$alpha_lambda) - parameters$alpha_lambda * log(parameters$beta_lambda))
      }
      # p(upsilon)
      lb <- lb + sum((parameters$alpha_upsilon - 1) * (digamma(upsilon$shape) + log(upsilon$scale)) - upsilon$shape * upsilon$scale / parameters$beta_upsilon - lgamma(parameters$alpha_upsilon) - parameters$alpha_upsilon * log(parameters$beta_upsilon))
      # p(a | lambda)
      for (o in 1:T) {
        lb <- lb - 0.5 * sum(diag(diag(as.vector(lambda[[o]]$shape * lambda[[o]]$scale), D[o], D[o]) %*% atimesaT.mean[[o]])) - 0.5 * (D[o] * log2pi - sum(log(lambda[[o]]$shape * lambda[[o]]$scale)))
      }
      # p(G | a, Km, upsilon)
      for (o in 1:T) {
        lb <- lb - 0.5 * sum(diag(GtimesGT.mean[[o]])) * upsilon$shape[o] * upsilon$scale[o] + crossprod(a[[o]]$mean, KmtimesGT.mean[[o]]) * upsilon$shape[o] * upsilon$scale[o] - 0.5 * sum(diag(KmKm[[o]] %*% atimesaT.mean[[o]])) * upsilon$shape[o] * upsilon$scale[o] - 0.5 * N[o] * P * (log2pi - log(upsilon$shape[o] * upsilon$scale[o]))
      }
      # p(gamma)
      lb <- lb + sum((parameters$alpha_gamma - 1) * (digamma(gamma$shape) + log(gamma$scale)) - gamma$shape * gamma$scale / parameters$beta_gamma - lgamma(parameters$alpha_gamma) - parameters$alpha_gamma * log(parameters$beta_gamma))
      # p(b | gamma)
      lb <- lb - 0.5 * sum(diag(diag(as.vector(gamma$shape * gamma$scale), T, T) %*% btimesbT.mean)) - 0.5 * (T * log2pi - sum(log(gamma$shape * gamma$scale)))
      # p(omega)
      lb <- lb + sum((parameters$alpha_omega - 1) * (digamma(omega$shape) + log(omega$scale)) - omega$shape * omega$scale / parameters$beta_omega - lgamma(parameters$alpha_omega) - parameters$alpha_omega * log(parameters$beta_omega))
      # p(e | omega)
      lb <- lb - 0.5 * sum(diag(diag(as.vector(omega$shape * omega$scale), P, P) %*% etimeseT.mean)) - 0.5 * (P * log2pi - sum(log(omega$shape * omega$scale)))
      # p(epsilon)
      lb <- lb + sum((parameters$alpha_epsilon - 1) * (digamma(epsilon$shape) + log(epsilon$scale)) - epsilon$shape * epsilon$scale / parameters$beta_epsilon - lgamma(parameters$alpha_epsilon) - parameters$alpha_epsilon * log(parameters$beta_epsilon))
      # p(y | b, e, G, epsilon)
      for (o in 1:T) {
        lb <- lb - 0.5 * crossprod(y[[o]], y[[o]]) * epsilon$shape[o] * epsilon$scale[o] + crossprod(y[[o]], crossprod(G[[o]]$mean, be$mean[(T + 1):(T + P)])) * epsilon$shape[o] * epsilon$scale[o] + sum(be$mean[o] * y[[o]]) * epsilon$shape[o] * epsilon$scale[o] - 0.5 * sum(diag(etimeseT.mean %*% GtimesGT.mean[[o]])) * epsilon$shape[o] * epsilon$scale[o] - sum(crossprod(G[[o]]$mean, etimesb.mean[,o])) * epsilon$shape[o] * epsilon$scale[o] - 0.5 * N[o] * btimesbT.mean[o,o] * epsilon$shape[o] * epsilon$scale[o] - 0.5 * N[o] * (log2pi - log(epsilon$shape[o] * epsilon$scale[o]))
      }

      # q(lambda)
      for (o in 1:T) {
        lb <- lb + sum(lambda[[o]]$shape + log(lambda[[o]]$scale) + lgamma(lambda[[o]]$shape) + (1 - lambda[[o]]$shape) * digamma(lambda[[o]]$shape))
      }
      # q(upsilon)
      lb <- lb + sum(upsilon$shape + log(upsilon$scale) + lgamma(upsilon$shape) + (1 - upsilon$shape) * digamma(upsilon$shape))
      # q(a)
      for (o in 1:T) {
        lb <- lb + 0.5 * (D[o] * (log2pi + 1) + logdet(a[[o]]$covariance))
      }
      # q(G)
      for (o in 1:T) {
        lb <- lb + 0.5 * N[o] * (P * (log2pi + 1) + logdet(G[[o]]$covariance))
      }
      # q(gamma)
      lb <- lb + sum(gamma$shape + log(gamma$scale) + lgamma(gamma$shape) + (1 - gamma$shape) * digamma(gamma$shape))
      # q(omega)
      lb <- lb + sum(omega$shape + log(omega$scale) + lgamma(omega$shape) + (1 - omega$shape) * digamma(omega$shape))
      # q(epsilon)
      lb <- lb + sum(epsilon$shape + log(epsilon$scale) + lgamma(epsilon$shape) + (1 - epsilon$shape) * digamma(epsilon$shape))
      # q(b, e)
      lb <- lb + 0.5 * ((T + P) * (log2pi + 1) + logdet(be$covariance))

      bounds[iter] <- lb
    }
  }
  
  if (parameters$progress == 1) {
    state <- list(lambda = lambda, upsilon = upsilon, a = a, gamma = gamma, omega = omega, epsilon = epsilon, be = be, bounds = bounds, parameters = parameters)
  }
  else {
    state <- list(lambda = lambda, upsilon = upsilon, a = a, gamma = gamma, omega = omega, epsilon = epsilon, be = be, parameters = parameters)
  }
}
