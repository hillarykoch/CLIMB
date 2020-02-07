# Functions to call anything from cgibbs.jl

prepare_julia <- function() {
  # Find julia v1.0.2 binary
  julia <- JuliaCall::julia_setup()
  ver <-
    as.numeric(stringr::str_split(string = julia$VERSION, pattern = "\\.")[[1]][1])
  if (ver < 1) {
    stop("Julia version > 1.0 required for this package to run.")
  }

  # load cgibbs module
  JuliaCall::julia_command("using cgibbs")
  JuliaCall::julia_command("using StatsBase")
  JuliaCall::julia_command("using DataFrames")
  JuliaCall::julia_command("using LinearAlgebra")
  JuliaCall::julia_command("using Distributions")
}


run_mcmc <- function(dat, hyp, nstep, retained_classes) {
  if (length(hyp$alpha) != nrow(retained_classes)) {
    stop(
      "length(hyp$alpha) must be the same as nrow(retained_classes). These both
      correspond to the cluster number."
    )
  }

  dat <- data.frame(dat)

  prepare_julia()

  # Compute some values from the inputs
  dm <- ncol(dat)
  n <- nrow(dat)
  nm <- length(hyp$alpha)
  nw <- 1

  # Assign inputs names in Julia
  JuliaCall::julia_assign("dat", dat)
  JuliaCall::julia_assign("kappa0", hyp$kappa0)
  JuliaCall::julia_assign("mu0", hyp$mu0)
  JuliaCall::julia_assign("Psi0", hyp$Psi0)
  JuliaCall::julia_command("hyp = (kappa0, mu0, Psi0);")
  JuliaCall::julia_assign("alph", hyp$alpha)
  JuliaCall::julia_assign("reduced_classes", retained_classes)

  # Conver floats to integers when appropriate
  JuliaCall::julia_assign("nw", 1)
  JuliaCall::julia_assign("nw", JuliaCall::julia_eval("Int64(nw)"))

  JuliaCall::julia_assign("dm", dm)
  JuliaCall::julia_assign("dm", JuliaCall::julia_eval("Int64(dm)"))

  JuliaCall::julia_assign("n", n)
  JuliaCall::julia_assign("n", JuliaCall::julia_eval("Int64(n)"))

  JuliaCall::julia_assign("nm", nm)
  JuliaCall::julia_assign("nm", JuliaCall::julia_eval("Int64(nm)"))

  JuliaCall::julia_assign("nstep", nstep)
  JuliaCall::julia_assign("nstep", JuliaCall::julia_eval("Int64(nstep)"))

  JuliaCall::julia_assign("tune_df", rep(1000.0, nm))

  # Generate starting values
  JuliaCall::julia_assign(
    "param",
    JuliaCall::julia_eval(
      "Array{Tuple{Array{Dict{String,Array{Float64,N} where N},1},Array{Float64,1},Array{Int64,1}}}(undef, (nw, 1));"
    )
  )
  for (i in 1:nw) {
    JuliaCall::julia_assign(
      "dictionary",
      JuliaCall::julia_eval("Dict{String,Array{Float64,N} where N}[];")
    )
    for (m in 1:nm) {
      JuliaCall::julia_assign("m", m)
      JuliaCall::julia_assign("m", JuliaCall::julia_eval("Int64(m)"))

      JuliaCall::julia_assign("Sigma",
                              JuliaCall::julia_eval("Matrix(Hermitian(Psi0[:,:,m]));"))
      JuliaCall::julia_assign("mu", JuliaCall::julia_eval("mu0[m,:];"))

      JuliaCall::julia_command("push!(dictionary, Dict(\"mu\" => mu, \"Sigma\" => Sigma));")
    }
    JuliaCall::julia_assign("z",
                            JuliaCall::julia_eval("wsample(1:1:nm, alph, n; replace = true);"))

    JuliaCall::julia_assign("i", i)
    JuliaCall::julia_command("param[i,1] = (dictionary, alph, z);")
  }

  labels <-
    apply(retained_classes + 1, 1, function(X)
      paste0(X, collapse = ""))
  JuliaCall::julia_assign("labels", labels)

  out <-
    JuliaCall::julia_eval("cgibbs.run_mcmc(dat, param, hyp, alph, nstep, labels, tune_df);")
  names(out) <-
    c("chain", "acceptance_rate_chain", "tune_df_chain")
  out
  }


extend_mcmc <- function(dat, hyp, nstep, retained_classes, mcmc) {
  prepare_julia()

  dat <- data.frame(dat)

  # Compute some values from the inputs
  dm <- ncol(dat)
  n <- nrow(dat)
  nm <- length(hyp$alpha)
  nw <- 1

  # Assign inputs names in Julia
  JuliaCall::julia_assign("dat", dat)
  JuliaCall::julia_assign("kappa0", hyp$kappa0)
  JuliaCall::julia_assign("mu0", hyp$mu0)
  JuliaCall::julia_assign("Psi0", hyp$Psi0)
  JuliaCall::julia_command("hyp = (kappa0, mu0, Psi0);")
  JuliaCall::julia_assign("alph", hyp$alpha)

  # Assign inputs names in Julia
  JuliaCall::julia_assign("nstep", nstep)
  JuliaCall::julia_assign("nstep", JuliaCall::julia_eval("Int64(nstep)"))

  JuliaCall::julia_assign("chain", mcmc$chain)
  JuliaCall::julia_assign("tune_df_chain", mcmc$tune_df_chain)
  JuliaCall::julia_assign("labels", apply(retained_classes + 1, 1, function(X)
    paste0(X, collapse = "")))

  out <-
    JuliaCall::julia_eval("cgibbs.extend_mcmc(dat, chain, tune_df_chain, hyp, alph, nstep, labels);")
  names(out) <- c("chain", "acceptance_rate_chain", "tune_df_chain")
  out
}

extract_chains <- function(mcmc) {
  prepare_julia()

  acceptance_rate_chain <- t(mcmc$acceptance_rate_chain)
  tune_df_chain <- t(mcmc$tune_df_chain)
  nm <- ncol(acceptance_rate_chain)

  JuliaCall::julia_assign("chain", mcmc$chain)
  JuliaCall::julia_assign("nm", nm)
  JuliaCall::julia_assign("nm", JuliaCall::julia_eval("Int64(nm)"))

  mu_chains <- Sigma_chains <- list()
  for (m in 1:nm) {
    JuliaCall::julia_assign("mu_chains",
                            JuliaCall::julia_eval("[cgibbs.get_mu_chain(chain, m) for m in 1:nm]"))
    JuliaCall::julia_assign("m", m)
    JuliaCall::julia_assign("m", JuliaCall::julia_eval("Int64(m)"))

    mu_chains[[m]] <- t(JuliaCall::julia_eval("mu_chains[m]"))
    Sigma_chains[[m]] <-
      JuliaCall::julia_eval("cgibbs.get_Sigma_chain(chain, m)")
  }

  z_chain <- JuliaCall::julia_eval("cgibbs.get_z_chain(chain);")
  prop_chain <- t(JuliaCall::julia_eval("cgibbs.get_prop_chain(chain);"))

  list(
    "mu_chains" = mu_chains,
    "Sigma_chains" = Sigma_chains,
    "prop_chain" = prop_chain,
    "z_chain" = z_chain
  )
}
