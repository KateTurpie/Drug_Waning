# ---------------------------------------------------------------------------- #

# ---- RSV Transmission Model ---- #

# Purpose : Defines the core age-structured RSV dynamic transmission model
# By : Kate Turpie
# Updated : 07-July-26


# ---------------------------------------------------------------------------- #

model_sirs <- odin2::odin(
  {
    # ---- ODE & Output ----
    
    # -- Population Level --
    update(S[]) <- S[i] - n_SIu[i] - n_SP[i] + n_PS[i]
    update(S2[]) <- S2[i] - n_S2Iu[i] + n_RS[i]
    update(P[]) <- P[i] + n_SP[i] - n_PS[i] - n_PIp[i]
    
    update(Iu[]) <- Iu[i] + n_SIu[i] + n_S2Iu[i] - n_IuR[i]
    update(Ip[]) <- Ip[i] + n_PIp[i] - n_IpR[i]
    update(I[]) <- Iu[i] + n_SIu[i] + n_S2Iu[i] - n_IuR[i] + Ip[i] + n_PIp[i] - n_IpR[i]
    
    update(R[]) <- R[i] + n_IuR[i] + n_IpR[i] - n_RS[i]
    
    # -- Trial Level --
    update(Sc[]) <- Sc[i] - n_ScIc[i]
    update(Stx[]) <- Stx[i] - n_StxItx[i]
    
    update(Ic[]) <- Ic[i] + n_ScIc[i]
    update(Itx[]) <- Itx[i] + n_StxItx[i]
    
    # -- Medically Attended Cases --
    update(Iu_ma[]) <- Iu[i] * p_ma[i]
    update(Ip_ma[]) <- Ip[i] * p_ma[i]
    update(I_ma[]) <- (Iu[i] + Ip[i]) * p_ma[i]
    update(Ic_ma[]) <- Ic[i] * p_ma[i]
    update(Itx_ma[]) <- Itx[i] * p_ma[i]
    #assumption: vaccination does not impact proportion of cases that are 
    #medically attended, just the number of cases overall (which in turn lowers 
    #number of medically attended cases in protected group)
    
    # -- Cumulative Cases --
    update(TotalCases_pop[]) <- TotalCases_pop[i] + n_SIu[i] + n_S2Iu[i] + n_PIp[i]
    update(TotalCases_trial[]) <- TotalCases_trial[i] + n_StxItx[i] + n_ScIc[i]
    update(TotalCases_popMA[]) <- TotalCases_pop[i] * p_ma[i]
    update(TotalCases_trialMA[]) <- TotalCases_trial[i] * p_ma[i]
    
    # -- Population Size Check --
    update(n_pop) <- sum(S) + sum(S2) + sum(P) +  sum(Iu) + sum(Ip) + sum(R)
    update(n_trial) <- sum(Stx) + sum(Sc) + sum(Itx) + sum(Ic)
    
    # -- FoI Check --
    update(lambda_check[]) <- lambda[i]
    update(beta_check) <- beta
    
    # -- Trial Efficacy --
    update(efficacy_true) <- (chi*exp(-omega*time))
    update(efficacy_measured) <-((Ic[1]/(n_trial/2))-(Itx[1]/(n_trial/2)))/(Ic[1]/(n_trial/2))
    
    # ---- Transitions & Transition Rates ----
    
    # -- Population-level Transitions --
    n_SP[] <- Binomial(S[i], p_SP[i])
    n_PS[] <- Binomial(P[i], p_PS[i])
    
    n_SIu[] <- Binomial(S[i], p_SIu[i])
    n_S2Iu[] <- Binomial(S2[i], p_S2Iu[i])
    n_PIp[] <- Binomial(P[i], p_PIp[i])
    
    n_IuR[] <- Binomial(Iu[i], p_IR[i])
    n_IpR[] <- Binomial(Ip[i], p_IR[i])
    
    n_RS[] <- Binomial(R[i], p_RS[i])
    
    # -- Trial-level Transitions --
    n_ScIc[] <- if (time >= tau_start && time < tau_start + tau_end) Binomial(Sc[i], p_SIu[i]) else 0
    
    n_StxItx[] <- if (time >= tau_start && time < tau_start + tau_end) Binomial(Stx[i], p_StxItx[i]) else 0
    
    # -- Transition Rates --
    p_SP[] <- 1 - exp(-rho[i] * dt)
    p_PS[] <- 1 - exp(-omega * dt)
    
    p_SIu[] <- 1 - exp(-lambda[i] * dt)
    p_S2Iu[] <- 1 - exp(-lambda2[i] * dt)
    p_PIp[] <- 1 - exp(-lambda_p[i] * dt)
    p_StxItx[] <- 1 - exp(-lambda_tx[i]*dt)
    
    p_IR[] <- 1 - exp(-gamma * dt)
    p_RS[] <- 1 - exp(-delta * dt)
    
    # ---- User Defined Parameters ----
    
    # -- Population Parameters --
    start_S <- parameter()
    start_Iu <- parameter()
    start_Ip <- parameter()
    start_R <- parameter()
    start_P <- parameter()
    
    # -- Trial Parameters -- 
    start_Sc <- parameter()
    start_Stx <- parameter()
    #tau = trial duration = tau_end - tau_start
    tau_start <- parameter() #trial start day (model burn in period)
    tau_end <- parameter() #trial end day (trial duration + model burn in period)
    
    # -- Contact Matrix --
    contacts <- parameter()
    
    # -- Drug Parameters -- 
    chi <- parameter() #drug efficacy (initial)
    omega <- parameter() #duration of protection of drug
    rho <- parameter() #prophylactic coverage rate
    
    # -- Infection Parameters -- 
    beta <- parameter() #infection rate (natural)
    gamma <- parameter() #recovery rate
    delta <- parameter() #re-susceptibility rate
    immunity <- parameter() #protective impact of a previous infection
    p_ma <- parameter() #proportion of infections that are medically attended
    
    # -- Seasonal Forcing --
    offset <- parameter() #offset for seasonal forcing
    alpha <- parameter() #amplitude for seasonal forcing
    
    # ---- Dimensions ----
    
    # -- Variables -- 
    n_age_groups <- parameter()
    n_age_groups_trial <- parameter()
    
    # -- Dimensions --
    dim(S, S2, I, Iu, Ip, R, P, n_age, age_groups) <- n_age_groups #for age stratification
    dim(start_S, start_I, start_Iu, start_Ip, start_R, start_P) <- n_age_groups #for age stratification
    dim(p_SIu, p_S2Iu, p_SP,p_PS, p_PIp, p_IR, p_RS) <- n_age_groups
    dim(lambda, lambda2, lambda_p, lambda_check, rho) <- n_age_groups
    dim(n_SIu, n_S2Iu, n_SP, n_PS, n_PIp, n_IuR, n_IpR, n_RS) <- n_age_groups
    dim(Iu_ma, Ip_ma, I_ma, p_ma) <- n_age_groups
    dim(contacts, contacts_ij) <- c(n_age_groups, n_age_groups)
    dim(Sc, Ic, Stx, Itx) <- n_age_groups_trial
    dim(start_Sc, start_Stx) <- n_age_groups_trial
    dim(n_ScIc, n_StxItx, p_StxItx, lambda_tx) <- n_age_groups_trial
    dim(Ic_ma, Itx_ma) <- n_age_groups_trial
    dim(TotalCases_pop, TotalCases_trial, TotalCases_popMA, TotalCases_trialMA) <- n_age_groups
    
    # ---- Initial Conditions ----
    
    # -- Initial Compartment Values --
    initial(S[]) <- start_S[i]
    initial(S2[]) <- 0
    initial(Iu[]) <- start_Iu[i]
    initial(Ip[]) <- start_Ip[i]
    initial(R[]) <- start_R[i]
    initial(P[]) <- start_P[i]
    initial(I[]) <- start_I[i]
    
    # -- Initial Trial Compartment Values -- 
    initial(Sc[]) <- start_Sc[i]
    initial(Stx[]) <- start_Stx[i]
    initial(Ic[]) <- 0
    initial(Itx[]) <- 0
    
    # -- Medically Attended Cases --
    initial(Iu_ma[]) <- 0
    initial(Ip_ma[]) <- 0
    initial(I_ma[]) <- 0
    initial(Ic_ma[]) <- 0
    initial(Itx_ma[]) <- 0
    
    # -- Cumulative Cases --
    initial(TotalCases_pop[]) <- start_I[i]
    initial(TotalCases_trial[]) <- 0
    initial(TotalCases_popMA[]) <- 0
    initial(TotalCases_trialMA[]) <- 0
    
    # -- Population Size --
    initial(n_pop) <- sum(n_age)
    initial(n_trial) <- sum(start_Sc) + sum(start_Stx)
    
    # -- FoI Check --
    initial(lambda_check[]) <- 0
    initial(beta_check) <- 0
    
    # -- Trial Efficacy --
    initial(efficacy_measured) <- 0
    initial(efficacy_true) <- chi
    
    # ---- Calculations ----
    
    ## ---- Contact Matrix ----
    start_I[] <- start_Ip[i] + start_Iu[i]
    n_age[] <- start_S[i] + start_I[i] + start_P[i] + start_R[i]
    
    contacts_ij[, ] <- contacts[i, j] * (I[j]/n_age[j]) #overall contacts by age group
    
    ## ---- Seasonal Forcing ----
    seasonal_forcing <- 1 + alpha * cos(2*pi*(time-offset)/365)
    
    ## ---- Force of Infection ----
    lambda[] <- beta * seasonal_forcing * sum(contacts_ij[i, ])
    lambda2[] <- beta * immunity * seasonal_forcing * sum(contacts_ij[i, ]) #FoI for those in S2 (previously infected but not protected by drug)
    
    beta_p <- beta*(1-chi) #infection rate (protected)
    lambda_p[] <- beta_p * seasonal_forcing * sum(contacts_ij[i,])
    
    lambda_tx[] <- lambda[i]*(1-(chi*exp((-omega)*time)))
    
  })
