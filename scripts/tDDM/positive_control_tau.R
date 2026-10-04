#!/usr/bin/env Rscript
# ==============================================================================
# POSITIVE-CONTROL τ RECOVERY  (JOB 1 — the missing half of the null-τ sim)
# ==============================================================================
# Question the null-τ sim did NOT answer: is the time-contrast estimator
#   SENSITIVE to a TRUE group onset gap? The null-sim proved BIAS (true Δτ=0 ->
#   Δτ̂≈-0.19). This proves/refutes SENSITIVITY (does Δτ̂ grow when true gap grows?).
#
# Design (mirrors nulltau_generative_sim.R exactly, engine verbatim):
#   Convention: gap = τ_DEP − τ_HC  (positive = DEP later, "depression delays onset").
#   HC subjects   -> tau_true = TAU_BASE                (level-independent)
#   DEP subjects  -> tau_true = TAU_BASE + gap          (level-dependent)
#   gap sweep ∈ {0, 0.085, 0.169, 0.34}  (0, ½×, 1×, 2× the observed 0.169 s gap)
#   Simulate trials with each subject's own signals -> refit 4j time-contrast ->
#   recovered_gap = mean(τ̂|DEP) − mean(τ̂|HC).  Repeat over reps.
#   gap=0 level reproduces the null-sim EXACTLY (sanity anchor; expect ≈ +0.19).
#
#   Efficiency: HC fits are identical across gap levels -> fit HC ONCE per rep,
#   reuse across levels; only DEP is re-fit per level.
#
# Decision rule (slope = d(recovered_gap)/d(true_gap)):
#   slope ≈ 0  -> estimator IGNORES true τ, recovered saturated at bias -> WORLD A
#                 (pure artifact) -> discard τ, robust multi-axis-slowing paper.
#   slope ≈ 1  -> recovered ≈ true + (+0.19 offset) -> SENSITIVE, bias is additive
#                 -> WORLD B possible (entangled-real) -> proceed to JOB 3.
#   0<slope<1  -> partial sensitivity -> cautious World B.
#
# Usage: Rscript scripts/tDDM/positive_control_tau.R --reps 20
# ==============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga2 <- function(f,d){ if(f %in% args){i<-which(args==f); if(i<length(args)) return(args[i+1])}; d }
ga <- function(f,d){ if(f %in% args){i<-which(args==f); if(i<length(args)) return(as.numeric(args[i+1]))}; d }
N_REP   <- ga("--reps", 20)
DE_ITER <- ga("--de-iter", 60)
DE_NP   <- ga("--de-np", 40)
# 2026-09-23: the canonical fit now keeps the BEST-NLL run of R repeats, because
# 4j does not converge in 200 generations (docs/PAPER/Memo_Analyses_0919.md S4).
# The control has to use the identical procedure or its bias is not comparable.
RESTARTS <- as.integer(ga("--restarts", 1))
sd <- tryCatch(dirname(sys.frame(1)$ofile), error=function(e){
  ca<-commandArgs(FALSE); fa<-sub("^--file=","",ca[grep("^--file=",ca)]); if(length(fa))dirname(fa) else getwd()})
source(file.path(sd,"..","..","config.R"))
NC<-ga("--cores", max(1,round((detectCores()-2)*2/3)))  # default 2/3 load to limit CPU heat
GAPS <- as.numeric(strsplit(ga2("--gaps","0,0.085,0.169,0.34"),",")[[1]])

# ---- engine (SHARED) --------------------------------------------------------
TRANS_MODE <- ga2("--trans", "image")  # default was "pipeline" until 2026-09-27
cat(sprintf("transition density: %s\n", TRANS_MODE))
# 2026-08-10: replaced this script's private copies of the Navarro-Fuss density,
# the split-at-tau transition density and the piecewise likelihood with the one
# shared definition, and moved the x-grid 25 -> 30 (TDDM_CANON) so every script
# integrates the split identically. Regression: scripts/tDDM/test_engine_equivalence.R
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
SIGMA <- TDDM_SIGMA; EPS <- TDDM_EPS; PEN <- TDDM_PENALTY
NG    <- as.integer(ga("--ng", TDDM_CANON$ngrid)); TDDM_NGRID <- NG
# mk_ll_tau() reads td$RT_ddm / td$choice_ddm, so fit1() below builds those
# column names (this script used to call them RT / ch).
ll4j <- mk_ll_tau("s1", "s2", trans = get_transition(TRANS_MODE))

# 2026-09-27: body replaced by the shared bridge-corrected simulator (models/tDDM/tddm_engine.R,
# tddm_simulate_trial). The old Euler copy lengthened mean DT by 17-39 ms vs the likelihood.
sim <- function(d1, d2, th, nDT, tau, s1, s2, sg = SIGMA, dt = .001, mx = 10)
  tddm_simulate_trial(d1, d2, th, nDT, tau, s1, s2, sg = sg, dt = dt, mx = mx)

# ---- data + real fits (group differences in d/thres/nDT PRESERVED) ----
cat("=== POSITIVE-CONTROL τ RECOVERY (sweep true gap) ===\n")
cat(sprintf("reps=%d/level | gaps={%s} | DE itermax=%d NP=%d | restarts=%d | cores=%d\n\n",
            N_REP,paste(GAPS,collapse=","),DE_ITER,DE_NP,RESTARTS,NC))
files<-list.files(DATA_DIR,pattern="^study4_.*\\.csv$",full.names=TRUE)
dat<-do.call(rbind,lapply(files,read.csv,stringsAsFactors=FALSE))
if("sub"%in%names(dat))dat$subject<-dat$sub
dat<-dat[is.finite(dat$choice_RT)&dat$choice_RT>0.15&dat$choice_RT<=10,]
dat$ms<-(dat$RPE_scale+dat$EPE_X1_scale)/2; dat$cs<-(dat$RPE_scale-dat$EPE_X1_scale)/2
# --fits: generating parameters.  The default reproduces the historical run; pass
# the current canonical fits when the question is about the CURRENT model.
FITS_PATH<-ga2("--fits","results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
fits<-read.csv(if(grepl("^/",FITS_PATH)) FITS_PATH else file.path(PROJECT_ROOT,FITS_PATH),stringsAsFactors=FALSE)
# The baseline onset must sit at the OPERATING POINT of the density being tested:
# the shipped fits give tau ~ 1.80, but the corrected density puts the real-data
# tau at ~0.72, and the induced bias depends on how many trials reach phase 2.
TAU_BASE<-as.numeric(ga2("--tau-base", NA_real_))
if(!is.finite(TAU_BASE)){
  if(TRANS_MODE=="image") stop("--tau-base is REQUIRED with --trans image.\n",
    "  The baseline tau must sit at the OPERATING POINT of the density under test.\n",
    "  This fits file was produced with the SHIPPED density (mean tau ",sprintf("%.3f",mean(fits$tau)),"),\n",
    "  but the corrected density puts the real-data tau near 0.72. Running at 1.80 gave\n",
    "  slope -0.408; rerunning at 0.72 gave +0.731. Pass e.g. --tau-base 0.720.")
  TAU_BASE<-mean(fits$tau)   # shipped-density default, only for --trans pipeline
}
cat(sprintf("observed Δτ (DEP-HC) = %+.3f | TAU_BASE (HC) = %.3f\n\n",
            mean(fits$tau[fits$group=="DEP"])-mean(fits$tau[fits$group=="HC"]),TAU_BASE))
subj<-list()
for(i in 1:nrow(fits)){sid<-fits$subject[i];tr<-dat[dat$subject==sid,c("ms","cs")];tr<-tr[complete.cases(tr),]
 if(nrow(tr)<5)next
 subj[[length(subj)+1]]<-list(d1=fits$d_mean[i],d2=fits$d_contrast[i],th=fits$thres[i],nDT=fits$nDT[i],
   grp=fits$group[i],s1=tr$ms,s2=tr$cs)}
NTRIALS<-as.integer(ga2("--ntrials", 0))
if(NTRIALS>0) subj<-lapply(subj,function(s){i<-sample(seq_along(s$s1),NTRIALS,replace=TRUE)
  s$s1<-s$s1[i]; s$s2<-s$s2[i]; s})
FRAC<-as.numeric(ga2("--nsub-frac", 1))
if(FRAC<1){set.seed(99); keep<-sample(seq_along(subj), round(FRAC*length(subj))); subj<-subj[keep]}
cat(sprintf("trials/subject = %s | subject fraction = %.2f | TAU_BASE = %.3f\n",
            if(NTRIALS>0) as.character(NTRIALS) else "real", FRAC, TAU_BASE))
subj_HC <-Filter(function(s)s$grp=="HC", subj)
subj_DEP<-Filter(function(s)s$grp=="DEP",subj)
cat(sprintf("subjects: HC=%d DEP=%d (total=%d)\n\n",length(subj_HC),length(subj_DEP),length(subj)))

# returns recovered tau_hat for one subject simulated at tau_true
fit1<-function(s,tau_true){n<-length(s$s1);rt<-numeric(n);ch<-integer(n)
 for(j in 1:n){r<-sim(s$d1,s$d2,s$th,s$nDT,tau_true,s$s1[j],s$s2[j]);rt[j]<-r$rt;ch[j]<-r$ch}
 td<-data.frame(s1=s$s1,s2=s$s2,RT_ddm=rt,choice_ddm=ch)
 td<-td[td$RT_ddm>0.15&td$RT_ddm<=10,];if(nrow(td)<5)return(NA_real_)
  f<-NULL
 for(k in seq_len(RESTARTS)){
   fk<-tryCatch(DEoptim(ll4j,c(-3,-3,0.5,0.01,-2),c(3,3,4,3,3),
     DEoptim.control(itermax=DE_ITER,NP=DE_NP,trace=FALSE),td=td),error=function(e)NULL)
   if(!is.null(fk)&&(is.null(f)||fk$optim$bestval<f$optim$bestval)) f<-fk }
 if(is.null(f))return(NA_real_);max(f$optim$bestmem[5],0)}  # 2026-09-28: tau<=0 clamped at 0 (= contrast model)

RNGkind("L'Ecuyer-CMRG");set.seed(20260715);ts<-format(Sys.time(),"%Y%m%d_%H%M%S")
od<-file.path(RESULTS_DIR,"tDDM",sprintf("positive_control_tau_%s",ts));dir.create(od,recursive=TRUE,showWarnings=FALSE)
pf<-file.path(od,"per_point.csv");t0<-Sys.time();rows<-list()
for(rep in 1:N_REP){
 # HC fit once per rep (level-independent: tau_true = TAU_BASE for all levels)
 # 2026-09-28: mclapply seeds each worker from the parent's L'Ecuyer stream WITHOUT advancing it,
 # so every call needs its own set.seed(); otherwise all reps reproduce the same simulated data.
 set.seed(20260715+1000*rep)
 hc<-unlist(mclapply(subj_HC,function(s)fit1(s,TAU_BASE),mc.cores=NC,mc.preschedule=FALSE))
 hc<-hc[is.finite(hc)];if(length(hc)<10){cat(sprintf("rep %d HC SKIP\n",rep));next}
 tHC<-mean(hc)
 for(gap in GAPS){
  set.seed(20260715+1000*rep+match(gap,GAPS))   # distinct stream per rep x gap, and from the HC call
  dep<-unlist(mclapply(subj_DEP,function(s)fit1(s,TAU_BASE+gap),mc.cores=NC,mc.preschedule=FALSE))
  dep<-dep[is.finite(dep)];if(length(dep)<10){cat(sprintf("rep %d gap %.3f DEP SKIP\n",rep,gap));next}
  tDEP<-mean(dep)
  rec<-tDEP-tHC   # recovered gap in DEP-HC convention
  row<-data.frame(rep=rep,true_gap=gap,recovered_gap=rec,tauhat_HC=tHC,tauhat_DEP=tDEP,
    n_HC=length(hc),n_DEP=length(dep))
  rows[[length(rows)+1]]<-row;h<-!file.exists(pf);write.table(row,pf,sep=",",row.names=FALSE,col.names=h,append=!h)
  el<-as.numeric(difftime(Sys.time(),t0,units="mins"))
  cat(sprintf("  rep %d/%d gap=%+.3f: recovered Δτ̂=%+.3f (τ̂_HC=%.3f τ̂_DEP=%.3f) [%.1f min]\n",
              rep,N_REP,gap,rec,tHC,tDEP,el));flush.console()}}
res<-do.call(rbind,rows)

# ---- slope / intercept + per-level summary ----
lm.fit<-lm(recovered_gap~true_gap,data=res)
co<-summary(lm.fit)$coefficients
slope<-co["true_gap","Estimate"];slope_se<-co["true_gap","Std. Error"]
intc<-co["(Intercept)","Estimate"];intc_se<-co["(Intercept)","Std. Error"]
slope_ci<-slope+c(-1.96,1.96)*slope_se
lev<-aggregate(recovered_gap~true_gap,res,function(x)c(m=mean(x),sd=sd(x),n=length(x)))
levdf<-data.frame(true_gap=lev$true_gap,mean_recovered=lev$recovered_gap[,"m"],
                  sd=lev$recovered_gap[,"sd"],n=lev$recovered_gap[,"n"])
cat("\n============================================================\n")
cat("PER-LEVEL recovered Δτ̂ (DEP-HC convention):\n")
for(i in 1:nrow(levdf)) cat(sprintf("  true gap %+.3f -> recovered %+.3f (sd=%.3f, n=%d)\n",
   levdf$true_gap[i],levdf$mean_recovered[i],levdf$sd[i],levdf$n[i]))
cat(sprintf("\nSLOPE d(Δτ̂)/d(true gap) = %.3f  [95%% CI %.3f, %.3f]  (SE=%.3f)\n",slope,slope_ci[1],slope_ci[2],slope_se))
cat(sprintf("INTERCEPT (bias at true gap=0) = %+.3f (SE=%.3f)  <- should ≈ null-sim +0.19\n",intc,intc_se))
cat("\nDECISION:\n")
# 2026-08-10: with a single rep per gap level the slope SE (and hence the CI) is
# NaN, and this branch died with "missing value where TRUE/FALSE needed" -- the
# whole run was lost at its very last step, after the fits had already finished.
if(!is.finite(slope_ci[2])){
  cat("  slope CI undefined (needs >= 2 reps per gap level) -- no verdict\n")
}else if(slope_ci[2] < 0.25){
 cat("  slope ≈ 0 (CI upper < 0.25) -> estimator INSENSITIVE to true τ\n")
 cat("  -> WORLD A (pure artifact). Discard τ; robust multi-axis-slowing paper.\n")
}else if(slope_ci[1] > 0.6){
 cat("  slope substantial (CI lower > 0.6) -> estimator SENSITIVE to true τ\n")
 cat("  -> WORLD B possible (entangled-real). Proceed to JOB 3 (latent-factor model).\n")
}else{
 cat("  slope INTERMEDIATE / uncertain -> partial sensitivity; interpret cautiously,\n")
 cat("     consider more reps or wider gap sweep before committing to World A vs B.\n")
}
write.csv(res,file.path(od,"per_point.csv"),row.names=FALSE)
write.csv(levdf,file.path(od,"summary_by_level.csv"),row.names=FALSE)
write.csv(data.frame(slope=slope,slope_se=slope_se,slope_lo=slope_ci[1],slope_hi=slope_ci[2],
  intercept=intc,intercept_se=intc_se,n_points=nrow(res),reps=N_REP),
  file.path(od,"slope_summary.csv"),row.names=FALSE)
cat(sprintf("\nSaved: %s\nTotal: %.1f min\n",od,as.numeric(difftime(Sys.time(),t0,units="mins"))))
