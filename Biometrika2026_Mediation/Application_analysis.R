suppressMessages(library(KernSmooth))
suppressMessages(library(riskCommunicator))

df       <- read.csv("framingham_sub.csv")

Y <- df$HEARTRTE; X <- df$CIGPDAY; M <- df$BMI; N <- nrow(df)
ord_lvls <- c(20, 30, 40)

ord_by   <- 1
ord_lvls <- seq(0, 90, by = ord_by)   # 0,5,10,...,90  -> includes 20,30,40

# ---------- model fitting ----------------------------------------------------
silverman <- function(z) 1.06*min(sd(z),IQR(z)/1.34)*length(z)^(-0.2)
ll1 <- function(Ztr,Ytr,z0,bw){
  d<-sweep(Ztr,2,z0,"-"); w<-exp(-0.5*rowSums(sweep(d,2,bw,"/")^2))
  if(sum(w)<1e-300) return(weighted.mean(Ytr,rep(1,length(Ytr))))
  X2<-cbind(1,d); Xw<-X2*w
  b<-tryCatch(solve(crossprod(Xw,X2),crossprod(Xw,Ytr)),
              error=function(e) c(weighted.mean(Ytr,w),rep(0,ncol(d))))
  b[1] }
llp <- function(Ztr,Ytr,Znew,bw){
  if(is.null(dim(Znew))) Znew<-matrix(Znew,ncol=length(bw))
  if(is.null(dim(Ztr)))  Ztr <-matrix(Ztr, ncol=length(bw))
  apply(Znew,1,function(z0) ll1(Ztr,Ytr,z0,bw)) }

bwX<-silverman(X)*4; bwM<-silverman(M)*4
bwY<-c(bwX,bwM); bwMm<-c(bwX)
ZYt<-cbind(X,M); ZMt<-matrix(X,ncol=1)
sigM<-sd(M-llp(ZMt,M,ZMt,bwMm))

mkmod <- function(ZY,Yv,ZM,Mv,sM) list(
  pY=function(x,m) llp(ZY,Yv,cbind(x,m),bwY),
  pM=function(x)   llp(ZM,Mv,matrix(x,ncol=1),bwMm), sM=sM)
FM <- mkmod(ZYt,Y,ZMt,M,sigM)

# ---------- identification ---------------------------------------------------
EYxMxs <- function(x,xs,md,nq=50){
  phi<-md$pM(xs); mg<-seq(phi-4*md$sM,phi+4*md$sM,length.out=nq)
  sum(md$pY(x,mg)*dnorm(mg,phi,md$sM))*(mg[2]-mg[1]) }
EYx <- function(x,md,nq=50) EYxMxs(x,x,md,nq)

dthdx<-function(x,m,md,h=0.5) (md$pY(x+h,m)-md$pY(x-h,m))/(2*h)
dthdm<-function(x,m,md,h=0.2) (md$pY(x,m+h)-md$pY(x,m-h))/(2*h)
dPdx <-function(x,m,md,h=0.5){
  phi<-md$pM(x); pp<-(md$pM(x+h)-md$pM(x-h))/(2*h)
  dnorm((m-phi)/md$sM)/md$sM*pp }
LNDE<-function(x,md,nq=40){
  phi<-md$pM(x); mg<-seq(phi-4*md$sM,phi+4*md$sM,length.out=nq)
  sum(dthdx(x,mg,md)*dnorm(mg,phi,md$sM))*(mg[2]-mg[1]) }
LNIE<-function(x,md,nq=40){
  phi<-md$pM(x); mg<-seq(phi-4*md$sM,phi+4*md$sM,length.out=nq)
  sum(dthdm(x,mg,md)*dPdx(x,mg,md))*(mg[2]-mg[1]) }
trap<-function(vs,x1,x2,K){ dx<-(x2-x1)/(K-1); (sum(vs)-0.5*(vs[1]+vs[K]))*dx }

# ---------- effect functions -------------------------------------------------
TE   <-function(x2,x1,md,nq=50) EYx(x2,md,nq)-EYx(x1,md,nq)

# NDE(x'',x') Def2 and NDE(x',x'') swapped — both needed (not skew-sym)
NDE  <-function(x2,x1,md,nq=50) EYxMxs(x2,x1,md,nq)-EYx(x1,md,nq)

# NIE(x'',x') Def2:    E[Y_{x',M_{x''}}] - E[Y_{x'}]
# NIE(x',x'') swapped: E[Y_{x'',M_{x'}}] - E[Y_{x''}]  — both needed (not skew-sym)
NIE  <-function(x2,x1,md,nq=50) EYxMxs(x1,x2,md,nq)-EYx(x1,md,nq)   # Def2
NIEsw<-function(x2,x1,md,nq=50) EYxMxs(x2,x1,md,nq)-EYx(x2,md,nq)   # swapped

# CNDE/CNIE — skew-symmetric, one form only (continuous treatment, via local effects)
CNDE <-function(x2,x1,md,K=20,nq=40){
  xg<-seq(x1,x2,length.out=K); trap(sapply(xg,LNDE,md=md,nq=nq),x1,x2,K) }
CNIE <-function(x2,x1,md,K=20,nq=40){
  xg<-seq(x1,x2,length.out=K); trap(sapply(xg,LNIE,md=md,nq=nq),x1,x2,K) }

# S-NDE/S-NIE — skew-symmetric, one form only
SNDE <-function(x2,x1,md,nq=50) 0.5*(NDE(x2,x1,md,nq)-NDE(x1,x2,md,nq))
SNIE <-function(x2,x1,md,nq=50) 0.5*(NIE(x2,x1,md,nq)-NIE(x1,x2,md,nq))

step_effects <- function(md, ls, nq=40){
  K <- length(ls) - 1
  NDE_up    <- numeric(K); NDE_down   <- numeric(K)
  NIE_up    <- numeric(K); NIE_down   <- numeric(K)
  NIEsw_up  <- numeric(K); NIEsw_down <- numeric(K)
  SNDE_up   <- numeric(K); SNIE_up    <- numeric(K)
  for(k in seq_len(K)){
    a <- ls[k]; b <- ls[k+1]              # step from a -> b  (a < b)
    NDE_up[k]    <- NDE(b,a,md,nq)        # NDE(b,a)
    NDE_down[k]  <- NDE(a,b,md,nq)        # NDE(a,b)
    NIE_up[k]    <- NIE(b,a,md,nq)        # NIE(b,a)
    NIE_down[k]  <- NIE(a,b,md,nq)        # NIE(a,b)
    NIEsw_up[k]  <- NIEsw(b,a,md,nq)      # NIEsw(b,a)
    NIEsw_down[k]<- NIEsw(a,b,md,nq)      # NIEsw(a,b)
    SNDE_up[k]   <- SNDE(b,a,md,nq)       # = -SNDE(a,b) by skew-symmetry
    SNIE_up[k]   <- SNIE(b,a,md,nq)       # = -SNIE(a,b) by skew-symmetry
  }
  # prefix sums, prefix[m] = sum_{k=1}^{m-1} step[k], length K+1, prefix[1]=0
  pre <- function(v) c(0, cumsum(v))
  list(
    pNDE_up    = pre(NDE_up),    pNDE_down   = pre(NDE_down),
    pNIE_up    = pre(NIE_up),    pNIE_down   = pre(NIE_down),
    pNIEsw_up  = pre(NIEsw_up),  pNIEsw_down = pre(NIEsw_down),
    pSNDE_up   = pre(SNDE_up),   pSNDE_down  = pre(-SNDE_up),
    pSNIE_up   = pre(SNIE_up),   pSNIE_down  = pre(-SNIE_up)
  )
}

rangeSumO <- function(a,b,ls,prefix_up,prefix_down){
  ia<-which(ls==a); ib<-which(ls==b)
  if(ia>ib) prefix_up[ia]-prefix_up[ib]      # a>b : sum over steps ib..ia-1 of "up"
  else      prefix_down[ib]-prefix_down[ia]  # a<b : sum over steps ia..ib-1 of "down"
}

cumNDE_O    <-function(x2,x1,se,ls) rangeSumO(x2,x1,ls,se$pNDE_up,  se$pNDE_down)  # CNDE-O(x2,x1)
cumNDE_O_sw <-function(x2,x1,se,ls) rangeSumO(x1,x2,ls,se$pNDE_up,  se$pNDE_down)  # CNDE-O(x1,x2)
cumNIE_O    <-function(x2,x1,se,ls) rangeSumO(x2,x1,ls,se$pNIE_up,  se$pNIE_down)  # CNIE-O(x2,x1)
cumNIE_O_sw <-function(x2,x1,se,ls) rangeSumO(x2,x1,ls,se$pNIEsw_up,se$pNIEsw_down) # CNIE-O(x1,x2)
cumSNDE_O   <-function(x2,x1,se,ls) rangeSumO(x2,x1,ls,se$pSNDE_up, se$pSNDE_down) # S-CNDE-O(x2,x1)
cumSNIE_O   <-function(x2,x1,se,ls) rangeSumO(x2,x1,ls,se$pSNIE_up, se$pSNIE_down) # S-CNIE-O(x2,x1)

# ---------- all measures for one contrast ------------------------------------
all_meas<-function(x2,x1,md,ls,se,K=20,nq=40){
  in_l<-(x2%in%ls)&&(x1%in%ls)
  c(TE       = TE(x2,x1,md,nq),
    NDE      = NDE(x2,x1,md,nq),       # Def2
    NDE_sw   = NDE(x1,x2,md,nq),       # swapped
    NIE      = NIE(x2,x1,md,nq),       # Def2
    NIE_sw   = NIEsw(x2,x1,md,nq),     # swapped
    CNDE     = CNDE(x2,x1,md,K,nq),
    CNIE     = CNIE(x2,x1,md,K,nq),
    SNDE     = SNDE(x2,x1,md,nq),
    SNIE     = SNIE(x2,x1,md,nq),
    CNDE_O   = if(in_l) cumNDE_O(x2,x1,se,ls)    else NA,
    CNDE_O_sw= if(in_l) cumNDE_O_sw(x2,x1,se,ls) else NA,
    CNIE_O   = if(in_l) cumNIE_O(x2,x1,se,ls)    else NA,
    CNIE_O_sw= if(in_l) cumNIE_O_sw(x2,x1,se,ls) else NA,
    SCNDE_O  = if(in_l) cumSNDE_O(x2,x1,se,ls)   else NA,
    SCNIE_O  = if(in_l) cumSNIE_O(x2,x1,se,ls)   else NA) }

# =============================================================================
contrasts <- list(
  list(x2=40,x1=20,lab="From 20 to 40"),
  list(x2=30,x1=20,lab="From 20 to 30"),
  list(x2=40,x1=30,lab="From 30 to 40"))

cat("Computing point-estimate step effects on ordinal grid...\n")
SE_pt <- step_effects(FM, ord_lvls, nq=40)

cat("Computing point estimates...\n")
pt <- lapply(contrasts,function(cc){
  cat(sprintf("  %s ...\n",cc$lab))
  c(list(lab=cc$lab,x2=cc$x2,x1=cc$x1),
    as.list(all_meas(cc$x2,cc$x1,FM,ord_lvls,SE_pt))) })

# =============================================================================
nms <- c("TE","NDE","NDE_sw","NIE","NIE_sw","CNDE","CNIE","SNDE","SNIE",
         "CNDE_O","CNDE_O_sw","CNIE_O","CNIE_O_sw","SCNDE_O","SCNIE_O")
B   <- 100
bmt <- lapply(contrasts,function(cc)
  matrix(NA,nrow=B,ncol=length(nms),dimnames=list(NULL,nms)))

cat(sprintf("\nBootstrapping (B=%d)...\n",B))
set.seed(1)
for(b in seq_len(B)){
  if(b%%10==0) cat(sprintf("  rep %d/%d\n",b,B))
  idx<-sample(N,N,replace=TRUE)
  sMb<-sd(M[idx]-llp(matrix(X[idx],ncol=1),M[idx],matrix(X[idx],ncol=1),bwMm))
  mb<-mkmod(cbind(X[idx],M[idx]),Y[idx],matrix(X[idx],ncol=1),M[idx],sMb)
  SE_b <- tryCatch(step_effects(mb, ord_lvls, nq=30), error=function(e) NULL)
  for(ci in seq_along(contrasts)){
    cc<-contrasts[[ci]]
    tryCatch({bmt[[ci]][b,]<-all_meas(cc$x2,cc$x1,mb,ord_lvls,SE_b,K=15,nq=30)},
             error=function(e) NULL) } }

# =============================================================================
# PLAIN TEXT TABLE
# =============================================================================
ci95<-function(bm,nm){v<-bm[!is.na(bm[,nm]),nm];c(quantile(v,.025),quantile(v,.975))}
pfmt<-function(mu,lo,hi) sprintf("%6.3f  [%6.3f, %6.3f]",mu,lo,hi)
SEP <-paste(rep("-",66),collapse="")
DOT <-paste(rep("\u00b7",66),collapse="")
DAS <-paste(rep(".",66),collapse="")

print_text_block<-function(ci,p,bm){
  cc<-contrasts[[ci]]; x2<-cc$x2; x1<-cc$x1
  bm<-bm[!is.na(bm[,"TE"]),,drop=FALSE]
  cat(sprintf("\n(%d). %s\n",ci,p$lab))
  cat(sprintf("  %-26s  %s\n","Estimator","Mean      95% CI"))
  cat(SEP,"\n")
  row<-function(lbl,nm){
    mu<-p[[nm]]; if(is.na(mu)){cat(sprintf("  %-26s  \u2014\n",lbl));return()}
    ci9<-ci95(bm,nm)
    cat(sprintf("  %-26s  %s\n",lbl,pfmt(mu,ci9[1],ci9[2]))) }
  
  row(sprintf("TE(%g,%g)",x2,x1),"TE"); cat(DOT,"\n")
  row(sprintf("NDE(%g,%g)",x2,x1),"NDE")
  row(sprintf("NDE(%g,%g)",x1,x2),"NDE_sw")
  cat(DOT,"\n")
  row(sprintf("NIE(%g,%g)",x2,x1),"NIE")
  row(sprintf("NIE(%g,%g)",x1,x2),"NIE_sw")
  cat(DOT,"\n")
  row(sprintf("CNDE(%g,%g)",x2,x1),"CNDE")
  row(sprintf("CNIE(%g,%g)",x2,x1),"CNIE")
  cat(DOT,"\n")
  row(sprintf("S-NDE(%g,%g)",x2,x1),"SNDE")
  row(sprintf("S-NIE(%g,%g)",x2,x1),"SNIE")
  if(!is.na(p$CNDE_O)){
    cat(DAS,"\n")
    row(sprintf("CNDE-O(%g,%g)",x2,x1),"CNDE_O")
    row(sprintf("CNDE-O(%g,%g)",x1,x2),"CNDE_O_sw")
    cat(DOT,"\n")
    row(sprintf("CNIE-O(%g,%g)",x2,x1),"CNIE_O")
    row(sprintf("CNIE-O(%g,%g)",x1,x2),"CNIE_O_sw")
    cat(DOT,"\n")
    row(sprintf("S-CNDE-O(%g,%g)",x2,x1),"SCNDE_O")
    row(sprintf("S-CNIE-O(%g,%g)",x2,x1),"SCNIE_O")
  }
  cat(SEP,"\n")
}

cat("\n"); cat(paste(rep("=",66),collapse=""),"\n")
cat("Table 2. Means and 95% CIs\n")
cat("X=cigarettes/day  M=BMI  Y=heart rate  (n=224)\n")
cat(sprintf("Ordinal grid for *-O measures: {%s} (by=%d)\n",
            paste(range(ord_lvls),collapse=","), ord_by))
cat(paste(rep("=",66),collapse=""),"\n")
for(ci in seq_along(contrasts)) print_text_block(ci,pt[[ci]],bmt[[ci]])

# =============================================================================
# LATEX TABLE
# =============================================================================
lfmt<-function(mu,lo,hi) sprintf("$%.3f\\;[%.3f,\\;%.3f]$",mu,lo,hi)

latex_block<-function(ci,p,bm){
  cc<-contrasts[[ci]]; x2<-cc$x2; x1<-cc$x1
  bm<-bm[!is.na(bm[,"TE"]),,drop=FALSE]
  lines<-character(0)
  add<-function(...) lines<<-c(lines,paste0(...))
  row<-function(lbl,nm){
    mu<-p[[nm]]; if(is.na(mu)) return()
    ci9<-ci95(bm,nm)
    add(lbl," & ",lfmt(mu,ci9[1],ci9[2])," \\\\") }
  
  add("\\multicolumn{2}{l}{\\textbf{(",ci,"). ",cc$lab,"}} \\\\")
  add("\\hline\\hline")
  
  row(sprintf("$\\mathrm{TE}(%g,%g)$",x2,x1),"TE")
  add("\\hline")
  
  row(sprintf("$\\mathrm{NDE}(%g,%g)$",x2,x1),"NDE")
  row(sprintf("$\\mathrm{NDE}(%g,%g)$",x1,x2),"NDE_sw")
  row(sprintf("$\\mathrm{NIE}(%g,%g)$",x2,x1),"NIE")
  row(sprintf("$\\mathrm{NIE}(%g,%g)$",x1,x2),"NIE_sw")
  add("\\hline")
  
  row(sprintf("$\\mathrm{CNDE}(%g,%g)$",x2,x1),"CNDE")
  row(sprintf("$\\mathrm{CNIE}(%g,%g)$",x2,x1),"CNIE")
  add("\\hline")
  
  row(sprintf("$\\mathrm{S\\text{-}NDE}(%g,%g)$",x2,x1),"SNDE")
  row(sprintf("$\\mathrm{S\\text{-}NIE}(%g,%g)$",x2,x1),"SNIE")
  add("\\hline")
  
  if(!is.na(p$CNDE_O)){
    row(sprintf("$\\mathrm{CNDE\\text{-}O}(%g,%g)$",x2,x1),"CNDE_O")
    row(sprintf("$\\mathrm{CNDE\\text{-}O}(%g,%g)$",x1,x2),"CNDE_O_sw")
    row(sprintf("$\\mathrm{CNIE\\text{-}O}(%g,%g)$",x2,x1),"CNIE_O")
    row(sprintf("$\\mathrm{CNIE\\text{-}O}(%g,%g)$",x1,x2),"CNIE_O_sw")
    add("\\hline")
    row(sprintf("$\\mathrm{S\\text{-}CNDE\\text{-}O}(%g,%g)$",x2,x1),"SCNDE_O")
    row(sprintf("$\\mathrm{S\\text{-}CNIE\\text{-}O}(%g,%g)$",x2,x1),"SCNIE_O")
    add("\\hline")
  }
  add("\\hline")
  lines }

cat("\n\n% =====================================================================\n")
cat("% LaTeX Table 2\n")
cat("% =====================================================================\n")
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Means and 95\\% CIs. ")
cat("Non-skew-symmetric measures (NDE, NIE, CNDE-O, CNIE-O) shown in both ")
cat("$(x'',x')$ and $(x',x'')$ forms. ")
cat("Skew-symmetric measures (CNDE, CNIE, S-NDE, S-NIE, S-CNDE-O, S-CNIE-O) ")
cat("shown in $(x'',x')$ form only. CNDE-O/CNIE-O/S-CNDE-O/S-CNIE-O are ")
cat(sprintf("accumulated over the ordinal grid \\{0,%d,...,90\\}.}\n", ord_by))
cat("\\label{tab:table2}\n")
cat("\\begin{tabular}{lc}\n\\hline\n")
cat("Estimator & Mean $[95\\%\\;\\mathrm{CI}]$ \\\\\n")
for(ci in seq_along(contrasts)){
  lns<-latex_block(ci,pt[[ci]],bmt[[ci]])
  for(l in lns) cat(l,"\n")
  cat("\n") }
cat("\\end{tabular}\n\\end{table}\n")