model_fit<-readRDS("./growth fits/vbgf-rpma-4-mec/rpma-4-vbgf-re.Rds")

#---VBGF model: extract parameter estimates

#---parameters on log scale ~MVN
#---Linf and k
exp(model_fit$report$mu) #
#---variance for Linf and k 
model_fit$report$sd 
#---correlation
model_fit$report$rho 
#---log normal observation error
model_fit$report$log_obs_er
#---t0 on real scale
model_fit$report$t0 


#---predictions for family lots fit as a random effect
#---linf
L_inf=exp(model_fit$report$ln_linf)
#---k
k=exp(model_fit$report$ln_k)
t0=model_fit$report$t0

#---data needed for model fits by family lot
tmp<-data.table(Age=model_fit$data$Age,
    ln_length=model_fit$data$ln_length,
    family_lot=model_fit$data$family_lot)
tmp[,Length:=exp(ln_length)]

library(lattice)
#---plot model fits to data for each family lot
xyplot(Length ~ Age | as.factor(family_lot),
       data = tmp,
       as.table=T,
       panel = function(x,y) {
           # level <- dimnames(trellis.last.object())[["family_lot"]][packet.number()]
           panel.xyplot(x,y)
           #panel.abline(h=L_inf[panel.number()])
           panel.curve(L_inf[panel.number()]* (1 - exp(-k[panel.number()] * (x-t0))),
            from=0, to=max(tmp$Age), n = 101)
           },
       main="von Bertalanffy growth model",
       xlab="Age (years)",
       ylab="Length (mm)")
       
 
