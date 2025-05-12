
library("R.matlab")
library("gamlss")

graphics.off()

# fi_name <- 'dev_homogeneity_homo_mean_network_weight'
string_list <- c('dev_fc_homogeneity_large')
#  
for (str in string_list) {
  fi_name <- str
for (neti in 2:18) {



testdata <- readMat(paste0("/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/curv_info/", fi_name, ".mat"))
data0 <- testdata[["d.curv"]]

idp_name <- paste0( "rs_global_fc_", neti )
data1 <- data0[[neti]]
colnames(data1) <- c("phenotype", "Age", "Sex", "meanFD", "ScannerSite")
plotdata = as.data.frame(data1)

plotdata$Sex <- as.factor(plotdata$Sex)
plotdata$ScannerSite <- as.factor(plotdata$ScannerSite)


con<-gamlss.control(n.cyc=100)

mdl_bs1<-gamlss(phenotype~bs(Age,df=3) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=3) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs2<-gamlss(phenotype~bs(Age,df=4) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=3) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs3<-gamlss(phenotype~bs(Age,df=5) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=3) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs4<-gamlss(phenotype~bs(Age,df=6) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=3) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs5<-gamlss(phenotype~bs(Age,df=4) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=4) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs6<-gamlss(phenotype~bs(Age,df=5) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=4) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs7<-gamlss(phenotype~bs(Age,df=6) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=4) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs8<-gamlss(phenotype~bs(Age,df=5) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=5) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs9<-gamlss(phenotype~bs(Age,df=6) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=5) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)
mdl_bs10<-gamlss(phenotype~bs(Age,df=6) + Sex + meanFD + random(ScannerSite),
                sigma.fo=~bs(Age,df=6) + Sex + meanFD + random(ScannerSite),
                nu.fo=~1,
                nu.tau=~1,
                family=JSU,
                data=plotdata,
                control=con)


conv<-matrix(data=0,nrow=1,ncol=10)
for (i in 1:10) {
  command <- paste0( "conv[1, i] <- mdl_bs", as.character(i), "$converged" )
  eval(parse(text = command))
}
conv
bic_val<-matrix(data=0,nrow=1,ncol=10)
for (i in 1:10) {
  command <- paste0( "bic_val[1, i] <- mdl_bs", as.character(i), "$sbc" )
  eval(parse(text = command))
}



bic_val
fit_best_order <- as.character(which.min(bic_val))
fit_best_order

################################################################################ step 2: get quantiles
quantiles<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
x<-matrix(data=NA,ncol=length(quantiles),nrow=8001)
y<-matrix(data=NA,ncol=length(quantiles),nrow=8001)

for (i in 1:length(quantiles)){
  command <- paste0( "Qua <- getQuantile(mdl_bs", fit_best_order, ", quantile=quantiles[i], term='Age', n.points = 8000)" )
  eval(parse(text = command))
  c<-curve(Qua, 0, 80,  lwd=2, lty=1, add=T,col="red", ylim=c (0,1), n = 8001)
  x[,i]<-c$x
  y[,i]<-c$y
}


################################################################################ step 3: get quantiles by sex
quantiles<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
x_male<-matrix(data=NA,ncol=length(quantiles),nrow=8001)
y_male<-matrix(data=NA,ncol=length(quantiles),nrow=8001)
x_female<-matrix(data=NA,ncol=length(quantiles),nrow=8001)
y_female<-matrix(data=NA,ncol=length(quantiles),nrow=8001)

for (i in 1:length(quantiles)){
  command <- paste0("Qua <- getQuantile(mdl_bs", fit_best_order ,", quantile=quantiles[i],term='Age',fixed.at=list(Sex=0), n.points = 8000)")
  eval(parse(text = command))
  c<-curve(Qua, 0, 80,  lwd=2, lty=1, add=T,col="red", ylim=c (0,1), n = 8001)
  x_female[,i]<-c$x
  y_female[,i]<-c$y
  
  command <- paste0("Qua <- getQuantile(mdl_bs", fit_best_order ,", quantile=quantiles[i],term='Age',fixed.at=list(Sex=1), n.points = 8000)")
  eval(parse(text = command))
  c<-curve(Qua, 0, 80,  lwd=2, lty=1, add=T,col="black", ylim=c (0,1), n = 8001)
  x_male[,i]<-c$x
  y_male[,i]<-c$y
}

################################################################################ step 4: get rate of change
age1=x[1:8000,4]
age2=x[2:8001,4]
v1=y[1:8000,4]
v2=y[2:8001,4]
a=age2-age1
v=v2-v1
velocity=v/a

v1=y_male[1:8000,4]
v2=y_male[2:8001,4]
v=v2-v1
velocity_male=v/a

v1=y_female[1:8000,4]
v2=y_female[2:8001,4]
v=v2-v1
velocity_female=v/a


################################################################################ step 5 :save model and other data
command <- paste0( "best_model <- mdl_bs", fit_best_order)
eval(parse(text = command))
#save(best_model, file = paste0("/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/gamlss_fit_allData_", idp_name ,".RData")) # save model

fname<-paste0("/home/jinlong/VDisk2/Jinlong/2024_IHFPS/IHFPs_result/prediction_ret/", fi_name, "/gamlss_fit_allData_", idp_name ,".mat")
command <- paste0("writeMat(fname,conv=conv,bic_val=bic_val,fit_best_order=fit_best_order,velocity_age=age1,all_velocity=velocity,velocity_male=velocity_male,velocity_female=velocity_female,all_cen_x=x,all_cen_y=y,male_cen_y=y_male,female_cen_y=y_female, quantiles=quantiles)")
eval(parse(text = command))


}
}