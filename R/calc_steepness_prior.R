#' Calculate prior for steepness
#'
#' \code{calc_steepness_prior} calculates a prior for steepness given likelihood profiles for individual species
#'
#' @param LogLike_zp matrix of log-likelihood values for each species \code{p} and likelihood value \code{z}
#' @param h_z Steepness value used to calculate each row of \code{LogLike_zp}

#' @export
calc_steepness_prior <-
function( LogLike_zp, h_z, Interpolate=TRUE, Degree=1, Ninterp=203, StockWeight=rep(1,ncol(LogLike_zp)),
    Downweight_outliers=1, Remove_at_bounds=FALSE, NburnIn=1e6, Nsamp=1e6, Nthin=1e2, RunDir=getwd(),
    AdmbDir=system.file("executables",package="SebastesSteepness") ){

  ################
  # Local functions
  ################

  calc_mode = function(Vec){
    Density = density(Vec)
    return( Density$x[which.max(Density$y)] )
  }
  par_fn=function(Npanel=1){
    Ncol = ceiling(sqrt(Npanel))
    Nrow = ceiling(Npanel/Ncol)
    Par = list(mgp=c(2,0.75,0),mar=c(3,3,2,0.05),tck=-0.02,mfrow=c(Nrow,Ncol))
    return(Par)
  }

  ################
  # Pre-process inputs
  ################

  # Possibly kick out species at bounds
  WhichMin = apply(LogLike_zp, MARGIN=2, FUN=which.max)
  WhichRemove = which(WhichMin==1 | WhichMin==nrow(LogLike_zp))
  if(Remove_at_bounds==TRUE & length(WhichRemove)>0){
    LogLike_zp = LogLike_zp[,-WhichRemove]
  }

  # Relax outliers
  WhichMin = apply(LogLike_zp, MARGIN=2, FUN=which.max)
  WhichOutliers = which(WhichMin==1 | WhichMin==nrow(LogLike_zp))
  if( Downweight_outliers!=1 & length(WhichOutliers)>0 ){
    LogLike_zp[,WhichOutliers] = LogLike_zp[,WhichOutliers] / Downweight_outliers
  }

  # Interpolate
  if(Interpolate==TRUE){
    hInterp_z = seq(0.2,1,length=Ninterp+2)[-c(1,Ninterp+2)]
    LogLikeInterp_zp = matrix(NA, ncol=ncol(LogLike_zp), nrow=Ninterp)

    for(StockI in 1:ncol(LogLikeInterp_zp)){
      Lm = lm(LogLike_zp[,StockI] ~ 0 + splines::bs(h_z, df=sum(!is.na(LogLike_zp[,StockI])), degree=Degree, intercept=TRUE) )
      LogLikeInterp_zp[,StockI] = predict( Lm, newdata=data.frame("h_z"=hInterp_z) )
    }
  }else{
    hInterp_z = h_z
    LogLikeInterp_zp = LogLike_zp[,,2]
  }

  ################
  # Run MCMC using ADMB
  ################

  # Files
  dir.create(RunDir)
  RunCommand = "profile_mc_JTT_v2"

  write(c("#Nstocks",ncol(LogLikeInterp_zp)), file=paste0(RunDir,"/",RunCommand,".dat"))
  write(c("#Nprofile",nrow(LogLikeInterp_zp)), file=paste0(RunDir,"/",RunCommand,".dat"), append=TRUE)
  write("#h_profX", file=paste0(RunDir,"/",RunCommand,".dat"), append=TRUE)
  write.table(matrix(hInterp_z,nrow=1), file=paste0(RunDir,"/",RunCommand,".dat"), row.names=FALSE, col.names=FALSE, append=TRUE)
  write("#StockWeight", file=paste0(RunDir,"/",RunCommand,".dat"), append=TRUE)
  write.table(matrix(StockWeight,nrow=1), file=paste0(RunDir,"/",RunCommand,".dat"), row.names=FALSE, col.names=FALSE, append=TRUE)
  write("#h_profY", file=paste0(RunDir,"/",RunCommand,".dat"), append=TRUE)  # NLL
  write.table(-1*t(LogLikeInterp_zp), file=paste0(RunDir,"/",RunCommand,".dat"), row.names=FALSE, col.names=FALSE, append=TRUE)
  if(Version=="Thorson") write(c("#TestVal",123456), file=paste0(RunDir,"/",RunCommand,".dat"), append=TRUE)

  # Run ADMB
  file.copy(from=paste0(AdmbDir,"/",RunCommand,".exe"), to=paste0(RunDir,"/",RunCommand,".exe"), overwrite=TRUE)
  setwd(RunDir)
  shell(RunCommand)

  # Read in diagnostics
  Par = scan(file=paste0(RunDir,"/",RunCommand,".par"),comment.char="#") # ,what="character"
  Beta = scan(file=paste0(RunDir,"/",RunCommand,".rep"),skip=9,nlines=1)
  Jacob = scan(file=paste0(RunDir,"/",RunCommand,".rep"),skip=11,nlines=1)
  h_profX = scan(file=paste0(RunDir,"/",RunCommand,".rep"),skip=13,nlines=1)
  h_profY = matrix(scan(file=paste0(RunDir,"/",RunCommand,".rep"),skip=15,nlines=ncol(LogLikeInterp_zp)), nrow=ncol(LogLikeInterp_zp), byrow=TRUE)
  LikeMat = matrix(scan(file=paste0(RunDir,"/",RunCommand,".rep"),skip=15+ncol(LogLikeInterp_zp)+1,nlines=ncol(LogLikeInterp_zp)), nrow=ncol(LogLikeInterp_zp), byrow=TRUE)

  # Run MCMC
  Files2Remove = c("mu.out","tau.out",paste0(RunCommand,".psv"))
  if( any(Files2Remove %in% list.files(RunDir)) ){
    file.remove( paste0(RunDir,"/",Files2Remove[Files2Remove %in% list.files(RunDir)]) )
  }
  rm(list=c("Mu","Sd_beta")[c("Mu","Sd_beta")%in%ls()])
  shell(paste0(RunCommand," -mcmc ",NburnIn," -mcscale ",NburnIn))
  shell(paste0(RunCommand," -mcr -mcmc ",Nsamp," -mcsave ",Nthin))
  shell(paste0(RunCommand," -mceval"))

  Mu = scan(paste0(RunDir,"/mu.out"))
  Mu_h = plogis(Mu)*0.8+0.2
  Sd_beta = scan(paste0(RunDir,"/tau.out"))
    mean(Sd_beta); calc_mode(Sd_beta)
  Pred_h = plogis(rnorm(length(Mu), mean=Mu, sd=Sd_beta))*0.8 + 0.2
    mean(Pred_h); calc_mode(Pred_h)
    sd(Pred_h)
  Pred_logitnorm = rnorm(length(Mu), mean=Mu, sd=Sd_beta)
    mean(Pred_logitnorm); calc_mode(Pred_logitnorm)
    sd(Pred_logitnorm)

  ################
  # Generate plots
  ################

  # Plot results
  png(file=paste0(RunDir,"/Results.png"), width=8,height=8,res=200,units="in")
    par( par_fn(4) )
    hist(Mu_h, breaks=seq(0.2,1,by=0.05),main="Mu_h_Post")
    hist(Sd_beta, breaks=seq(0,max(Sd_beta)+0.05,by=0.05),main="Sd_beta_Post")
    plot(x=Mu_h,y=Sd_beta,main="Both_Post")
    points(x=mean(Mu_h), y=mean(Sd_beta), col="red", cex=2)
    hist(Pred_h, breaks=seq(0.2,1,by=0.05),main="h_Pred")
    text(x=mean(axTicks(1)),y=max(axTicks(2)),labels=paste("Mean=",format(mean(Pred_h),digits=3),"\n","SD=",format(sd(Pred_h),digits=3),sep=""), adj=c(0.5,1))  # , font=2
  dev.off()

  ##### Convergence diagnostics
  # Save Trace
  png(paste0(RunDir,"/Trace.png"),width=2*2,height=2*1,units="in",res=200)
    par(mfrow=c(1,2), mar=c(0,1,1,0),mgp=c(1.25,0.25,0),tck=-0.02)
    for(ParI in 1:2){
      matplot(cbind(Mu_h,Sd_beta)[,ParI],type="l",xaxt="n",main=c("Mu_h","Sd_beta")[ParI]) # , ylim=c(0,min(5,max(Derived[seq(1,nrow(Derived),length=1000),ParI,])))
    }
  dev.off()   #-0.295287  , 0.247201

  # Save Autocorrelation
  png(paste0(RunDir,"/ACF.png"),width=2*2,height=2*1,units="in",res=200)
    par(mfrow=c(1,2), mar=c(0,1,1,0),mgp=c(1.25,0.25,0),tck=-0.02)
    for(ParI in 1:2){
      Acf = apply(cbind(Mu_h,Sd_beta)[,ParI,drop=FALSE],MARGIN=2,FUN=function(Vec){acf(Vec,plot=FALSE)$acf})
      if(!any(is.na(Acf))){
        matplot(Acf, type="h", lty="solid", xaxt="n", xlab="",ylab="", ylim=c(0,1), lwd=2, main=c("Mu_h","Sd_beta")[ParI])
      }else{
        plot.new()
      }
    }
  dev.off()   #-0.295287  , 0.247201
}
