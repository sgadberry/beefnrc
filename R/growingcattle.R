#' @title Nutrient Requirements for All Growing Cattle.
#'
#' @description \code{growing.cattle} Nutrient requirements for growing cattle.
#'
#' @param breed.code Default is 1 = Angus.  Breed code must match a code listed
#' in the table 19_1 data file.
#' @param Bull Default is false. Adjustment for NEm applied when TRUE.
#' @param mbw.lbs Default is 1200.  This value represents the mature, unshrunk bw
#' (lbs) for growing cattle at finishing or at maturity.
#' @param yearling Default is FALSE. This value adjusts DMI for growing cattle > 11
#' months of age.
#' @param target.adg Default is 2. This value should represent realistic rates of
#' gain such as 0.5 as a low and 3.2 as a high.
#' @param from.bw.lbs Avoid using values < 400. This sets the minimum of the range
#' for prediction.
#' @param to.bw.lbs This value sets the maximum of the range for prediction.  Maximum
#' should be realistic for age (calf vs yearling).
#' @param increment.bw.lbs This value sets the incremental change between
#'  \code{from.bw.lbs} to \code{to.bw.lbs}.
#' @param implant Default is TRUE.  When \code{Bull} is TRUE, implant is overridden to non-implant
#' factor of 0.94.
#' @param ionophore Default is TRUE. When \code{Bull} is TRUE, ionophore is overridden to
#' no ionophore.
#' @param monensin Default is TRUE.  When \code{ionophore} is FALSE, monensin is overridden.
#'
#' @return Two tables are printed including nutrient requirement totals and nutrient
#' requirement as a proportion of dry matter intake.
#'
#' @examples
#' growing.cattle()
#' growing.cattle(implant=FALSE,ionophore=FALSE)
#'
#' @export
growing.cattle<-function(
  breed.code=1,
  Bull=FALSE,
  mbw.lbs=1200,
  yearling=FALSE,
  target.adg=2,
  from.bw.lbs=400,
  to.bw.lbs=800,
  increment.bw.lbs=100,
  implant=TRUE,
  ionophore=TRUE,
  monensin=TRUE){

  implant.dmi<-ifelse(isTRUE(Bull),0.94,
                      ifelse(isTRUE(implant),1,0.94))
  ionophore.dmi<-ifelse(isTRUE(Bull),1,
                        ifelse(isTRUE(ionophore)&isTRUE(monensin),0.97,1))
  ionophore.nema<-ifelse(isTRUE(Bull),1,
                         ifelse(isTRUE(ionophore)&isTRUE(monensin),1.023,
                                ifelse(ionophore==TRUE&monensin != TRUE,1.015,1)))

  #DATE INPUTS
  #1 starts as 1 month since conception,requirement at 30d increments
  code<-breed.code
  bw.range<-seq(from=from.bw.lbs,to=to.bw.lbs,by=increment.bw.lbs)
  bw.length<-length(bw.range)

  implant.dmi2<-rep(implant.dmi,length=bw.length)
  ionophore.dmi2<-rep(ionophore.dmi,length=bw.length)
  ionophore.nema2<-rep(ionophore.nema,length=bw.length)

  #BCS INPUTS
  current.bcs<-6

  #INPUTS to METRIC
  mature.bw.kg<-mbw.lbs/2.2
  mature.sbw.kg<-mature.bw.kg*0.96
  current.bw.kg<-bw.range/2.2
  current.sbw.kg<-current.bw.kg*0.96
  current.mbbw.kg<-current.sbw.kg^0.75
  adg.kg<-target.adg/2.2



  #MAINTENANCE CALC comp creates vector
  be<-ifelse(code==1,1,table19_1[table19_1[,1]==code,3])
  comp<-0.8+(current.bcs-1)*0.05
  sex<-ifelse(isTRUE(Bull),1.15,1)
  NEm.a1<-0.077
  NEm<-current.mbbw.kg*(be*comp*sex*NEm.a1)  #used for total NEm
  MPm<-current.mbbw.kg*3.8 #used for total MP


  #GROWTH
  srw<-478
  eqsbw<-(current.sbw.kg)*(srw/mature.sbw.kg)
  eqebw<-0.891*eqsbw
  bwg<-adg.kg*.96
  ewg<-bwg*0.96
  re<-(ewg/(12.341*eqebw^-0.6837))^(1/0.9116)

  #Total NEm
  total.NEm<-round(NEm,1)

  #DMI CALF PREDICTION 1996 NRC

  #DMI NEma estimation
  #create empty cells and starting tdn value
  tdn.C<-rep(45,each=bw.length)
  de.C<-tdn.C/100*4.409
  me.C<-de.C*.82
  nema.C<-1.37*me.C-0.138*me.C^2+0.0105*me.C^3-1.12
  start.x.C<-nema.C
  dmi.C<-rep(current.sbw.kg*0.015)
  nega.C<-1.42*me.C-0.174*me.C^2+0.0122*me.C^3-1.65
  ffm.C<-total.NEm/nema.C
  ffg.C<-dmi.C-ffm.C
  total.NEg.C<-ffg.C*nega.C

  #loop solve for calves dmi and nema
  for(i in 1:bw.length){
    while(round(total.NEg.C[i],2) < re[i]){
      tdn.C[i]<-tdn.C[i]+0.01
      de.C[i]<-tdn.C[i]/100*4.409
      me.C[i]<-de.C[i]*0.82
      nema.C[i]<-1.37*me.C[i]-0.138*me.C[i]^2+0.0105*me.C[i]^3-1.12
      nega.C[i]<-1.42*me.C[i]-0.174*me.C[i]^2+0.0122*me.C[i]^3-1.65
      start.x.C[i]<-nema.C[i]*ionophore.nema
      dmi.C[i]<-(current.mbbw.kg[i]*(0.2435*start.x.C[i]-0.0466*start.x.C[i]^2-0.1128)/start.x.C[i])*implant.dmi*ionophore.dmi
      ffm.C[i]<-total.NEm[i]/start.x.C[i]
      ffg.C[i]<-dmi.C[i]-ffm.C[i]
      total.NEg.C[i]<-ffg.C[i]*nega.C[i]
    }
  }

  #DMI YEARLING PREDICTION 1996 NRC

  #DMI NEma estimation
  #create empty cells and starting tdn value
  tdn.Y<-rep(45,each=bw.length)
  de.Y<-tdn.Y/100*4.409
  me.Y<-de.Y*.82
  nema.Y<-1.37*me.Y-0.138*me.Y^2+0.0105*me.Y^3-1.12
  start.x.Y<-nema.Y
  dmi.Y<-rep(current.sbw.kg*0.015)
  nega.Y<-1.42*me.Y-0.174*me.Y^2+0.0122*me.Y^3-1.65
  ffm.Y<-total.NEm/nema.Y
  ffg.Y<-dmi.Y-ffm.Y
  total.NEg.Y<-ffg.Y*nega.Y

  #loop solve for calves dmi and nema
  for(i in 1:bw.length){
    while(round(total.NEg.Y[i],2) < re[i]){
      tdn.Y[i]<-tdn.Y[i]+0.01
      de.Y[i]<-tdn.Y[i]/100*4.409
      me.Y[i]<-de.Y[i]*0.82
      nema.Y[i]<-1.37*me.Y[i]-0.138*me.Y[i]^2+0.0105*me.Y[i]^3-1.12
      nega.Y[i]<-1.42*me.Y[i]-0.174*me.Y[i]^2+0.0122*me.Y[i]^3-1.65
      start.x.Y[i]<-nema.Y[i]*ionophore.nema
      dmi.Y[i]<-(current.mbbw.kg[i]*(0.2435*start.x.Y[i]-0.0466*start.x.Y[i]^2-0.0869)/start.x.Y[i])*implant.dmi*ionophore.dmi
      ffm.Y[i]<-total.NEm[i]/start.x.Y[i]
      ffg.Y[i]<-dmi.Y[i]-ffm.Y[i]
      total.NEg.Y[i]<-ffg.Y[i]*nega.Y[i]
    }
  }

  #DMI ALTERNATE WITH NEma

  #DMI NEma estimation
  #create empty cells and starting tdn value
  #tdn.ALT<-rep(45,each=bw.length)
  #de.ALT<-tdn.ALT/100*4.409
  #me.ALT<-de.ALT*.82
  #nema.ALT<-1.37*me.ALT-0.138*me.ALT^2+0.0105*me.ALT^3-1.12
  #start.x.ALT<-nema.ALT
  #dmi.ALT<-rep(current.sbw.kg*0.015)
  #nega.ALT<-1.42*me.ALT-0.174*me.ALT^2+0.0122*me.ALT^3-1.65
  #ffm.ALT<-total.NEm/nema.ALT
  #ffg.ALT<-dmi.ALT-ffm.ALT
  #total.NEg.ALT<-ffg.ALT*nega.ALT

  #loop solve for calves dmi and nema
  #for(i in 1:bw.length){
  #    while(round(total.NEg.ALT[i],2) < re[i]){
  #		tdn.ALT[i]<-tdn.ALT[i]+0.01
  #		de.ALT[i]<-tdn.ALT[i]/100*4.409
  #		me.ALT[i]<-de.ALT[i]*0.82
  #		nema.ALT[i]<-1.37*me.ALT[i]-0.138*me.ALT[i]^2+0.0105*me.ALT[i]^3-1.12
  #		nega.ALT[i]<-1.42*me.ALT[i]-0.174*me.ALT[i]^2+0.0122*me.ALT[i]^3-1.65
  #		start.x.ALT[i]<-nema.ALT[i]*ionophore.nema
  #		dmi.ALT[i]<-(current.sbw.kg[i]*(0.012425+0.019218*start.x.ALT[i]-0.007259*start.x.ALT[i]^2))*implant.dmi*ionophore.dmi
  #		ffm.ALT[i]<-total.NEm[i]/start.x.ALT[i]
  #		ffg.ALT[i]<-dmi.ALT[i]-ffm.ALT[i]
  #		total.NEg.ALT[i]<-ffg.ALT[i]*nega.ALT[i]
  #	}
  #}

  #protein totals
  NPg<-bwg*(268-29.4*re/bwg)
  MPg<-NPg/max(c(0.492,0.834-0.00114*eqsbw))
  total.MP<-round(MPm+MPg,1)
  total.CP<-round(total.MP/0.67,1)


  #calcium requirement
  Ca.m<-0.0154*current.sbw.kg/0.5
  Ca.g<-NPg*0.071/0.5
  ca.total<-Ca.m+Ca.g

  #phosphorus requirement
  P.m<-0.016*current.sbw.kg/0.68
  P.g<-NPg*0.039/0.68
  p.total<-P.m+P.g

  #1996 DMI Equivalent Tables

  dmi.96<-ifelse(rep(isTRUE(yearling),times=bw.length),dmi.Y,dmi.C)
  tdn.96<-ifelse(rep(isTRUE(yearling),times=bw.length),tdn.Y,tdn.C)
  NEma.96<-ifelse(rep(isTRUE(yearling),times=bw.length),nema.Y,nema.C)
  NEga.96<-ifelse(rep(isTRUE(yearling),times=bw.length),nega.Y,nega.C)
  NE.intake.96<-dmi.96*NEma.96

  cat('Breed:',ifelse(code==1,'Angus',as.character(table19_1[table19_1[,1]==code,2])),"\n")
  cat('Mature BW, lbs:',mbw.lbs,"\n")
  cat('Age Group:',ifelse(isTRUE(yearling),'Yearling','Calf'),"\n")
  cat('ADG, lbs:',target.adg,"\n")
  implant.use<-ifelse(implant==TRUE,"GROWTH IMPLANT ADJUSTMENT","NO GROWTH IMPLANT ADJUSTMENT")
  ionophore.use<-ifelse(ionophore==TRUE,"IONOPHORE ADJUSTMENT","NO IONOPHORE ADJUSTMENT")
  cat(implant.use,"\n")
  cat("implant dmi adjustment:",implant.dmi,"\n")
  cat(ionophore.use,"\n")
  cat("ionophore dmi adjustment:",ionophore.dmi,"\n")
  cat("ionophore nema adjustment:",ionophore.nema,"\n")
  cat("\n")
  cat('1996 DMI Equation Daily Total Calculations',"\n")
  print(data.frame(Unshrunk.BW.lb=bw.range,
                   TDN.lb=round(tdn.96*dmi.96*2.2/100,2),
                   NEm.Mcal=total.NEm,
                   NEg.Mcal=round(ifelse(isTRUE(yearling),total.NEg.Y,total.NEg.C),1),
                   CP.lb=round(total.CP/454.54,2),
                   Ca.lb=round(ca.total/454.54,3),
                   P.lb=round(p.total/454.54,3)))
  cat("\n")
  cat('1996 DMI Equation Daily Proportion',"\n")
  print(data.frame(Unshrunk.BW.lb=bw.range,
                   TDN.pct=round(tdn.96,1),
                   DMI.lb=round(dmi.96*2.2,1),
                   NEma.Mcal.lb=round(NEma.96/2.2,2),
                   NEg.Mcal.lb=round(NEga.96/2.2,2),
                   CP.pct=round(total.CP/1000/dmi.96*100,1),
                   Ca.pct=round(ca.total/1000/dmi.96*100,2),
                   P.pct=round(p.total/1000/dmi.96*100,2)))
  cat("\n")
  #cat('2016 ALT NEma DMI Equation Daily Total Calculations',"\n")
  #print(data.frame(Unshrunk.BW.lb=bw.range,
  #				TDN.lb=round(tdn.ALT*dmi.ALT*2.2/100,2),
  #				NEm.Mcal=total.NEm,
  #				NEg.Mcal=round(total.NEg.ALT,1),
  #				CP.lb=round(total.CP/454.54,2),
  #				Ca.lb=round(ca.total/454.54,3),
  #				P.lb=round(p.total/454.54,3)))

  #cat("\n")
  #cat('2016 ALT NEma DMI Equation Daily Proportion',"\n")
  #print(data.frame(Unshrunk.BW.lb=bw.range,
  #				TDN.pct=tdn.ALT,
  #				DMI.lb=round(dmi.ALT*2.2,1),
  #				NEma.Mcal.lb=round(nema.ALT/2.2,2),
  #				NEg.Mcal.lb=round(nega.ALT/2.2,2),
  #				CP.pct=round(total.CP/1000/dmi.ALT*100,1),
  #				Ca.pct=round(ca.total/1000/dmi.ALT*100,2),
  #				P.pct=round(p.total/1000/dmi.ALT*100,2)))

}
