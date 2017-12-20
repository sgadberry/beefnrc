#' @title Nutrient Requirements for Bred Heifers Calving as 2 Yr Old.
#'
#' @description \code{preg.heifer} Nutrient requirements for bred heifers from
#' conception to calving.
#'
#' @param breed.code Default is 1 = Angus.  Breed code must match a code listed
#' in the table 19_1 data file.
#' @param mbw.lbs Default is 1200.  This value represents the mature, unshrunk bw
#' (lbs) for cows > 3 yrs old.  For first calf heifers, calving at 2 yr is assumed.
#' mbw.lbs is the estimated mature weight for first calf heifers, and equations are
#' based on #' a shrunk bw of 0.85(mbw) at first calving and 0.9(mbw) at second
#' calving.
#' @param breeding.pct.bw Default is 60. This value represents the percentage of
#' mature body weight at heifer breeding.
#' @param calving.pct.bw Default is 85. This value represents the target percentage
#' of mature body weight obtained by first calving.  This weight does not reflect the
#' additional weight of the conceptus.
#' @return Two tables are printed including nutrient requirement totals and nutrient
#' requirement as a proportion of dry matter intake.
#'
#' @examples
#' preg.heifer()
#' preg.heifer(breed.code=4,breeding.pct.bw=65)
#'
#' @export
#'
preg.heifer<-function(
    breed.code=1,
    mbw.lbs=1200,
    breeding.pct.bw=60,
    calving.pct.bw=85){

  #DATE INPUTS
  #1 starts as 1 month since conception,requirement at 30d increments
  code<-breed.code
  mnth<-rep(1:9)
  days.pregnant<-(mnth*30)

  #COW INPUTS
  current.bcs<-6

  #INPUTS to METRIC
  bw.kg<-mbw.lbs/2.2
  mbw.breeding<-bw.kg*breeding.pct.bw/100
  mbw.calving<-bw.kg*calving.pct.bw/100
  mbw.gain<-(mbw.calving-mbw.breeding)/285
  mbw.mnth.preg<-mbw.breeding+(mbw.gain*days.pregnant)
  sbw<-mbw.mnth.preg*0.96
  mbbw<-sbw^0.75
  calf.bw.kg<-(bw.kg*0.96*(calving.pct.bw/100))*0.066

  #MAINTENANCE CALC comp creates vector
  be<-ifelse(code==1,1,table19_1[table19_1[,1]==code,3])
  comp<-0.8+(current.bcs-1)*0.05
  sex<-1
  NEm.a1<-0.077
  NEm<-mbbw*(be*comp*sex*NEm.a1)  #used for total NEm
  MPm<-mbbw*3.8 #used for total MP


  #############STOPPED HERE#################

  #PREGNANCY CALC days.pregnant creates vector
  NEy<-calf.bw.kg*(0.05855-0.0000996*days.pregnant)*exp(0.03233*days.pregnant-0.0000275*days.pregnant^2)/1000
  MEy<-NEy/0.13
  Km.l<-0.6    #partial efficiency for ME to NE
  NEp<-MEy*Km.l    #used for total NEm
  Ypn<-calf.bw.kg*(0.001669-0.00000211*days.pregnant)*exp(0.0278*days.pregnant-0.0000176*days.pregnant^2)*6.25
  MPp<-Ypn/0.65   #used for total MP


  #GROWTH
  cw<-calf.bw.kg*0.01828*exp(0.02*days.pregnant-0.0000143*days.pregnant*2)
  srw<-478
  msbw<-bw.kg*0.96
  eqsbw<-(sbw-cw)*(srw/msbw)
  eqebw<-0.891*eqsbw
  bwg<-mbw.gain*0.96
  ewg<-bwg*0.96
  re<-(ewg/(12.341*eqebw^-0.6837))^(1/0.9116)

  adg.preg<-calf.bw.kg*(18.28*(0.02-0.0000286*days.pregnant)*exp(0.02*days.pregnant-0.0000143*days.pregnant^2))
  amnthg<-adg.preg*30/1000
  cum.amnthg<-rep(0,times=9)
  for(i in 1:length(amnthg)){
    if(i==1){
      cum.amnthg[i]<-amnthg[i]
    }else{
      cum.amnthg[i]<-cum.amnthg[i-1]+amnthg[i]
    }
  }






  #Total NEm
  total.NEm<-round(NEm+NEp,1)


  #DMI NEma estimation
  #create empty cells and starting tdn value
  tdn94<-rep(48,each=length(days.pregnant))
  de94<-tdn94/100*4.409
  me94<-de94*.82
  nema94<-1.37*me94-0.138*me94^2+0.0105*me94^3-1.12
  start.x94<-nema94
  dmi94<-rep(8,times=length(days.pregnant))
  nega94<-1.42*me94-0.174*me94^2+0.0122*me94^3-1.65
  ffm94<-total.NEm/nema94
  ffg94<-dmi94-ffm94
  total.NEg<-ffg94*nega94

  #loop solve for dmi and nema
  for(i in 1:length(days.pregnant)){
    while(round(total.NEg[i],2) < re[i]){
      tdn94[i]<-tdn94[i]+0.01
      de94[i]<-tdn94[i]/100*4.409
      me94[i]<-de94[i]*0.82
      nema94[i]<-1.37*me94[i]-0.138*me94[i]^2+0.0105*me94[i]^3-1.12
      nega94[i]<-1.42*me94[i]-0.174*me94[i]^2+0.0122*me94[i]^3-1.65
      start.x94[i]<-nema94[i]
      dmi94[i]<-mbbw[i]*((0.0384+0.04997*start.x94[i]^2)/start.x94[i])
      ffm94[i]<-total.NEm[i]/nema94[i]
      ffg94[i]<-dmi94[i]-ffm94[i]
      total.NEg[i]<-ffg94[i]*nega94[i]
    }
  }

  #days pregnant < 94 formula

  #create empty cells and starting tdn value
  #create empty cells and starting tdn value
  tdn093<-rep(48,each=length(days.pregnant))
  de093<-tdn093/100*4.409
  me093<-de093*.82
  nema093<-1.37*me093-0.138*me093^2+0.0105*me093^3-1.12
  start.x093<-nema093
  dmi093<-rep(8,times=length(days.pregnant))
  nega093<-1.42*me093-0.174*me093^2+0.0122*me093^3-1.65
  ffm093<-total.NEm/nema093
  ffg093<-dmi093-ffm093
  total.NEg<-ffg093*nega093

  #loop solve for dmi and nema
  for(i in 1:length(days.pregnant)){
    while(round(total.NEg[i],2) < re[i]){
      tdn093[i]<-tdn093[i]+0.01
      de093[i]<-tdn093[i]/100*4.409
      me093[i]<-de093[i]*0.82
      nema093[i]<-1.37*me093[i]-0.138*me093[i]^2+0.0105*me093[i]^3-1.12
      nega093[i]<-1.42*me093[i]-0.174*me093[i]^2+0.0122*me093[i]^3-1.65
      start.x093[i]<-nema093[i]
      dmi093[i]<-mbbw[i]*((0.0384+0.04997*start.x093[i]^2)/start.x093[i])
      ffm093[i]<-total.NEm[i]/nema093[i]
      ffg093[i]<-dmi093[i]-ffm093[i]
      total.NEg[i]<-ffg093[i]*nega093[i]
    }
  }

  dmi<-ifelse(days.pregnant < 94, dmi093,dmi94)
  tdn<-ifelse(days.pregnant < 94, tdn093,tdn94)
  NEma<-ifelse(days.pregnant < 94, start.x093,start.x94)
  NEga<-ifelse(days.pregnant < 94, nega093,nega94)
  NE.intake<-dmi*NEma
  dmi.NEm<-total.NEm/NEma
  total.NEg<-(dmi-dmi.NEm)*NEga
  rbind(dmi,tdn,NEma,NE.intake,total.NEm)

  #protein totals
  NPg<-bwg*(268-29.4*re/bwg)
  MPg<-NPg/max(c(0.492,0.834-0.00114*eqsbw))
  total.MP<-round(MPm+MPp+MPg,1)
  total.CP<-round(total.MP/0.67,1)


  #calcium requirement
  Ca.m<-0.0154*sbw/0.5
  Ca.p<-ifelse(mnth > 6,calf.bw.kg*(13.7/90)/0.5,0)
  Ca.g<-NPg*0.071/0.5
  ca.total<-Ca.m+Ca.p+Ca.g

  #phosphorus requirement
  P.m<-0.016*sbw/0.68
  P.p<-ifelse(mnth > 6,calf.bw.kg*(7.6/90)/0.68,0)
  P.g<-NPg*0.039/0.68
  p.total<-P.m+P.g+P.p

  cat('Breed:',ifelse(code==1,'Angus',as.character(table19_1[table19_1[,1]==code,2])),"\n")
  cat('Mature BW, lbs:',mbw.lbs,"\n")
  cat('Calf birth weight set at 0.066(dams shrunk BW):',calf.bw.kg*2.2,"\n")
  cat("\n")
  cat('Daily total Calculations',"\n")
  print(data.frame(Mnth.Since.Conception=mnth,
                   TDN.lb=round(tdn*dmi*2.2/100,2),
                   NEm.Mcal=total.NEm,
                   NEg.Mcal=round(total.NEg,1),
                   CP.lb=round(total.CP/454.54,2),
                   Ca.lb=round(ca.total/454.54,3),
                   P.lb=round(p.total/454.54,3)))
  cat("\n")
  cat('Daily Proportion',"\n")
  print(data.frame(Mnth.Since.Conception=mnth,
                   Unshrunk.BW.lb=round(mbw.mnth.preg*2.2),
                   BWwPreg=round(mbw.mnth.preg*2.2+(cum.amnthg*2.2)),
                   TDN.pct=round(tdn,1),
                   DMI.lb=round(dmi*2.2,1),
                   NEma.Mcal.lb=round(NEma/2.2,2),
                   NEg.Mcal.lb=round(NEga/2.2,2),
                   CP.pct=round(total.CP/1000/dmi*100,1),
                   Ca.pct=round(ca.total/1000/dmi*100,2),
                   P.pct=round(p.total/1000/dmi*100,2)))
}
#'
