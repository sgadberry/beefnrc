#' @title Calculate requirements for mature cows and first calf heifers
#'
#' @description
#' \code{mature.cow} returns calculated requirements for mature cows and first calf heifers as
#' daily totals and proportion of dry matter intake.
#'
#' @param breed.code Default is 1 = Angus.  Breed code must match a code listed
#' in the table 19_1 data file.
#' @param mbw.lbs Default is 1200.  This value represents the mature, unshrunk bw
#' (lbs) for cows > 3 yrs old.  For first calf heifers, calving at 2 yr is assumed.
#' mbw.lbs is the estimated mature weight for first calf heifers, and equations are based on
#' a shrunk bw of 0.85(mbw) at first calving and 0.9(mbw) at second calving.
#' @param milk.pk.lb Default is 20.  Users can modify peak milk (lbs). Peak milks is estimated
#' to occur at 8.5 weeks after calving.
#' @param RMY Default is 5.  Relative Milk Yield (RMY) is a convenient method of adjusting
#' milk production.
#' @param first.calf Default is FALSE.  Change to TRUE to adjust requirements for first
#' calf heifers.
#' @param first.calf.NEmadj Default is 1.25.  This value can be changed to NRC default 1.
#' @return Two tables are printed including nutrient requirement totals and nutrient
#' requirement as a proportion of dry matter intake.
#'
#' @examples
#' mature.cow()
#' mature.cow(breed.code=4)
#' mature.cow(breed.code=4,first.calf=TRUE)
#' @export

mature.cow<-function(
  breed.code=1,
  mbw.lbs=1200,
  milk.pk.lb=20,
  RMY=5,
  first.calf=FALSE,
  first.calf.NEmadj=1.25){

  #DATE INPUTS
  #mnth 1 starts at first calving;requirement at 30d
  code<-breed.code
  mnth<-rep(1:12)
  wean.mnth<-7
  days.in.lactation<-ifelse(mnth > wean.mnth,mnth*0,mnth*30)
  days.pregnant<-ifelse(mnth <= 3,mnth*0,mnth*30-90)

  #COW INPUTS
  age.yr<-ifelse(isTRUE(first.calf),2,4)
  age.yr2<-rep(age.yr,times=length(mnth))
  BCS.by.mnth<-c(5,5,5,5,5,5,5,5,5,5,5,5)
  milk.fat<-ifelse(code==1,4.0,table19_1[table19_1[,1]==code,7])
  milk.prot<-ifelse(code==1,3.8,table19_1[table19_1[,1]==code,8])
  milk.solid<-ifelse(code==1,8.3,table19_1[table19_1[,1]==code,9])


  #INPUTS to METRIC
  bw.kg<-mbw.lbs/2.2
  first.calf.bw<-bw.kg*0.85
  second.calf.bw<-bw.kg*0.9
  first.calf.bw.gain<-(second.calf.bw-first.calf.bw)/12
  cumulative.bw.gain<-mnth*first.calf.bw.gain
  sbw.first.calf<-(first.calf.bw+cumulative.bw.gain)*0.96
  sbw.cow<-rep(bw.kg*0.96,times=12)
  mbbw.first.calf<-sbw.first.calf^0.75
  mbbw.cow<-sbw.cow^0.75
  mbbw<-ifelse(age.yr2==2,mbbw.first.calf,mbbw.cow)
  milk.pk.kg<-milk.pk.lb/2.2
  calf.bw.kg<-ifelse(isTRUE(first.calf),0.9*(bw.kg)*0.96*0.066,bw.kg*0.96*0.066)

  #MAINTENANCE CALC comp creates vector
  be<-ifelse(code==1,1,table19_1[table19_1[,1]==code,3])
  l.fact<-ifelse(code==1,1.2,table19_1[table19_1[,1]==code,4])
  l.fact2<-ifelse(days.in.lactation > 0,l.fact,1)
  comp<-0.8+(BCS.by.mnth-1)*0.05
  sex<-1
  NEm.a1<-0.077*ifelse(isTRUE(first.calf) & days.in.lactation >0,first.calf.NEmadj,1)
  NEm<-mbbw*(be*l.fact2*comp*sex*NEm.a1)  #used for total NEm
  MPm<-mbbw*3.8 #used for total MP

  #LACTATION CALC lact.n creates vector
  t<-8.5
  k<-1/t
  #RMY<-5
  aPKYD<-milk.pk.kg*(0.125*RMY+0.375)
  a<-1/(aPKYD*k*exp(1))
  lact.n<-days.in.lactation/7
  age.fact<-ifelse(isTRUE(first.calf),0.74,1)#adjust milk for first calf heifers
  E<-0.092*milk.fat+0.049*milk.solid-0.0569
  Yn<-lact.n/(a*exp(k*lact.n))*age.fact #lact.n creates vector
  NEl<-Yn*E                 #used for total NEm
  YProtn<-Yn*milk.prot/100
  MPl<-(YProtn/0.65)*1000   #used for total MP

  #PREGNANCY CALC days.pregnant creates vector
  NEy<-calf.bw.kg*(0.05855-0.0000996*days.pregnant)*exp(0.03233*days.pregnant-0.0000275*days.pregnant^2)/1000
  MEy<-NEy/0.13
  Km.l<-0.6    #partial efficiency for ME to NE
  NEp<-MEy*Km.l    #used for total NEm
  Ypn<-calf.bw.kg*(0.001669-0.00000211*days.pregnant)*exp(0.0278*days.pregnant-0.0000176*days.pregnant^2)*6.25
  MPp<-Ypn/0.65   #used for total MP


  #FIRST.CALF HEIFER GROWTH
  cw<-calf.bw.kg*0.01828*exp(0.02*days.pregnant-0.0000143*days.pregnant*2)
  srw<-478
  msbw<-bw.kg*0.96
  eqsbw<-(sbw.first.calf-cw)*(srw/msbw)
  eqebw<-0.891*eqsbw
  bwg<-first.calf.bw.gain/30
  ewg<-bwg*0.96
  re<-(ewg/(12.341*eqebw^-0.6837))^(1/0.9116)
  adg.preg<-calf.bw.kg*(18.28*(0.02-0.0000286*days.pregnant)*exp(0.02*days.pregnant-0.0000143*days.pregnant^2))
  amnthg<-adg.preg*30/1000
  cum.amnthg<-rep(0,times=12)
  for(i in 1:length(amnthg)){
    if(i==1){
      cum.amnthg[i]<-amnthg[i]
    }else{
      cum.amnthg[i]<-cum.amnthg[i-1]+amnthg[i]
    }
  }

  #Total NEm
  total.NEm<-round(NEm+NEl+NEp,1)


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
    while(ifelse(isTRUE(first.calf),round(total.NEg[i],2) < re[i],start.x94[i]*dmi94[i] <= total.NEm[i])){
      tdn94[i]<-tdn94[i]+0.01
      de94[i]<-tdn94[i]/100*4.409
      me94[i]<-de94[i]*0.82
      nema94[i]<-1.37*me94[i]-0.138*me94[i]^2+0.0105*me94[i]^3-1.12
      nega94[i]<-1.42*me94[i]-0.174*me94[i]^2+0.0122*me94[i]^3-1.65
      start.x94[i]<-nema94[i]
      dmi94[i]<-mbbw[i]*((0.0384+0.04997*start.x94[i]^2)/start.x94[i])+0.2*Yn[i]
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
    while(ifelse(isTRUE(first.calf),round(total.NEg[i],2) < re[i],start.x093[i]*dmi093[i] <= total.NEm[i])){
      tdn093[i]<-tdn093[i]+0.01
      de093[i]<-tdn093[i]/100*4.409
      me093[i]<-de093[i]*0.82
      nema093[i]<-1.37*me093[i]-0.138*me093[i]^2+0.0105*me093[i]^3-1.12
      nega093[i]<-1.42*me093[i]-0.174*me093[i]^2+0.0122*me093[i]^3-1.65
      start.x093[i]<-nema093[i]
      dmi093[i]<-mbbw[i]*((0.0384+0.04997*start.x093[i]^2)/start.x093[i])+0.2*Yn[i]
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
  rbind(dmi,tdn,NEma,NE.intake,total.NEm)

  #protein totals
  NPg<-ifelse(isTRUE(first.calf),bwg*(268-29.4*re/bwg),0)
  MPg<-NPg/max(c(0.492,0.834-0.00114*eqsbw))
  total.MP<-round(MPm+MPl+MPp+MPg,1)
  total.CP<-round(total.MP/0.67,1)


  #calcium requirement
  Ca.m<-0.0154*ifelse(age.yr==2,sbw.first.calf,sbw.cow)/0.5
  Ca.l<-Yn*1.23/0.5
  Ca.p<-ifelse(mnth > 9,calf.bw.kg*(13.7/90)/0.5,0)
  Ca.g<-ifelse(isTRUE(first.calf),NPg*0.071/0.5,0)
  ca.total<-Ca.m+Ca.l+Ca.p+Ca.g

  #phosphorus requirement
  P.m<-0.016*ifelse(age.yr==2,sbw.first.calf,sbw.cow)/0.68
  P.l<-Yn*0.95/0.68
  P.p<-ifelse(mnth > 9,calf.bw.kg*(7.6/90)/0.68,0)
  P.g<-ifelse(isTRUE(first.calf),NPg*0.039/0.68,0)
  p.total<-P.m+P.l+P.p+P.g

  breed<-as.character(table19_1[table19_1$code==breed.code,2])

  cat('Breed:',breed,"\n")
  cat('Mature BW, lbs:',mbw.lbs,"\n")
  cat('User Input Peak Milk Production, lbs:',milk.pk.lb,"\n")
  cat('NRC Default Peak Milk Production, lbs:',ifelse(code==1,8*2.2,table19_1[table19_1[,1]==code,6]*2.2),"\n")
  cat('Relative Milk Yield:',RMY,"\n")
  cat('Calf birth weight set at 0.066(dams shrunk BW):',calf.bw.kg*2.2,"\n")
  cat('Adjusted for First Calf Heifer Calving at 2 yr:',first.calf,"\n")
  cat('Daily total Calculations',"\n")
  cat("\n")
  print(data.frame(TDN.lb=round(tdn*dmi*2.2/100,2),
                   NEm.Mcal=total.NEm,
                   NEg.Mcal=round(ifelse(isTRUE(first.calf),total.NEg,NA),1),
                   CP.lb=round(total.CP/454.54,2),
                   Ca.lb=round(ca.total/454.54,3),
                   P.lb=round(p.total/454.54,3)))
  cat("\n")
  cat('Daily Proportion',"\n")
  print(data.frame(BW=round(ifelse(age.yr2==2,sbw.first.calf/0.96*2.2,sbw.cow/0.96*2.2)),
                   BWwPreg=round(ifelse(age.yr2==2,sbw.first.calf/0.96*2.2+(cum.amnthg*2.2),sbw.cow/0.96*2.2+(cum.amnthg*2.2))),
                   TDN.pct=round(tdn,1),
                   DMI.lb=round(dmi*2.2,1),
                   NEma.Mcal.lb=round(NEma/2.2,2),
                   NEg.Mcal.lb=round(ifelse(isTRUE(first.calf),NEga/2.2,NA),2),
                   CP.pct=round(total.CP/1000/dmi*100,1),
                   Ca.pct=round(ca.total/1000/dmi*100,2),
                   P.pct=round(p.total/1000/dmi*100,2)))
}
