#' @title Nutrient Requirements to Increase Cow Body Condition Score
#'
#' @description \code{bcs.gain.cow} Nutrient requirements to increase body condition score of mature cows
#' or previously Calved first-calf heifers during the 3rd trimester following weaning
#'
#' @param breed.code Default is 1 = Angus.  Breed code must match a code listed
#' in the table 19_1 data file.
#' @param mbw.lbs Default is 1200.  This value represents the mature, unshrunk bw
#' (lbs) for cows > 3 yrs old.  For first calf heifers, calving at 2 yr is assumed.
#' mbw.lbs is the estimated mature weight for first calf heifers, and equations are based on
#' a shrunk bw of 0.85(mbw) at first calving and 0.9(mbw) at second calving.
#' @param initial.bcs Default is 4. This value represents the body condition score at
#' the end of the 2nd trimester or at calf weaning.
#' @param target.bcs Default is 6. This value represents the target body condition score
#' at calving.
#' @param first.calf Default is FALSE.  Change to TRUE to adjust requirements for first
#' calf heifers.
#' @return Three tables are printed including nutrient requirement totals, nutrient
#' requirement as a proportion of dry matter intake, and a third table whereby crude
#' protein is adjusted to a ratio of 6:1 for tdn:cp.
#'
#' @examples
#' bcs.gain.cow()
#' bcs.gain.cow(target.bcs=5)
#' bcs.gain.cow(first.calf=TRUE)
#'
#' @export
#'
bcs.gain.cow<-function(
                  breed.code=1,
                  mbw.lbs=1200,
                  initial.bcs=4,
                  target.bcs=6,
                  first.calf=FALSE){

  if(target.bcs <=initial.bcs){
    cat('ERROR: initial BCS >= target BCS')
  }

  #DATE INPUTS
  #mnth 1 starts at first calving;requirement at 30d
  code<-breed.code
  mnth<-rep(1:3)
  mnth.value<-c(7,8,9)
  days.pregnant<-(mnth.value*30)

  #COW INPUTS
  age.yr<-ifelse(isTRUE(first.calf),2,4)
  bcs.chg<-(target.bcs-initial.bcs)/3
  previous.bcs<-c(initial.bcs,initial.bcs+bcs.chg,initial.bcs+bcs.chg*2)
  current.bcs<-c(target.bcs-bcs.chg*2,target.bcs-bcs.chg,target.bcs)

  #INPUTS to METRIC
  bw.kg<-mbw.lbs/2.2
  sbw<-ifelse(isTRUE(first.calf),bw.kg*0.85*0.96,bw.kg*0.96)#first calf heifers 85% at calving
  mbbw<-sbw^0.75
  calf.bw.kg<-ifelse(isTRUE(first.calf),0.9*(bw.kg)*0.96*0.066,sbw*0.066)

  #MAINTENANCE CALC comp creates vector
  be<-ifelse(code==1,1,table19_1[table19_1[,1]==code,3])
  l.fact<-1
  comp<-0.8+(current.bcs-1)*0.05
  sex<-1
  NEm.a1<-0.077
  NEm<-mbbw*(be*l.fact*comp*sex*NEm.a1)  #used for total NEm
  MPm<-mbbw*3.8 #used for total MP


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
  eqsbw<-(sbw-cw)*(srw/msbw)
  eqebw<-0.891*eqsbw
  bw.yr3<-0.9*msbw
  bwg<-(bw.yr3-sbw)/(12*30)
  ewg<-bwg*0.96
  re<-(ewg/(12.341*eqebw^-0.6837))^(1/0.9116)

  #energy and protein reserves
  previous.afat<-previous.bcs*0.037683
  previous.aprot<-0.200886-previous.bcs*0.0066762
  ebw<-0.851*sbw
  waf.adj<-1
  previous.bcs.waf<-1-0.07105*(5-previous.bcs)
  previous.waf.ebw<-ebw/(previous.bcs.waf*waf.adj)
  previous.ebw.bcs<-previous.waf.ebw*previous.bcs.waf
  previous.tfat<-previous.afat*previous.ebw.bcs
  previous.tprot<-previous.aprot*previous.ebw.bcs
  previous.te<-9.4*previous.tfat+5.7*previous.tprot

  current.afat<-current.bcs*0.037683
  current.aprot<-0.200886-current.bcs*0.0066762
  ebw<-0.851*sbw
  waf.adj<-1
  current.bcs.waf<-1-0.07105*(5-current.bcs)
  current.waf.ebw<-ebw/(current.bcs.waf*waf.adj)
  current.ebw.bcs<-current.waf.ebw*current.bcs.waf
  current.tfat<-current.afat*current.ebw.bcs
  current.tprot<-current.aprot*current.ebw.bcs
  current.te<-9.4*current.tfat+5.7*current.tprot

  change.te<-current.te-previous.te
  change.te.d<-change.te/30
  NE.bcgain<-change.te.d


  #Total NEm
  total.NEm<-round(NEm+NEp+NE.bcgain,1)


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
      dmi94[i]<-mbbw*((0.0384+0.04997*start.x94[i]^2)/start.x94[i])
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
      dmi093[i]<-mbbw*((0.0384+0.04997*start.x093[i]^2)/start.x093[i])
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
  NPg<-ifelse(isTRUE(first.calf),bwg*(268-29.4*re/bwg),0)
  MPg<-NPg/max(c(0.492,0.834-0.00114*eqsbw))
  total.MP<-round(MPm+MPp+MPg,1)
  total.CP<-round(total.MP/0.67,1)


  #calcium requirement
  Ca.m<-0.0154*sbw/0.5
  Ca.p<-calf.bw.kg*(13.7/90)/0.5
  Ca.g<-ifelse(isTRUE(first.calf),NPg*0.071/0.5,0)
  ca.total<-Ca.m+Ca.p+Ca.g

  #phosphorus requirement
  P.m<-0.016*sbw/0.68
  P.p<-calf.bw.kg*(7.6/90)/0.68
  P.g<-ifelse(isTRUE(first.calf),NPg*0.039/0.68,0)
  p.total<-P.m+P.g+P.p

  cat('Breed:',ifelse(code==1,'Angus',as.character(table19_1[table19_1[,1]==code,2])),"\n")
  cat('Mature BW, lbs:',mbw.lbs,"\n")
  cat('Calf birth weight set at 0.066(dams shrunk BW):',calf.bw.kg*2.2,"\n")
  cat('Adjusted for First Calf Heifer Calving at 2 yr:',first.calf,"\n")
  cat("\n")
  cat('Daily total Calculations',"\n")
  print(data.frame(Mnth.Since.Calving=mnth.value,
                   TDN.lb=round(tdn*dmi*2.2/100,2),
                   NEm.Mcal=total.NEm,
                   NEbcgain.Mcal=NE.bcgain,
                   NEg.Mcal=round(ifelse(isTRUE(first.calf),total.NEg,NA),1),
                   CP.lb=round(total.CP/454.54,2),
                   Ca.lb=round(ca.total/454.54,3),
                   P.lb=round(p.total/454.54,3)))

  cat("\n")
  cat('Daily Proportion',"\n")
  print(data.frame(Mnth.Since.Calving=mnth.value,
                   Previous.bcs=round(previous.bcs,1),
                   Current.bcs=round(current.bcs,1),
                   TDN.pct=round(tdn,1),
                   DMI.lb=round(dmi*2.2,1),
                   NEma.Mcal.lb=round(NEma/2.2,2),
                   NEg.Mcal.lb=round(ifelse(isTRUE(first.calf),NEga/2.2,NA),2),
                   CP.pct=round(total.CP/1000/dmi*100,1),
                   Ca.pct=round(ca.total/1000/dmi*100,2),
                   P.pct=round(p.total/1000/dmi*100,2)))
  cat("\n")
  cat('Alternate protein as 6:1 TDN:CP',"\n")
  print(data.frame(CP.pct=round(tdn/6,1),
                   CP.lb=round(dmi*2.2*(tdn/6)/100,1)))


}

#'
