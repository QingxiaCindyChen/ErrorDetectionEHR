data.all=data.table(
  line,     # Ordering of data measurements 
  subjid,   # ID of the participant
  param,    # Height (cm) or weight (kg)
  meas_id,  # Measurement ID 
  agedays,  # Age (days)
  ageyears, # Age (years)
  v,        # Measurement
  sex,      # Sex at birth (0 if male, 1 if female, 2 if another sex) 
  concept,  # Concept ID of height/weight used
  date,     # Measurement Date
  z.orig,   # Z-score
  sd.orig)  # MSDS score

# Assign z-scores and MSDS scores
if (param == "HEIGHTCM") {
  anthro = as.data.table(rbind(
    with(CDCdat, data.frame(src='CDC', param='HEIGHTCM', sex=sex, mu=ht_mu, sigma=ht_sigma, nu=ht_nu, tau=ht_tau, gampos=ht_gampos, gamneg=ht_gamneg))
  ))
}
if (param == "WEIGHTKG") {
  anthro = as.data.table(rbind(
    with(CDCdat, data.frame(src='CDC', param='WEIGHTKG', sex=sex, mu=ht_mu, sigma=ht_sigma, nu=ht_nu, tau=ht_tau, gampos=ht_gampos, gamneg=ht_gamneg))
  ))
}
measurement.to.z <- function(param, sex, measurement, csd=F){
  src = 'CDC'
  dt = data.table(src,param,sex,measurement)
  dt$sex <- factor(dt$sex, levels=c(0, 1, 2), labels=c("m", "f", "o"))
  dt = anthro[dt]
  
  dt[, ret := as.double(NA)]
  dt[, tmp := as.double(NA)]
  dt[measurement>0, tmp:=pBCT(measurement, mu=mu, sigma=sigma, nu=nu, tau=tau)]
  dt[measurement>0, ret:=qnorm(tmp)]
  dt[measurement>0 & is.infinite(ret), ret := qnorm(1-1e-16)]
  dt[measurement>0 & is.na(ret) & tmp>0, ret := qnorm(1-1e-16)]
  
  if(csd) {
    dt[nu<1 & measurement<mu & measurement>0, ret:= (measurement - mu)/gamneg]
    dt[nu<1 & measurement>=mu, ret:= (measurement - mu)/gampos]
  }
  
  return(dt$ret)
}

# Exponentially Weighted Moving Average (EWMA)
ewma = function(agedays, z, ewma.exp, ewma.adjacent=T) {
  n = length(agedays)
  ewma.all = ewma.before = ewma.after = vector('numeric',0)
  if(n>0) {
    index = order(agedays)
    data = data.table(agedays, z, key='agedays')
    
    delta = as.matrix.delta(data$agedays)
    delta = ifelse(delta == 0, 0, (delta + 5)^ewma.exp)
    
    ewma.all[index] = delta %*% data$z / apply(delta,1,sum)
    
    if(ewma.adjacent) {
      if (n>2) {
        delta2 = delta
        delta2[col(delta2)==row(delta2)-1]=0
        ewma.before[index] = delta2 %*% data$z / apply(delta2,1,sum)
        delta3 = delta
        delta3[col(delta3)==row(delta3)+1]=0
        ewma.after[index] = delta3 %*% data$z / apply(delta3,1,sum)
      } else {
        ewma.before = ewma.after = ewma.all
      }
    }
  }
  return(if(ewma.adjacent) data.frame(ewma.all,ewma.before,ewma.after) else data.frame(ewma.all))
}
as.matrix.delta = function(agedays) {
  n = length(agedays)
  delta = matrix(0,n,n)
  for(i in 1:n) {
    delta[i,] = abs(agedays - agedays[i])
  }
  return(delta)
}
sd.median = function(param, sex, agedays, sd.orig) {
  dt = data.table(param, sex, agedays, ageyears = floor(agedays/365.25), sd.orig)
  setkey(dt, param, sex, agedays)
  agedays.to.cover = with(dt, floor(min(ageyears)*365.25):floor(max(ageyears+1)*365.25))
  #dt[ageyears > 19, ageyears := 19]
  dt.median = dt[!is.na(sd.orig), list(sex=c(0,1,2),sd.median=rep(median(sd.orig),3)), by=.(param,ageyears)]
  dt.median=dt.median[, list(agedays=agedays.to.cover, sd.median=approx(floor((ageyears + 0.5)*365.25), sd.median, xout=agedays.to.cover, rule=2)$y), by=.(param,sex)]
  setkey(dt.median, param, sex, agedays)
  return(dt.median)
}
if(!is.data.table(sd.recenter)) {
  sd.recenter = data.all[exclude < 'Exclude', sd.median(param, sex, agedays, sd.orig)]
}
setkey(data.all, param, sex, agedays)
data.all = sd.recenter[data.all]
setkey(data.all, subjid, param, agedays)
data.all[, tbc.sd := sd.orig - sd.median]

# Error Detection Algorithm
classify <- function(data) {
  data.df = data.table(data, key=c('subjid','param','agedays','index'))
  if(!quietly & .parallel) sink(sprintf("cleangrowth-%s-batch-%s.log", Sys.Date(), data.df$batch[1]))
  if(!quietly) cat(sprintf("[%s] Processing Batch #%s...\n", Sys.time(), data.df$batch[1]))
  
  # Extra Functions (these functions assign temporary duplicates and subset later error detections to only those still eligible for error assignment)
  na.as.false = function(v) {
    v[is.na(v)] = F
    v
  }
  temporary.duplicates = function(df = data.df) {
    if(is.null(df$subjid)) df[, subjid := NA]
    if(is.null(df$param)) df[, param := NA]
    valid.rows = valid(df, include.temporary.duplicates=T)
    df = df[j=.(tbc.sd),keyby=.(subjid,param,agedays,index)]
    df[, `:=`(median.sd=as.double(NA), delta.median.sd=as.double(NA), duplicates.this.day=F, duplicate=F)]
    df[valid.rows, duplicates.this.day := (.N > 1), by=.(subjid,param,agedays)]
    df[valid.rows & !duplicates.this.day, median.sd := median(tbc.sd), by=.(subjid,param)]
    df[valid.rows, median.sd := sort(median.sd)[1], by=.(subjid,param)]
    df[valid.rows & duplicates.this.day, delta.median.sd := abs(tbc.sd - median.sd), by=.(subjid, param)]
    subj.all.dups = df[valid.rows & duplicates.this.day & is.na(delta.median.sd), unique(subjid)]
    df[valid.rows & subjid %in% subj.all.dups, delta.median.sd := (function(subj.df) {
      for(p in subj.df[is.na(delta.median.sd) & duplicates.this.day, unique(param)]) {
        median.other.sd = subj.df[param != p & !duplicates.this.day, median(tbc.sd)]
        subj.df[param == p, delta.median.sd := abs(tbc.sd - median.other.sd)]
      }
      return(subj.df$delta.median.sd)
    })(copy(.SD)), by=.(subjid)]
    df[valid.rows & duplicates.this.day & is.na(delta.median.sd), delta.median.sd := abs(tbc.sd)]
    df[valid.rows & duplicates.this.day, duplicate := seq_along(delta.median.sd) != which.min(delta.median.sd), by=.(subjid, param, agedays)]
    return(df$duplicate & valid.rows)
  }
  valid = function(df=data.df, include.temporary.duplicates=F, include.duplicates=F, include.carryforward=F) {
    if(is.data.frame(df)){
      exclude <- df$exclude
      flag <- df$flag
      out <- (exclude < 'Exclude' & flag=='None')
      return(out 
             | include.temporary.duplicates & flag == 'Warning-Temporary-Duplicate'
             | include.duplicates & flag == 'Warning-Repeated-Measure'
             | include.carryforward & flag == 'Warning-Carried-Forward')
    }else {
      exclude <- df
      out <- (exclude < 'Exclude')
      return(out)
    }
  }
  
  
  #Step 1 - Assign temporary duplicates
  data.df[, v.orig := v]
  if(!quietly) cat(sprintf("[%s] Preliminarily identify potential duplicates...\n", Sys.time()))
  data.df$flag[temporary.duplicates()] = 'Warning-Temporary-Duplicate'
  subj.dup = data.df[flag == 'Warning-Temporary-Duplicate', unique(subjid)]
  if(!quietly) cat(sprintf("[%s] End Step 1\n", Sys.time()))
  
  #Step 2 - Assign measurements carried forward
  if(!quietly) cat(sprintf("[%s] Flag measurements carried forward...\n", Sys.time()))
  data.df[, prev.v := as.double(NaN)]
  data.df[valid(), prev.v := c(NA, v.orig[-.N]), by=.(subjid,param)]
  data.df[!(subjid %in% subj.dup) & v.orig == prev.v, flag := 'Warning-Carried-Forward']
  data.df[subjid %in% subj.dup & valid(include.temporary.duplicates=T), flag := (function(df) {
    setkey(df,agedays)
    ages = unique(agedays)
    if (length(ages) > 1) {
      for(i in 2:length(ages)) {
        all.prev.v = df[agedays == ages[i-1], v.orig]
        df[agedays == ages[i] & v.orig %in% all.prev.v, flag := 'Warning-Carried-Forward']
      }
    }
    return(df$flag)
  })(copy(.SD)), .SDcols=c('agedays','flag','v.orig'), by=.(subjid,param)]
  data.df[flag == 'Warning-Carried-Forward', carried.forward := 1]
  
  
  data.df[flag == 'Warning-Temporary-Duplicate', flag := 'None']
  data.df[temporary.duplicates(), flag := 'Warning-Temporary-Duplicate']
  if(!quietly) cat(sprintf("[%s] End Step 2\n", Sys.time()))
  
  #Step 3 - Assign extreme measures based on MSDS score
  if(!quietly) cat(sprintf("[%s] Exclude extreme measurements based on SD...\n", Sys.time()))
  data.df[na.as.false(valid(include.temporary.duplicates=T) & abs(tbc.sd) > sd.extreme 
                      | flag %in% c('None', 'Warning-Temporary-Duplicate', 'Warning-Carried-Forward') & abs(z.orig) > z.extreme),
          exclude := 'Exclude-SD-Cutoff']
  data.df[flag == 'Warning-Temporary-Duplicate', flag := 'None']
  data.df[temporary.duplicates(), flag := 'Warning-Temporary-Duplicate']    
  if(!quietly) cat(sprintf("[%s] End Step 3\n", Sys.time()))
  
  #Step 4 - Assign extreme measures based on EWMA
  if(!quietly) cat(sprintf("[%s] Exclude extreme measurements based on EWMA...\n", Sys.time()))
  data.df = data.df[, exclude := (function(df) {
    num.ewma.excluded = 0
    has.duplicates = subjid %in% subj.dup
    while(T) {
      df[, (ewma.fields) := as.double(NaN)]
      df[valid(exclude), (ewma.fields) := ewma(agedays, tbc.sd, ewma.exp, T)]
      if(has.duplicates) {
        df[, `:=`(ewma.all=sort(ewma.all)[1], ewma.before=sort(ewma.before)[1], ewma.after=sort(ewma.after)[1]), by=.(agedays)]
      }
      df[, `:=`(dewma.all=tbc.sd - ewma.all, dewma.before=tbc.sd - ewma.before, dewma.after=tbc.sd - ewma.after)]
      num.valid = sum(valid(df))
      rep = na.as.false(with(df, num.valid >= 3 & valid(df)
                             & (dewma.all > 3.5 & dewma.before > 3 & dewma.after > 3 & tbc.sd > 3.5
                                | dewma.all < -3.5 & dewma.before < -3 & dewma.after < -3 & tbc.sd < -3.5)))
      num.exclude = sum(rep)
      if(num.exclude == 1) df[rep, exclude := 'Exclude-EWMA-Extreme']
      if(num.exclude > 1) {
        worst.row = with(df, order(rep, abs(tbc.sd + (tbc.sd - ewma.all)), decreasing=T))[1]
        df[worst.row, exclude := 'Exclude-EWMA-Extreme']
      }
      rep = na.as.false(with(df, num.valid == 2 & (tbc.sd - ewma.all > 3.5 & tbc.sd > 3.5 | tbc.sd - ewma.all < -3.5 & tbc.sd < -3.5)))
      num.exclude = sum(rep)
      if(num.exclude == 1) df[rep, exclude := 'Exclude-EWMA-Extreme-Pair']
      if(num.exclude > 1) {
        worst.row = with(df, order(rep, abs(tbc.sd), decreasing=T))[1]
        df[worst.row, exclude := 'Exclude-EWMA-Extreme-Pair']
      }
      if(has.duplicates) {
        df[flag == 'Warning-Temporary-Duplicate', flag := 'None']
        df[temporary.duplicates(df), flag := 'Warning-Temporary-Duplicate']
      }
      newly.excluded = sum(df$exclude %in% c('Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair'))
      if(newly.excluded > num.ewma.excluded) {
        num.ewma.excluded = newly.excluded
      } else {
        break
      }
    }  
    return(df$exclude)
  })(copy(.SD)), by=.(subjid, param), .SDcols=c('index', 'sex', 'agedays', 'tbc.sd', 'exclude', 'flag')]
  if(!quietly) cat(sprintf("[%s] End Step 4\n", Sys.time()))
  
  #Step 5 - Assign official repeated measures
  if(!quietly) cat(sprintf("[%s] Exclude duplicates based on EWMA...\n", Sys.time()))
  data.df[flag == 'Warning-Temporary-Duplicate', flag := 'None']
  temp.dups = temporary.duplicates()
  data.df[temp.dups, flag := 'Warning-Temporary-Duplicate'] 
  valid.rows = valid()
  data.df[, `:=`(ewma.all=as.double(NaN), abssum2=as.double(NaN), median.other.sd=as.double(NaN), duplicate=F)]
  dup.ratio.df = data.df[subjid %in% subj.dup & (valid.rows | temp.dups), list(dup=(.N>1)), by=.(subjid, param, agedays)][j=list(dup.ratio=mean(dup)),keyby=.(subjid, param)]
  subj.param.not.all.dups = dup.ratio.df[dup.ratio<1.0,list(subjid,param)]
  subj.param.all.dups = dup.ratio.df[dup.ratio == 1, list(subjid,param)]
  subj.all.dups = subj.param.all.dups[,unique(subjid)]
  data.df[subjid %in% subj.dup & valid.rows, ewma.all := ewma(agedays, tbc.sd, ewma.exp, ewma.adjacent=F), by=.(subjid, param)]
  data.df[subjid %in% subj.dup, ewma.all := sort(ewma.all)[1], by=.(subjid, param, agedays)]
  data.df[, abssum2 := 2*abs(tbc.sd - ewma.all) + abs(tbc.sd)]
  data.df[J(subj.param.not.all.dups), duplicate := seq_along(abssum2) != which.min(abssum2), by=.(subjid, param, agedays)]
  data.df[temp.dups, flag := 'None']
  data.df[(valid.rows | temp.dups) & duplicate, flag := 'Warning-Repeated-Measure']
  data.df[J(dup.ratio.df[dup.ratio > 1/2,list(subjid,param)]), flag := (function(df) {
    df[,`:=`(tbc.sd.min=as.double(NaN), tbc.sd.max=as.double(NaN))]
    df[valid(data.frame(exclude,flag), include.duplicates=T, include.temporary.duplicates=T), `:=`(tbc.sd.min=min(tbc.sd), tbc.sd.max=max(tbc.sd))]
    df[tbc.sd.max - tbc.sd.min > ifelse(param=='HEIGHTCM', 3, ifelse(param=='WEIGHTKG', ifelse(tbc.sd.min < 10, 0.25, ifelse(tbc.sd.min < 30, 0.5, 1)), NA)),
       flag := 'Warning-Repeated-Measure']
    return(df$flag)
  })(copy(.SD)), .SDcols=c('exclude','flag','tbc.sd'), by=.(subjid,param,agedays)]
  
  data.df[subjid %in% subj.all.dups, flag := (function(subj.df) { 
    subj.df[,`:=`(duplicates.this.day=F, duplicate=F, tbc.sd.min=as.double(NaN), tbc.sd.max=as.double(NaN))]
    valid.rows = valid(df=subj.df, include.duplicates=T, include.temporary.duplicates=T, include.carryforward=F)
    subj.df[valid.rows, duplicates.this.day := (.N > 1), by=.(param,agedays)]
    for(p in subj.df[j=unique(param)]) {
      median.sd = subj.df[param != p & !duplicates.this.day, median(tbc.sd)]
      subj.df[param == p, median.other.sd := median.sd]
    }
    subj.df[is.na(median.other.sd), median.other.sd := 0]
    subj.df[duplicates.this.day==T, duplicate := (seq_along(median.other.sd) != which.min(abs(tbc.sd - median.other.sd))), by=.(param,agedays)]
    subj.df[duplicates.this.day & !duplicate, flag := 'None']
    subj.df[duplicates.this.day &  duplicate, flag := 'Warning-Repeated-Measure']
    subj.df[duplicates.this.day==T, `:=`(tbc.sd.min=min(tbc.sd), tbc.sd.max=max(tbc.sd)), by=.(param,agedays)]
    subj.df[tbc.sd.max - tbc.sd.min > ifelse(param=='HEIGHTCM', 3, ifelse(param=='WEIGHTKG', ifelse(tbc.sd.min < 10, 0.25, ifelse(tbc.sd.min < 30, 0.5, 1)), NA)),
            flag := 'Warning-Repeated-Measure']
    subj.df[, duplicates.this.day := F]
    subj.df[exclude != 'Missing', duplicates.this.day := (.N > 1), by=.(param,agedays)]
    subj.df[duplicates.this.day & exclude %in% c('Exclude-SD-Cutoff', 'Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair'), flag := 'Warning-Repeated-Measure']
    return(subj.df$flag)
  })(copy(.SD)), .SDcols=c('param', 'agedays','exclude','flag','tbc.sd'), by=.(subjid)]
  
  data.df[subjid %in% subj.dup, flag := (function(subj.df) {
    if(.N > 1) {
      subj.df[exclude %in% c('Exclude-SD-Cutoff', 'Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair'), flag := 'Warning-Repeated-Measure']
    }
    return(subj.df$flag)
  })(copy(.SD)), .SDcols=c('exclude','flag'), by=.(subjid, param, agedays)]
  if(!quietly) cat(sprintf("[%s] End Step 5\n", Sys.time()))
  
  
  valid = function(df=data.df, include.temporary.duplicates=F, include.duplicates=F, include.carryforward=F) {
    exclude = if(is.data.frame(df)) df$exclude else df
    return(exclude < 'Exclude' 
           | include.temporary.duplicates & exclude == 'Exclude-Temporary-Duplicate'
           | include.duplicates & exclude == 'Exclude-Duplicate'
           | include.carryforward & exclude == 'Exclude-Carried-Forward')
  }
  
  #Step 6
  data.df[, delta := ifelse(param == 'WEIGHTKG', .05*v, 1)]
  data.df[, `:=`(v.minus=v - delta, v.plus=v + delta)]
  if (all(data.df$param == "HEIGHTCM")) {
    data.df[, `:=`(tbc.sd.minus=measurement.to.z(param, sex, v.minus, T),
                   tbc.sd.plus=measurement.to.z(param, sex, v.plus, T))]
  }
  
  if (all(data.df$param == "WEIGHTKG")) {
    data.df[, `:=`(tbc.sd.minus=measurement.to.z(param, sex, v.minus, T),
                   tbc.sd.plus=measurement.to.z(param, sex, v.plus, T))]
  }
  
  if(!quietly) cat(sprintf("[%s] End Step 6\n", Sys.time()))
  
  #Step 7 - Assign moderate, inlier EWMA errors based on individual trajectories
  if(!quietly) cat(sprintf("[%s] Exclude moderate errors based on EWMA...\n", Sys.time()))
  
  
  data.df[, exclude := (function(subj.df) {
    num.ewma.excluded = 0
    while(T) {
      valid.rows = valid(subj.df)
      subj.df[, `:=`(ewma.all=as.double(NaN), ewma.before=as.double(NaN), ewma.after=as.double(NaN), tbc.sd.prev=as.double(NaN), tbc.sd.next=as.double(NaN),
                     delta.agedays.prev=as.integer(NaN), delta.agedays.next=as.integer(NaN), abs.2ndlast.sd=as.double(NaN),
                     tbc.other.sd=as.double(NaN))]
      subj.df[valid.rows, (ewma.fields) := ewma(agedays, tbc.sd, ewma.exp, T), by=param]
      subj.df[, `:=`(dewma.all=tbc.sd - ewma.all, dewma.before=tbc.sd - ewma.before, dewma.after=tbc.sd - ewma.after)]
      subj.df[valid.rows, `:=`(tbc.sd.prev=c(NA, tbc.sd[-.N]), tbc.sd.next=c(tbc.sd[-1], NA)), by=param]
      subj.df[, `:=`(dprev.sd=tbc.sd - tbc.sd.prev, dprev.sd.minus=tbc.sd.minus - tbc.sd.prev, dprev.sd.plus=tbc.sd.plus - tbc.sd.prev,
                     dnext.sd=tbc.sd - tbc.sd.next, dnext.sd.minus=tbc.sd.minus - tbc.sd.next, dnext.sd.plus=tbc.sd.plus - tbc.sd.next)]
      subj.df[valid.rows, `:=`(delta.agedays.prev=as.integer(agedays - c(NA, agedays[-.N])),
                               delta.agedays.next=as.integer(agedays - c(agedays[-1], NA))), by=param]
      subj.df[valid.rows, abs.2ndlast.sd := ifelse(.N >= 2, abs(tbc.sd[.N-1]), as.double(NaN)), by=param]
      valid.other=subj.df[valid.rows,list(param.other=ifelse(param=='WEIGHTKG', 'HEIGHTCM', 'WEIGHTKG'),agedays.other=agedays,tbc.sd)]
      setkey(valid.other,param.other,agedays.other)
      subj.df[valid.rows, tbc.other.sd := valid.other[list(param,agedays),tbc.sd]]
      subj.df[, exclude := (function(df) {
        median.tbc.other.sd=df[j=median(tbc.other.sd, na.rm=T)]
        num.valid=df[j=sum(valid(exclude))]
        df[, temp.exclude := factor(NA, levels=exclude.levels, ordered=T)]
        df[dewma.all > 1 & dewma.before > 1 & dewma.after > 1 & dprev.sd > 1 & dprev.sd.plus > 1 & dprev.sd.minus > 1 & dnext.sd > 1 &  dnext.sd.plus > 1 & dnext.sd.minus > 1
           | dewma.all < -1 & dewma.before < -1 & dewma.after < -1 & dprev.sd < -1 & dprev.sd.plus < -1 & dprev.sd.minus < -1 & dnext.sd < -1 &  dnext.sd.plus < -1 & dnext.sd.minus < -1,
           temp.exclude := 'Exclude-EWMA-8']
        df[, `:=`(first.of.three.or.more=F, last.of.three.or.more=F)]
        df[is.na(delta.agedays.prev) & num.valid >= 3, first.of.three.or.more := T]
        df[is.na(delta.agedays.next) & num.valid >= 3, last.of.three.or.more := T]
        df[first.of.three.or.more & delta.agedays.next < 365.25
           & (dewma.all > 2 & dewma.after > 1 & dnext.sd > 1 &  dnext.sd.plus > 1 & dnext.sd.minus > 1
              | dewma.all < -2 & dewma.after < -1 & dnext.sd < -1 &  dnext.sd.plus < -1 & dnext.sd.minus < -1),
           temp.exclude := 'Exclude-EWMA-9']
        df[first.of.three.or.more & delta.agedays.next >= 365.25
           & (dewma.all > 3 & dewma.after > 1 & dnext.sd > 1 &  dnext.sd.plus > 1 & dnext.sd.minus > 1
              | dewma.all < -3 & dewma.after < -1 & dnext.sd < -1 &  dnext.sd.plus < -1 & dnext.sd.minus < -1),
           temp.exclude := 'Exclude-EWMA-10']
        df[last.of.three.or.more & delta.agedays.prev < 730.5 & abs.2ndlast.sd < 2
           & (dewma.all > 2 & dewma.before > 1 & dprev.sd > 1 &  dprev.sd.plus > 1 & dprev.sd.minus > 1
              | dewma.all < -2 & dewma.before < -1 & dprev.sd < -1 &  dprev.sd.plus < -1 & dprev.sd.minus < -1),
           temp.exclude := 'Exclude-EWMA-11']
        df[last.of.three.or.more & delta.agedays.prev < 730.5 & abs.2ndlast.sd >= 2
           & (dewma.all > abs.2ndlast.sd & dewma.before > 1 & dprev.sd > 1 &  dprev.sd.plus > 1 & dprev.sd.minus > 1
              | dewma.all < -abs.2ndlast.sd & dewma.before < -1 & dprev.sd < -1 &  dprev.sd.plus < -1 & dprev.sd.minus < -1),
           temp.exclude := 'Exclude-EWMA-12']
        df[last.of.three.or.more & delta.agedays.prev >= 730.5 & abs.2ndlast.sd < 2
           & (dewma.all > 3 & dewma.before > 1 & dprev.sd > 1 &  dprev.sd.plus > 1 & dprev.sd.minus > 1 
              & (tbc.sd - tbc.other.sd > 4 | is.na(tbc.other.sd) & tbc.sd - median.tbc.other.sd > 4 | is.na(median.tbc.other.sd))
              | dewma.all < -3 & dewma.before < -1 & dprev.sd < -1 &  dprev.sd.plus < -1 & dprev.sd.minus < -1 
              & (tbc.sd - tbc.other.sd < -4 | is.na(tbc.other.sd) & tbc.sd - median.tbc.other.sd < -4 | is.na(median.tbc.other.sd))),
           temp.exclude := 'Exclude-EWMA-13']
        df[last.of.three.or.more & delta.agedays.prev >= 730.5 & abs.2ndlast.sd >= 2
           & (dewma.all > 1 + abs.2ndlast.sd & dewma.before > 1 & dprev.sd > 1 &  dprev.sd.plus > 1 & dprev.sd.minus > 1
              & (tbc.sd - tbc.other.sd > 4 | is.na(tbc.other.sd) & tbc.sd - median.tbc.other.sd > 4 | is.na(median.tbc.other.sd))
              | dewma.all < -1 - abs.2ndlast.sd & dewma.before < -1 & dprev.sd < -1 &  dprev.sd.plus < -1 & dprev.sd.minus < -1 
              & (tbc.sd - tbc.other.sd < -4 | is.na(tbc.other.sd) & tbc.sd - median.tbc.other.sd < -4 | is.na(median.tbc.other.sd))),
           temp.exclude := 'Exclude-EWMA-14']
        rep = !is.na(df$temp.exclude)
        num.exclude = sum(rep)
        if(num.exclude == 1) df[rep, exclude := temp.exclude]
        if(num.exclude > 1) {
          worst.row = with(df, order(rep, abs(tbc.sd + dewma.all), decreasing=T))[1]
          df[worst.row, exclude := temp.exclude]
        }
        return(df$exclude)
      })(copy(.SD)),by=param]
      newly.excluded = sum(subj.df$exclude >= 'Exclude-EWMA-8' & subj.df$exclude <= 'Exclude-EWMA-14')
      if(newly.excluded > num.ewma.excluded) {
        num.ewma.excluded = newly.excluded
      } else {
        break
      }
    }
    return(subj.df$exclude)
  })(copy(.SD)),by=subjid, .SDcols=c('index', 'sex', 'param', 'agedays', 'v', 'tbc.sd', 'v.minus', 'v.plus', 'tbc.sd.minus', 'tbc.sd.plus', 'exclude')]
  if(!quietly) cat(sprintf("[%s] End Step 7\n", Sys.time()))
  
  #Step 8 - Exclude heights based on large changes between consecutive measures (skip section if measures are weights)
  if(!quietly) cat(sprintf("[%s] Exclude heights based on growth velocity...\n", Sys.time()))
  data.df[param == 'HEIGHTCM', exclude := (function(subj.df) {
    subj.df[, `:=` (index=1:.N, sd.orig=sd.orig, sd.median=sd.median)]
    num.height.excluded = 0
    while(T) {
      subj.df[valid(exclude), exclude := (function (df) {
        
        df[, (ewma.fields) := as.double(NaN)]
        df[, `:=`(v.prev=as.double(NaN), v.next=as.double(NaN), dewma.after.prev=as.double(NaN), dewma.before.next=as.double(NaN),
                  abs.tbc.sd.prev=as.double(NaN), abs.tbc.sd.next=as.double(NaN), agedays.next=as.integer(NaN), abs.2ndlast.sd=as.double(NaN),
                  mindiff.prev.ht=as.double(NaN), mindiff.next.ht=as.double(NaN), maxdiff.prev.ht=as.double(NaN), maxdiff.next.ht=as.double(NaN),
                  pair.prev=F, pair.next=F)]
        df[, (ewma.fields) := ewma(agedays, tbc.sd, ewma.exp, T)]# calculate some usefule values (e.g. dewma values and tbc.sd) for use in later steps
        df[, `:=`(dewma.all=tbc.sd - ewma.all, dewma.before=tbc.sd - ewma.before, dewma.after=tbc.sd - ewma.after, abs.tbc.sd=abs(tbc.sd))]
        df[, `:=`(agedays.next=c(agedays[-1],NA), v.next=c(v[-1],NA), dewma.before.next=c(dewma.before[-1],NA), abs.tbc.sd.next=c(abs.tbc.sd[-1],NA))]
        df[, delta.agedays.next := agedays.next - agedays]
        df[, mid.agedays := 0.5*(agedays.next + agedays)]
        
        df[, mindiff.next.ht := -8]
        df[, maxdiff.next.ht := 8]
        df[, maxdiff.prev.ht := 8]
        df[, `:=`(v.prev=c(NA, v[-.N]), dewma.after.prev=c(NA, dewma.after[-.N]), abs.tbc.sd.prev=c(NA, abs.tbc.sd[-.N]),
                  mindiff.prev.ht=c(NA, mindiff.next.ht[-.N]), maxdiff.prev.ht=c(NA, maxdiff.next.ht[-.N]))]
        df[, `:=`(delta.prev.ht= v - v.prev,
                  delta.next.ht=v.next - v)]
        df[, pair := na.as.false(delta.prev.ht < mindiff.prev.ht | delta.next.ht < mindiff.next.ht | delta.prev.ht > maxdiff.prev.ht | delta.next.ht > maxdiff.next.ht)]
        df[, `:=`(pair.prev=c(F, pair[-.N]), pair.next=c(pair[-1], F))]
        df[, `:=`(bef.g.aftm1=na.as.false(abs(dewma.before) > abs(dewma.after.prev)  & pair & pair.prev),
                  aft.g.befp1=na.as.false(abs(dewma.after)  > abs(dewma.before.next) & pair & pair.next))]
        num.valid = .N
        df[, temp.diff := as.double(NaN)]
        df[, temp.exclude := factor(NA, levels=exclude.levels, ordered=T)]
        
        df[delta.prev.ht < mindiff.prev.ht & bef.g.aftm1,
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Height-Change')]
        df[delta.next.ht < mindiff.next.ht & aft.g.befp1,
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Max-Height-Change')]
        
        df[delta.prev.ht > maxdiff.prev.ht & bef.g.aftm1,
           `:=`(temp.diff=abs(dewma.before), temp.exclude='Exclude-Max-Height-Change')]
        df[delta.next.ht > maxdiff.next.ht & aft.g.befp1,
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Min-Height-Change')]
        df[delta.next.ht > 8 & is.na(delta.prev.ht),
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Height-Change')]
        
        df[delta.next.ht > 8 & agedays == agedays[1],
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Height-Change')]
        
        df[delta.prev.ht > 8 & is.na(delta.next.ht),
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Max-Height-Change')]
        
        df[delta.prev.ht > 8 & agedays == agedays[.N],
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Max-Height-Change')]
        
        df[delta.prev.ht < mindiff.prev.ht & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev
           | delta.next.ht < mindiff.next.ht & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next,
           temp.exclude := 'Exclude-Min-Height-Change']
        df[delta.prev.ht > maxdiff.prev.ht & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev
           | delta.next.ht > maxdiff.next.ht & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next,
           temp.exclude := 'Exclude-Max-Height-Change']         
        rep = !is.na(df$temp.exclude)
        num.exclude = sum(rep)
        if(num.exclude == 1) df[rep, exclude := temp.exclude]
        if(num.exclude > 1) {
          worst.row = order(rep, df$temp.diff, decreasing=T)[1]
          df[worst.row, exclude := temp.exclude]
        }
        
        return(df$exclude)
      })(copy(.SD))]
      newly.excluded = sum(subj.df$exclude %in% c('Exclude-Min-Height-Change', 'Exclude-Max-Height-Change'))
      if(newly.excluded > num.height.excluded) {
        num.height.excluded = newly.excluded
      } else {
        break
      }
    }  
    setkey(subj.df,index)
    return(subj.df$exclude)
  })(copy(.SD)), by=.(subjid),.SDcols=c('sex', 'agedays', 'v', 'tbc.sd', 'exclude')]
  if(!quietly) cat(sprintf("[%s] End Step 8\n", Sys.time()))
  
  #Step 9 - Exclude weights based on large changes between consecutive measures (skip section if measures are heights)
  if(!quietly) cat(sprintf("[%s] Exclude weights based on implausible change within 180 days...\n", Sys.time()))
  data.df[param == 'WEIGHTKG', exclude := (function(subj.df) {
    subj.df[, index:=1:.N]
    num.weight.excluded = 0
    while(T) {
      subj.df[valid(exclude), exclude := (function (df) {
        
        df[, (ewma.fields) := as.double(NaN)]
        df[, `:=`(v.prev=as.double(NaN), v.next=as.double(NaN), dewma.after.prev=as.double(NaN), dewma.before.next=as.double(NaN),
                  abs.tbc.sd.prev=as.double(NaN), abs.tbc.sd.next=as.double(NaN), agedays.next=as.integer(NaN), abs.2ndlast.sd=as.double(NaN),
                  mindiff.prev.wt=as.double(NaN), mindiff.next.wt=as.double(NaN), maxdiff.prev.wt=as.double(NaN), maxdiff.next.wt=as.double(NaN),
                  pair.prev=F, pair.next=F)]
        df[, (ewma.fields) := ewma(agedays, tbc.sd, ewma.exp, T)]
        df[, `:=`(dewma.all=tbc.sd - ewma.all, dewma.before=tbc.sd - ewma.before, dewma.after=tbc.sd - ewma.after, abs.tbc.sd=abs(tbc.sd))]
        
        df[, `:=`(agedays.next=c(agedays[-1],NA), agedays.prev=c(NA,agedays[-.N]), v.next=c(v[-1],NA), dewma.before.next=c(dewma.before[-1],NA), abs.tbc.sd.next=c(abs.tbc.sd[-1],NA))]
        df[, delta.agedays.next := agedays.next - agedays]
        df[, delta.agedays.prev := agedays - agedays.prev]
        
        df[, mindiff.next.wt := -15]
        
        df[, maxdiff.next.wt := 15]
        df[, maxdiff.prev.wt := 15]
        
        df[, `:=`(v.prev=c(NA, v[-.N]), dewma.after.prev=c(NA, dewma.after[-.N]), abs.tbc.sd.prev=c(NA, abs.tbc.sd[-.N]),
                  mindiff.prev.wt=c(NA, mindiff.next.wt[-.N]), maxdiff.prev.wt=c(NA, maxdiff.next.wt[-.N]))]
        df[, `:=`(delta.prev.wt= v - v.prev,
                  delta.next.wt=v.next - v)]
        df[, pair := na.as.false(delta.prev.wt < mindiff.prev.wt | delta.next.wt < mindiff.next.wt | delta.prev.wt > maxdiff.prev.wt | delta.next.wt > maxdiff.next.wt)]
        df[, `:=`(pair.prev=c(F, pair[-.N]), pair.next=c(pair[-1], F))]
        df[, `:=`(bef.g.aftm1=na.as.false(abs(dewma.before) > abs(dewma.after.prev)  & pair & pair.prev),
                  aft.g.befp1=na.as.false(abs(dewma.after)  > abs(dewma.before.next) & pair & pair.next))]
        num.valid = .N
        df[, temp.diff := as.double(NaN)]
        df[, temp.exclude := factor(NA, levels=exclude.levels, ordered=T)]
        
        df[delta.prev.wt < mindiff.prev.wt & bef.g.aftm1 & (delta.agedays.prev < 180),
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Weight-Change')]
        df[delta.next.wt < mindiff.next.wt & aft.g.befp1 & (delta.agedays.next < 180),
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Min-Weight-Change')]
        df[delta.next.wt > 15 & is.na(delta.prev.wt) & (delta.agedays.next < 180),
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Weight-Change')]
        
        df[delta.next.wt > 15 & agedays == agedays[1] & (delta.agedays.next < 180),
           `:=`(temp.diff=abs(dewma.before), temp.exclude = 'Exclude-Min-Weight-Change')]
        
        df[delta.prev.wt > 15 & is.na(delta.next.wt) & (delta.agedays.prev < 180),
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Min-Weight-Change')]
        
        df[delta.prev.wt > 15 & agedays == agedays[.N] & (delta.agedays.prev < 180),
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Min-Weight-Change')]
        
        df[delta.prev.wt > maxdiff.prev.wt & bef.g.aftm1 & (delta.agedays.prev < 180),
           `:=`(temp.diff=abs(dewma.before), temp.exclude='Exclude-Max-Weight-Change')]
        df[delta.next.wt > maxdiff.next.wt & aft.g.befp1 & (delta.agedays.next < 180),
           `:=`(temp.diff=abs(dewma.after), temp.exclude = 'Exclude-Max-Weight-Change')]
        df[delta.prev.wt < mindiff.prev.wt & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev & delta.agedays.prev < 180
           | delta.next.wt < mindiff.next.wt & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next & delta.agedays.next < 180,
           temp.exclude := 'Exclude-Min-Weight-Change']
        df[delta.prev.wt > maxdiff.prev.wt & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev & delta.agedays.prev < 180
           | delta.next.wt > maxdiff.next.wt & num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next & delta.agedays.next < 180,
           temp.exclude := 'Exclude-Max-Weight-Change']            
        rep = !is.na(df$temp.exclude)
        num.exclude = sum(rep)
        if(num.exclude == 1) df[rep, exclude := temp.exclude]
        if(num.exclude > 1) {
          worst.row = order(rep, df$temp.diff, decreasing=T)[1]
          df[worst.row, exclude := temp.exclude]
        }
        
        return(df$exclude)
      })(copy(.SD))]
      newly.excluded = sum(subj.df$exclude %in% c('Exclude-Min-Weight-Change', 'Exclude-Max-Weight-Change'))
      if(newly.excluded > num.weight.excluded) {
        num.weight.excluded = newly.excluded
      } else {
        break
      }
    }  
    setkey(subj.df,index)
    
    return(subj.df$exclude)
  })(copy(.SD)), by=.(subjid),.SDcols=c('sex', 'agedays', 'v', 'tbc.sd', 'exclude')]
  if(!quietly) cat(sprintf("[%s] End Step 9\n", Sys.time()))
  
  if(!quietly) cat(sprintf("[%s] Completed Batch #%s...\n", Sys.time(), data.df$batch[1]))
  if(!quietly & .parallel) sink()
  return(data.df)
}