---
title: "algorithm2"
output: html_document
---

## R package
```{r message=FALSE,warning=FALSE}
library(MASS)
library(ggplot2)
```

## Preparation for Missing Data
```{r message=FALSE, warning=FALSE}
##initial steps
a=read.csv(file="mr_sites_in_regions.csv",header=T,sep=",")
set.seed(1010)

##remove first two chromosome columns
data1=a[,-c(1:2)]

##Preparation of twins for absolute difference data
k=matrix(0,nrow=nrow(data1),ncol=ncol(data1)/2)
k=as.data.frame(k)

##keep intermediate data
data2=data1

##get twins for absolute difference data
for(i in 1:130){
  k[,i]=abs(data2[,2*i-1]-data2[,2*i])
}
b=k


##Filter the data according to the number of 0 or 1 in each row of the data
k=rep(0,nrow(b))
for(i in 1:nrow(b)){
  k[i]=length(which(b[i,]==1| b[i,]==0))
}
k=as.data.frame(k)


##Select the 15 rows with the least 0 or 1
k$ID=1:nrow(a)
colnames(k)=c("Count","ID")
k=k[order(k$Count,decreasing=F),]
id=head(k$ID,15)
data=b[id,]
dat=a[id,]
dat=dat[,-c(1:2)]

##Select the region where these 15 rows are located
#id_1:2 dim region  id_2:3 dim region  id_new=id_1+id_2
id_new=c(7,8,11,12,17,18,23,24,27,28,45,46,53,54,    95,96,97,   98,99,100,   
         107,108,109,   110,111,112,   122,123,124,   131,132,133    )
id_1=c(7,8,11,12,17,18,23,24,27,28,45,46,53,54) #id_1=c(1:14)
id_2=c(95,96,97    ,98,99,100,   107,108,109,   110,111,112,   122,123,124,   
       131,132,133)
data1=b[id_1,]
data2=b[id_2,]
dat1=a[id_1,]
dat2=a[id_2,]
dat1=dat1[,-c(1:2)]
dat2=dat2[,-c(1:2)]

##2d region data frame NA preparation mainly for the row of id_1
fill_NA_2d=function(x,data,i){
  missing_id1=sample(x=1:ncol(data),floor(ncol(data)*x))
  id2=setdiff(1:ncol(data),missing_id1)
  missing_id2=sample(x=id2,floor(ncol(data)*x))
  c=data[(2*i-1):(2*i),]
  c[1,missing_id1]=NA
  c[2,missing_id2]=NA
  return(c)
}
a=matrix(0,14,130)
a=as.data.frame(a)
m_df_2d=function(x,data){
  for(i in 1:7){
    a[(2*i-1):(2*i),]=fill_NA_2d(x,data,i)
  }
  return(a)
}
dat21.=m_df_2d(0.05,data1)
dat22.=m_df_2d(0.1,data1)
dat23.=m_df_2d(0.15,data1)
dat24.=m_df_2d(0.2,data1)
dat25.=m_df_2d(0.3,data1)

##3d region data frame NA preparation mainly for the row of id_2
fill_NA_3d=function(x,data,i){
  missing_id1=sample(x=1:ncol(data),floor(ncol(data)*x))
  id2=setdiff(1:ncol(data),missing_id1)
  missing_id2=sample(x=id2,floor(ncol(data)*x))
  id3=setdiff(1:ncol(data),c(missing_id1,missing_id2))
  missing_id3=sample(x=id3,floor(ncol(data)*x))
  c=data[(3*i-2):(3*i),]
  c[1,missing_id1]=NA
  c[2,missing_id2]=NA
  c[3,missing_id3]=NA
  return(c)
}

b=matrix(0,18,130)
b=as.data.frame(b)
m_df_3d=function(x,data){
  for(i in 1:6){
    b[(3*i-2):(3*i),]=fill_NA_3d(x,data,i)
  }
  return(b)
}

dat31.=m_df_3d(0.05,data2)
dat32.=m_df_3d(0.1,data2)
dat33.=m_df_3d(0.15,data2)
dat34.=m_df_3d(0.2,data2)
dat35.=m_df_3d(0.3,data2)


##Mainly to prepare the missing data frame required for the beta distribution(260 cols),
##because the FN distribution data involves the subtraction of the absolute value of the 
##twins,restore it to the data frame of the absolute value subtraction and missing at the
##corresponding position
get_NA=function(x){
  rd=sample(1:3,1,prob=c(1,1,1))
  if(rd==1){a=c(NA,NA)
  }else if(rd==2){a=c(x[1,1],NA)
  }else if(rd==3){a=c(NA,x[1,2])}
  return(a)
}

##prepartion for the missing data frame for beta distribution
a=dat1;data=dat21.
for(i in 1:14){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat21=a

p=rep(0,14)
for(j in 1:14){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat1;data=dat22.
for(i in 1:14){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat22=a


p=rep(0,14)
for(j in 1:14){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}


a=dat1;data=dat23.
for(i in 1:14){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat23=a

p=rep(0,14)
for(j in 1:14){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat1;data=dat24.
for(i in 1:14){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat24=a

p=rep(0,14)
for(j in 1:14){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat1;data=dat25.
for(i in 1:14){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat25=a


p=rep(0,14)
for(j in 1:14){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat2;data=dat31.
for(i in 1:18){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat31=a
p=rep(0,18)
for(j in 1:18){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat2;data=dat32.
for(i in 1:18){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat32=a

p=rep(0,18)
for(j in 1:18){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat2;data=dat33.
for(i in 1:18){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat33=a

p=rep(0,18)
for(j in 1:18){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat2;data=dat34.
for(i in 1:18){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat34=a

p=rep(0,18)
for(j in 1:18){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}



a=dat2;data=dat35.
for(i in 1:18){
  k=which(is.na(data[i,]))
  n=length(k)
  for(j in 1:n){
    a[i,(2*k[j]-1):(2*k[j])]=get_NA(a[i,(2*k[j]-1):(2*k[j])])
  }
}
dat35=a

p=rep(0,18)
for(j in 1:18){
  p[j]=length(which(is.na(a[j,])))/length(a[j,])
}


##df1:FN distribtion missing data frame; df2:beta distribution missing data frame
df1=list(dat21,dat22,dat23,dat24,dat25,dat31,dat32,dat33,dat34,dat35)
df2=list(dat21.,dat22.,dat23.,dat24.,dat25.,dat31.,dat32.,dat33.,dat34.,dat35.)

```


## Applying beta to imputate

```{r}
## estimate the mle of beta 
##2d region
mean1=apply(dat1[,],1,mean)
var1=apply(dat1[,],1,var)
par1=data.frame(mean=mean1,var=var1)

##3d region 
mean2=apply(dat2[,],1,mean)
var2=apply(dat2[,],1,var)
par2=data.frame(mean=mean2,var=var2)

##2d mle of beta
alpha1=rep(0,nrow(par1))
beta1=rep(0,nrow(par1))
for(i in 1:nrow(par1)){
  alpha1[i]=(mean1[i]*(1-mean1[i])/(var1[i]^2)-1)*mean1[i]
  beta1[i]=(1-mean1[i])*(mean1[i]*(1-mean1[i])/(var1[i]^2)-1)
}
par1=data.frame(mean=mean1,var=var1,alpha=alpha1,beta=beta1)

##3d mle of beta
alpha2=rep(0,nrow(par2))
beta2=rep(0,nrow(par2))
for(i in 1:nrow(par2)){
  alpha2[i]=(mean2[i]*(1-mean2[i])/(var2[i]^2)-1)*mean2[i]
  beta2[i]=(1-mean2[i])*(mean2[i]*(1-mean2[i])/(var2[i]^2)-1)
}
par2=data.frame(mean=mean2,var=var2,alpha=alpha2,beta=beta2)

f1=function(x){
  n=length(x);
  a=rep(0,n);
  for(i in 1:n){
    if(x[i]%%2==1){a[i]=x[i]+1}else{a[i]=x[i]-1}
  }
  return(a)}

##2d region
mean1=apply(dat1[,],1,mean)
var1=apply(dat1[,],1,var)
par1=data.frame(mean=mean1,var=var1)

##3d region
mean2=apply(dat2[,],1,mean)
var2=apply(dat2[,],1,var)
par2=data.frame(mean=mean2,var=var2)

##2d beta sampling
beta_2_sample_REM=function(data,i,seed){
  b=data[(2*i-1):(2*i),]
  true_value1=dat1[2*i-1,which(is.na(data[2*i-1,]))]
  true_value2=dat1[2*i,which(is.na(data[2*i,]))]
  true_value1_p=dat1[2*i-1,f1(which(is.na(data[2*i-1,])))]
  true_value2_p=dat1[2*i,f1(which(is.na(data[2*i,])))]
  true_value1=abs(true_value1-true_value1_p)
  true_value2=abs(true_value2-true_value2_p)
  est_value1=rep(0,length(which(is.na(data[1,]==T))))
  for(j in 1:length(est_value1)){
    est_value1[j]=rbeta(1,alpha1[2*i-1],beta1[2*i-1])
  }
  est_value2=rep(0,length(which(is.na(data[2,]==T)))   )
  for(j in 1:length(est_value2)){
    est_value2[j]=rbeta(1,alpha1[2*i],beta1[2*i])
  }
  est_value1=abs(est_value1-true_value1_p)
  est_value2=abs(est_value2-true_value2_p)
  REM=sum(abs(true_value1 - est_value1))+ sum(abs(true_value2 - est_value2))
  k=length(which(is.na(data[2*i-1,])))/length(data[2*i-1,])
  a=data.frame("SAD"=REM,"method"="beta","regionDim"=2,"regionNumber"=i,"percentage"=k,"repeatNumber"=seed)
  return(a)
}


##3d beta sampling
beta_3_sample_REM=function(data,i,seed){
  b=data[(3*i-2):(3*i),]
  true_value1=dat2[3*i-2,which(is.na(data[3*i-2,]))]
  true_value2=dat2[3*i-1,which(is.na(data[3*i-1,]))]
  true_value3=dat2[3*i,which(is.na(data[3*i,]))]
  true_value1_p=dat2[3*i-2,f1(which(is.na(data[3*i-2,])))]
  true_value2_p=dat2[3*i-1,f1(which(is.na(data[3*i-1,])))]
  true_value3_p=dat2[3*i,f1(which(is.na(data[3*i,])))]
  true_value1=abs(true_value1-true_value1_p)
  true_value2=abs(true_value2-true_value2_p)
  true_value3=abs(true_value3-true_value3_p)
  est_value1=rep(0,length(which(is.na(data[1,]==T))))
  for(j in 1:length(est_value1)){
    est_value1[j]=rbeta(1,alpha2[3*i-2],beta2[3*i-2])
  }
  est_value2=rep(0,length(which(is.na(data[2,]==T))))
  for(j in 1:length(est_value2)){
    est_value2[j]=rbeta(1,alpha2[3*i-1],beta2[3*i-1])
  }
  est_value3=rep(0,length(which(is.na(data[3,]==T))))
  for(j in 1:length(est_value3)){
    est_value3[j]=rbeta(1,alpha2[3*i],beta2[3*i])
  }
  est_value1=abs(est_value1-true_value1_p)
  est_value2=abs(est_value2-true_value2_p)
  est_value3=abs(est_value3-true_value3_p)
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))+sum(abs(true_value3 - est_value3))
  k=length(which(is.na(data[3*i-2,])))/length(data[3*i-2,])
  a=data.frame("SAD"=REM,"method"="beta","regionDim"=3,"regionNumber"=i+7,"percentage"=k,"repeatNumber"=seed)
  return(a)
}


B2=matrix(0,700,6)
B2=as.data.frame(B2)
for(k in 1:5){
  for(j in 1:20){
    for(i in 1:7){
      B2[140*(k-1)+20*(i-1)+j,]=beta_2_sample_REM(df1[[k]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(B2)=colnames(b)


B3=matrix(0,600,6)
B3=as.data.frame(B3)
for(k in 1:5){
  for(j in 1:20){
    for(i in 1:6){
      B3[120*(k-1)+20*(i-1)+j,]=beta_3_sample_REM(df1[[k+5]],i,j)}
  }
}
b=beta_3_sample_REM(dat31,2,1)
colnames(B3)=colnames(b)
```



## Applying rnorm to imputate
```{r message=FALSE, warning=FALSE}
##estimate mle of rnorm(1d FN)
##numerical method nlm to get mle of FN distribution
f=function(rr,m,s) dnorm(rr,mean=m,sd=s)+dnorm(-rr,mean=m,sd=s)
minusl=function(params,r){
  mu=params[1]
  sigma=params[2]
  -sum(log(f(r,mu,sigma)))
}

##2d region
a1=matrix(0,14,2)
colnames(a1)=c("mean","sd")
for(i in 1:14){
  start=c(mean(as.numeric(data1[i,])),sd(as.numeric(data1[i,])))
  nlm(minusl,start, r=as.numeric(data1[i,]))
  a1[i,]=nlm(minusl,start, r=as.numeric(data1[i,]))$estimate
}

##3d region
a2=matrix(0,18,2)
colnames(a2)=c("mean","sd")
for(i in 1:18){
  start=c(mean(as.numeric(data2[i,])),sd(as.numeric(data2[,i])))
  nlm(minusl,start, r=as.numeric(data2[i,]))
  a2[i,]=nlm(minusl,start, r=as.numeric(data2[i,]))$estimate
}
a1=as.data.frame(a1)
a2=as.data.frame(a2)
write.csv(a1,"1d FN_2region.csv")
write.csv(a2,"1d FN_3region.csv")


##2d region rnorm sampling
FN_2_sample_REM=function(data,i,seed){
  b=data[(2*i-1):(2*i),]
  true_value1=data1[2*i-1,which(is.na(data[2*i-1,]))]
  true_value2=data1[2*i,which(is.na(data[2*i,]))]
  est_value1=abs(rnorm(length(true_value1),mean=a1$mean[2*i-1],sd=a1$sd[2*i-1]))
  est_value2=abs(rnorm(length(true_value2),mean=a1$mean[2*i],sd=a1$sd[2*i]))
  REM=sum(abs(true_value1 - est_value1))+  sum(abs(true_value2 - est_value2))
  k=mean(length(which(is.na(data[2*i-1,])))/length(data[2*i-1,]),length(which(is.na(data[2*i,])))/length(data[2*i,]))
  a=data.frame("SAD"=REM,"method"="rnorm","regionDim"=2,"regionNumber"=i,"percentage"=k,"repeatNumber"=seed)
  return(a)
}



##3d region rnorm sampling 
FN_3_sample_REM=function(data,i,seed){
  b=data[(3*i-2):(3*i),]
  true_value1=data2[3*i-2,which(is.na(data[3*i-2,]))]
  true_value2=data2[3*i-1,which(is.na(data[3*i-1,]))]
  true_value3=data2[3*i,which(is.na(data[3*i,]))]
  est_value1=abs(rnorm(length(true_value1),mean=a2$mean[3*i-2],sd=a2$sd[3*i-2]))
  est_value2=abs(rnorm(length(true_value2),mean=a2$mean[3*i-1],sd=a2$sd[3*i-1]))
  est_value3=abs(rnorm(length(true_value3),mean=a2$mean[3*i],sd=a2$sd[3*i]))
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))+sum(abs(true_value3 - est_value3))
  k=mean(length(which(is.na(data[3*i-2,])))/length(data[3*i-2,]),length(which(is.na(data[3*i-1,])))/length(data[3*i-1,]),length(which(is.na(data[3*i,])))/length(data[3*i,]))
  a=data.frame("SAD"=REM,"method"="rnorm","regionDim"=3,"regionNumber"=i+7,"percentage"=k,"repeatNumber"=seed)
  return(a)
}


R2=matrix(0,700,6)
R2=as.data.frame(R2)
for(k in 1:5){
  for(i in 1:7){
    for(j in 1:20){
      R2[140*(k-1)+20*(i-1)+j,]=FN_2_sample_REM(df2[[k]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(R2)=colnames(b)

R3=matrix(0,600,6)
R3=as.data.frame(R3)
for(k in 1:5){
  for(i in 1:6){
    for(j in 1:20){
      R3[120*(k-1)+20*(i-1)+j,]=FN_3_sample_REM(df2[[k+5]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(R3)=colnames(b)
```

## Applying mvtnorm to imputate
```{r message=FALSE,warning=FALSE}
##estimate mle of multi-dimensional FN distribution
##mainly by numerical method 
a1=read.csv(file="a1_par.csv")
a2=read.csv(file="a2_par.csv")
a1=a1[,-1]
a2=a2[,-1]
colnames(a1)=c("u1","u2","s1","s2","rho")
colnames(a2)=c("u1","u2","u3","delta11","delta12","delta13","delta22","delta23","delta33")


##2d mvtnorm sampling
mvrnorm_2_sample_order=function(data,i,seed){
  b=data[(2*i-1):(2*i),]
  true_value12=data1[(2*i-1):(2*i),which(is.na(data[2*i-1,]))]
  true_value12_order=true_value12[,order(true_value12[2,])]
  true_value1=true_value12_order[1,]
  true_value21=data1[(2*i-1):(2*i),which(is.na(data[2*i,]))]
  true_value21_order=true_value21[,order(true_value21[1,])]
  true_value2=true_value21_order[2,]
  u=c(a1[i,1],a1[i,2]);sigma=matrix(c(a1[i,3]^2,a1[i,5]*(a1[i,3]*a1[i,4]),a1[i,5]*
                                        (a1[i,3]*a1[i,4]),a1[i,4]^2),2,2)
  est_value=abs(mvrnorm(length(true_value1),u,sigma))
  est_value1=est_value[1,]
  est_value1=sort(est_value1)
  est_value2=est_value[2,]
  est_value2=sort(est_value2)
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))
  k=mean(length(which(is.na(data[2*i-1,])))/length(data[2*i-1,]),length(which(is.na(data[2*i,])))/length(data[2*i,]))
  a=data.frame("SAD"=REM,"method"="mvtnorm","regionDim"=2,"regionNumber"=i,"percentage"=k,"repeatNumber"=seed)
  return(a)
}



##3d mvtnorm sampling
mvrnorm_3_sample_order=function(data,i,seed){
  b=data[(3*i-2):(3*i),]
  true_value12=data2[(3*i-2):(3*i),which(is.na(data[3*i-2,]))]
  true_value12_order=true_value12[,order(true_value12[2,]+true_value12[3,])]
  true_value1=true_value12_order[1,]
  true_value21=data2[(3*i-2):(3*i),which(is.na(data[3*i-1,]))]
  true_value21_order=true_value21[,order(true_value21[1,]+true_value21[3,])]
  true_value2=true_value21_order[2,]
  true_value13=data2[(3*i-2):(3*i),which(is.na(data[3*i,]))]
  true_value13_order=true_value13[,order(true_value13[1,]+true_value13[2,])]
  true_value3=true_value13_order[3,]
  u=c(a2[i,1],a2[i,2],a2[i,3]);sigma=matrix(c(a2[i,4],a2[i,5],a2[i,6],a2[i,5],a2[i,7],a2[i,8],a2[i,6],a2[i,8],a2[i,9]),3,3,byrow=T)
  est_value=abs(mvrnorm(length(true_value1),u,sigma))
  est_value=t(est_value)
  est_value1=est_value[1,]
  est_value1=sort(est_value1)
  est_value2=est_value[2,]
  est_value2=sort(est_value2)
  est_value3=est_value[3,]
  est_value3=sort(est_value3)
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))+sum(abs(true_value3 - est_value3))
  k=mean(length(which(is.na(data[3*i-2,])))/length(data[3*i-2,]),length(which(is.na(data[3*i-1,])))/length(data[3*i-1,]),length(which(is.na(data[3*i,])))/length(data[3*i,]))
  a=data.frame("SAD"=REM,"method"="mvtnorm","regionDim"=3,"regionNumber"=i+7,"percentage"=k,"repeatNumber"=seed)
  return(a)
}


T2=matrix(0,700,6)
T2=as.data.frame(T2)
for(k in 1:5){
  for(i in 1:7){
    for(j in 1:20){
      T2[140*(k-1)+20*(i-1)+j,]=mvrnorm_2_sample_order(df2[[k]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(T2)=colnames(b)


T3=matrix(0,600,6)
T3=as.data.frame(T3)
for(k in 1:5){
  for(i in 1:6){
    for(j in 1:20){
      T3[120*(k-1)+20*(i-1)+j,]=mvrnorm_3_sample_order(df2[[k+5]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(T3)=colnames(b)

```

## Applying Gibbs to imputate
```{r message=FALSE,warning=FALSE}
##2d Gibbs sampling
###initial parameter
g2=function(mu1,mu2,sigma1,sigma2,rho){
  ##Gibbs 2 dim FN sampling with mean not equal to 0 
  N <- 10000 #length of chain
  burn<- 2500 #burn-in length
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  ###### generate the chain #####
  X[1, ] <- c(mu1, mu2) #initialize
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <-abs( rnorm(1, m1, s1))
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <-abs(rnorm(1, m2, s2))
  }
  b <- burn + 1
  sample <- X[b:N, ]
}

Gibbs_2_sample_order=function(data,i,seed){
  b=data[(2*i-1):(2*i),]
  true_value12=data1[(2*i-1):(2*i),which(is.na(data[2*i-1,]))]
  true_value12_order=true_value12[,order(true_value12[2,])]
  true_value1=true_value12_order[1,]
  true_value21=data1[(2*i-1):(2*i),which(is.na(data[2*i,]))]
  true_value21_order=true_value21[,order(true_value21[1,])]
  true_value2=true_value21_order[2,]
  sample=g2(mu1=a1$u1[i],mu2=a1$u2[i],sigma1=a1$s1[i],sigma2=a1$s2[i],rho=a1$rho[i])
  est_value=tail(sample,length(true_value1))
  est_value=t(est_value)
  est_value1=est_value[1,]
  est_value1=sort(est_value1)
  est_value2=est_value[2,]
  est_value2=sort(est_value2)
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))
  k=mean(length(which(is.na(data[2*i-1,])))/length(data[2*i-1,]),length(which(is.na(data[2*i,])))/length(data[2*i,]))
  a=data.frame("SAD"=REM,"method"="Gibbs","regionDim"=2,"regionNumber"=i,"percentage"=k,"repeatNumber"=seed)
  return(a)
}

##3d Gibbs sampling
g3=function(u,sig){
  ##3 dim Gibbs of FN distribution
  N <- 10000 #length of chain
  burn<- 2500 #burn-in length
  X <- matrix(0, N, 3) #the chain, a bivariate sample
  s1=sqrt(sig[1,1]-sig[1,-1]%*%solve(sig[-1,-1])%*%sig[-1,1])
  s2=sqrt(sig[2,2]-sig[2,-2]%*%solve(sig[-2,-2])%*%sig[-2,2])
  s3=sqrt(sig[3,3]-sig[3,-3]%*%solve(sig[-3,-3])%*%sig[-3,3])
  X[1, ] <-u #initialize
  for (i in 2:N) {
    x1=X[i-1,-1]
    m1=u[1]+sig[1,-1]%*%solve(sig[-1,-1])%*%(x1-u[-1])
    X[i,1]=abs(rnorm(1,m1,s1))
    
    x2=c(X[i,1],X[i-1,3])
    m2 <-u[2]+sig[2,-2]%*%solve(sig[-2,-2])%*%(x2-u[-2])
    X[i,2]=abs(rnorm(1,m2,s2))
    
    x3=X[i,-3]
    m3  =u[3]+sig[3,-3]%*%solve(sig[-3,-3])%*%(x3-u[-3])
    X[i,3]=abs(rnorm(1,m3,s3))
    
  }
  b <- burn + 1
  sample <- X[b:N, ]
}

Gibbs_3_sample_order=function(data,i,seed){
  b=data[(3*i-2):(3*i),]
  true_value12=data2[(3*i-2):(3*i),which(is.na(data[3*i-2,]))]
  true_value12_order=true_value12[,order(true_value12[2,]+true_value12[3,])]
  true_value1=true_value12_order[1,]
  true_value21=data2[(3*i-2):(3*i),which(is.na(data[3*i-1,]))]
  true_value21_order=true_value21[,order(true_value21[1,]+true_value21[3,])]
  true_value2=true_value21_order[2,]
  true_value13=data2[(3*i-2):(3*i),which(is.na(data[3*i,]))]
  true_value13_order=true_value13[,order(true_value13[1,]+true_value13[2,])]
  true_value3=true_value13_order[3,]
  sample=g3(u=c(a2[i,1],a2[i,2],a2[i,3]),sig=matrix(c(a2[i,4],a2[i,5],a2[i,6],a2[i,5],a2[i,7],a2[i,8],a2[i,6],a2[i,8],a2[i,9]),3,3,byrow=T))
  est_value=tail(sample,length(true_value1))
  est_value=t(est_value)
  est_value1=est_value[1,]
  est_value1=sort(est_value1)
  est_value2=est_value[2,]
  est_value2=sort(est_value2)
  est_value3=est_value[3,]
  est_value3=sort(est_value3)
  REM=sum(abs(true_value1 - est_value1))+sum(abs(true_value2 - est_value2))+sum(abs(true_value3 - est_value3))
  k=mean(length(which(is.na(data[3*i-2,])))/length(data[3*i-2,]),length(which(is.na(data[3*i-1,])))/length(data[3*i-1,]),length(which(is.na(data[3*i,])))/length(data[3*i,]))
  a=data.frame("SAD"=REM,"method"="Gibbs","regionDim"=3,"regionNumber"=i+7,"percentage"=k,"repeatNumber"=seed)
  return(a)
}


G2=matrix(0,700,6)
G2=as.data.frame(G2)
for(k in 1:5){
  for(i in 1:7){
    for(j in 1:20){
      G2[140*(k-1)+20*(i-1)+j,]=Gibbs_2_sample_order(df2[[k]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(G2)=colnames(b)

G3=matrix(0,600,6)
G3=as.data.frame(G3)
for(k in 1:5){
  for(i in 1:6){
    for(j in 1:20){
      G3[120*(k-1)+20*(i-1)+j,]=Gibbs_3_sample_order(df2[[k+5]],i,j)}
  }
}
b=beta_2_sample_REM(dat21,2,1)
colnames(G3)=colnames(b)
```


## Visualization
```{r message=FALSE,warning=FALSE}
B2$pc.=round(R2$percent,2)
B3$pc.=round(R3$percent,2)
R2$pc.=round(R2$percent,2)
R3$pc.=round(R3$percent,2)
T2$pc.=round(R2$percent,2)
T3$pc.=round(R3$percent,2)
G2$pc.=round(R2$percent,2)
G3$pc.=round(R3$percent,2)

dat=rbind(B2,B3,R2,R3,T2,T3,G2,G3)
dat$pc.=paste(dat$pc.*100,"%",sep="")
library(ggplot2)
dat$pc.=as.factor(dat$pc.)
dat$pc.=ordered(dat$pc.,levels=c("5%","10%","15%","20%","30%"))

dat$method=as.factor(dat$method)
dat$method<- factor(dat$method, levels = c("Gibbs", "mvtnorm", "beta","rnorm"))
dat$pc.=as.factor(dat$pc.)


##fig1
p1=ggplot(dat,aes(x=pc.,y=SAD))+geom_boxplot(aes(fill=method),outlier.colour=NA)+labs(x="The percentage of missing data",y="SAD")+theme(axis.title.x =element_text(size=20,vjust=-0.25), axis.title.y=element_text(size=20))+ theme(axis.text.x = element_text(size = 15,vjust=1))+theme(axis.text.y = element_text(size = 15))+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ theme(
  legend.key.size = unit(20, "pt")
)+theme(legend.text=element_text(size=12),legend.title = element_text(size=14))
ggsave("REMtable2.png",height=6,width=10,dpi=300,plot=p1)

```

