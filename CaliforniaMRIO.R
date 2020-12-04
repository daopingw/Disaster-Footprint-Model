
rm(list = ls())
cat("\014")

library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gdxrrw)
igdx(gamsSysDir = "D:/GAMS/win64/24.8")

library(data.table)
library(readxl)

library(Rcpp)
sourceCpp("matrix_multiplication.cpp")


R = 59
N = 80
U = 4
Set = as.matrix(read_excel(path = "Data/IOT for MRIO_test.xlsx", sheet = "SetCalifornia", col_names = T))
Region = list(c(Set[1:R,1]))
Sector = list(c(Set[1:N,2]))
Factor = list(c(Set[1:U,3]))
Regionlst = list(name = "r", type = "set", uels = Region, ts = "Region")
Sectorlst = list(name = "i", type = "set", uels = Sector, ts = "Sector")
Factorlst = list(name = "u", type = "set", uels = Factor, ts = "Factor")

# IOT data
#IOT = as.matrix(fread(input = "Data/IO/MRIO_final.csv",header = F))
load("MRIO_f.RData")
IOT = MRIO
IOT = IOT/365
#IOT[1:(N*R),N*R+R] = IOT[1:(N*R),N*R+R] + apply(IOT[,1:(N*R)],2,sum) - apply(IOT[1:(N*R),],1,sum)

a = sum((apply(IOT[,1:(N*R)],2,sum) - apply(IOT[1:(N*R),],1,sum))*(apply(IOT[,1:(N*R)],2,sum) - apply(IOT[1:(N*R),],1,sum)))

Z_0_raw = IOT[1:(R*N),1:(R*N)]
Z_0 = array(0,dim = c(R,N,R,N))
for(i in 1:R)
{
    for(j in 1:N)
    {
        Z_0[i,j,,] = matrix(Z_0_raw[N*(i-1)+j,],nrow = R,byrow = T)
    }
}
Z_0lst = list(name = "Z_0", type = "parameter", form = "full", ts = "Z0", val = Z_0, uels = c(Region,Sector,Region,Sector))

V_0_raw = IOT[(R*N+1):(R*N+U),1:(R*N)]
V_0 = array(0,dim = c(U,R,N))
for(i in 1:R)
{
    V_0[,i,] = V_0_raw[,(N*(i-1)+1):(N*i)]
}
V_0lst = list(name = "V_0", type = "parameter", form = "full", ts = "V0", val = V_0, uels = c(Factor,Region,Sector))

F_0_raw = IOT[1:(R*N),(R*N+1):(R*N+R)]
F_0 = array(0,dim = c(R,N,R))
for(i in 1:R)
{
    F_0[i,,] = F_0_raw[(N*(i-1)+1):(N*i),]
}
F_0lst = list(name = "F_0", type = "parameter", form = "full", ts = "F0", val = F_0, uels = c(Region,Sector,Region))

X_0_raw = apply(IOT[,1:(R*N)],2,sum)

A = eigenMapMatMult(Z_0_raw,diag(1/X_0_raw))
temp = X_0_raw - (eigenMapMatMult(A,X_0_raw) + apply(F_0_raw,1,sum))
for(i in 1:(N*R))
{
    if(temp[i]<0)
    {
        a = which(F_0_raw[i,]>abs(temp[i]))[1]
        if(is.na(a))
        {print("bad");print(i)}
        else
            F_0_raw[i,a] = F_0_raw[i,a] + temp[i]
    }
}


A_0 = array(0,dim = c(R,N,R,N))
for(i in 1:R)
{
    for(j in 1:N)
    {
        A_0[i,j,,] = matrix(A[N*(i-1)+j,],nrow = R,byrow = T)
    }
}
A_0lst = list(name = "A", type = "parameter", form = "full", ts = "A", val = A_0, uels = c(Region,Sector,Region,Sector))

X_0 = matrix(0,nrow = R,ncol = N)
for(i in 1:R)
{
    X_0[i,] = X_0_raw[(N*(i-1)+1):(N*i)]
}
X_0lst = list(name = "X_0", type = "parameter", form = "full", ts = "X_0", val = X_0, uels = c(Region,Sector))

F_02 = matrix(0,nrow = R,ncol = N)
F_02_raw = apply(F_0_raw,1,sum)
for(i in 1:R)
{
    F_02[i,] = F_02_raw[(N*(i-1)+1):(N*i)]
}
F_02lst = list(name = "F_0", type = "parameter", form = "full", ts = "F_0", val = F_02, uels = c(Region,Sector))


F_up_raw = apply(F_0_raw,1,sum)
X_up_raw = X_0_raw

F_up = matrix(0,nrow = R,ncol = N)
X_up = matrix(0,nrow = R,ncol = N)


# Date
Ts = as.POSIXct(strptime("2018-06-01", "%Y-%m-%d"), tz = "UTC") #Staring date
Te = as.POSIXct(strptime("2020-12-31", "%Y-%m-%d"), tz = "UTC") #End date
Td = as.numeric(Te - Ts) + 1
K = Td

# Key Simulation Results
X_mat = matrix(0,nrow = K,ncol = N*R)
X_mat[1,] = X_0_raw


IOT_0 = matrix(0,R+1,R+1)
for(i in 1:R)
{
    for (j in 1:R) {
        IOT_0[i,j] = sum(IOT[(N*(i-1)+1):(N*i),(N*(j-1)+1):(N*j)])
    }
    IOT_0[R+1,i] = sum(IOT[,(N*(i-1)+1):(N*i)])
    IOT_0[i,R+1] = sum(IOT[(N*(i-1)+1):(N*i),(N*R+1):(N*R+R)])
}







######################

FinalMatrix = matrix(0,nrow = R*N,ncol = R+1)
FinalMatrix[,1:R] = F_0_raw

CapacityMatrix = matrix(1,nrow = 3,ncol = N*R)


# Transpotation Constraint
#####
data = as.data.frame(read_excel(path = "Direct/Wildfires.xlsx",sheet = "Domain",col_names = TRUE))
Sn = 80
Cn = 58

result_trans_c = array(0,dim = c(Sn,Td,Cn)) # dimnames = (Sector,Period,County)

for(i in 1:Td)
{
    for(j in 1:length(data[,1]))
    {
        Tc = Ts + (i-1)*86400
        if((Tc>=data$StartDate[j])&(Tc<=data$EndingDate[j]))
        {
            Cc = data$CountyCode[j]
            result_trans_c[,i,Cc] = result_trans_c[,i,Cc] + (data$Size[j]/data$Area[j])/(as.numeric(data$EndingDate[j]-data$StartDate[j]) +1)
        }
    }
}

lambda_trans = 7  # recovery time

result_trans = array(0,dim = c(Sn,Td,Cn))
for(i in 1:Td)
{
    for(j in 1:(min(lambda_trans,i)))
    {
        result_trans[,i,] = result_trans[,i,] + ((lambda_trans-j+1)/lambda_trans)*result_trans_c[,i-j+1,]
    }
}
#####


# Capital
#####
Capital_mat = matrix(0,nrow = K,ncol = R*N)
Capital_0 = 365*5*apply(V_0_raw[1:2,],2,sum)
Capital_mat[1,] = Capital_0

Capital_0_h = rep(10000000,R)

Capital_mat_recovery = matrix(0,nrow = K,ncol = R*N)
Recovery = matrix(0,nrow = R*N,ncol = R*N)
Recovery_h = matrix(0,nrow = R*N,ncol = R)
Capital_mat_demand = matrix(0,nrow = R*N,ncol = R*N)
Capital_mat_demand_h = matrix(0,nrow = R*N,ncol = R)


Damage_cap = matrix(0,nrow = K,ncol = N*R)
Damage_cap_h = matrix(0,nrow = K,ncol = R)


data = as.data.frame(read_excel(path = "Direct/Wildfires.xlsx",sheet = "Domain Capital",col_names = TRUE))
data[is.na(data)] = 0

Sn = 80
Cn = 58

#确定性数据
p1 = 3.9/16344*1000  #Camp fire mean price (lo)
p2 = 4.66/6700*1000   #Woolsey fire mean price (up)
p3 = (p1+p2)/2   #mean price


#####居民建筑损失
#####

#不确定性
alpha = 0.5    #(0.05---0.95)  #损毁程度的不确定性

result_capital_h = matrix(0,nrow = Td,ncol = Cn)
for(i in 1:Td)
{
    for(j in 1:length(data[,1]))
    {
        Price = p3   #(p1---p2)  #其他火灾区域房屋价格的不确定性
        Tc = Ts + (i-1)*86400
        if((Tc>=data$StartDate[j])&(Tc<=data$EndingDate[j]))
        {
            if(data$IncidentNumber[j] == "CA-BTU-016737")
            {
                Price = p1
            }
            if(data$IncidentNumber[j] == "CA-VNC-091023")
            {
                Price = p2
            }
            Cc = data$CountyCode[j]
            result_capital_h[i,Cc] = result_capital_h[i,Cc] + ((data$ResidencesDestroyed[j]+alpha*data$ResidencesDamaged[j])*Price)/(as.numeric(data$EndingDate[j]-data$StartDate[j]) +1)
        }
    }
}
sum(result_capital_h)

Damage_cap_h[,1:(R-1)] = result_capital_h/10000000


#####生产性资本损失
#####
p4 = (16.5*1000 -(13972+462*0.5)*p1)/(528+102*0.5+4293) 
p5 = (5.2*1000 -(1115+280*0.5)*p2)/(42+62*0.5+343)
p6 = (p4+p5)/2

VAK = matrix(V_0_raw[2,],nrow = N)

I = matrix(1,nrow = 1,ncol = N)

VAK = VAK%*%diag(as.vector((1/(I%*%VAK))))

#不确定性
beta = 0.5   #(0.2---1)#损毁程度不确定性
otherp = 1   

result_capital = array(0,dim = c(Sn,Td,Cn))   # dimnames = (Sector,Period,County)
for(i in 1:Td)
{
    for(j in 1:length(data[,1]))
    {
        Price = p6  #(p4---p5)
        Tc = Ts + (i-1)*86400
        if((Tc>=data$StartDate[j])&(Tc<=data$EndingDate[j]))
        {
            if(data$IncidentNumber[j] == "CA-BTU-016737")
            {
                Price = p4
            }
            if(data$IncidentNumber[j] == "CA-VNC-091023")
            {
                Price = p5
            }
            Cc = data$CountyCode[j]
            result_capital[,i,Cc] = result_capital[,i,Cc] + (((data$CommercialDestroyed[j]+beta*data$CommercialDamaged[j])*Price + (data$OtherDestroyed[j]+0.5*data$OtherDamaged[j])*Price*otherp)/(as.numeric(data$EndingDate[j]-data$StartDate[j]) +1))*VAK[,Cc]  
        }
    }
}
sum(result_capital)


for(i in 1:K)
{
    Damage_cap[i,1:(N*(R-1))] = result_capital[,i,]
    Damage_cap[i,] = Damage_cap[i,]/Capital_0
}

T_lag = 60

#####
Distribute_calc = function(R,N,n,proportion,Final)
{
    temp = Final[seq(n,R*(N-1)+n,N),]
    temp = temp%*%diag(proportion/apply(temp,2,sum))
    result = matrix(0,R*N,R)
    result[seq(n,R*(N-1)+n,N),] = temp
    return(result)
}


Distribute_mat = matrix(0,nrow = R*N,ncol = R*N)
Dis1 = Distribute_calc(R,N,7,0.95,F_0_raw)
Dis2 = Distribute_calc(R,N,9,0.05,F_0_raw)
Dis = Dis1 + Dis2
Distribute_mat = Dis[,sort(rep(1:R,N))]

Dis1h = Distribute_calc(R,N,8,0.95,F_0_raw)
Dis2h = Distribute_calc(R,N,10,0.05,F_0_raw)
Dish = Dis1h + Dis2h
Distribute_mat_h = Dish


#####
Ti = as.POSIXct(strptime("2018-12-31", "%Y-%m-%d"), tz = "UTC")


data = as.data.frame(read_excel(path = "Direct/Health.xlsx",
                                sheet = "Sheet4",        #不确定性Sheet2--损失大,shee3--损失小
                                col_names = TRUE))


a = apply(matrix(V_0_raw[1,1:(N*(R-1))],nrow = N),2,sum)
b = data[,2]/1000000/365
result_labor = matrix(t(matrix(rep(b/a,N),nrow = R-1)),nrow = 1)



a = data[,2]
b = t(matrix(rep(a,80),58))


temp = MRIO[4721,1:4720]
temp = temp[1:4640]
temp = matrix(temp,80)

temp = temp%*%diag(1/apply(temp,2,sum))
L_mat_a = b*temp

write.table(x = L_mat_a,file = "L_mat_aa.csv",sep = ",",row.names = F,col.names = F)




##############################

# Simulation
for(k in 2:365) #最大K-1
{
    print(k)
    
    print(Sys.time())
    
        
    Capital_mat_recovery[k,] = apply(Recovery,2,sum)
    Capital_mat_demand = Capital_mat_demand + eigenMapMatMult(Distribute_mat,diag(Damage_cap[k,]*Capital_0)) - Recovery
    if(k>T_lag)
    {
        Capital_mat[k,] = Capital_mat[k-1,] - Damage_cap[k,]*Capital_0 + Capital_mat_recovery[k-T_lag,]
    }    else    {
        Capital_mat[k,] = Capital_mat[k-1,] - Damage_cap[k,]*Capital_0
    }
    Capital_mat[k,] = apply(Capital_mat[k:(k+1),],2,max)
    Capital_mat_demand_h = Capital_mat_demand_h + eigenMapMatMult(Distribute_mat_h,diag(Damage_cap_h[k,]*Capital_0_h)) - Recovery_h


    CapacityMatrix[2,] = Capital_mat[k,]/Capital_0
    
    CapacityMatrix[3,1:(N*(R-1))] = 1 - matrix(result_trans[,k,],nrow = 1)
    
    CapacityMatrix[1,1:(N*(R-1))] = 1 - result_labor
    
    FinalMatrix[,(R+1)] = apply(Capital_mat_demand,1,sum) + apply(Capital_mat_demand_h,1,sum)

    
    F_up_raw = apply(FinalMatrix,1,sum)
    
    temp = apply(CapacityMatrix,2,min)
    temp[temp<0.4] = 0.4 + runif(1,0,0.01)
    X_up_raw = temp*X_0_raw
    
    t1t = F_up_raw/F_02_raw
    min(t1t)
    max(t1t)
    min(t1t[4641:4720])
    max(t1t[4641:4720])
    t2t = X_up_raw/X_0_raw
    min(t2t)
    if(max(t2t)>1)
    {
        print("maxt2t")
        break
    }
    if(min(t2t[4641:4720])!=1)
    {
        print("mint2tw")
        break
    }
    if(max(t2t[4641:4720])!=1)
    {
        print("maxt2tw")
        break
    }
    
    
    # Upper list
    for(i in 1:R)
    {
        F_up[i,] = F_up_raw[(N*(i-1)+1):(N*i)]
    }
    F_uplst = list(name = "F_up", type = "parameter", form = "full", ts = "F_up", val = F_up, uels = c(Region,Sector))
    
    #X_up_raw = X_0_raw
    #X_up_raw[241:320] = 0.9*X_up_raw[241:320]
    for(i in 1:R)
    {
        X_up[i,] = X_up_raw[(N*(i-1)+1):(N*i)]
    }
    X_uplst = list(name = "X_up", type = "parameter", form = "full", ts = "X_up", val = X_up, uels = c(Region,Sector))
    
    # Data Transfer    
    wgdx("Data/Database.GDX", Regionlst,Sectorlst,Factorlst,A_0lst,X_0lst,F_02lst,X_uplst,F_uplst)
    
    # Solve
    a = gams("Program/CaliforniaMRIO02.gms")
    if(a!=0){break}
    
    # Result Transfer
    X_A = matrix(t(rgdx("Result/Result.GDX",list(name = "X_0", uels = c(Region,Sector), form = "full"))$val),nrow = 1)
    X_mat[k,] = X_A
    print(sum(X_mat[k,]))
    print(sum(X_0_raw[1:4640])*k-sum(X_mat[1:k,1:4640]))
    print(sum(X_0_raw)*k-sum(X_mat[1:k,]))

    if(0)
    {
    IOTn = IOT
    IOTn[1:(N*R),1:(N*R)] = eigenMapMatMult(A,diag(X_mat[k,]))
    IOTn[1:(N*R),(N*R+1)] = t(X_A)-eigenMapMatMult(A,t(X_A))
    IOT_1 = matrix(0,R+1,R+1)
    for(i in 1:R)
    {
        for (j in 1:R) {
            IOT_1[i,j] = sum(IOTn[(N*(i-1)+1):(N*i),(N*(j-1)+1):(N*j)])
        }
        IOT_1[R+1,i] = sum(X_A[(N*(i-1)+1):(N*i)])
        IOT_1[i,R+1] = sum(IOTn[(N*(i-1)+1):(N*i),(N*R+1)])
    }
    
    a = IOT_1/IOT_0-1
    
    write.csv(a,"temp_4f.csv")
    }
    
    ratio = eigenMapMatMult((diag(1/apply(FinalMatrix,1,sum))),FinalMatrix)
    FinalMatrixA = eigenMapMatMult(diag(as.vector((t(X_A)-eigenMapMatMult(A,t(X_A))))),ratio)
    
    RecoveryDemand = cbind(Capital_mat_demand,Capital_mat_demand_h)
    RecoveryDemand[which(RecoveryDemand<0)] = 0
    temp = diag(1/apply(RecoveryDemand,1,sum))
    temp[temp == Inf] = 0
    temp[is.na(temp)] = 0
    temp[is.nan(temp)] = 0
    ratio = eigenMapMatMult(temp,RecoveryDemand)
    RecoveryDemandA = eigenMapMatMult(diag(as.vector(FinalMatrixA[,R+1])),ratio)
    
    Recovery = RecoveryDemandA[,1:(N*R)]
    Recovery_h = RecoveryDemandA[,(N*R+1):(N*R+R)]
    print(Sys.time())
}










