import numpy as np
from numpy.linalg import inv
from sympy import *
import math
import matplotlib.pyplot as plt
N=6;gamma=0.5;alpha=0.5;theta_cho=0

arr = np.loadtxt('/Users/islet/Desktop/虚拟网络/DATA/connections.txt',delimiter=',')
connect=arr.astype(np.int64)
vertex=[2,5];vertexResource=[5,10];bandwith=[1,5]
calres=[100]*N#calculate recource

def degree(con):
    num=[];bor=[]
    for i in list(range(N)):
        numof=0;borof=[0]*N
        xline = con[i,:]
        yline = con[:,i]
        for j in list(range(N)):
            if xline[j] != 0:
                numof+=1
                borof[j] = 1
            if yline[j] != 0:
                numof+=1
                borof[j] = 1
        num.append(numof)
        bor.append(borof)

    return num,bor
def clustering(num,bor):
    num=np.array(num)
    cluster=[0]*N
    for i in list(range(N)):
        borof=np.array(bor[i])
        deg_bor=sum(num*borof)/num[i]
        cluster[i]=deg_bor
    return cluster
def resource(con,calres,bor):
    res=np.zeros(N)
    RE_res=np.zeros(N)
    P=np.zeros((N,N))

    for i in list(range(N)):
        resof=((sum(con[i,:])+sum(con[:,i]))*calres[i])/2
        res[i]=resof
    resout=res
    res=np.array(res)

    for i in list(range(N)):
        borof=np.array(bor[i])
        RE_resof=res[i]/sum(res*borof)
        RE_res[i]=RE_resof

    bor=np.array(bor)
    for i in list(range(N)):
        Pof=(RE_res*bor[i])/sum(RE_res*bor[i])
        P[i]=Pof

    return resout,RE_res,P
def Markov(RE_res,P):
    I=np.zeros((N,N))
    V=(1-gamma)*np.dot(RE_res,inv(I-gamma*P))

    return V
def normalize(num):
    ans=(num-min(num)*np.ones(N))/(max(num)-min(num))
    return ans
def Normal_Ogive(RE_a,RE_b,theta):
    Q=np.zeros(N)
    P=np.zeros(N)
    a_fin=sum(RE_a)/N
    b_fin=sum(RE_b)/N
    theta = 4*(theta-min(theta)*np.ones(N))/(max(theta)-min(theta))-2
    # theta = 4*theta/(max(theta)-min(theta))
    for i in list(range(N)):
        x = symbols('x')
        expr1 = exp(-x**2/2)/sqrt(2*math.pi)
        P[i] = integrate(expr1, (x,a_fin,a_fin*(theta[i]-b_fin)))
        if theta[i]>theta_cho:
            Q[i]=P[i]
        else:Q[i]=1-P[i]

    return Q,P

[num,bor]=degree(connect)
a=clustering(num,bor)
[res,RE_res,P]=resource(connect,calres,bor)
b=Markov(RE_res,P)
RE_a = normalize(a)
RE_b = normalize(b)
theta = RE_a + alpha*RE_b
print(theta)
[Q,P] = Normal_Ogive(RE_a,RE_b,theta)
fig = plt.figure(figsize=(12, 8))
x=theta;y=Q
plt.scatter(x,y,color='r', marker='x')
plt.show()
print(Q)
# print(P)
