import numpy as np
from numpy.linalg import inv
from sympy import *
import math
import csv
import matplotlib.pyplot as plt


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
    count=0;perc=0;
    a_fin=sum(RE_a)/N
    b_fin=sum(RE_b)/N
    theta = 4*(theta-min(theta)*np.ones(N))/(max(theta)-min(theta))-2
    # theta = 4*theta/(max(theta)-min(theta))
    for i in list(range(N)):
        x = symbols('x')
        expr1 = exp(-x**2/2)/sqrt(2*math.pi)
        ans_con = integrate(expr1, (x,a_fin,a_fin*(theta[i]-b_fin)))
        if theta[i]>theta_cho:
            count+=1
            Q[i]=ans_con
        else:
            Q[i]=1-ans_con
    perc=count/N

    return Q,count,perc


if __name__ == '__main__':

    N = 24
    arr = np.loadtxt('/Users/islet/Desktop/虚拟网络/DATA/Topology/USA network.txt', delimiter=',')
    connect = arr.astype(np.int64)
    vertex = [2, 5]
    vertexResource = [5, 10]
    bandwith = [1, 5]
    calres = [100] * N  # calculate recource

    gamma_list=np.linspace(0.2,0.5,4)
    alpha_list=np.linspace(0.3,0.7,5)
    theta_cho_list = np.linspace(-2,2,11)

    name_str = '/Users/islet/Desktop/虚拟网络/DATA/Result/RES_USA network.csv'
    f = open(name_str, 'w', encoding='utf-8')
    csv_writer = csv.writer(f)
    csv_writer.writerow(['count', 'perc','gamma','alpha','theta_cho'])

    for gamma in gamma_list:
        for alpha in alpha_list:
            for theta_cho in theta_cho_list:

                [num,bor]=degree(connect)
                a=clustering(num,bor)
                [res,RE_res,P]=resource(connect,calres,bor)
                b=Markov(RE_res,P)
                RE_a = normalize(a)
                RE_b = normalize(b)
                theta = RE_a + alpha*RE_b
                # print(theta)
                Q,count,perc = Normal_Ogive(RE_a,RE_b,theta)

                csv_writer.writerow([count, perc,gamma,alpha,theta_cho])

                # fig = plt.figure(figsize=(12, 8))
                # x=theta;y=Q
                # plt.scatter(x,y,color='r', marker='x')
                # plt.show()

    f.close()

