import initParticlePos as initPos
import numpy as np
import math
import WallProcess as wall
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
import time
import Acou_Temp as process

def initParticle(Cmp,N,main):
    theta = np.zeros(N);
    fine = np.zeros(N);
    Vx = np.zeros(N);
    Vy = np.zeros(N);
    Vz = np.zeros(N);
    if main == 1:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]);
            Vy[i] = -Cmp  * math.cos(theta[i]);
            Vx[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);
    
    if main ==2:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]);
            Vy[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);   
            Vx[i] = Cmp  * math.cos(theta[i]); 
             
    if main ==3:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]);
            Vy[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);   
            Vx[i] = -Cmp  * math.cos(theta[i]);  
                  
    if main ==4:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = -Cmp  * math.cos(theta[i]);
            Vy[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);   
            Vx[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]);        
        
    if main ==5:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]);  
            Vy[i] = Cmp  * math.cos(theta[i]);
            Vx[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);     

    if main ==6:
        for i in range(N):
            theta[i] = math.asin(math.pow(np.random.rand(),1/2));
            fine[i] = 2 * math.pi * np.random.rand();
            Vz[i] = Cmp  * math.cos(theta[i]);
            Vy[i] = Cmp  * math.sin(theta[i]) * math.cos(fine[i]); 
            Vx[i] = Cmp  * math.sin(theta[i]) * math.sin(fine[i]);  
    
    return Vx,Vy,Vz
##############################################################    
def reflect(main, Num, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G):
    Ts = 1e-4
    T = 15  # ℃
    Cmp = 331.4 + 0.6 * T
    N = np.arange(Num)
    #
    Vx, Vy, Vz = initParticle(Cmp, Num, main)
    maxtime = 2000
    #
    ForeX, ForeY = [], []
    RsideY, RsideZ = [], []
    LsideY, LsideZ = [], []
    TopX, TopZ = [], []
    BotX, BotZ = [], []
    BackX, BackY = [], []
    #
    for i in N:
        for t in range(maxtime):
            X[i] += Vx[i] * Ts
            Y[i] += Vy[i] * Ts
            Z[i] += Vz[i] * Ts

            if Z[i] > Lz:
                Vx[i] = Vy[i] = Vz[i] = 0
                Z[i] = Lz
                ForeX.append(X[i])
                ForeY.append(Y[i])
                break
            if X[i] > Lx / 2:
                Vx[i] = Vy[i] = Vz[i] = 0
                X[i] = Lx / 2
                RsideZ.append(Z[i])
                RsideY.append(Y[i])
                break
            if X[i] < -Lx / 2:
                Vx[i] = Vy[i] = Vz[i] = 0
                X[i] = -Lx / 2
                LsideZ.append(Z[i])
                LsideY.append(Y[i])
                break
            if Y[i] > Ly / 2:
                Vx[i] = Vy[i] = Vz[i] = 0
                Y[i] = Ly / 2
                TopX.append(X[i])
                TopZ.append(Z[i])
                break
            if Y[i] < -Ly / 2:
                Vx[i] = Vy[i] = Vz[i] = 0
                Y[i] = -Ly / 2
                BotZ.append(Z[i])
                BotX.append(X[i])
                break
            if Z[i] < 0:
                Vx[i] = Vy[i] = Vz[i] = 0
                Z[i] = 0
                BackY.append(Y[i])
                BackX.append(X[i])
                break
            ###################################################   
            ##################################################
    ptop, top = wall.Top(TopX, TopZ, len(N), Glength, Gwide ,G)
    plside, lside = wall.Lside(LsideZ, LsideY, len(N), Gheight, Gwide ,G)
    prside, rside = wall.Rside(RsideZ, RsideY, len(N), Gheight, Gwide ,G)
    pfore, fore = wall.PG(ForeX, ForeY, len(N), Gheight, Glength ,G)
    pbot, bot = wall.Bottom(BotX, BotZ, len(N), Glength, Gwide ,G)
    pback, back = wall.Back(BackX, BackY, len(N), Gheight, Glength ,G)

    O = np.concatenate((top, lside, rside, fore, bot, back), axis=0)
    return O

######################################################
def TOP(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C1 = np.empty((0,NUM));
    for i in range(Gwide):
        for j in range(Glength):
            X1,Y1,Z1 = initPos.initParticlePos(1,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);  #生成点的坐标
            O = reflect(1,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X1,Y1,Z1,G);
            C1 = np.vstack((C1,O))
    return C1

def LSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C2 = np.empty((0, NUM))
    for i in range(Gheight):
        for j in range(Gwide):
            X2,Y2,Z2 = initPos.initParticlePos(2,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);
            O = reflect(2,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X2,Y2,Z2,G);
            C2 = np.vstack((C2,O))  
    return C2

def RSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C3 = np.empty((0, NUM))
    for i in range(Gheight):
        for j in range(Gwide):
            X3,Y3,Z3 = initPos.initParticlePos(3,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);
            O = reflect(3,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X3,Y3,Z3,G);
            C3 = np.vstack((C3,O))      
    return C3

def PG(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C4 = np.empty((0,NUM));
    for i in range(Gheight):
        for j in range(Glength):
            X4,Y4,Z4 = initPos.initParticlePos(4,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);
            O = reflect(4,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X4,Y4,Z4,G);
            C4 = np.vstack((C4,O)) 
    return C4

def BOT(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C5 = np.empty((0,NUM));
    for i in range(Gwide):
        for j in range(Glength):
            X5,Y5,Z5 = initPos.initParticlePos(5,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);
            O = reflect(5,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X5,Y5,Z5,G);
            C5 = np.vstack((C5,O)) 
    return C5

def BAC(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    C6 = np.empty((0,NUM));
    for i in range(Gheight):
        for j in range(Glength):
            X6,Y6,Z6 = initPos.initParticlePos(6,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G);
            O = reflect(6,N,Lx,Ly,Lz,Glength,Gwide,Gheight,X6,Y6,Z6,G);
            C6 = np.vstack((C6,O)) 
    return C6
######################################################
def main(N, length, wide, height, G):
    
    Glength = int(length / G)
    Gwide = int(wide / G)
    Gheight = int(height / G)

    NUM = Glength * Gwide * 2 + Glength * Gheight * 2 + Gwide * Gheight * 2
    Lx = length * 0.01
    Ly = height * 0.01
    Lz = wide * 0.01

    RG = np.empty((0,NUM));
    TOP_res=TOP(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)
    LSIDE_res=LSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)
    RSIDE_res=RSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)
    PG_res=PG(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)
    BOT_res=BOT(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)
    BAC_res=BAC(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)

    RG = np.vstack((RG,TOP_res));
    RG = np.vstack((RG,LSIDE_res))  
    RG = np.vstack((RG,RSIDE_res))  
    RG = np.vstack((RG,PG_res))  
    RG = np.vstack((RG,BOT_res))  
    RG = np.vstack((RG,BAC_res))
    ########################
    data = RG
    np.savez_compressed('./Data/VacRG.npz', data=data)

if __name__ == '__main__':
    start_time = time.time() 
    length = 210
    wide = 1900
    height = 300
    G = 15
    main(50, length, wide, height, G)
    end_time = time.time() 
    run_time = end_time - start_time
    print(f"code computer：{run_time}second")