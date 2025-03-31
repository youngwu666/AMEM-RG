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
    ptop, top = wall.Top(TopX, TopZ, len(N), Glength, Gwide ,G)
    plside, lside = wall.Lside(LsideZ, LsideY, len(N), Gheight, Gwide ,G)
    prside, rside = wall.Rside(RsideZ, RsideY, len(N), Gheight, Gwide ,G)
    pfore, fore = wall.PG(ForeX, ForeY, len(N), Gheight, Glength ,G)
    pbot, bot = wall.Bottom(BotX, BotZ, len(N), Glength, Gwide ,G)
    pback, back = wall.Back(BackX, BackY, len(N), Gheight, Glength ,G)

    O = np.concatenate((top, lside, rside, fore, bot, back), axis=0)
    return O

def process_grid1(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(1, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(1, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O

def process_grid2(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(2, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(2, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O

def process_grid3(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(3, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(3, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O

def process_grid4(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(4, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(4, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O

def process_grid5(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(5, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(5, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O

def process_grid6(i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    X, Y, Z = initPos.initParticlePos(6, N, Glength, Gwide, Gheight, Lx, Ly, Lz, i, j, G)
    O = reflect(6, N, Lx, Ly, Lz, Glength, Gwide, Gheight, X, Y, Z ,G)
    return O
######################################################
def TOP(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer top')
    C1 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid1, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gwide) for j in range(Glength)]
        for future in futures:
            C1 = np.vstack((C1, future.result()))
    print('top complete')
    return C1

def LSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer lside')
    C2 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid2, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gheight) for j in range(Gwide)]
        for future in futures:
            C2 = np.vstack((C2, future.result()))
    print('lside complete')
    return C2

def RSIDE(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer rside')
    C3 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid3, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gheight) for j in range(Gwide)]
        for future in futures:
            C3 = np.vstack((C3, future.result()))
    print('rside complete')
    return C3

def PG(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer pg')
    C4 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid4, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gheight) for j in range(Glength)]
        for future in futures:
            C4 = np.vstack((C4, future.result()))
    print('pg complete')
    return C4

def BOT(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer bot')
    C5 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid5, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gwide) for j in range(Glength)]
        for future in futures:
            C5 = np.vstack((C5, future.result()))
    print('bot complete')
    return C5

def BAC(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G):
    print('computer bac')
    C6 = np.empty((0, NUM))
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_grid6, i, j, N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G) for i in range(Gheight) for j in range(Glength)]
        for future in futures:
            C6 = np.vstack((C6, future.result()))
    print('bac complete')
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

    num_processes = 6
    if __name__ == '__main__':
        multiprocessing.freeze_support()
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = [
                pool.apply_async(TOP, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)),
                pool.apply_async(LSIDE, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)),
                pool.apply_async(RSIDE, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)),
                pool.apply_async(PG, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)),
                pool.apply_async(BOT, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G)),
                pool.apply_async(BAC, args=(N, Glength, Gwide, Gheight, Lx, Ly, Lz, NUM, G))
            ]
            results = [result.get() for result in results]
        RG = np.vstack(results)
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