import matplotlib.pyplot as plt
import numpy as np
import math 
def initParticlePos(main,N,Glength,Gwide,Gheight,Lx,Ly,Lz,i,j,G):
    X = np.zeros(N); 
    Y = np.zeros(N);
    Z = np.zeros(N);
    g = 0.01 * G 
    Half = 0;
###TOP   
    if main ==1:
        # for i in range(wide):
        #     for j in range(length):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            X[k] = j * g + g/2 - Lx/2 + Half;
            Z[k] = (Gwide - i - 1) * g + g/2 + Half;
            Y[k] = Ly/2;
###LSIDE   
    if main ==2: 
        # for i in range(height):
        #     for j in range(wide):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            Z[k] = j * g + g/2  + Half;
            Y[k] = (Gheight - i - 1) * g + g/2 -Ly/2 + Half;   
            X[k] = -Lx/2;        
###RSIDE
    if main ==3:
        # for i in range(height):
        #     for j in range(wide):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            Z[k] = j * g + g/2  + Half;
            Y[k] = (Gheight - i - 1) * g + g/2 -Ly/2 + Half;   
            X[k] = Lx/2;  
###PG 
    if main ==4:
        # for i in range(height):
        #     for j in range(length):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            X[k] = j * g + g/2 - Lx/2 + Half;
            Y[k] = (Gheight - i - 1) * g + g/2 -Ly/2 + Half;
            Z[k] = Lz;
###BOT
    if main ==5:
        # for i in range(wide):
        #     for j in range(length):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            X[k] = j * g + g/2 - Lx/2 + Half;
            Z[k] = (Gwide - i - 1) * g + g/2 + Half;
            Y[k] = -Ly/2;
###BAC
    if main ==6:
        # for i in range(height):
        #     for j in range(length):
        for k in range(N):
            # Half = (2 * np.random.rand() - 1) * g/2;
            X[k] = j * g + g/2 - Lx/2 + Half;
            Y[k] = (Gheight - i - 1) * g + g/2 -Ly/2 + Half;
            Z[k] = 0;
            
    return X,Y,Z
