import numpy as np
import math 
def Top(TopX,TopZ,SUM,Glength,Gwide,G):
    g = 0.01 * G
    A = Glength;
    B = Gwide;
    Ptop = np.zeros((B,A));
    ptop = np.zeros((B,A));
    for i in range(len(TopX)):
        
        if TopX[i] <= (-Glength * g/2):
            TopX[i] = -Glength * g/2 + 0.001;
        if TopX[i] >= Glength * g/2:
            TopX[i] = Glength * g/2 - 0.001;
        if TopZ[i] <= 0:
            TopZ[i] = 0 + 0.001;
        if TopZ[i] >= Gwide * g:
            TopZ[i] = Gwide * g - 0.001;

        m =  math.floor( (Gwide * g - TopZ[i])/g );
        n = math.floor( (TopX[i] + Glength * g/2)/g );
        Ptop[m,n] = Ptop[m,n] + 1;
    ###################################################
    top = np.zeros(A * B);
    ptop = Ptop/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            top[t] = ptop[i,j];
    return Ptop,top

def Lside(LZ,LY,SUM,Gheight,Gwide,G):

    g = 0.01 * G
    A = Gwide;
    B = Gheight;
    Plside = np.zeros((B,A));
    plside = np.zeros((B,A));
    
    for i in range(len(LZ)):
        if LY[i] <= (-Gheight * g/2):
            LY[i] = -Gheight * g/2 + 0.001;
        if LY[i] >= Gheight * g/2:
            LY[i] = Gheight * g/2 - 0.001;
        if LZ[i] <= 0:
            LZ[i] = 0 + 0.001;
        if LZ[i] >= Gwide * g:
            LZ[i] = Gwide * g - 0.001;

        m =  math.floor( (Gheight * g/2 - LY[i])/g );
        n = math.floor( LZ[i]/g );
        Plside[m,n] = Plside[m,n] + 1;
    ###################################################
    lside = np.zeros(B * A);
    plside = Plside/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            lside[t] = plside[i,j];
    return Plside,lside
            
def Rside(RZ,RY,SUM,Gheight,Gwide,G):

    g = 0.01 * G    
    A = Gwide;
    B = Gheight;
    Prside = np.zeros((B,A));
    prside = np.zeros((B,A));
    
    for i in range(len(RZ)):
        if RY[i] <= (-Gheight * g/2):
            RY[i] = -Gheight * g/2 + 0.001;
        if RY[i] >= Gheight * g/2:
            RY[i] = Gheight * g/2 - 0.001;
        if RZ[i] <= 0:
            RZ[i] = 0 + 0.001;
        if RZ[i] >= Gwide * g:
            RZ[i] = Gwide * g - 0.001;

        m =  math.floor( (Gheight * g/2 - RY[i])/g );
        n = math.floor( RZ[i]/g );
        Prside[m,n] = Prside[m,n] + 1;
    ###################################################
    rside = np.zeros(B * A);
    prside = Prside/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            rside[t] = prside[i,j];
    return Prside,rside

def PG(PGX,PGY,SUM,Gheight,Glength,G):

    g = 0.01 * G    
    A = Glength;
    B = Gheight;
    Ppg = np.zeros((B,A));
    ppg = np.zeros((B,A));
    for i in range(len(PGX)):
        
        if PGY[i] <= (-Gheight * g/2):
            PGY[i] = -Gheight * g/2 + 0.001;
        if PGY[i] >= Gheight * g/2:
            PGY[i] = Gheight *g/2 - 0.001;
        if PGX[i] <= (-Glength * g/2):
            PGX[i] = -Glength * g/2 + 0.001;
        if PGX[i] >= Glength * g/2:
            PGX[i] = Glength * g/2 - 0.001;

        m =  math.floor( (Gheight * g/2 - PGY[i])/g );
        n = math.floor( (Glength * g/2 + PGX[i])/g );
        Ppg[m,n] = Ppg[m,n] + 1;
    ###################################################
    pg = np.zeros(B * A);
    ppg = Ppg/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            pg[t] = ppg[i,j];
    return Ppg,pg

def Back(BackX,BackY,SUM,Gheight,Glength,G):

    g = 0.01 * G    
    A = Glength;
    B = Gheight;
    Pback = np.zeros((B,A));
    pback = np.zeros((B,A));
    for i in range(len(BackX)):

        if BackY[i] <= (-Gheight * g/2):
            BackY[i] = -Gheight * g/2 + 0.001;         
        if BackY[i] >= Gheight * g/2:
            BackY[i] = Gheight * g/2 - 0.001;           
        if BackX[i] <= (-Glength * g/2):
            BackX[i] = -Glength * g/2 + 0.001;           
        if BackX[i] >= Glength * g/2:
            BackX[i] = Glength * g/2 - 0.001;
            
        m =  math.floor( (Gheight * g/2 - BackY[i])/g );
        n = math.floor( (Glength * g/2 + BackX[i])/g );

        Pback[m,n] = Pback[m,n] + 1;
    ###################################################
    back = np.zeros(B * A);
    pback = Pback/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            back[t] = pback[i,j];
    return Pback,back

def Bottom(BottomX,BottomZ,SUM,Glength,Gwide,G):

    g = 0.01 * G    
    A = Glength;
    B = Gwide;
    Pbot = np.zeros((B,A));
    pbot = np.zeros((B,A));
    for i in range(len(BottomX)):
        
        if BottomX[i] <= (-Glength * g/2):
            BottomX[i] = -Glength * g/2 + 0.001;
        if BottomX[i] >= Glength * g/2:
            BottomX[i] = Glength * g/2-0.001;
        if BottomZ[i] <= 0:
            BottomZ[i] = 0 + 0.001;
        if BottomZ[i] >= Gwide * g:
            BottomZ[i] = Gwide * g - 0.001;

        m =  math.floor( (Gwide * g - BottomZ[i])/g );
        n = math.floor( (BottomX[i] + Glength * g/2)/g );
        Pbot[m,n] = Pbot[m,n] + 1;
    ###################################################
    bot = np.zeros(A * B);
    pbot = Pbot/SUM;
    for i in range(B):
        for j in range(A):           
            t = i * A + j;
            bot[t] = pbot[i,j];
    return Pbot,bot