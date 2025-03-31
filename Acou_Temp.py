import numpy as np
import math
def Col(Vx, Vy, Vz, grad_x, grad_y, grad_z, temperature):
    Cmp = 331.3 + 0.6 * temperature;
    theta_e = np.random.rand() * 2 * np.pi;
    fine_e = math.acos( 2 * np.random.rand() - 1 );
    k = 10
    V_x = Vx - k * 0.6 * grad_x;
    V_y = Vy - k * 0.6 * grad_y;
    V_z = Vz + k * 0.6 * grad_z;
    return V_x, V_y, V_z