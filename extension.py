from matplotlib import pyplot as plt
import numpy as np
import scipy.integrate as scp

#variables and constants
G = 9.81
k = 1000
lengthx = 16
lengthy = 25
m = 0.01
initial = [0/lengthx,1.6/lengthx,1/lengthy,0/lengthy]

def func(t, inp, lx, ly, k, m):
    """
    Assumptions made:
        Extensible string with elastic force directly proportional to extension
        Size of bottle small enough relative to string length to be considered a material point
        No air resistance
        Motion in one axis is independent of the motion in the other (no coupling)
        The angle the string makes with the vertical at any point in either plane is small enough that the small angle approximation applies
    """

    x1 = inp[0]
    x2 = inp[1]
    y1 = inp[2]
    y2 = inp[3]
    if k == -1:
        dx1 = x2
        dx2 = -(G/lx)*x1
        dy1 = y2
        dy2 = -(G/ly)*y1
    elif k != 0: 
        dx1 = x2
        dx2 = -(k/(m*(lx**2))) * (np.sqrt(((x1**2)*(lx + (m*G/k)*(1-(x1**2/2)))**2) + ((y1**2)*(ly + (m*G/k)*(1-(y1**2/2)))**2) + (ly + (m*G/k)*(1-(y1**2/2)))**2)-np.sqrt((lx*x1)**2 + (ly*y1)**2 + (ly)**2)) * ((x1*(lx + (m*G/k)*(1-(x1**2/2)))**2 + ((m*G*x1**3)/k)*(lx + (m*G/k)*(1-(x1**2/2))))/np.sqrt(((x1**2)*(lx + (m*G/k)*(1-(x1**2/2)))**2) + ((y1**2)*(ly + (m*G/k)*(1-(y1**2/2)))**2) + ((ly + (m*G/k)*(1-(y1**2/2)))**2)) - ((x1 * lx**2)/(np.sqrt((lx*x1)**2 + (ly*y1)**2 + (ly)**2))))
        dy1 = y2 
        dy2 = -(k/(m*(ly**2))) * (np.sqrt(((x1**2)*(lx + (m*G/k)*(1-(x1**2/2)))**2) + ((y1**2)*(ly + (m*G/k)*(1-(y1**2/2)))**2) + (ly + (m*G/k)*(1-(y1**2/2)))**2)-np.sqrt((lx*x1)**2 + (ly*y1)**2 + (ly)**2)) * ((y1*(ly + (m*G/k)*(1-(y1**2/2)))**2 + ((m*G)/k)*(y1 + y1**3)*(ly + (m*G/k)*(1-(y1**2/2))))/np.sqrt(((x1**2)*(lx + (m*G/k)*(1-(x1**2/2)))**2) + ((y1**2)*(ly + (m*G/k)*(1-(y1**2/2)))**2) + ((ly + (m*G/k)*(1-(y1**2/2)))**2)) - ((y1 * ly**2)/(np.sqrt((lx*x1)**2 + (ly*y1)**2 + (ly)**2)))) 

    return[dx1, dx2, dy1, dy2]

time = (0, 60)
time_eval = np.linspace(time[0], time[1], 500)

sol0 = scp.solve_ivp(func, time, initial, args=(lengthx, lengthy, -1 , m), t_eval = time_eval)
sol1 = scp.solve_ivp(func, time, initial, args=(lengthx, lengthy, k, m), t_eval = time_eval, method="LSODA")

xplot0 = lengthx * sol0.y[0]
yplot0 = lengthy * sol0.y[2]  
xplot1 = lengthx * sol1.y[0]
yplot1 = lengthy * sol1.y[2]   

plt.plot(xplot0, yplot0, label = "inextensible")
plt.plot(xplot1, yplot1, label = "elastic coefficient =" + str(k))
plt.xlabel("x axis (primary string)")
plt.ylabel("y axis (primary string and supporting strings)")
plt.legend()
plt.show()