from matplotlib import pyplot as plt
import numpy as np
import scipy.integrate as scp

#variables and constants
G = 9.81
k = 5000000
k1 = 2500000
lengthx = 25
lengthy = 36
m = 0.1
initial = [1/lengthx,0/lengthx,1/lengthy,0/lengthy]

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

    delta1 = np.sqrt(((x1**2)*(lx + (m*G/k1)*(1-(x1**2/2)))**2) + ((y1**2)*(ly + (m*G/k)*(1-(y1**2/2)))**2) + (1-(y1**2)/2)**2 * (ly + (m*G/k)*(1-(y1**2/2)))**2)
    delta2 = np.sqrt((lx*x1)**2 + (ly*y1)**2 + (ly*(1 - (y1**2)/2))**2)
    f1 = lx + (m*G/k1)*(1-(x1**2/2))
    f2 = ly + (m*G/k)*(1-(y1**2/2))
    f3 = (1-(y1**2)/2)
    f4 = (m*G/k)

    dx1 = x2
    dx2 = (1/(m*(f1 - f4*x1**2)**2))*(3*m*f4*(x2**2)*x1*(f1 - f4*x1) - k1*(delta1 - delta2)*(((x1*f1**2 - f4*f1*x1**3)/(delta1))-((x1*lx**2)/(delta2))))

    dy1 = y2 
    dy2 = (1/(m*(f2 - f4*y1**2)**2 + m*(y1*(ly + (2*m*G/k)*f3))**2))*(-2*m*y1*y2**2 * (ly + (2*m*G/k)*f3)**2 +4*f4*y2*y1**3*(ly + (2*m*G/k)*f3) + 3*m*f4*y1*y2**2 * (f2 - f4*y1**2) + 2*m*y1*y2**2* (ly + (2*m*G/k)*f3)**2 - 4*m*f4*(y1**3)*(y2**2)*(ly + (2*m*G/k)*f3) - m*G*y1*f2 - (((m*G)**2)/k)*y1*f3 - k*(delta1-delta2)*((((1/2 * y1**3)*f2**2 - f4*y1*f2*(y1 + f3))/(delta1)) - ((lx**2 * y1**2)/(2*delta2))))

    return[dx1, dx2, dy1, dy2]

time = (0, 300)
time_eval = np.linspace(time[0], time[1], 20000)

#sol0 = scp.solve_ivp(func, time, initial, args=(lengthx, lengthy, 999999999999 , m), t_eval = time_eval)
sol1 = scp.solve_ivp(func, time, initial, args=(lengthx, lengthy, k, m), t_eval = time_eval, method="BDF")

#xplot0 = lengthx * sol0.y[0]
#yplot0 = lengthy * sol0.y[2]  
xplot1 = lengthx*sol1.y[0]
yplot1 = lengthy*sol1.y[2]   

#plt.plot(xplot0, yplot0, label = "inextensible")
plt.plot(xplot1, yplot1)
plt.xlabel("x axis (primary string)")
plt.ylabel("y axis (primary string and supporting strings)")
plt.show()