import scipy.integrate as scp
from matplotlib import pyplot as plt
import numpy as np

G = 9.81

def plane_func(t, x, l, b, m):     #general form for the equation of motion for an individual plane. 
    """
    Assumptions made:
        Inextensible string
        Size of bottle small enough relative to string length to be considered a material point
        Turbulent flow of air around the bottle
        Motion in one axis is independent of the motion in the other (no coupling)
        The angle the string makes with the vertical at any point in either plane is small enough that the small angle approximation applies
    """
    x1 = x[0]
    x2 = x[1]
    dx1 = x2
    dx2 = -(G/l) * x1 -(b*l/m)*x2*abs(x2)                   #drag expressed for turbulent flow
    return [dx1, dx2]

time = (0, 60)
time_eval = np.linspace(time[0], time[1], 500)
lengthx = 16
lengthy = 25
initialx = [0/lengthx,1.6/lengthx]             #first is position, second is velocity
initialy = [1/lengthy,0/lengthy]             #first is position, second is velocity
rad = 0.03
height = 0.2
dense = 1.2
cd = 0.6   #coefficient of drag
b = 0.5 * ((rad ** 2) *np.pi + rad*height*np.pi)*dense*cd

solx0 = scp.solve_ivp(plane_func, time, initialx, args=(lengthx, 0, 0.1), t_eval=time_eval)           #general Lissajous curves
soly0 = scp.solve_ivp(plane_func, time, initialy, args=(lengthy, 0, 0.1), t_eval=time_eval)

solx = scp.solve_ivp(plane_func, time, initialx, args=(lengthx, b, 0.1), t_eval=time_eval)
soly = scp.solve_ivp(plane_func, time, initialy, args=(lengthy, b, 0.1), t_eval=time_eval)


xplot0 = lengthx * solx0.y[0]
yplot0 = lengthy * soly0.y[0]
xplot = lengthx * solx.y[0]
yplot = lengthy * soly.y[0]

plt.plot(xplot0, yplot0, label = "no drag")
plt.plot(xplot, yplot, label = "drag coefficient =" + str(b))
plt.xlabel("x axis (primary string)")
plt.ylabel("y axis (primary string and supporting strings)")
plt.legend()
plt.show()