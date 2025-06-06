import numpy as np
import matplotlib.pyplot as plt

def lissajous_curve(length = 6,period_x = 1, amplitude_x = 1, period_y = 1, amplitude_y = 1):
    """
    Returns a lissajous curve with the specified period and amplitude in the x-direction and y-direction
    """
    time = np.linspace(0, length, 400)
    xs = [np.sin(t * (2*np.pi) / period_x) * amplitude_x for t in time]
    ys = [np.sin(t * (2*np.pi) / period_y) * amplitude_y for t in time]
    plt.plot(xs, ys)


def pendulum_1D(string_length = 1, initial_angular_velocity = 1, initial_angle_with_downward_vertical = 0, duration = 3):
    """
    Returns a graph of the 1D pendulum with time as the x-axis and horizontal displacement as the y-axis
    Assumptions made includes:
    1. the bottle is a particle with constant mass so friction is negligible;
    2. the angle between the string and the downward vertical is small such that the SHM model is used;
    3. the string is inextensible and taut throughout the motion;
    """
    g = 9.81
    time = np.linspace(0, duration, 400)
    theta = [initial_angular_velocity * np.sin(g / string_length * t) + initial_angle_with_downward_vertical * np.cos(g / string_length * t) for t in time]
    xs = [angle * string_length for angle in theta]
    plt.plot(time,xs)


def pendulum_2D(string_length_x = 1, string_length_y = 1, initial_angular_velocity_x = 1, initial_angular_velocity_y = 1, initial_position = (0,0), duration = 3):
    """
    Returns a graph of the projection of the 2D pendulum on the xy-plane.
    Assumptions made includes:
    1. the bottle is a particle with constant mass so friction is negligible;
    2. the angle between the string and the downward vertical is small such that the SHM model is used;
    3. the string is inextensible and taut throughout the motion;
    4. the bottle's motion in the x-direction will not affect its motion in the y-direction
    """
    g = 9.81
    time = np.linspace(0, duration, 400)
    initial_x = initial_position[0]
    initial_y = initial_position[1]
    angle_x = [initial_angular_velocity_x * ((string_length_x / g) ** (1/2)) * np.sin(((g / string_length_x) ** (1/2)) * t)
               + initial_x / string_length_x * np.cos(((g / string_length_x) ** (1/2)) * t) for t in time]
    xs = [angle * string_length_x for angle in angle_x]
    angle_y = [initial_angular_velocity_y * ((string_length_y / g) ** (1/2)) * np.sin(((g / string_length_y) ** (1/2)) * t)
               + initial_y / string_length_y * np.cos(((g / string_length_y) ** (1/2)) * t) for t in time]
    ys = [angle * string_length_y for angle in angle_y]
    plt.plot(xs, ys)


def pendulum_2D_with_friction(string_length_x = 1, string_length_y = 1, initial_angular_velocity_x = 1, initial_angular_velocity_y = 1, initial_position = (0,0), duration = 3, factor_of_friction = 0.1, mass_of_the_particle = 1):
    """
    Returns a graph of the projection of the 2D pendulum on the xy-plane.
    Assumptions made includes:
    1. friction is proportional to the speed of the particle(damping force in SHM);
    2. the bottle has constant mass
    3. the angle between the string and the downward vertical is small such that the SHM model is used;
    4. the string is inextensible and taut throughout the motion;
    5. the bottle's motion in the x-direction will not affect its motion in the y-direction
    """
    g = 9.81
    time = np.linspace(0, duration, 400)
    initial_x = initial_position[0]
    initial_y = initial_position[1]
    c = factor_of_friction / mass_of_the_particle
    initial_angle_x = initial_x / string_length_x
    initial_angle_y = initial_y / string_length_y

    delta_x = c ** 2 - 4 * g / string_length_x
    delta_y = c ** 2 - 4 * g / string_length_y

    if delta_x < 0:
        delta_x_sqrt = (- delta_x) ** (1/2)
        angle_x = [np.exp(- c / 2 * t) * ((initial_angular_velocity_x + c / 2 * initial_angle_x) * 2 / delta_x_sqrt * np.sin(delta_x_sqrt / 2 * t) + initial_angle_x * np.cos(delta_x_sqrt / 2 * t)) for t in time]
    elif delta_x == 0:
        angle_x = [(c / 2 * initial_angle_x + initial_angular_velocity_x) * t * np.exp(- c / 2 * t) + initial_angle_x * np.exp(- c / 2 * t) for t in time]
    else:
        k1 = - c / 2 + delta_x ** (1/2) / 2
        k2 = - c / 2 - delta_x ** (1/2) / 2
        A = (initial_angle_x * k2 - initial_angular_velocity_x) / (k2 - k1)
        B = (initial_angle_x * k1 - initial_angular_velocity_x) / (k1 - k2)
        angle_x = [A * np.exp(k1 * t) + B * np.exp(k2 * t) for t in time]
    xs = [angle * string_length_x for angle in angle_x]

    if delta_y < 0:
        delta_y_sqrt = (- delta_y) ** (1/2)
        angle_y = [np.exp(- c / 2 * t) * ((initial_angular_velocity_y + c / 2 * initial_angle_y) * 2 / delta_y_sqrt * np.sin(delta_y_sqrt / 2 * t) + initial_angle_y * np.cos(delta_y_sqrt / 2 * t)) for t in time]
    elif delta_y == 0:
        angle_y = [(c / 2 * initial_angle_y + initial_angular_velocity_y) * t * np.exp(- c / 2 * t) + initial_angle_y * np.exp(- c / 2 * t) for t in time]
    else:
        k1 = - c / 2 + delta_y ** (1/2) / 2
        k2 = - c / 2 - delta_y ** (1/2) / 2
        A = (initial_angle_y * k2 - initial_angular_velocity_y) / (k2 - k1)
        B = (initial_angle_y * k1 - initial_angular_velocity_y) / (k1 - k2)
        angle_y = [A * np.exp(k1 * t) + B * np.exp(k2 * t) for t in time]
    ys = [angle * string_length_y for angle in angle_y]

    plt.plot(xs, ys)

#Changing the initial position of the particle stretches  corresponding direction if the particle starts from rest
#pendulum_2D(40,90,0,0,(1,1),60)
#pendulum_2D(40,90,0,0,(1,2),60)
#plt.show()


#pendulum_2D(10,40,0,0,(1,1),60)
#pendulum_2D(10,40,0.01,0.01,(1,1),60)
#plt.show()


#Changing the initial angular velocity of the particle changes the trajectory completely, but the trajectory is still similar to the original lissajous curve
#See Dual Axis Rotation
#Projection of a 3D object rotating along its axis???
#pendulum_2D(40,90,0,0,(1,1),60)
#lissajous_curve(10,2,1,3,1)
#pendulum_2D(40,90,0.01,0,(1,1),60)
#lissajous_curve(10,2,1,3,1)
#pendulum_2D(40,90,0,0.01,(1,1),60)
#lissajous_curve(10,2,1,3,3)
#plt.show()

#pendulum_2D(90,160,0,0,(1,1),60)
#lissajous_curve(30,3,1,4,1)
#pendulum_2D(40,90,0.01,0,(1,1),60)
#lissajous_curve(10,2,1,3,1)
#pendulum_2D(40,90,0,0.01,(1,1),60)
#lissajous_curve(10,2,1,3,3)
#plt.show()

pendulum_2D_with_friction(4,9,0,0,(1,1),60,0.01,3)
plt.show()
