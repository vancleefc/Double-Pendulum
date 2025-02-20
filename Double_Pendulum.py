import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, Slider

# Constants
g = 9.81  # Acceleration due to gravity (m/s^2)

# Initial conditions [theta1, omega1, theta2, omega2]
initial_conditions = [np.pi / 2, 0, np.pi / 2, 0]  # Start with both pendulums at 90 degrees

# Time parameters
dt = 0.01  # Time step (s)
t_max = 10  # Maximum time for each ODE solve (s)
t_eval = np.arange(0, t_max, dt)  # Time array for ODE solver

# Initial parameters
L1_init = 1.0  # Initial length of the first pendulum (m)
L2_init = 1.0  # Initial length of the second pendulum (m)
m1_init = 1.0  # Initial mass of the first pendulum (kg)
m2_init = 1.0  # Initial mass of the second pendulum (kg)

# Function to compute derivatives of the state vector
def derivatives(state, t, L1, L2, m1, m2):
    theta1, omega1, theta2, omega2 = state

    # Equations of motion for a double pendulum
    delta_theta = theta2 - theta1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta_theta) ** 2
    den2 = (L2 / L1) * den1

    dtheta1_dt = omega1
    domega1_dt = (
        m2 * L1 * omega1**2 * np.sin(delta_theta) * np.cos(delta_theta)
        + m2 * g * np.sin(theta2) * np.cos(delta_theta)
        + m2 * L2 * omega2**2 * np.sin(delta_theta)
        - (m1 + m2) * g * np.sin(theta1)
    ) / den1

    dtheta2_dt = omega2
    domega2_dt = (
        -L1 / L2 * omega1**2 * np.sin(delta_theta)
        - g * np.sin(theta2)
        + (m1 + m2) / m2 * g * np.sin(theta1) * np.cos(delta_theta)
        - L1 / L2 * omega1**2 * np.sin(delta_theta) * np.cos(delta_theta)
    ) / den2

    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt]

# Solve the ODE for the initial time period
solution = odeint(derivatives, initial_conditions, t_eval, args=(L1_init, L2_init, m1_init, m2_init))

# Extract angles and compute positions
theta1 = solution[:, 0]
theta2 = solution[:, 2]

x1 = L1_init * np.sin(theta1)
y1 = -L1_init * np.cos(theta1)

x2 = x1 + L2_init * np.sin(theta2)
y2 = y1 - L2_init * np.cos(theta2)

# Set up the plot
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.4)  # Adjust layout for sliders and reset button
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-2.5, 2.5)
ax.set_aspect("equal")

# Plot the pendulum rods
pendulum_line, = ax.plot([], [], "-", lw=2, color="black")

# Plot the balls (m1 and m2)
m1_dot, = ax.plot([], [], 'o', markersize=10, color='blue')  # Blue ball for m1
m2_dot, = ax.plot([], [], 'o', markersize=10, color='blue')  # Blue ball for m2

# Plot the trace of the second pendulum's tip
trace_line, = ax.plot([], [], "-", lw=1, color="red")

# Initialize the trace
trace_x, trace_y = [], []

# Function to reset the simulation
def reset(event):
    global initial_conditions, solution, theta1, theta2, x1, y1, x2, y2, trace_x, trace_y
    
    # Reset initial conditions
    initial_conditions = [np.pi / 2, 0, np.pi / 2, 0]  # Reset to initial state
    
    # Solve the ODE for the initial time period
    solution = odeint(derivatives, initial_conditions, t_eval, args=(L1_slider.val, L2_slider.val, m1_slider.val, m2_slider.val))
    
    # Extract angles and compute positions
    theta1 = solution[:, 0]
    theta2 = solution[:, 2]
    x1 = L1_slider.val * np.sin(theta1)
    y1 = -L1_slider.val * np.cos(theta1)
    x2 = x1 + L2_slider.val * np.sin(theta2)
    y2 = y1 - L2_slider.val * np.cos(theta2)
    
    # Clear the trace
    trace_x.clear()
    trace_y.clear()
    trace_line.set_data([], [])
    
    # Update the pendulum rods and balls
    pendulum_line.set_data([0, x1[0], x2[0]], [0, y1[0], y2[0]])
    m1_dot.set_data([x1[0]], [y1[0]])  # Pass as sequences
    m2_dot.set_data([x2[0]], [y2[0]])  # Pass as sequences
    
    # Update marker sizes based on mass
    m1_dot.set_markersize(5 + 5 * m1_slider.val)  # Scale m1 marker size
    m2_dot.set_markersize(5 + 5 * m2_slider.val)  # Scale m2 marker size
    
    plt.draw()

# Add a reset button
ax_reset = plt.axes([0.8, 0.05, 0.1, 0.04])
reset_button = Button(ax_reset, 'Reset')
reset_button.on_clicked(reset)

# Add sliders for L1, L2, m1, and m2
ax_L1 = plt.axes([0.2, 0.25, 0.6, 0.03])
L1_slider = Slider(ax_L1, 'L1 (m)', 0.1, 2.0, valinit=L1_init)

ax_L2 = plt.axes([0.2, 0.20, 0.6, 0.03])
L2_slider = Slider(ax_L2, 'L2 (m)', 0.1, 2.0, valinit=L2_init)

ax_m1 = plt.axes([0.2, 0.15, 0.6, 0.03])
m1_slider = Slider(ax_m1, 'm1 (kg)', 0.1, 2.0, valinit=m1_init)

ax_m2 = plt.axes([0.2, 0.10, 0.6, 0.03])
m2_slider = Slider(ax_m2, 'm2 (kg)', 0.1, 2.0, valinit=m2_init)

# Function to update the simulation when sliders are changed
def update(val):
    global solution, theta1, theta2, x1, y1, x2, y2, trace_x, trace_y
    
    # Solve the ODE with the new parameters
    solution = odeint(derivatives, initial_conditions, t_eval, args=(L1_slider.val, L2_slider.val, m1_slider.val, m2_slider.val))
    
    # Extract angles and compute positions
    theta1 = solution[:, 0]
    theta2 = solution[:, 2]
    x1 = L1_slider.val * np.sin(theta1)
    y1 = -L1_slider.val * np.cos(theta1)
    x2 = x1 + L2_slider.val * np.sin(theta2)
    y2 = y1 - L2_slider.val * np.cos(theta2)
    
    # Clear the trace
    trace_x.clear()
    trace_y.clear()
    trace_line.set_data([], [])
    
    # Update the pendulum rods and balls
    pendulum_line.set_data([0, x1[0], x2[0]], [0, y1[0], y2[0]])
    m1_dot.set_data([x1[0]], [y1[0]])  # Pass as sequences
    m2_dot.set_data([x2[0]], [y2[0]])  # Pass as sequences
    
    # Update marker sizes based on mass
    m1_dot.set_markersize(5 + 5 * m1_slider.val)  # Scale m1 marker size
    m2_dot.set_markersize(5 + 5 * m2_slider.val)  # Scale m2 marker size
    
    plt.draw()

# Attach the update function to the sliders
L1_slider.on_changed(update)
L2_slider.on_changed(update)
m1_slider.on_changed(update)
m2_slider.on_changed(update)

# Function to update the animation
def animate(i):
    global solution, theta1, theta2, x1, y1, x2, y2, initial_conditions

    # Update pendulum rods and balls
    pendulum_line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    m1_dot.set_data([x1[i]], [y1[i]])  # Pass as sequences
    m2_dot.set_data([x2[i]], [y2[i]])  # Pass as sequences
    
    # Update trace
    trace_x.append(x2[i])
    trace_y.append(y2[i])
    trace_line.set_data(trace_x, trace_y)
    
    # Check if we need to solve the ODE for the next time period
    if i == len(t_eval) - 1:
        # Update initial conditions for the next solve
        initial_conditions = solution[-1]
        
        # Solve the ODE for the next time period
        solution = odeint(derivatives, initial_conditions, t_eval, args=(L1_slider.val, L2_slider.val, m1_slider.val, m2_slider.val))
        
        # Extract angles and compute positions
        theta1 = solution[:, 0]
        theta2 = solution[:, 2]
        x1 = L1_slider.val * np.sin(theta1)
        y1 = -L1_slider.val * np.cos(theta1)
        x2 = x1 + L2_slider.val * np.sin(theta2)
        y2 = y1 - L2_slider.val * np.cos(theta2)
    
    return pendulum_line, m1_dot, m2_dot, trace_line

# Create animation
ani = FuncAnimation(fig, animate, frames=len(t_eval), interval=10, blit=True, cache_frame_data=False)

plt.show()