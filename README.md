Double Pendulum Simulation

Overview:
This Python script simulates the motion of a double pendulum using numerical integration. It provides an interactive visualization where users can adjust the lengths and masses of the pendulums and observe chaotic motion dynamics.

Features:
Real-time simulation of a double pendulum using scipy.integrate.odeint.
Adjustable pendulum lengths (L1, L2) and masses (m1, m2) via sliders.
Animated visualization of pendulum motion with a trace of the second mass.
Reset functionality to restart the simulation from the initial state.

Requirements:
Python 3.x
NumPy
Matplotlib
SciPy

Installation:
Install required dependencies using:
pip install numpy matplotlib scipy
Usage
Run the script using:
python Double_Pendulum.py

Controls:
L1 Slider: Adjusts the length of the first pendulum.
L2 Slider: Adjusts the length of the second pendulum.
m1 Slider: Adjusts the mass of the first pendulum.
m2 Slider: Adjusts the mass of the second pendulum.
Reset Button: Resets the simulation to initial conditions.

How It Works:
The script initializes the double pendulum system and solves the equations of motion using odeint.
It computes the motion of both pendulums and updates their positions at each time step.
The visualization updates dynamically, showing the pendulum movement and trace path.
Users can modify parameters in real-time to explore different behaviors.
