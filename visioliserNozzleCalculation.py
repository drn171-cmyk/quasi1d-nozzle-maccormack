# This code for visualizing the results of quasi 1D nozzle flow

# Importing libraries

import pandas as pd
import matplotlib.pyplot as plt

# Reading data

df = pd.read_csv('nozzleResults.csv') # Reading the calculated data from C++

# Defining Vectors

x = df['x'] # Position vector
Mach = df['Mach'] # Mach number vector
p = df['p'] # Static pressure vector
rho = df['rho'] # Density vector
T = df['T'] # Static temperature vector

# Creating Figure

fig, axs = plt.subplots(2, 2, figsize=(12, 8)) # 2x2 grid for 4 plots
fig.suptitle('Quasi-1D Flow Results') # Main title of the figure

# Plotting Mach Number

axs[0, 0].plot(x, Mach, color='blue') # Plotting x vs Mach
axs[0, 0].set_title('Mach Number') # Subplot title
axs[0, 0].set_xlabel('x') # x axis label
axs[0, 0].set_ylabel('Mach') # y axis label
axs[0, 0].grid(True) # Adding grid

# Plotting Static Pressure

axs[0, 1].plot(x, p, color='red') # Plotting x vs pressure
axs[0, 1].set_title('Static Pressure') # Subplot title
axs[0, 1].set_xlabel('x') # x axis label
axs[0, 1].set_ylabel('p') # y axis label
axs[0, 1].grid(True) # Adding grid

# Plotting Density

axs[1, 0].plot(x, rho, color='green') # Plotting x vs density
axs[1, 0].set_title('Density') # Subplot title
axs[1, 0].set_xlabel('x') # x axis label
axs[1, 0].set_ylabel('rho') # y axis label
axs[1, 0].grid(True) # Adding grid

# Plotting Static Temperature

axs[1, 1].plot(x, T, color='purple') # Plotting x vs temperature
axs[1, 1].set_title('Static Temperature') # Subplot title
axs[1, 1].set_xlabel('x') # x axis label
axs[1, 1].set_ylabel('T') # y axis label
axs[1, 1].grid(True) # Adding grid

# Showing the plots

plt.tight_layout() # Adjusting spaces between subplots
plt.show() # Display the figure