"""Generate plots"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt

n = 8
generations = 30000
# load file
save_location = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
arr = np.load(os.path.join(save_location,"single_env_8_30000_20220508-125943.npy"))
# generate colours for plots
color = plt.cm.rainbow(np.linspace(0,1,16))
print(color)

gen_array = np.arange(0,generations)

# Plot
fig, ax = plt.subplots(figsize=(10,6))
color_idx = 0
for i in range(n):
    for j in range(n):
        ax.plot(gen_array, arr[:,i,j])

# ax.grid(True,ls="--")
ax.set_title("Single Env Interaction Matrix Coefs", fontsize=14)
ax.set_xlabel("Generations", fontsize=12)
ax.set_ylabel("Interaction Coef Values", fontsize=12)
plt.show()