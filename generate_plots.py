"""Generate plots"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt

n = 8
gens = 200000
# load file
save_location = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
B_arr = np.load(os.path.join(save_location,"single_env_8_200000_20220508-133614.npy"))
B_final = np.load(os.path.join(save_location,"final_B_8_200000_20220508-133614.npy"))
# # generate colours for plots
# color = plt.cm.rainbow(np.linspace(0,1,16))


def plot_interaction_coefficients(n_indiv, generations, interaction_array,show_save_ind=0):
    """Plot the interaction matrix coefficents froma given npy filename
    representating changes across t_evo generations."""
    gen_array = np.arange(0,generations)

    fig, ax = plt.subplots(figsize=(10,6))
    color_idx = 0
    for i in range(n_indiv):
        for j in range(n_indiv):
            ax.plot(gen_array, interaction_array[:,i,j])

    # ax.grid(True,ls="--")
    ax.set_title("Single Env Interaction Matrix Coefs", fontsize=14)
    ax.set_xlabel("Generations", fontsize=12)
    ax.set_ylabel("Interaction Coef Values", fontsize=12)
    if show_save_ind==0:
        plt.show()
    elif show_save_ind==1:
        plt.savefig(f"interaction_coeffs_{n_indiv}_{generations}.pdf")

def plot_B_matrix_instance(B_matrix, show_save_ind=0):
    """Plot an instance of the B interaction coefficients matrix."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mat_plot = ax.matshow(B_final)
    fig.colorbar(mat_plot)
    ax.set_xlabel("gene i")
    ax.set_ylabel("gene j")
    if show_save_ind==0:
        plt.show()
    elif show_save_ind==1:
        plt.savefig(f"b_final_{B_matrix.shape[0]}.pdf")

# plot_B_matrix_instance(B_final)
# plot_interaction_coefficients(n, gens, B_arr)

