"""Generate plots"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt

n = 8
gens = 200000

# load file
save_location = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
# B_arr = np.load(os.path.join(save_location,"single_env_8_800000.npy"))
# B_final = np.load(os.path.join(save_location,"final_B_8_800000.npy"))


def plot_interaction_coefficients(n_indiv, generations, title_text, interaction_array,show_save_ind=0):
    """Plot the interaction matrix coefficents from a given npy filename
    representating changes across t_evo generations."""
    gen_array = np.arange(0,generations)

    fig, ax = plt.subplots(figsize=(12,6))
    color_idx = 0
    for i in range(n_indiv):
        for j in range(n_indiv):
            ax.plot(gen_array, interaction_array[:,i,j])

    # ax.grid(True,ls="--")
    ax.set_title(title_text, fontsize=14)
    ax.set_xlabel("Generations", fontsize=12)
    ax.set_ylabel("Interaction Coef Values", fontsize=12)
    if show_save_ind==0:
        plt.show()
    elif show_save_ind==1:
        plt.savefig(f"interaction_coeffs_{n_indiv}_{generations}.pdf")

def plot_B_matrix_instance(B_matrix, title_text,show_save_ind=0):
    """Plot an instance of the B interaction coefficients matrix."""
    plot_time_stamp = time.strftime("%Y%m%d-%H%M%S")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title_text)
    mat_plot = ax.matshow(B_matrix)
    fig.colorbar(mat_plot)
    ax.set_xlabel("gene i")
    ax.set_ylabel("gene j")
    if show_save_ind==0:
        plt.show()
    elif show_save_ind==1:
        plt.savefig(f"b_final_{plot_time_stamp}.pdf")
