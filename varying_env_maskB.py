"""This scripts

Perform an extension to the varying environment experiment by
gradually masking the interaction matrix and evaluating changes in the 
capacity to generate results"""


import numpy as np
import random
import sys,os, time
from single_env import GPMappingGRN
from generate_plots import plot_B_matrix_instance, plot_interaction_coefficients
from tqdm import tqdm
import matplotlib.pyplot as plt

SAVE_PATH = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
TIME_STAMP = time.strftime("%Y%m%d-%H%M%S")
TARGET_1 = np.array([1,1,-1,-1,-1,1,-1,1])
TARGET_2 = np.array([1,-1,1,-1,1,-1,-1,-1])
N = 8
T_DEV = 10
T_EVO = 400000
T_EVO_SUBGEN = 2000
IDX = 0
FLAG = 0
varyingEnvGrn = GPMappingGRN(N, T_DEV, T_EVO)
# placeholder for interaction matrix plot shape NxNxT_evo
plot_placeholder = np.zeros((T_EVO,N,N))

np.random.seed(7)

with tqdm(T_EVO) as pbar:
    while(IDX<T_EVO):
        if FLAG==0:
            for i in range(T_EVO_SUBGEN):
                # produce phenotype
                # if IDX in DEV_STEP_FLAGS: ### EXP2 bracket D####
                #     p_star,temp_dev_arr = varyingEnvGrn.update_phenotype3(varyingEnvGrn.geno_vec, varyingEnvGrn.interaction_matrix)
                #     DEV_STEP_LIST.append(temp_dev_arr)
                #     p_star_ftns = varyingEnvGrn.calculate_fitness(p_star, TARGET_1)
                # else:
                p_star = varyingEnvGrn.update_phenotype2(varyingEnvGrn.geno_vec, varyingEnvGrn.interaction_matrix)
                p_star_ftns = varyingEnvGrn.calculate_fitness(p_star, TARGET_1)
                # mutate G and B
                mut_geno_vec = varyingEnvGrn.mutate_G()
                if random.uniform(0, 1) <= 0.067:
                    mut_B_matrix = varyingEnvGrn.mutate_B_destruct(col_start=0, col_finish=7,row_start=0, row_finish=7, single_cell=1)
                else:
                    mut_B_matrix = varyingEnvGrn.interaction_matrix
                # run evo loop on mutant genotype
                p_star_mutant = varyingEnvGrn.update_phenotype2(mut_geno_vec, mut_B_matrix)
                p_star_mutant_ft = varyingEnvGrn.calculate_fitness(p_star_mutant, TARGET_1)
                if p_star_mutant_ft>p_star_ftns:
                    varyingEnvGrn.geno_vec = mut_geno_vec
                    varyingEnvGrn.interaction_matrix = mut_B_matrix
                # save interaction matrix instance
                plot_placeholder[IDX,:,:] = varyingEnvGrn.interaction_matrix
                IDX+=1
                pbar.update()
            FLAG=1
        #### Switch target ####
        elif FLAG==1:
            for j in range(T_EVO_SUBGEN):
                # produce phenotype
                p_star2 = varyingEnvGrn.update_phenotype2(varyingEnvGrn.geno_vec, varyingEnvGrn.interaction_matrix)
                p_star_ftns2 = varyingEnvGrn.calculate_fitness(p_star2, TARGET_2)
                # mutate G and B
                mut_geno_vec2 = varyingEnvGrn.mutate_G()
                if random.uniform(0, 1) <= 0.067:
                    mut_B_matrix2 = varyingEnvGrn.mutate_B_destruct(col_start=0, col_finish=7,row_start=0, row_finish=7, single_cell=1)
                else:
                    mut_B_matrix2 = varyingEnvGrn.interaction_matrix
                # run evo loop on mutant genotype
                p_star_mutant2 = varyingEnvGrn.update_phenotype2(mut_geno_vec2, mut_B_matrix2)
                p_star_mutant_ft2 = varyingEnvGrn.calculate_fitness(p_star_mutant2, TARGET_2)
                if p_star_mutant_ft2>p_star_ftns2:
                    varyingEnvGrn.geno_vec = mut_geno_vec2
                    varyingEnvGrn.interaction_matrix = mut_B_matrix2
                # save interaction matrix instance
                plot_placeholder[IDX,:,:] = varyingEnvGrn.interaction_matrix
                IDX+=1
                pbar.update()
            FLAG=0

plot_B_matrix_instance(varyingEnvGrn.interaction_matrix,"Interaction matrix without 49 genes expressing")

bracket_e_matrix = np.zeros((30,8))
for i in range(30):
    random_input = np.random.randn(8)
    p_star_rand = varyingEnvGrn.update_phenotype2(random_input, varyingEnvGrn.interaction_matrix)
    bracket_e_matrix[i,:] = p_star_rand
    print(f"input:{random_input} ---- output: {p_star_rand}")
    print("\n")

fig= plt.figure()
fig.set_figheight(8)
fig.set_figwidth(12)
ax = fig.add_subplot()
ax.set_title("Evolved Network - Phenotypic samples - \n 49 genes not expressing", fontsize=20)
mat_plot = plt.imshow(bracket_e_matrix,aspect="auto")
plt.colorbar(mat_plot,ax=ax)
ax.set_xlabel("genes", fontsize=16)
ax.set_ylabel("samples", fontsize=16)
# plt.savefig("evolved_samples_single33notexp.pdf")
plt.show()