import numpy as np
import random
import sys,os, time
from single_env import GPMappingGRN
from generate_plots import plot_B_matrix_instance, plot_interaction_coefficients
from tqdm import tqdm

SAVE_PATH = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
TIME_STAMP = time.strftime("%Y%m%d-%H%M%S")
TARGET_1 = np.array([1,1,-1,-1,-1,1,-1,1])
TARGET_2 = np.array([1,-1,1,-1,1,-1,-1,-1])
N = 8
T_DEV = 10
T_EVO = 800000
T_EVO_SUBGEN = 2000
IDX = 0
FLAG = 0

varyingEnvGrn = GPMappingGRN(N, T_DEV, T_EVO)
# placeholder for interaction matrix plot shape NxNxT_evo
plot_placeholder = np.zeros((T_EVO,N,N))

with tqdm(T_EVO) as pbar:
    while(IDX<T_EVO):
        if FLAG==0:
            for i in range(T_EVO_SUBGEN):
                # produce phenotype
                p_star = varyingEnvGrn.update_phenotype2(varyingEnvGrn.geno_vec, varyingEnvGrn.interaction_matrix)
                p_star_ftns = varyingEnvGrn.calculate_fitness(p_star, TARGET_1)
                # mutate G and B
                mut_geno_vec = varyingEnvGrn.mutate_G()
                if random.uniform(0, 1) <= 0.067:
                    mut_B_matrix = varyingEnvGrn.mutate_B()
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
                    mut_B_matrix2 = varyingEnvGrn.mutate_B()
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

plot_interaction_coefficients(N, T_EVO, plot_placeholder)
plot_B_matrix_instance(varyingEnvGrn.interaction_matrix)