import numpy as np
import random
import sys,os, time

save_location = os.path.join(os.getcwd(), "git_projects\evo-phenotypic-mem\plot_source")
time_stamp = time.strftime("%Y%m%d-%H%M%S")

class GPMappingGRN():

    def __init__(self, indiv_size, dev_steps, evo_steps):
        self.vec_size = indiv_size
        self.geno_vec = np.zeros(indiv_size)
        self.interaction_matrix = np.zeros((indiv_size, indiv_size))
        # self.pheno_vec = np.zeros(indiv_size)
        self.t_dev = dev_steps
        self.t_evo = evo_steps

    def mutate_G(self):
        """Mutation of the genotype vector. We sample randomly an index of 
        the genotype vector and then add a value from (-0.1,0.1) to the sampled 
        index. However, direct effects are capped on +/-1 so we run a check before
        updating the mutant genotype vector."""
        
        geno_vec_mutant = self.geno_vec.copy()
        #print(f"Initial geno vec: {geno_vec_mutant}")
        update_idx = random.randint(0, self.vec_size - 1) # select a random index uniformly
        #print(f"Index selected: {update_idx}")
        mu_1 = np.random.uniform(-0.1, 0.1) # sample mu_1 from (-0.1,0.1)
        #print(f"Value sampled: {mu_1}")
        # the magnitude of direct effects is capped on +/-1
        temp_val = geno_vec_mutant[update_idx]+mu_1
        #print(f"Temp val before assignment {temp_val}")
        if temp_val>=-1 and temp_val<=1:
            geno_vec_mutant[update_idx] = geno_vec_mutant[update_idx]+mu_1
        
        return geno_vec_mutant

    def mutate_B(self):
        """Mutation on B matrix - sample random indices i,j and add
            mu_2 parameter with probability 0.0067"""

        mut_inter_matrix = self.interaction_matrix.copy()
        idx_i = random.randint(0, self.vec_size - 1)
        idx_j = random.randint(0, self.vec_size - 1)
        mu_2 = np.random.uniform(-0.0067, 0.0067)

        mut_inter_matrix[idx_i,idx_j] = mut_inter_matrix[idx_i, idx_j]+mu_2

        return mut_inter_matrix


    def mutate_B_destruct(self, col_start=0, col_finish=0, row_start=0, row_finish=0, single_cell=0):
        """Mutation on B matrix - sample random indices i,j and add
            mu_2 parameter with probability 0.0067
            Dictate a part of the matrix that won't be mutated."""
        mut_inter_matrix = self.interaction_matrix.copy()
        if single_cell==0:
            masked_i = col_start
            masked_j = row_start
            idx_i = random.randint(0, self.vec_size - 1)
            idx_j = random.randint(0, self.vec_size - 1)
            if masked_i == idx_i and masked_j==idx_j:
                return mut_inter_matrix
            else:
                mu_2 = np.random.uniform(-0.0067, 0.0067)
                mut_inter_matrix[idx_i,idx_j] = mut_inter_matrix[idx_i, idx_j]+mu_2
                return mut_inter_matrix
                
        elif single_cell==1:
            masked_i_list = np.arange(row_start, row_finish)
            masked_j_list = np.arange(col_start, col_finish)
            idx_i = random.randint(0, self.vec_size - 1)
            idx_j = random.randint(0, self.vec_size - 1)
            if (idx_i in masked_i_list) and (idx_j in masked_j_list):
                return mut_inter_matrix
            else:
                mu_2 = np.random.uniform(-0.0067, 0.0067)
                mut_inter_matrix[idx_i,idx_j] = mut_inter_matrix[idx_i, idx_j]+mu_2
                return mut_inter_matrix

    def calculate_fitness(self, pheno_vec, target_vec):
        fitness = 1 + np.dot(pheno_vec, target_vec)

        return fitness
    

    def update_phenotype(self, t1_rate=1, t2_rate=0.2):
        self.pheno_vec = self.pheno_vec + t1_rate*np.tanh(np.dot(self.interaction_matrix, self.pheno_vec) -
                                                          t2_rate*self.pheno_vec)

    def update_phenotype2(self, genotype_vector, interaction_matrix, t1_rate=1, t2_rate=0.2):
        phenotype_vec = genotype_vector
        for i in range(self.t_dev):
            phenotype_vec = phenotype_vec + t1_rate*np.tanh(np.dot(interaction_matrix, phenotype_vec) -
                                                          t2_rate*phenotype_vec)
        
        return phenotype_vec

    def update_phenotype3(self, genotype_vector, interaction_matrix, t1_rate=1, t2_rate=0.2):
        """To generate developmental time step plots"""
        developmental_arr = np.zeros((10,8))
        phenotype_vec = genotype_vector
        for i in range(self.t_dev):
            phenotype_vec = phenotype_vec + t1_rate*np.tanh(np.dot(interaction_matrix, phenotype_vec) -
                                                          t2_rate*phenotype_vec)
            developmental_arr[i,:]=phenotype_vec
        
        return phenotype_vec, developmental_arr
                                                                                                        


n = 8
t_dev = 10
t_evo = 20

demotarget = np.array([1, 1,-1,-1,-1, 1,-1, 1])
demoGRN = GPMappingGRN(n, t_dev, 1)

# placeholder for interaction matrix plot shape NxNxT_evo
plot_placeholder = np.zeros((t_evo,n,n))

# flags for extraction of dev arrays
flags_list = [2000, 4000, 5000, 6000]
dev_mats_list = []
# evolution loop
for j in range(t_evo):
    # run one evo loop on the genotype
    if j in flags_list:
        p_star, temp_dev_arr = demoGRN.update_phenotype3(demoGRN.geno_vec, demoGRN.interaction_matrix)
        p_star_fitness = demoGRN.calculate_fitness(p_star, demotarget)
        dev_mats_list.append(temp_dev_arr)
    else:
        p_star = demoGRN.update_phenotype2(demoGRN.geno_vec, demoGRN.interaction_matrix)
        p_star_fitness = demoGRN.calculate_fitness(p_star, demotarget)
    # print(f"Step {i}")
    if j%10==0:
        print(f"Fitness at at evo {j} - ",p_star_fitness)
        print(f"Phenotype vector at evo {j}", p_star)
    # mutate G and B
    mut_geno_vec = demoGRN.mutate_G()
    if random.uniform(0, 1) <= 0.067:
        mut_B_matrix = demoGRN.mutate_B()
    else:
        mut_B_matrix = demoGRN.interaction_matrix
    
    # run evo loop on mutant genotype
    p_star_mutant = demoGRN.update_phenotype2(mut_geno_vec, mut_B_matrix)
    p_star_mutant_ft = demoGRN.calculate_fitness(p_star_mutant, demotarget)
    if j%10==0:
        print(f"Mutant phenotype vector at evo {j}", p_star_mutant)
    if p_star_mutant_ft>p_star_fitness:
        demoGRN.geno_vec = mut_geno_vec
        demoGRN.interaction_matrix = mut_B_matrix
    
    # save interaction matrix instance
    plot_placeholder[j,:,:] = demoGRN.interaction_matrix

# save interaction matrices to system 
#np.save(os.path.join(save_location, f"single_env_{n}_{t_evo}_{time_stamp}"),plot_placeholder)
#np.save(os.path.join(save_location, f"final_B_{n}_{t_evo}_{time_stamp}"),demoGRN.interaction_matrix)
