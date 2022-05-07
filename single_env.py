import numpy as np
import random


n = 8
t_dev = 10
t_evo = 2000


class GPSingleEnvMappingGRN():

    def __init__(self, indiv_size, dev_steps, evo_steps, s_target, pop_size=100):
        self.vec_size = indiv_size
        self.geno_vec = np.zeros(indiv_size)
        self.interaction_matrix = np.zeros((indiv_size, indiv_size))
        # self.pheno_vec = np.zeros(indiv_size)
        self.t_dev = dev_steps
        self.t_evo = evo_steps
        self.target_vec = s_target
        self.population_size = pop_size

    # def generate_population(self):
    #     """Generate a population of 'pop_size' random strings."""
    #     pop_size = self.population_size
    #     pop_list = []
    #     for i in range(pop_size):
    #         pop_list.append(np.zeros(self.vec_size))

    #     return pop_list

    def mutate_G(self):
        # select a random index uniformly
        update_idx = random.randint(0, self.vec_size - 1)
        # sample mu_1 from (-0.1,0.1)
        mu_1 = np.random.uniform(-0.1, 0.1)
        self.geno_vec[update_idx] = self.geno_vec[update_idx]+mu_1

    def mutate_B(self):
        """Mutation on B matrix - sample random indices i,j and add
            mu_2 parameter with probability 0.0067"""
        idx_i = random.randint(0, self.vec_size - 1)
        idx_j = random.randint(0, self.vec_size - 1)
        mu_2 = np.random.uniform(-0.0067, 0.0067)

        self.interaction_matrix[idx_i,
            idx_j] = self.interaction_matrix[idx_i, idx_j]+mu_2

    def calculate_fitness(self):
        fitness = 1 + np.dot(self.pheno_vec, self.target_vec)

        return fitness

    def update_phenotype(self, t1_rate=1, t2_rate=0.2):
        self.pheno_vec = self.pheno_vec + t1_rate*np.tanh(np.dot(self.interaction_matrix, self.pheno_vec) -
                                                          t2_rate*self.pheno_vec)

    def evo_loop(self):
        # run a sequence of developmental steps creating one developmental loop

        pass


demotarget = np.array([1, -1, -1, 1])
demoGRN = GPSingleEnvMappingGRN(4, 10, 1, demotarget)

population = demoGRN.generate_population()
print(population)
# fitness_init = demoGRN.calculate_fitness()
# # evolution loop
# for j in range(10):
#     # developmental loop
#     for i in range(10):
#         demoGRN.update_phenotype()
#         #temp_fitness = demoGRN.calculate_fitness()
#         # print(f"Step {i}")
#         # print(f"Fitness at {i} - ",temp_fitness)
#         print(f"Phenotype vector at dev{i} - evo {j}", demoGRN.pheno_vec)
#     demoGRN.mutate_G()
#     demoGRN.pheno_vec = demoGRN.geno_vec
#     if random.uniform(0, 1) <= 0.067:
#         demoGRN.mutate_B()
