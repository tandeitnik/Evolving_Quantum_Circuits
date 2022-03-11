# -*- coding: utf-8 -*-

import stabilizer_functions as sf
import utility_functions as uf
import evolution_functions as ef
import numpy as np

loops = 100

M = 10 #initial population size
N = 5  #number of qubits
T = 3  #N*T gates will be draw for random individuals
t = 1  #weight of errors to be corrected
adj_mat = uf.global_adj_mat(N) #topology of lattice
mutation_density = [1,1,1,1,1,1,1,2,2,2,3,3,4,5]
progenitors = 2 #number of progenitors to breed
death_rate = 0 #percentage of population that is replaced by a new random set of individuals
population_surplus = 3 #when the size of population is M*population_surplus, it is reduced to M(1-death_rate) (keeping the best chromossomes) and M*death_rate are added 
acq_time = 10  #at multiples of acq_time data recorded
max_generations = 5000 #if the max_generations is reached the loop ends
max_fitness_target = 1027.2

########################################################

errors_literal, affected_qubits = uf.list_of_errors(N,t)
stab_group, gp_group = sf.stabilizer_group_computational(N)
stab_group_X, gp_group_X = sf.stabilizer_group_X(N)

for loop in range(loops):
    
    population = []
    
    for i in range(M):
        
        rand_circuit = sf.randcirc(N,T,adj_mat)
        fitness = ef.fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)

        individual = [rand_circuit,fitness]
        population.append(individual)
      
    mean_evolution = []
    max_evolution  = []
    std_evolution  = []
      
    for i in range(max_generations):
        
        if (i % acq_time) == 0: #gattering data
            
            fitness_array = np.zeros([len(population),2])
            
            for j in range(len(population)):
                
                fitness_array[j,0] = j
                fitness_array[j,1] = population[j][1]
                
            fitness_array = fitness_array[np.argsort(fitness_array[:, 1])]
            mean_fitness = np.mean(fitness_array[:,1])
            max_fitness  = np.max(fitness_array[:,1])
            std_fitness  = np.std(fitness_array[:,1])
            mean_evolution.append(mean_fitness)
            max_evolution.append(max_fitness)
            std_evolution.append(std_fitness)
            
            mutation_rate   = np.exp(-1*std_fitness)
            crossover_rate  = 1-np.exp(-1*std_fitness)
            
            #printing new data
            # print("Loop = "+str(loop))
            # print("Generation "+str(i))
            # print("Mean fitness = "+str(mean_fitness))
            # print("Best fitness = "+str(max_fitness))
            # print("Std fitness = "+str(std_fitness))
            # print("Mut. rate   = "+str(mutation_rate))
            # print("Cross. rate = "+str(crossover_rate))
            # print("#####################################")
            # print("#####################################")
        
        population = ef.circle_of_life(progenitors,population,N,adj_mat,mutation_rate,crossover_rate,t,errors_literal,affected_qubits,mutation_density,stab_group, gp_group,stab_group_X, gp_group_X)
    
        #refreshing gene pool if population exceeds (1+death_rate)*M
        if len(population) >= population_surplus*M:

            population = ef.refresh_population(population,M,N,T,adj_mat,t,errors_literal,affected_qubits,death_rate,stab_group, gp_group,stab_group_X, gp_group_X)
            
        if max_fitness == max_fitness_target:
            break
            
    np.save("mean_loop_"+str(loop)+".npy",np.array(mean_evolution))
    np.save("max_loop_"+str(loop)+".npy",np.array(max_evolution))
    np.save("std_loop_"+str(loop)+".npy",np.array(std_evolution))