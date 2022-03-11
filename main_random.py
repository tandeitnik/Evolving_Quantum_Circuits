# -*- coding: utf-8 -*-

import stabilizer_functions as sf
import utility_functions as uf
import evolution_functions as ef
import numpy as np

loops = 100

N = 5  #number of qubits
T = 3  #N*T gates will be draw for random individuals
t = 1  #weight of errors to be corrected
adj_mat = uf.global_adj_mat(N) #topology of lattice
max_generations = 5000 #if the max_generations is reached the loop ends
max_fitness_target = 1027.2

########################################################

errors_literal, affected_qubits = uf.list_of_errors(N,t)
stab_group, gp_group = sf.stabilizer_group_computational(N)
stab_group_X, gp_group_X = sf.stabilizer_group_X(N)

for loop in range(loops):
    
    rand_circuit = sf.randcirc(N,T,adj_mat)
    fitness = ef.fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)
    max_fitness  = [fitness]
    
    for i in range(max_generations-1):
        
        print(str(max_fitness[-1]))
        
        rand_circuit = sf.randcirc(N,T,adj_mat)
        fitness = ef.fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)
        
        if fitness > max_fitness[-1]:
            
            max_fitness.append(fitness)
            
        else:
            
            max_fitness.append(max_fitness[-1])
            
        if fitness == max_fitness_target or i == max_generations-1:
            break
        
    np.save("max_fitness_loop_"+str(loop)+".npy",np.array(max_fitness))
        