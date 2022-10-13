# -*- coding: utf-8 -*-

import numpy as np
import stabilizer_functions as sf
import qecc_functions as qf
import utility_functions as uf
import random

#######################
## FITNESS FUNCTIONS ##
#######################


def weight(string):
    
    w = 0
    
    for i in range(len(string)):
        
        if string[i] != 'I':
            
            w += 1
            
    return w

#versão alternativa onde 2**N pares são testados
def fitness_qecc(N,circuit,t,errors_literal,affected_qubits,stab_group,gp_group,errors_literal_total,initialOrtCircuits):
    
    D = sf.depth(N,circuit)
        
    common_stabilizer_storage = []
    cdArray = np.zeros(2**N-1)
    
    stab_group_circuit, gp_group_circuit, stab_lit_circuit = sf.stabilizer_group(N,circuit,stab_group, gp_group)
    stab_lit_circuit_set = {tuple(stab_lit_circuit[i,:]) for i in range(2**N-1)}
    
    for k in range(2**N-1):
    #for k in range(10):
    
        #building the second codeword
        circuit1 = initialOrtCircuits[k]
        circuit1 = np.vstack([circuit1, circuit])
        circuit1 = circuit1.astype(int)
        
        #evaluating the stabilizer group of both codewords
        
        stab_group_circuit1, gp_group_circuit1, stab_lit_circuit1 = sf.stabilizer_group(N,circuit1,stab_group, gp_group)
        
        #building the common stabilizer list of the codewords
        
        stab_lit_ort_set = {tuple(stab_lit_circuit1[i,:]) for i in range(2**N-1)}
        stab_set = [stab_lit_circuit_set,stab_lit_ort_set]
        common_stabilizers = stab_lit_circuit_set.intersection(stab_lit_ort_set)
        common_stabilizers = list(common_stabilizers)
        common_stabilizers_set = {common_stabilizers[i][:N] for i in range(len(common_stabilizers))}
        common_stabilizers_list = [list(common_stabilizers[i][:N]) for i in range(len(common_stabilizers))]
        common_stabilizer_storage.append([common_stabilizers_set,common_stabilizers_list])
        #calculo do cd
        cd, t_d, css_group_d, css_d = qf.correctability_degree_color(N,t,errors_literal,affected_qubits,stab_set)
        #calculo do depth
        
        cdArray[k] = cd
    
        if cd == 1:
            break
        
    # #evaluating the distance
    # for k in range(2**N-1):
        
    #     if cdArray[k] == max(cdArray):
            
    #         common_stabilizers_set = common_stabilizer_storage[k][0]
    #         common_stabilizers_list = common_stabilizer_storage[k][1]
            
    #         for i in range(len(errors_literal_total)):
                
    #             if (tuple(errors_literal_total[i]) in common_stabilizers_set) == False:
                    
    #                 commutationTest = 1
                    
    #                 j = 0
                    
    #                 while commutationTest == 1 and j < len(common_stabilizers_list):
                        
    #                     commutationTest = uf.pauli_commutation(errors_literal_total[i],common_stabilizers_list[j])
    #                     j += 1
                        
    #                 if commutationTest == 1: #it is a logical operator
                    
    #                     distance = weight(errors_literal_total[i])
    #                     fitnessArray.append(cdArray[k]+1/D+distance)
    #                     break

    fitness = 1000*max(cdArray) - D
                    
    return fitness

#versão simplificado onde calculo a distância diretamente pelo grupo
# def fitness_qecc(N,circuit,t,errors_literal,affected_qubits,stab_group,gp_group,errors_literal_total):
    
#     #building the second codeword
#     circuit1 = np.array([[4,0,0]])
#     circuit1 = np.vstack([circuit1, circuit])
    
#     #evaluating the stabilizer group of both codewords
#     stab_group_circuit, gp_group_circuit, stab_lit_circuit = sf.stabilizer_group(N,circuit,stab_group, gp_group)
#     stab_group_circuit1, gp_group_circuit1, stab_lit_circuit1 = sf.stabilizer_group(N,circuit1,stab_group, gp_group)
    
#     #building the common stabilizer list of the codewords
#     stab_lit_circuit_set = {tuple(stab_lit_circuit[i,:]) for i in range(2**N-1)}
#     stab_lit_ort_set = {tuple(stab_lit_circuit1[i,:]) for i in range(2**N-1)}
#     stab_set = [stab_lit_circuit_set,stab_lit_ort_set]
#     common_stabilizers = stab_lit_circuit_set.intersection(stab_lit_ort_set)
#     common_stabilizers = list(common_stabilizers)
#     common_stabilizers_set = {common_stabilizers[i][:N] for i in range(len(common_stabilizers))}
#     common_stabilizers_list = [list(common_stabilizers[i][:N]) for i in range(len(common_stabilizers))]
    
#     #evaluating the distance
#     for i in range(len(errors_literal_total)):
        
#         if (tuple(errors_literal_total[i]) in common_stabilizers_set) == False:
            
#             commutationTest = 1
            
#             j = 0
            
#             while commutationTest == 1 and j < len(common_stabilizers_list):
                
#                 commutationTest = uf.pauli_commutation(errors_literal_total[i],common_stabilizers_list[j])
#                 j += 1
                
#             if commutationTest == 1: #it is a logical operator
            
#                 distance = weight(errors_literal_total[i])
#                 break
                
#     #calculo do cd
#     cd = qf.correctability_degree(N,t,errors_literal,affected_qubits,stab_set)
#     #calculo do depth
#     D = sf.depth(N,circuit)
    
#     #calculo final da fitness
    
#     if D == 0:
        
#         fitness = 0
        
#     else:
        
#         #fitness = (1+cd)**2+1/D+distance #standard qecc fitness
#         fitness = cd + 1/D + distance #standard qecc fitness

#     return fitness

def fitness_color(N,circuit,t,errors_literal,affected_qubits,stab_group,gp_group,errors_literal_total,stab_group_X):
    
    fitnessArray = []
    
    for k in range(2**N-1):
    
        #building the second codeword
        circuit1 = np.zeros([N,3])
        for i in range(N):
            circuit1[i,:] = np.array([stab_group_X[k,2*i]*4,i,0])

        circuit1 = np.vstack([circuit1, circuit])
        circuit1 = circuit1.astype(int)
        
        #evaluating the stabilizer group of both codewords
        stab_group_circuit, gp_group_circuit, stab_lit_circuit = sf.stabilizer_group(N,circuit,stab_group, gp_group)
        stab_group_circuit1, gp_group_circuit1, stab_lit_circuit1 = sf.stabilizer_group(N,circuit1,stab_group, gp_group)
        
        #building the common stabilizer list of the codewords
        stab_lit_circuit_set = {tuple(stab_lit_circuit[i,:]) for i in range(2**N-1)}
        stab_lit_ort_set = {tuple(stab_lit_circuit1[i,:]) for i in range(2**N-1)}
        stab_set = [stab_lit_circuit_set,stab_lit_ort_set]
        common_stabilizers = stab_lit_circuit_set.intersection(stab_lit_ort_set)
        common_stabilizers = list(common_stabilizers)
        common_stabilizers_set = {common_stabilizers[i][:N] for i in range(len(common_stabilizers))}
        common_stabilizers_list = [list(common_stabilizers[i][:N]) for i in range(len(common_stabilizers))]
        
        #evaluating the distance
        for i in range(len(errors_literal_total)):
            
            if (tuple(errors_literal_total[i]) in common_stabilizers_set) == False:
                
                commutationTest = 1
                
                j = 0
                
                while commutationTest == 1 and j < len(common_stabilizers_list):
                    
                    commutationTest = uf.pauli_commutation(errors_literal_total[i],common_stabilizers_list[j])
                    j += 1
                    
                if commutationTest == 1: #it is a logical operator
                
                    distance = weight(errors_literal_total[i])
                    break
                    
        #calculo do cd
        cd, t_d, css_group_d, css_d = qf.correctability_degree_color(N,t,errors_literal,affected_qubits,stab_set)
        #calculo do depth
        D = sf.depth(N,circuit)
        
        #calculo final da fitness
        fitnessArray.append((1+cd)**10+(1+css_d)**10+1/D+distance) #color fitness
        
        if cd == 1 and css_d == 1:
            break
        
    return max(fitnessArray)

def fitness_toy(N,circ):
    D = sf.depth(N,circ)
    mu = sf.mean_entropy(N,circ)
    fitness = float(mu)+(1/D)
    return fitness

####################################
## CROSSOVER & MUTATION FUNCTIONS ##
####################################


def mutation(circ,N,adj_mat):
    
    if len(circ) == 0:
        new_gene_test = 1
    else:
        new_gene_test = np.random.randint(0,2) #if 0 an old gene is select to be mutated, if 1 a new gene is inserted at a random point
    
    if new_gene_test == 0:
    
        gene = np.random.randint(0,len(circ)) #selects gene to be mutated
    
        gate = np.random.randint(0,4) #selects new gate
        
        new_circ = np.copy(circ) #copies original circuit
        
        if gate == 3: #if it is a CNOT, choose qubits
        
            qubit_A, qubit_B = uf.choose_qubits(adj_mat)
        
            new_circ[gene,:] = np.array([gate,qubit_A,qubit_B])
        
        elif gate == 0: #the gane is deleted
        
            new_circ = np.delete(new_circ, gene, 0)
        
        else:
            
            qubit_A = np.random.randint(0,N)
            
            new_circ[gene,:] = np.array([gate,qubit_A,0])

    else:
        
        insertion_location = np.random.randint(0,len(circ)+1) #selects insertion location of new gene
    
        gate = np.random.randint(1,4) #selects new gate, identity is ignored
        
        new_circ = np.zeros([len(circ)+1,3])
        
        #making room to new gene and copying the rest
        if insertion_location == 0:
            
            new_circ[1:len(circ)+1,:] = circ
            
        elif insertion_location == len(circ):
            
            new_circ[:len(circ),:] = circ
            
        else:
            
            new_circ[:insertion_location,:] = circ[:insertion_location,:]
            new_circ[insertion_location+1:,:]  = circ[insertion_location:,:]


        if gate == 3: #if it is a CNOT, choose qubits
    
            qubit_A, qubit_B = uf.choose_qubits(adj_mat)
        
            new_circ[insertion_location,:] = np.array([gate,qubit_A,qubit_B])
            
        else:
            
            qubit_A = np.random.randint(0,N)
            
            new_circ[insertion_location,:] = np.array([gate,qubit_A,0])
    
    
    #if the fitness function can spawn zero fitness-valued individuals, the above code may fail, i.e.,
    #it generates an empty selected_idx list. If this happens, random individuals are sampled from the
    #population.
    # if len(selected_idx) == 0:

    #     selected_idx = random.sample(range(len(population)),n_ind)                
        
    return new_circ.astype(int)

def crossover(population,selected_idx):
    
    if len(population[selected_idx[0]][0]) == 0 or len(population[selected_idx[1]][0]) == 0:
        
        offspring_A = population[selected_idx[0]][0]
        offspring_B = population[selected_idx[1]][0]
                       
    else:
        
        cut_point_A = np.random.randint(0,len(population[selected_idx[0]][0]))
        cut_point_B = np.random.randint(0,len(population[selected_idx[1]][0]))
        
        offspring_A = np.vstack([ population[selected_idx[0]][0][:cut_point_A,:]  , population[selected_idx[1]][0][cut_point_B:,:]  ])
        offspring_B = np.vstack([ population[selected_idx[1]][0][:cut_point_B,:]  , population[selected_idx[0]][0][cut_point_A:,:]  ])
                     
    
    return offspring_A, offspring_B


#########################
## SELECTION FUNCTIONS ##
#########################

def refresh_population(population,M,N,T,adj_mat,t,errors_literal,affected_qubits,death_rate,stab_group, gp_group):
        
    
    fitness_array = np.zeros([len(population),2])
                
    for i in range(len(population)):
        
        fitness_array[i,0] = i
        fitness_array[i,1] = population[i][1]
        
    fitness_array = fitness_array[np.argsort(fitness_array[:, 1])]
    
    n_del = int(len(population)-M+M*death_rate)
    
    population.sort(key=lambda x: x[1])
    
    del population[0:n_del]
    
    for i in range(int(M*death_rate)):
        
        rand_circuit = sf.randcirc(N,T,adj_mat)
        fitness = fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits,stab_group, gp_group)
    
        individual = [rand_circuit,fitness]
        population.append(individual)
        
    return population

def ind_selection(n_ind,population):

    selected_idx = []
    
    for j in range(n_ind):

        if j == 0:
    
            fitness_array = np.zeros([len(population),2])
            
            for i in range(len(population)):
                
                fitness_array[i,0] = i
                fitness_array[i,1] = population[i][1]
                
            sum_fitness = np.sum(fitness_array[:,1])
            fitness_array = fitness_array[np.argsort(fitness_array[:, 1])]
            fitness_array_cumul = np.copy(fitness_array)
            
            for i in range(1,len(population)):
                
                fitness_array_cumul[i,1] = np.sum(fitness_array[:i+1,1])
            
            draw = random.uniform(0.0,sum_fitness)
            
            if draw < fitness_array_cumul[0,1]:
                
                selected_idx.append(int(fitness_array_cumul[0,0]))
                fitness_array = np.delete(fitness_array, 0, 0)
            else:
                
                for i in range(len(population)):
                        
                    if fitness_array_cumul[-1-i,1] > draw and fitness_array_cumul[-1-i-1,1] <= draw:
                        
                        selected_idx.append(int(fitness_array_cumul[-1-i,0]))
                        fitness_array = np.delete(fitness_array, -1-i, 0)
                        break
                
        else:
            
            sum_fitness = np.sum(fitness_array[:,1])
            
            fitness_array_cumul = np.copy(fitness_array)
            
            for i in range(1,len(population)-j):
                
                fitness_array_cumul[i,1] = np.sum(fitness_array[:i+1,1])
            
            draw = random.uniform(0.0,sum_fitness)
            
            if draw < fitness_array_cumul[0,1]:
                
                selected_idx.append(int(fitness_array_cumul[0,0]))
                fitness_array = np.delete(fitness_array, 0, 0)
            
            else:
                
                for i in range(len(population)-j):
                        
                    if fitness_array_cumul[-1-i,1] > draw and fitness_array_cumul[-1-i-1,1] <= draw:
                        
                        selected_idx.append(int(fitness_array_cumul[-1-i,0]))
                        fitness_array = np.delete(fitness_array, -1-i, 0)
                        break
                
    return selected_idx

###################
## CoL FUNCTIONS ##
###################


def circle_of_life(progenitors,population,N,adj_mat,mutation_rate,crossover_rate,t,errors_literal,affected_qubits,mutation_density,stab_group, gp_group,errors_literal_total):    

    #select progenitors
    selected_idx = ind_selection(progenitors,population)
    
    #breed
    crossover_test = random.uniform(0.0,1)
    if crossover_test <= crossover_rate:
    
        offspring_A, offspring_B = crossover(population,selected_idx)
    
    else:
        
        offspring_A = np.copy(population[selected_idx[0]][0])
        offspring_B = np.copy(population[selected_idx[1]][0])
    
    #mutate
    mutation_test_A = random.uniform(0.0,1)
    mutation_test_B = random.uniform(0.0,1)
    
    if mutation_test_A <= mutation_rate:
        
        n_mut = random.sample(mutation_density, 1)[0]
        
        for i in range(n_mut):
        
            offspring_A = mutation(offspring_A,N,adj_mat)
        
    if mutation_test_B <= mutation_rate:
        
        n_mut = random.sample(mutation_density, 1)[0]
        
        for i in range(n_mut):
        
            offspring_B = mutation(offspring_B,N,adj_mat)
        
    #calculate fitness and add offspring to population
    fitness_A = fitness_qecc(N,offspring_A,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total)
    fitness_B = fitness_qecc(N,offspring_B,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total)
    
    population.append([offspring_A,fitness_A])
    population.append([offspring_B,fitness_B])

    return population

def circle_of_life_color(progenitors,population,N,adj_mat,mutation_rate,crossover_rate,t,errors_literal,affected_qubits,mutation_density,stab_group, gp_group,errors_literal_total,initialOrtCircuits):    

    #select progenitors
    selected_idx = ind_selection(progenitors,population)
    
    #breed
    crossover_test = random.uniform(0.0,1)
    if crossover_test <= crossover_rate:
    
        offspring_A, offspring_B = crossover(population,selected_idx)
    
    else:
        
        offspring_A = np.copy(population[selected_idx[0]][0])
        offspring_B = np.copy(population[selected_idx[1]][0])
    
    #mutate
    mutation_test_A = random.uniform(0.0,1)
    mutation_test_B = random.uniform(0.0,1)
    
    if mutation_test_A <= mutation_rate:
        
        n_mut = random.sample(mutation_density, 1)[0]
        
        for i in range(n_mut):
        
            offspring_A = mutation(offspring_A,N,adj_mat)
        
    if mutation_test_B <= mutation_rate:
        
        n_mut = random.sample(mutation_density, 1)[0]
        
        for i in range(n_mut):
        
            offspring_B = mutation(offspring_B,N,adj_mat)
        
    #calculate fitness and add offspring to population
    #fitness_A = fitness_color(N,offspring_A,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total,initialOrtCircuits)
    #fitness_B = fitness_color(N,offspring_B,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total,initialOrtCircuits)
    fitness_A = fitness_qecc(N,offspring_A,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total,initialOrtCircuits)
    fitness_B = fitness_qecc(N,offspring_B,t,errors_literal,affected_qubits,stab_group, gp_group,errors_literal_total,initialOrtCircuits)
    
    
    
    population.append([offspring_A,fitness_A])
    population.append([offspring_B,fitness_B])

    return population
