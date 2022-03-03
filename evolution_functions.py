# -*- coding: utf-8 -*-

import numpy as np
import stabilizer_functions as sf
import qecc_functions as qf
import utility_functions as uf
import random

#######################
## FITNESS FUNCTIONS ##
#######################

def fitness_qecc(N,circuit,t,errors_literal,affected_qubits):

    #in computational basis
    codeword_0 = np.copy(circuit)
    codeword_1 = np.zeros([len(circuit)+1,3])
    codeword_1[0,:] = np.array([4,0,0])
    codeword_1[1:,:] = circuit
    codeword_1 = codeword_1.astype(int)

    codewords = [codeword_0,codeword_1]

    #in diagonal basis
    codeword_0 = np.zeros([len(circuit)+1,3])
    codeword_0[0,:] = np.array([1,0,0])
    codeword_0[1:,:] = circuit
    codeword_0 = codeword_0.astype(int)

    codeword_1 = np.zeros([len(circuit)+2,3])
    codeword_1[0,:] = np.array([4,0,0])
    codeword_1[1,:] = np.array([1,0,0])
    codeword_1[2:,:] = circuit
    codeword_1 = codeword_1.astype(int)

    codewords_diagonal = [codeword_0,codeword_1]
    
    #this part of the code is to speed up
    stab_lit  = []
 
    for i in range(len(codewords)):
        
        M , gp = sf.stabilizer(N, codewords[i])
        stab_group, gp_group, stab_lit_ = sf.stabilizer_group(M,gp)

        stab_lit.append(stab_lit_)
        
    stab_set = [{tuple(stab_lit[j][i,:]) for i in range(len(stab_lit[j]))} for j in range(len(codewords))]
            

    c_d, t_d, css_group_d, css_d = qf.correctability_degree(N,codewords,t,errors_literal,affected_qubits,stab_set)
    f_d = qf.phase_distance(N,codewords_diagonal,stab_set)
    p_d = qf.phase_distance(N,codewords,stab_set)
    #D   = sf.depth(N,circuit)
    
    fitness = (1+c_d)**10  + ( (f_d*p_d) /(np.sqrt(abs(f_d-p_d))+1))/81  + (1+css_d)**10 #+ t_d + 1/D
    
    #fitness = (c_d/D)#*( (f_d*p_d)**2 /(np.sqrt(abs(f_d-p_d))+1))
        
    return fitness


####################################
## CROSSOVER & MUTATION FUNCTIONS ##
####################################


def mutation(circ,N,adj_mat):
    
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
        
        
        
    return new_circ.astype(int)

def crossover(population,selected_idx):

    min_len = min(len(population[selected_idx[0]][0]),len(population[selected_idx[1]][0]))
    
    cut_point = np.random.randint(0,min_len)
    
    
    offspring_A = np.zeros([len(population[selected_idx[1]][0]),3])
    offspring_B = np.zeros([len(population[selected_idx[0]][0]),3])
    
    offspring_A[:cut_point,:] = population[selected_idx[0]][0][:cut_point,:]
    offspring_A[cut_point:,:] = population[selected_idx[1]][0][cut_point:,:]
    
    offspring_B[:cut_point,:] = population[selected_idx[1]][0][:cut_point,:]
    offspring_B[cut_point:,:] = population[selected_idx[0]][0][cut_point:,:]
    
    offspring_A = offspring_A.astype(int)
    offspring_B = offspring_B.astype(int)
    
    return offspring_A, offspring_B


#########################
## SELECTION FUNCTIONS ##
#########################

def refresh_population(population,M,N,T,adj_mat,t,errors_literal,affected_qubits,death_rate):
        
    
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
        fitness = fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits)
    
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


def circle_of_life(progenitors,population,N,adj_mat,mutation_rate,crossover_rate,t,errors_literal,affected_qubits,mutation_density):    

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
    fitness_A = fitness_qecc(N,offspring_A,t,errors_literal,affected_qubits)
    fitness_B = fitness_qecc(N,offspring_B,t,errors_literal,affected_qubits)
    
    population.append([offspring_A,fitness_A])
    population.append([offspring_B,fitness_B])

    return population