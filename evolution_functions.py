# -*- coding: utf-8 -*-

import numpy as np
import stabilizer_functions as sf
import qecc_functions as qf
import utility_functions as uf
import random

#######################
## FITNESS FUNCTIONS ##
#######################

def fitness_qecc(N,circuit,t,errors_literal,affected_qubits,stab_group,gp_group,stab_group_X, gp_group_X):

    #primeiro eu gero o grupo do circuito, vou precisar desse grupo para calcular
    #o phase flip
    stab_group_circuit, gp_group_circuit, stab_lit_circuit = sf.stabilizer_group(N,circuit,stab_group, gp_group)

    #calculo do menor número de literais num stabilizer - faço esse calculo para
    #ter um low bound nos literais 
    literals = np.zeros(2**N-1)

    for k in range(2**N-1):
        
        literals[k] = N - len(uf.find_literal('I',stab_lit_circuit[k,:]))
        if literals[k] == 1:
            min_lit = 1
            break
            
    if k == 2**N-2:
        min_lit = min(literals)
        
    #agora vou descobrir os operadores de bit-flip.

    stab_group_temp_X = np.copy(stab_group_X)
    gp_group_temp_X = np.copy(gp_group_X)

    ###################################################
    #applying the circuit to the logical operators  X##
    ###################################################

    if len(circuit) != 0:

        for i in range(len(circuit)): #running through circuit
        
            if circuit[i,0] == 1: #Hadamard
            
                for j in range(2**N-1):
                    
                    if  stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1:
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)
                        
                    temp_x = stab_group_temp_X[j,2*circuit[i,1]]
                    stab_group_temp_X[j,2*circuit[i,1]] = stab_group_temp_X[j,2*circuit[i,1]+1]
                    stab_group_temp_X[j,2*circuit[i,1]+1] = temp_x
                    
            elif circuit[i,0] == 2: #Phase
                
                for j in range(2**N-1):
                    
                    if  stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1:
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)
                
                    stab_group_temp_X[j,2*circuit[i,1]+1] = np.mod(stab_group_temp_X[j,2*circuit[i,1]]+stab_group_temp_X[j,2*circuit[i,1]+1],2)
            
            elif circuit[i,0] == 3: #CNOT (left)
            
                 for j in range(2**N-1):
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1 and stab_group_temp_X[j,2*circuit[i,2]] == 1 and stab_group_temp_X[j,2*circuit[i,2]+1] == 1) or (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 0 and stab_group_temp_X[j,2*circuit[i,2]] == 0 and stab_group_temp_X[j,2*circuit[i,2]+1] == 1):
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)
                        
            
                    stab_group_temp_X[j,2*circuit[i,2]] = np.mod(stab_group_temp_X[j,2*circuit[i,2]] + stab_group_temp_X[j,2*circuit[i,1]] ,2)
                    stab_group_temp_X[j,2*circuit[i,1]+1] = np.mod(stab_group_temp_X[j,2*circuit[i,1]+1] + stab_group_temp_X[j,2*circuit[i,2]+1],2)
            
            elif circuit[i,0] == 4: #X
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp_X[j,2*circuit[i,1]] == 0 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1) or (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)
                        
            elif circuit[i,0] == 5: #Y
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 0) or (stab_group_temp_X[j,2*circuit[i,1]] == 0 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)
                        
            elif circuit[i,0] == 6: #Z
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 0) or (stab_group_temp_X[j,2*circuit[i,1]] == 1 and stab_group_temp_X[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp_X[j] = np.mod(gp_group_temp_X[j]+2,4)

    stab_lit_comp = [['I' for columns in range(N+1)] for rows in range(len(gp_group_temp_X))]
        
    for i in range(len(gp_group_temp_X)):
        for j in range(N):
            if stab_group_temp_X[i,2*j] == 0 and stab_group_temp_X[i,2*j+1] == 1:
                stab_lit_comp[i][j] = 'Z'
            elif stab_group_temp_X[i,2*j] == 1 and stab_group_temp_X[i,2*j+1] == 0:
                stab_lit_comp[i][j] = 'X'
            elif stab_group_temp_X[i,2*j] == 1 and stab_group_temp_X[i,2*j+1] == 1:
                stab_lit_comp[i][j] = 'Y'
                
        if gp_group_temp_X[i] == 0:
            stab_lit_comp[i][-1] = '1'
        elif gp_group_temp_X[i] == 1:
            stab_lit_comp[i][-1] = 'i'
        elif gp_group_temp_X[i] == 2:
            stab_lit_comp[i][-1] = '-1'
        else:
            stab_lit_comp[i][-1] = '-i'
            
    stab_lit_comp = np.array(stab_lit_comp)

    #o bit flip distance é calculado
    bit_flip_distance = np.array([N - len(uf.find_literal('I',stab_lit_comp[i,:])) for i in range(len(stab_lit_comp))])

    #Agora vem o calculo do phase-flip. A ideia é a seguinte: eu sei que os 
    #operadores de phase-flip são os stabilizers com a fase trocada. Eu tenho
    #a lista dos stabilziers em stab_lit_circuit e eu tenho a lista dos operadores
    #de bit-flip lógico em stab_lit_comp. Agora, seja um dos operadores de bit-flip
    #U que será aplicado ao circuito para transformar ele no estado ortogonal. Para
    #cada S_i do grupo de stabilizer do circuito, será feito US_iU onde o efeito será
    #trocar ou não a fase global já que U é um operador de Pauli geral. Para saber se
    #a fase foi trocada, podemos testar a comutação entre U e cada S_i: se eles comutam
    #a fase não foi trocada, se eles não comutam a fase foi trocada, pois não comutar
    #significa que U e S_i se interceptam em um número ímpar de letras de Pauli, portanto
    #o resultado final de US_iU terá um fator -1 multiplicando.

    #O procedimento então é o seguinte: para cada U_i pertecente a stab_lit_comp, eu
    #testo sua comutatividade para cada S_i pertecente a stab_lit_circuit. Se o
    #não comutar, eu registro o peso de S_i. No final será o menor peso de todos re-
    #gistrados. Já que eu calculei o low bound do peso em min_lit, se eu encontrar
    #um S_i que não comuta com peso igual a min_lit, eu posso parar a busca e já registrar
    #o peso mínimo

    phase_flip_distance = np.ones(2**N-1)*min_lit

    for i in range(len(stab_lit_comp)):
        
        phase_temp = []
        
        for j in range(len(stab_lit_circuit)):
            
            if uf.pauli_commutation(stab_lit_comp[i,:N],stab_lit_circuit[j,:N]) == -1:
                
                phase_temp.append(N - len(uf.find_literal('I',stab_lit_circuit[j,:])))
                
                if phase_temp[-1] == min_lit:
                    break
                
        if j == len(stab_lit_circuit)-1:
            
            phase_flip_distance[i] = min(phase_temp)
            
    f_d = [ 1/(abs(bit_flip_distance[i]-phase_flip_distance[i])+1) for i in range(2**N-1)]        

    D = sf.depth(N,circuit)

    stab_lit_circuit_set = {tuple(stab_lit_circuit[i,:]) for i in range(2**N-1)}

    fitness = []

    for k in range(2**N-1):
        
        if f_d[k] == max(f_d):

            stab_group_ort = np.copy(stab_group_circuit)
            gp_group_ort = np.copy(gp_group_circuit)
            
            bit_flip_operator = stab_lit_comp[k,:N]
            
            circuit_ort = np.zeros([N - len(uf.find_literal('I',bit_flip_operator)),3])

            count = 0
            
            for i in range(N):
                
                if bit_flip_operator[i] == 'X':
                    circuit_ort[count,:] = np.array([4,i,0])
                    count += 1
                    
                elif bit_flip_operator[i] == 'Y':
                    circuit_ort[count,:] = np.array([5,i,0])
                    count += 1
                    
                elif bit_flip_operator[i] == 'Z':
                    circuit_ort[count,:] = np.array([6,i,0])
                    count += 1
            
            circuit_ort = circuit_ort.astype(int)
            
            ###################################################
            #applying the bit-flip operator to the circuit   ##
            ###################################################
            
            if len(circuit_ort) != 0:
            
                for i in range(len(circuit_ort)): #running through circuit_ort
                
                    if circuit_ort[i,0] == 1: #Hadamard
                    
                        for j in range(2**N-1):
                            
                            if  stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1:
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
                                
                            temp_x = stab_group_ort[j,2*circuit_ort[i,1]]
                            stab_group_ort[j,2*circuit_ort[i,1]] = stab_group_ort[j,2*circuit_ort[i,1]+1]
                            stab_group_ort[j,2*circuit_ort[i,1]+1] = temp_x
                            
                    elif circuit_ort[i,0] == 2: #Phase
                        
                        for j in range(2**N-1):
                            
                            if  stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1:
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
                        
                            stab_group_ort[j,2*circuit_ort[i,1]+1] = np.mod(stab_group_ort[j,2*circuit_ort[i,1]]+stab_group_ort[j,2*circuit_ort[i,1]+1],2)
                    
                    elif circuit_ort[i,0] == 3: #CNOT (left)
                    
                         for j in range(2**N-1):
                             
                            #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                            if (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1 and stab_group_ort[j,2*circuit_ort[i,2]] == 1 and stab_group_ort[j,2*circuit_ort[i,2]+1] == 1) or (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 0 and stab_group_ort[j,2*circuit_ort[i,2]] == 0 and stab_group_ort[j,2*circuit_ort[i,2]+1] == 1):
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
                                
                    
                            stab_group_ort[j,2*circuit_ort[i,2]] = np.mod(stab_group_ort[j,2*circuit_ort[i,2]] + stab_group_ort[j,2*circuit_ort[i,1]] ,2)
                            stab_group_ort[j,2*circuit_ort[i,1]+1] = np.mod(stab_group_ort[j,2*circuit_ort[i,1]+1] + stab_group_ort[j,2*circuit_ort[i,2]+1],2)
                    
                    elif circuit_ort[i,0] == 4: #X
                    
                        for j in range(2**N-1):
                            
                             if (stab_group_ort[j,2*circuit_ort[i,1]] == 0 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1) or (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1):
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
                                
                    elif circuit_ort[i,0] == 5: #Y
                    
                        for j in range(2**N-1):
                            
                             if (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 0) or (stab_group_ort[j,2*circuit_ort[i,1]] == 0 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1):
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
                                
                    elif circuit_ort[i,0] == 6: #Z
                    
                        for j in range(2**N-1):
                            
                             if (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 0) or (stab_group_ort[j,2*circuit_ort[i,1]] == 1 and stab_group_ort[j,2*circuit_ort[i,1]+1] == 1):
                                gp_group_ort[j] = np.mod(gp_group_ort[j]+2,4)
            
            stab_lit_ort = [['I' for columns in range(N+1)] for rows in range(2**N-1)]
                
            for i in range(2**N-1):
                for j in range(N):
                    if stab_group_ort[i,2*j] == 0 and stab_group_ort[i,2*j+1] == 1:
                        stab_lit_ort[i][j] = 'Z'
                    elif stab_group_ort[i,2*j] == 1 and stab_group_ort[i,2*j+1] == 0:
                        stab_lit_ort[i][j] = 'X'
                    elif stab_group_ort[i,2*j] == 1 and stab_group_ort[i,2*j+1] == 1:
                        stab_lit_ort[i][j] = 'Y'
                        
                if gp_group_ort[i] == 0:
                    stab_lit_ort[i][-1] = '1'
                elif gp_group_ort[i] == 1:
                    stab_lit_ort[i][-1] = 'i'
                elif gp_group_ort[i] == 2:
                    stab_lit_ort[i][-1] = '-1'
                else:
                    stab_lit_ort[i][-1] = '-i'
                    
            stab_lit_ort = np.array(stab_lit_ort)

            stab_lit_ort_set = {tuple(stab_lit_ort[i,:]) for i in range(2**N-1)}

            stab_set = [stab_lit_circuit_set,stab_lit_ort_set]
            
            #c_d, t_d, css_group_d, css_d = qf.correctability_degree(N,t,errors_literal,affected_qubits,stab_set)
            c_d = qf.correctability_degree(N,t,errors_literal,affected_qubits,stab_set)
            #appending the overall fitness score 
            fitness.append((1+c_d)**10+1/D+min(bit_flip_distance[k],phase_flip_distance[k]))
            
            if c_d == 1:
                break
            
    return(max(fitness))


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

    cut_point_A = np.random.randint(0,len(population[selected_idx[0]][0]))
    cut_point_B = np.random.randint(0,len(population[selected_idx[1]][0]))
    
    offspring_A = np.vstack([ population[selected_idx[0]][0][:cut_point_A,:]  , population[selected_idx[1]][0][cut_point_B:,:]  ])
    offspring_B = np.vstack([ population[selected_idx[1]][0][:cut_point_B,:]  , population[selected_idx[0]][0][cut_point_A:,:]  ])
                 
    
    return offspring_A, offspring_B


#########################
## SELECTION FUNCTIONS ##
#########################

def refresh_population(population,M,N,T,adj_mat,t,errors_literal,affected_qubits,death_rate,stab_group, gp_group,stab_group_X, gp_group_X):
        
    
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
        fitness = fitness_qecc(N,rand_circuit,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)
    
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


def circle_of_life(progenitors,population,N,adj_mat,mutation_rate,crossover_rate,t,errors_literal,affected_qubits,mutation_density,stab_group, gp_group,stab_group_X, gp_group_X):    

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
    fitness_A = fitness_qecc(N,offspring_A,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)
    fitness_B = fitness_qecc(N,offspring_B,t,errors_literal,affected_qubits,stab_group, gp_group,stab_group_X, gp_group_X)
    
    population.append([offspring_A,fitness_A])
    population.append([offspring_B,fitness_B])

    return population
