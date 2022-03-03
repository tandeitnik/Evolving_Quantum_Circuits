# -*- coding: utf-8 -*-

import numpy as np
import stabilizer_functions as sf
#import qecc_functions as qf
import evolution_functions as ef
import random
#import networkx as nx
import matplotlib.pyplot as plt
import itertools

### LIST OF FUNTIONS ###
#

# group_generator(generators,signals) = group, signal
# group2gen(group,N_gen) = list_gen
# pauli_mult(operator_1,operator_2) = ans, signal
# pauli_multi_list(operator_1,operator_2) = ans, signal
# pauli_decoder(pauli_codeword) = literal
# pauli_decoder_list(pauli_bit) = pauli_list
# find_literal(literal,array) = positions
# update_msg(fit_mean,mu_mean,D_mean,generated_circuits,running_time,fit_register,mu_register,D_register)
# gen_zero_toy = circ_set, fit_register, mu_register, D_register
# suplementary_variables(M,fit_register,mu_register,D_register) = fit_mean, fit_std, mu_mean, mu_std, D_mean, D_std, running_time, generated_circuits
# donut_adj_mat(a,b) = adj_mat
# ring_adj_mat(N) = adj_mat
# square_adj_mat(N,M) = adj_mat
# circ_adj_mat(N,circ) = adj_mat
# global_adj_mat(N) = adj_mat
# draw_graph(adj_mat,file_name,label,close) = 
# draw_circ_graph(N,circ,file_name,label,close)
# choose_qubits(adj_mat) = qubit_A, qubit_B
# arb_adj_mat(N,neighbors) = adj_mat
# pauli_commutation(pauli_1,pauli_2) = (-1)**n
# rank_mod2(matrix) = rank


def group_generator(generators,signals):
    
    """
    Input: a list of generators and their oveall signal in the literal form
    Output: a lisf of all combinations of operators and the signal
    """
    
    N = len(generators[0])
    
    lst = list(itertools.product([0, 1], repeat=N))
    
    group = []
    signal = []
    
    for i in range(1,len(lst)):
        
        if sum(lst[i]) == 1:
            
            for j in range(N):
                
                if lst[i][j] == 1:
                    
                    group.append(generators[j])
                    signal.append(1)
                    break
                    
        else:
            
            gen_idx = np.nonzero(lst[i])[0]
            
            for j in range(len(gen_idx) - 1):
                
                if j == 0:
                    
                    ans_temp ,signal_temp  = pauli_multi_list(generators[gen_idx[0]],generators[gen_idx[1]])
                    
                else:
                    
                    ans_temp ,signal_temp_2  = pauli_multi_list(ans_temp,generators[gen_idx[1+j]])
                    signal_temp *= signal_temp_2
                    
            group.append(ans_temp)
            signal.append(signal_temp)
                
    return group, np.real(signal)
                

def group2gen(group,N_gen):
    
    """
    given a list of operators (in literal format), and how many generators
    should be found (N_gen), the function sorts all combinations and returns
    a list with sets of N_gen LI generators
    
    note that this function can be extremelly slow since the number of combinations
    grows exponentially
    """
    
    N = len(group[0])
    comb_gen = list(itertools.combinations(range(len(group)), N_gen))
    list_gen = []
    
    for i in range(len(comb_gen)):
        
        print((i/len(comb_gen))*100)
        
        matrix = np.zeros([N_gen,2*N])
        
        for j in range(N_gen):
            
            for k in range(N):
                
                if group[comb_gen[i][j]][k] == 'X':
                
                    matrix[j,2*k] = 1
                    
                if group[comb_gen[i][j]][k] == 'Z':
                
                    matrix[j,2*k+1] = 1
                    
                if group[comb_gen[i][j]][k] == 'Y':
                    
                    matrix[j,2*k] = 1
                    matrix[j,2*k+1] = 1
        
        if rank_mod2(matrix) == N_gen:
        
            gen_temp = []
            
            for j in range(N_gen):
                
                gen_temp.append(group[comb_gen[i][j]])
                
            list_gen.append(gen_temp)
            
    return list_gen
        

def pauli_mult(operator_1,operator_2):
    
    """
    this function return the product of operator_1*operator_2 and the global phase
    
    0 -> 1
    1 -> i
    2 -> -1
    3 -> -i
    
    00 -> I
    01 -> Z
    10 -> X
    11 -> Y
    """
    
    op_1_lit = pauli_decoder(operator_1)
    op_2_lit = pauli_decoder(operator_2)
        
    if op_1_lit == 'I' or op_2_lit == 'I' or op_1_lit == op_2_lit:
        signal = 0
    elif (op_1_lit == 'X' and op_2_lit == 'Y') or (op_1_lit == 'Y' and op_2_lit == 'Z') or (op_1_lit == 'Z' and op_2_lit == 'X'):
        signal = 1
    elif (op_1_lit == 'X' and op_2_lit == 'Z') or (op_1_lit == 'Y' and op_2_lit == 'X') or (op_1_lit == 'Z' and op_2_lit == 'Y'):
        signal = 3
    
    
    ans = np.mod(operator_1 + operator_2,2)
    
    return ans, signal

def pauli_multi_list(operator_1,operator_2):
    
    """
    Operator 1 and 2 are expected to be in literal form as ['Z',...]
    """
    
    ans = []
    signal = 1
    N = len(operator_1)
    
    for i in range(N):
        
        if (operator_1[i] == 'I' and operator_2[i] == 'I') or (operator_1[i] == operator_2[i]):
            
            ans.append('I')
            
        elif operator_1[i] == 'I' and operator_2[i] != 'I':
            
            ans.append(operator_2[i])
            
        elif operator_1[i] != 'I' and operator_2[i] == 'I':
        
             ans.append(operator_1[i])
             
        elif operator_1[i] == 'X' and operator_2[i] == 'Y':
            
            ans.append('Z')
            signal *= 1j
            
        elif operator_1[i] == 'X' and operator_2[i] == 'Z':
            
            ans.append('Y')
            signal *= -1j
            
        elif operator_1[i] == 'Y' and operator_2[i] == 'Z':
            
            ans.append('X')
            signal *= 1j
            
        elif operator_1[i] == 'Y' and operator_2[i] == 'X':
            
            ans.append('Z')
            signal *= -1j
        
        elif operator_1[i] == 'Z' and operator_2[i] == 'X':
            
            ans.append('Y')
            signal *= 1j
            
        elif operator_1[i] == 'Z' and operator_2[i] == 'Y':
            
            ans.append('X')
            signal *= -1j
            
    return ans,signal
        

def pauli_decoder(pauli_codeword):
    
    if pauli_codeword[0] == 0 and pauli_codeword[1] == 0:
        literal = 'I'
    if pauli_codeword[0] == 1 and pauli_codeword[1] == 0:
        literal = 'X'
    if pauli_codeword[0] == 1 and pauli_codeword[1] == 1:
        literal = 'Y'
    if pauli_codeword[0] == 0 and pauli_codeword[1] == 1:
        literal = 'Z'
        
    return literal

def pauli_decoder_list(pauli_bit):
    
    pauli_list = []
    
    for i in range(int(len(pauli_bit)/2)):
    
        if pauli_bit[2*i] == 0 and pauli_bit[2*i+1] == 0:
            literal = 'I'
        if pauli_bit[2*i] == 1 and pauli_bit[2*i+1] == 0:
            literal = 'X'
        if pauli_bit[2*i] == 1 and pauli_bit[2*i+1] == 1:
            literal = 'Y'
        if pauli_bit[2*i] == 0 and pauli_bit[2*i+1] == 1:
            literal = 'Z'
            
        pauli_list.append(literal)
        
    return pauli_list

def find_literal(literal,array):
    
    positions = np.where(array == literal)[0]
        
    return positions

def update_msg(fit_mean,mu_mean,D_mean,generated_circuits,running_time,fit_register,mu_register,D_register):
    
    best_cicuit_pos = np.argsort(-1*fit_register)[0]
    
    print("mean fitness = "+str(fit_mean[-1]))
    print("mean entropy = "+str(mu_mean[-1]))
    print("mean depth = "+str(D_mean[-1]))
    print("Best circuit has fitness = "+str(fit_register[best_cicuit_pos])+", entropy = "+str(mu_register[best_cicuit_pos])+" and depth = "+str(D_register[best_cicuit_pos]))
    print("generated circuits = "+str(generated_circuits[-1]))
    print("the program has been running for "+str(running_time[-1]/3600.)+" hours.")
    print("#####################################")
    print("#####################################")
    return

def gen_zero_toy(M,N,T,adj_mat,possible_gates):

    #initializing generation 0
    circ_set= [] #stores the current competing circuits
    fit_register = np.zeros(M) #stores the fitness of each circuit
    mu_register = np.zeros(M)  #stores the entropy of each circuit
    D_register = np.zeros(M)   #stores the depth of each circuit
    
    #creating generation 0
    for i in range(M):
        circ = sf.randcirc(N,T,adj_mat,possible_gates)
        circ_set.append(circ)
        fit_register[i], mu_register[i], D_register[i] = ef.fitness_toy(N, circ)
    
    return circ_set, fit_register, mu_register, D_register

def suplementary_variables(M,fit_register,mu_register,D_register):

    fit_mean = np.array([np.mean(fit_register)]) #stores the mean fitness @ multiples of acq_time
    fit_std = np.array([np.std(fit_register)]) #stores the std fitness @ multiples of acq_time
    mu_mean = np.array([np.mean(mu_register)]) #stores the mean entropy @ multiples of acq_time
    mu_std = np.array([np.std(mu_register)]) #stores the std entropy @ multiples of acq_time
    D_mean = np.array([np.mean(D_register)]) #stores the mean depth @ multiples of acq_time
    D_std = np.array([np.std(D_register)]) #stores the std depth @ multiples of acq_time
    running_time = np.array([0])  #stores the elapsed time @ multiples of acq_time
    generated_circuits = np.array([M]) #stores how many unique circuits were generated (including dead circuits)
    
    return fit_mean, fit_std, mu_mean, mu_std, D_mean, D_std, running_time, generated_circuits

def donut_adj_mat(a,b):

    N = a*b
    
    adj_mat= np.zeros([N,N])
    
    for i in range(N):
        #linking to left and right neighbours
        if (i % b) == 0:
            adj_mat[i,i+b-1] = 1
            adj_mat[i+b-1,i] = 1
            adj_mat[i,i+1] = 1
            adj_mat[i+1,i] = 1
        elif (i % b) == b-1:
            adj_mat[i,i-1] = 1
            adj_mat[i-1,i] = 1
        else:
            adj_mat[i,i+1] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i,i-1] = 1
            adj_mat[i-1,i] = 1
            
        #linking to upper and lower neighbours
        if i < b:
            adj_mat[i,i+b*(a-1)] = 1
            adj_mat[i+b*(a-1),i] = 1
        else:
            adj_mat[i,i-b] = 1
            adj_mat[i-b,i] = 1
    return adj_mat.astype(int)

def ring_adj_mat(N):

    adj_mat = np.zeros([N,N])
    
    for i in range(N-1):
        adj_mat[i,np.mod(i+1,N)] = 1
        adj_mat[i,np.mod(i-1,N)] = 1
        adj_mat[np.mod(i+1,N),i] = 1
        adj_mat[np.mod(i-1,N),i] = 1
    
    return adj_mat.astype(int)

def square_adj_mat(N,M):

    #adj matrix for a square lattice of NxM qubits    

    L = N*M
    
    adj_mat = np.zeros([L,L])
    
    for i in range(L):
        if i == 0:
            adj_mat[i,i+1] = 1
            adj_mat[i,i+M] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i+M,i] = 1
        elif i == M-1:
            adj_mat[i,i-1] = 1
            adj_mat[i,i+M] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i+M,i] = 1
        elif i == M*(N-1):
            adj_mat[i,i+1] = 1
            adj_mat[i,i-M] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i-M,i] = 1
        elif i == L -1:
            adj_mat[i,i-1] = 1
            adj_mat[i,i-M] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i-M,i] = 1
        elif (i % M) == 0 and i != M*(N-1) and i != 0:
            adj_mat[i,i-M] = 1
            adj_mat[i,i+1] = 1
            adj_mat[i,i+M] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i+M,i] = 1
            adj_mat[i-M,i] = 1
        elif (i % M) == M -1 and i != M-1 and i != L-1:
            adj_mat[i,i-M] = 1
            adj_mat[i,i-1] = 1
            adj_mat[i,i+M] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i+M,i] = 1
            adj_mat[i-M,i] = 1
        elif i < M-1 and i != 0:
            adj_mat[i,i+M] = 1
            adj_mat[i,i+1] = 1
            adj_mat[i,i-1] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i+M,i] = 1
        elif i > M*(N-1) and i != L-1:
            adj_mat[i,i-M] = 1
            adj_mat[i,i+1] = 1
            adj_mat[i,i-1] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i-M,i] = 1
        else:
            adj_mat[i,i+M] = 1
            adj_mat[i,i-M] = 1
            adj_mat[i,i+1] = 1
            adj_mat[i,i-1] = 1
            adj_mat[i-1,i] = 1
            adj_mat[i+1,i] = 1
            adj_mat[i-M,i] = 1
            adj_mat[i+M,i] = 1
            
    return adj_mat.astype(int)

def circ_adj_mat(N,circ):

    circuit = np.copy(circ)
    
    neighbor = []
    
    for i in range(len(circuit)):
        if circuit[i,0] == 3:
            neighbor_a = int(min(circuit[i,1:]))
            neighbor_b = int(max(circuit[i,1:]))
            neighbor.append([neighbor_a,neighbor_b])
            
    neighbor_unique = []
    
    for i in neighbor:
        if i not in neighbor_unique:
            neighbor_unique.append(i)
            
    adj_mat = np.zeros([N,N])
    
    for i in range(len(neighbor_unique)):
        adj_mat[neighbor_unique[i][0],neighbor_unique[i][1]] = 1
        adj_mat[neighbor_unique[i][1],neighbor_unique[i][0]] = 1
        
    return adj_mat.astype(int)

def global_adj_mat(N):

    adj_mat = np.ones([N,N])
    
    for i in range(N):
        adj_mat[i,i] = 0
        
    return adj_mat.astype(int)

def draw_graph(adj_mat,file_name,label,close):
    Graph = nx.from_numpy_matrix(adj_mat)
    nx.draw(Graph,node_size=300, with_labels=label)
    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    plt.savefig(file_name+".png", format="PNG")
    if close == True:
        plt.close()
    return

def draw_circ_graph(N,circ,file_name,label,close):
    adj_mat = circ_adj_mat(N,circ)
    draw_graph(adj_mat,file_name,label,close)
    return

def choose_qubits(adj_mat):

    candidates = range(len(adj_mat))
    qubit_A = random.sample(candidates,  1)[0]
    
    controle = 0
    
    while controle == 0:
        qubit_B = random.sample(candidates,  1)[0]
        if adj_mat[qubit_A,qubit_B] == 1:
            controle = 1
            
    return qubit_A, qubit_B

def choose_qubit(adj_mat,qubit_A):

    candidates = range(len(adj_mat))
    
    controle = 0
    
    while controle == 0:
        qubit_B = random.sample(candidates,  1)[0]
        if adj_mat[qubit_A,qubit_B] == 1:
            controle = 1
            
    return qubit_B

def arb_adj_mat(N,neighbors):
    #neighbors should be a list of lists of the neighbors, e.g., neighbors = [ [0,1], [5,2], [5,1] ]
    adj_mat = np.zeros([N,N])

    for i in range(len(neighbors)):
        adj_mat[neighbors[i][0],neighbors[i][1]] = 1
        adj_mat[neighbors[i][1],neighbors[i][0]] = 1
        
    return adj_mat.astype(int)

def pauli_commutation(pauli_1,pauli_2):
    
    #Given two n-qubit pauli operators (as n-arrays where the operators are literals)
    #this function evaluates in they commute or anti-commute.
    #If output = 1, commute, if output = -1, anti-commute
    
    n = 0
    
    for i in range(len(pauli_1)):
        if (pauli_1[i] != pauli_2[i]) and (pauli_1[i] != 'I' and pauli_2[i] != 'I'):
            n += 1
            
    return (-1)**n

def rank_mod2(matrix):
    """
    https://stackoverflow.com/questions/56856378/fast-computation-of-matrix-rank-over-gf2
    """
    
    a,b = np.shape(matrix)
    
    rows = []
    
    for i in range(a):
        string = ''
        for j in range(b):
            string += str(int(matrix[i,j]))
        rows.append(int(string, 2))
    
    rank = 0
    while rows:
        pivot_row = rows.pop()
        if pivot_row:
            rank += 1
            lsb = pivot_row & -pivot_row
            for index, row in enumerate(rows):
                if row & lsb:
                    rows[index] = row ^ pivot_row
    return rank

def list_of_errors(N,t):

    errors_literal = []
    affected_qubits = []
    
    for i in range(1,t+1):
        
        combinations = list(itertools.combinations(range(N),i))
        
        e_list = []
            
        for j in range(len(combinations)*3**i):
            e_list.append(['I']*N)
        
        for j in range(len(combinations)):  
    
            for k in range(3**i):
                affected_qubits.append(combinations[j])
    
            for k in range(i):
    
                for l in range(3**i):
                    
                    if (l // 3**k) % 3 == 0:
                    
                        e_list[j*3**i+l][combinations[j][k]] = 'X'
                        
                    if (l // 3**k) % 3 == 1:
                    
                        e_list[j*3**i + l][combinations[j][k]] = 'Y'
                        
                    if (l // 3**k) % 3 == 2:
                    
                        e_list[j*3**i + l][combinations[j][k]] = 'Z'
                    
                    
                    
        for k in range(len(e_list)):
            errors_literal.append(e_list[k])
            
    return errors_literal, affected_qubits