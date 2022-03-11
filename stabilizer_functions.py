# -*- coding: utf-8 -*-

import numpy as np
import random
#from qiskit import QuantumCircuit, assemble, Aer
#from qiskit.tools.visualization import circuit_drawer
import copy
from math import ceil
import itertools
#import qecc_functions as qf
#import evolution_functions as ef
import utility_functions as uf

### LIST OF FUNTIONS ###
#
# randcirc(N,T,adj_mat) = circuit
# gene_express(N,circ) = circuit
# final_state(N,circ) = out_state
# stabilizer(N,circuit) = stabilizer, global_phase
# measurement(M,gp,operator,qubit) = M_post_measurement, gp_post_measurement
# stabilizer_literal(stab_bit, global_phase) = stab_lit
# stabilizer_bit(stab_lit) = M, gp
# stabilizer_group(M,gp) = stab_group, gp_group, stab_lit
# stabilizer_canonical(stabilizer,global_phase) = stab_can, gp_can
# basisnormcirc(M, gp) = C, M_c_l
# check_orthogonality(circ_1,circ_2,N) = orthogonality
# inner_product(circ_1,circ_2,N) = inner_product
# perpendicular_circuits(N,circ) = perp_circuits
# circuit_reverser(C,M_c_l) = C_reversed
# generators2circuit(M,gp) = C_reversed
# entropy_region(M,region) = entropy
# mutual_information(N,circuit,region_A,region_B) = m_info
# topological_entropy(stab, N, M, draws) = S
# error_operator(adj_mat,size) = error_string
# mean_entropy(N,circ) = mean_entropy
# depth(N,circ) = D
# remove_id(circ) = new_circ
# simplify_circuit(N,circ,loops) = simplified_circuit
# draw_circuit(N,circ, file_name, scale, fold)
# random_stabilizer(N,T) = circuit

def randcirc(N,T,adj_mat):
    
    circuit = []
    
    for i in range(T):
        
        for j in range(N):
            
            operation = np.random.randint(0,4)
            
            if operation == 3:
                qubit_B = uf.choose_qubit(adj_mat,j)
                
                circuit.append([operation,j,qubit_B])
                
            elif operation != 0:
                
                circuit.append([operation,j,0])
                
    circuit = np.array(circuit)
    
    return circuit


def gene_express(N,circ):
    
    ''' 
    creates a qiskit circuit from an array circuit
    '''
    
    circuit = QuantumCircuit(int(N))
    # initial_state = [1,0]
    
    # for i in range(N):
    #     circuit.initialize(initial_state, i) 
    
    
    for i in range(len(circ)):
        gate = circ[i][0]
        
        if gate == 1:
            v0 =  int(circ[i][1]) # select qubit for H
            circuit.h(v0)
        elif gate == 2:
            v0 =  int(circ[i][1]) # select qubit for Phase gate
            circuit.s(v0)
        elif gate == 3:
            v0 =  int(circ[i][1]) # select control qubit for CNOT
            v1 = int(circ[i][2]) # select target qubit for CNOT
            circuit.cx(v0, v1)
        elif gate == 4:
            v0 =  int(circ[i][1]) # select qubit for X gate
            circuit.x(v0)
        elif gate == 5:
            v0 =  int(circ[i][1]) # select qubit for Y gate
            circuit.y(v0)
        elif gate == 6:
            v0 =  int(circ[i][1]) # select qubit for Z gate
            circuit.z(v0)
            
    return circuit

def final_state(N,circ):
    
    q_circ = gene_express(N,circ)
    qobj = assemble(q_circ)
    sim = Aer.get_backend('statevector_simulator')
    result = sim.run(qobj).result()
    out_state = result.get_statevector()
    
    return out_state

def stabilizer(N,circuit):
    
    ''' 
    evaluates the stabilizer matrix of a given circuit
    
    00 -> I
    01 -> Z
    10 -> X
    11 -> Y
    
    '''

    stabilizer = np.zeros([N,2*N])
    global_phase = np.zeros(N)
    
    #initializing generators for initial |0>^N state
    for i in range(N):
        stabilizer[i,2*i+1] = 1

    if len(circuit) != 0:
        for i in range(len(circuit)): #running through circuit
        
            if circuit[i,0] == 1: #Hadamard
            
                for j in range(N):
                    
                    if  stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1:
                        global_phase[j] = np.mod(global_phase[j]+2,4)
                        
                    temp_x = stabilizer[j,2*circuit[i,1]]
                    stabilizer[j,2*circuit[i,1]] = stabilizer[j,2*circuit[i,1]+1]
                    stabilizer[j,2*circuit[i,1]+1] = temp_x
                    
            elif circuit[i,0] == 2: #Phase
                
                for j in range(N):
                    
                    if  stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1:
                        global_phase[j] = np.mod(global_phase[j]+2,4)
                
                    stabilizer[j,2*circuit[i,1]+1] = np.mod(stabilizer[j,2*circuit[i,1]]+stabilizer[j,2*circuit[i,1]+1],2)
            
            elif circuit[i,0] == 3: #CNOT (left)
            
                 for j in range(N):
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1 and stabilizer[j,2*circuit[i,2]] == 1 and stabilizer[j,2*circuit[i,2]+1] == 1) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0 and stabilizer[j,2*circuit[i,2]] == 0 and stabilizer[j,2*circuit[i,2]+1] == 1):
                        global_phase[j] = np.mod(global_phase[j]+2,4)
                        
        
                    stabilizer[j,2*circuit[i,2]] = np.mod(stabilizer[j,2*circuit[i,2]] + stabilizer[j,2*circuit[i,1]] ,2)
                    stabilizer[j,2*circuit[i,1]+1] = np.mod(stabilizer[j,2*circuit[i,1]+1] + stabilizer[j,2*circuit[i,2]+1],2)
            
            elif circuit[i,0] == 4: #X
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 0 and stabilizer[j,2*circuit[i,1]+1] == 1) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[j] = np.mod(global_phase[j]+2,4)
                        
            elif circuit[i,0] == 5: #Y
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0) or (stabilizer[j,2*circuit[i,1]] == 0 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[j] = np.mod(global_phase[j]+2,4)
                        
            elif circuit[i,0] == 6: #Z
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[j] = np.mod(global_phase[j]+2,4)

    return stabilizer.astype(int), global_phase.astype(int)

def measurement(M,gp,operator,qubit):
       #INPUT: M and gp are the usual objects. Operador is 'X', or 'Y' or 'Z',
       #it indicates which operator will be measured. qubit is the number of the
       #qubit measured. 
       #OUTPUT: M and gp post-measurement
       
       #measurement_lit and _bit are measurements of the N+i qubit of the 
       #Z operator
       
       N_qubits = len(gp)
       M_lit = stabilizer_literal(M,gp)
       
       measurement_lit = ['I']*(N_qubits)
       measurement_lit[qubit] = operator
       measurement_bit = [0]*(2*N_qubits)
       if operator == 'Z':
           measurement_bit[2*qubit+1] = 1
       elif operator == 'X':
           measurement_bit[2*qubit] = 1
       elif operator == 'Y':
           measurement_bit[2*qubit] = 1
           measurement_bit[2*qubit+1] = 1
   
       M_post_measurement = []
       gp_post_measurement = []
       
       for j in range(N_qubits):
               
           if uf.pauli_commutation(measurement_lit,M_lit[j,:-1]) == 1:
               M_post_measurement.append(list(M[j,:]))
               gp_post_measurement.append(gp[j])
           else:
                M_post_measurement.append(measurement_bit)
                gp_post_measurement.append(0)
                
       M_post_measurement = np.array(M_post_measurement)
       gp_post_measurement = np.array(gp_post_measurement)
       
       return M_post_measurement, gp_post_measurement

def stabilizer_literal(stab_bit, global_phase):
    
    ''' 
    return the literal stabilizer matrix for a given circuit
    '''
    
    N = int(np.shape(stab_bit)[1]/2.)
    
    stab_lit = [['I' for columns in range(N+1)] for rows in range(len(global_phase))]
    
    for i in range(len(global_phase)):
        for j in range(N):
            if stab_bit[i,2*j] == 0 and stab_bit[i,2*j+1] == 1:
                stab_lit[i][j] = 'Z'
            elif stab_bit[i,2*j] == 1 and stab_bit[i,2*j+1] == 0:
                stab_lit[i][j] = 'X'
            elif stab_bit[i,2*j] == 1 and stab_bit[i,2*j+1] == 1:
                stab_lit[i][j] = 'Y'
                
        if global_phase[i] == 0:
            stab_lit[i][-1] = '1'
        elif global_phase[i] == 1:
            stab_lit[i][-1] = 'i'
        elif global_phase[i] == 2:
            stab_lit[i][-1] = '-1'
        else:
            stab_lit[i][-1] = '-i'
            
    stab_lit = np.array(stab_lit)
                
    return stab_lit

def stabilizer_bit(stab_lit):
    
    
    N = np.shape(stab_lit)[1]-1
    n_op = np.shape(stab_lit)[0]
    
    M = np.zeros([n_op,2*N])
    gp = np.zeros(n_op)
    
    for i in range(n_op):
        for j in range(N):
            if stab_lit[i,j] == 'X':
                M[i,2*j] = 1
            elif stab_lit[i,j] == 'Z':
                M[i,2*j+1] = 1
            elif stab_lit[i,j] == 'Y':
                M[i,2*j] = 1
                M[i,2*j+1] = 1
        if stab_lit[i,-1] == '-1':
            gp[i] = 2
        elif stab_lit[i,-1] == 'i':
             gp[i] = 1
        elif stab_lit[i,-1] == '-i':
             gp[i] = 3
                
    return M, gp

def stabilizer_group(N,circuit, stab_group = None, gp_group = None):
    

    #Calculates the initial stabilizer tableau if not provided. The user
    #may provide the initial tableau if it is going to be evaluated many
    #times, because the next lines of code, inside the if, does not scale
    #well with N (it grows as 2**N), therefore avoiding it can save a lot
    #of time.
    if str(type(stab_group)) == "<class 'NoneType'>" or str(type(gp_group)) == "<class 'NoneType'>":

        M , gp = stabilizer(N, np.array([[0,0,0],[0,0,0]])) #calculates the tableau for |0>**N
        
        
        stab_group = np.zeros([2**N-1,2*N])
        gp_group = np.zeros(2**N-1)
        
        lst = list(itertools.product([0, 1], repeat=N))[1:]
        
            
        for i in range(len(lst)):
        
             R = np.zeros(2*N)
             R_phase = 0 
        
             for j in range(N):
                if lst[i][j] == 1:
                    R_phase = np.mod(R_phase+gp[j],4)
                    for n in range(N):
                        ans, signal = uf.pauli_mult(M[j,2*n:2*n+2],R[2*n:2*n+2])
                        R[2*n:2*n+2] = ans
                        R_phase = np.mod(R_phase+signal,4)
                        
             stab_group[i,:] = R
             gp_group[i] = R_phase
         

    stab_group_temp = np.copy(stab_group)
    gp_group_temp = np.copy(gp_group)
    
    
    #applying the circuit to the stabilizer group
    if len(circuit) != 0:
    
        for i in range(len(circuit)): #running through circuit
        
            if circuit[i,0] == 1: #Hadamard
            
                for j in range(2**N-1):
                    
                    if  stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 1:
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
                        
                    temp_x = stab_group_temp[j,2*circuit[i,1]]
                    stab_group_temp[j,2*circuit[i,1]] = stab_group_temp[j,2*circuit[i,1]+1]
                    stab_group_temp[j,2*circuit[i,1]+1] = temp_x
                    
            elif circuit[i,0] == 2: #Phase
                
                for j in range(2**N-1):
                    
                    if  stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 1:
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
                
                    stab_group_temp[j,2*circuit[i,1]+1] = np.mod(stab_group_temp[j,2*circuit[i,1]]+stab_group_temp[j,2*circuit[i,1]+1],2)
            
            elif circuit[i,0] == 3: #CNOT (left)
            
                 for j in range(2**N-1):
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 1 and stab_group_temp[j,2*circuit[i,2]] == 1 and stab_group_temp[j,2*circuit[i,2]+1] == 1) or (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 0 and stab_group_temp[j,2*circuit[i,2]] == 0 and stab_group_temp[j,2*circuit[i,2]+1] == 1):
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
                        
            
                    stab_group_temp[j,2*circuit[i,2]] = np.mod(stab_group_temp[j,2*circuit[i,2]] + stab_group_temp[j,2*circuit[i,1]] ,2)
                    stab_group_temp[j,2*circuit[i,1]+1] = np.mod(stab_group_temp[j,2*circuit[i,1]+1] + stab_group_temp[j,2*circuit[i,2]+1],2)
            
            elif circuit[i,0] == 4: #X
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp[j,2*circuit[i,1]] == 0 and stab_group_temp[j,2*circuit[i,1]+1] == 1) or (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
                        
            elif circuit[i,0] == 5: #Y
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 0) or (stab_group_temp[j,2*circuit[i,1]] == 0 and stab_group_temp[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
                        
            elif circuit[i,0] == 6: #Z
            
                for j in range(2**N-1):
                    
                     if (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 0) or (stab_group_temp[j,2*circuit[i,1]] == 1 and stab_group_temp[j,2*circuit[i,1]+1] == 1):
                        gp_group_temp[j] = np.mod(gp_group_temp[j]+2,4)
    
    #making literal version
    
    stab_lit = [['I' for columns in range(N+1)] for rows in range(len(gp_group_temp))]
        
    for i in range(len(gp_group_temp)):
        for j in range(N):
            if stab_group_temp[i,2*j] == 0 and stab_group_temp[i,2*j+1] == 1:
                stab_lit[i][j] = 'Z'
            elif stab_group_temp[i,2*j] == 1 and stab_group_temp[i,2*j+1] == 0:
                stab_lit[i][j] = 'X'
            elif stab_group_temp[i,2*j] == 1 and stab_group_temp[i,2*j+1] == 1:
                stab_lit[i][j] = 'Y'
                
        if gp_group_temp[i] == 0:
            stab_lit[i][-1] = '1'
        elif gp_group_temp[i] == 1:
            stab_lit[i][-1] = 'i'
        elif gp_group_temp[i] == 2:
            stab_lit[i][-1] = '-1'
        else:
            stab_lit[i][-1] = '-i'
            
    stab_lit = np.array(stab_lit)

    return stab_group_temp, gp_group_temp, stab_lit

def stabilizer_group_computational(N):
    
    #This function calculates the stabilizer grupo for the |0>**N state
    #The output is use extensively throughout the script and its evaluation
    #is time consuming, therefore a lot of time is gained by evaluating one
    #time and reusing it
    
    M , gp = stabilizer(N, np.array([[0,0,0],[0,0,0]]))
    stab_group = np.zeros([2**N-1,2*N])
    gp_group = np.zeros(2**N-1)
    
    lst = list(itertools.product([0, 1], repeat=N))[1:]
    
    for i in range(len(lst)):
    
          R = np.zeros(2*N)
          R_phase = 0 
    
          for j in range(N):
            if lst[i][j] == 1:
                R_phase = np.mod(R_phase+gp[j],4)
                for n in range(N):
                    ans, signal = uf.pauli_mult(M[j,2*n:2*n+2],R[2*n:2*n+2])
                    R[2*n:2*n+2] = ans
                    R_phase = np.mod(R_phase+signal,4)
                    
          stab_group[i,:] = R
          gp_group[i] = R_phase
          
    return stab_group, gp_group

def stabilizer_group_X(N):
    
    circuit_temp = np.zeros([N,3])
    
    for i in range(N):
        
        circuit_temp[i,:] = np.array([1,i,0])
        
    circuit_temp = circuit_temp.astype(int)
    
    M , gp = stabilizer(N, circuit_temp)
    stab_group_X = np.zeros([2**N-1,2*N])
    gp_group_X = np.zeros(2**N-1)
    
    lst = list(itertools.product([0, 1], repeat=N))[1:]
    
    for i in range(len(lst)):
    
          R = np.zeros(2*N)
          R_phase = 0 
    
          for j in range(N):
            if lst[i][j] == 1:
                R_phase = np.mod(R_phase+gp[j],4)
                for n in range(N):
                    ans, signal = uf.pauli_mult(M[j,2*n:2*n+2],R[2*n:2*n+2])
                    R[2*n:2*n+2] = ans
                    R_phase = np.mod(R_phase+signal,4)
                    
          stab_group_X[i,:] = R
          gp_group_X[i] = R_phase
          
    return stab_group_X, gp_group_X

def stabilizer_canonical(stabilizer,global_phase):
    
    stab_can = np.copy(stabilizer)
    gp_can = np.copy(global_phase)
    
    N = len(gp_can)
    
    #X BLOCK
    
    i = 0 #initialize first row of submatrix
    
    j = 0 #initialize first column of submatrix
    
    search = False
    
    k = i #reset k
    
    while search == False: 
        
        literal = uf.pauli_decoder(stab_can[k,2*j:2*j+2])
    
        if literal == 'X' or literal == 'Y': #if True seach succeded
            stab_can[[i,k]] = stab_can[[k,i]] #row swap
            gp_can[[i,k]] = gp_can[[k,i]] #row swap
            
            #for m in range(i+1,N): #row multiplication ; 1st case, multiplication only occurs in submatrix
            for m in range(N): #row multiplication ; 2nd case, multiplication occurs in whole matrix (probably the right case!!!)
                
                literal = uf.pauli_decoder(stab_can[m,2*j:2*j+2])
                
                if m != i and (literal == 'X' or literal == 'Y'):
                    
                    gp_can[m] = np.mod(gp_can[m]+gp_can[i],4)
                    
                    for n in range(N):
                            ans, signal = uf.pauli_mult(stab_can[i,2*n:2*n+2],stab_can[m,2*n:2*n+2])
                            stab_can[m,2*n:2*n+2] = ans
                            gp_can[m] = np.mod(gp_can[m]+signal,4)
                            
            i += 1 #reduce submatrix height
            j += 1 #reduce submatrix width
            k = i #reset k
            
        else:
            
            k += 1 #if search failed, look for next row
            
            if k == N: #search failed for all rows for column j
                j += 1 #go to next column
                k = i #reset k
            
            
        if i == N or j == N: # at this point the X block is finished
            search = True
            
    #Z BLOCK
    
    j = 0 #initilize first column of submatrix, notice that the first row remains the same (i is not changed)
    
    if i != N:
        search = False
    
    k = i #reset k
    
    while search == False: 
        
        literal = uf.pauli_decoder(stab_can[k,2*j:2*j+2])
    
        if literal == 'Z': #if True seach succeded
            stab_can[[i,k]] = stab_can[[k,i]] #row swap
            gp_can[[i,k]] = gp_can[[k,i]] #row swap
            
            #for m in range(i+1,N): #row multiplication ; #row multiplication ; 1st case, multiplication only occurs in submatrix
            for m in range(N): #row multiplication ; 2nd case, multiplication occurs in whole matrix (probably the right case!!!)
                
                literal = uf.pauli_decoder(stab_can[m,2*j:2*j+2])
                
                if m != i and (literal == 'Z' or literal == 'Y'):
                    
                    gp_can[m] = np.mod(gp_can[m]+gp_can[i],4)
                    
                    for n in range(N):
                            ans, signal = uf.pauli_mult(stab_can[i,2*n:2*n+2],stab_can[m,2*n:2*n+2])
                            stab_can[m,2*n:2*n+2] = ans
                            gp_can[m] = np.mod(gp_can[m]+signal,4)
                            
            i += 1 #reduce submatrix height
            j += 1 #reduce submatrix width
            k = i #reset k
            
        else:
            
            k += 1 #if search failed
            
            if k == N: #search failed for all rows for column j
                j += 1 #width of submatrix decreased
                k = i #reset k
            
            
        if i == N or j == N: # at this point the Z block is finished
            search = True
            
    return stab_can, gp_can

def basisnormcirc(M, gp):

    M_c, gp_c = stabilizer_canonical(M,gp) #canonical stabilizer matrix
    
    N = len(gp_c)
    
    hadamard_block_1 = []
    CNOT_block_2 = []
    CZ_block_3 = []
    phase_block_4 = []
    hadamard_block_5 = []
    
    
    #Hadamard blocks
    
    i = 0
    
    M_c_l = stabilizer_literal(M_c, gp_c)
    
    
    for j in range(N):
        
        k = uf.find_literal('X',M_c_l[i:,j]) #indices of the rows where the j'th entry is X, if exists
        
        if len(k) != 0: #X was found
            k = k[0]+i #grab first index
            M_c[[i,k]] = M_c[[k,i]] #row swap
            gp_c[[i,k]] = gp_c[[i,k]] #row swap
            M_c_l = stabilizer_literal(M_c, gp_c) #update M_c_l
            
        elif len(k) == 0:
            
            k = uf.find_literal('Y',M_c_l[i:,j]) #indice of the rows where the j'th entry is X, if exists
            
            if len(k) != 0: #Y was found
                k = k[0]+i #grab first index
                M_c[[i,k]] = M_c[[k,i]] #row swap
                gp_c[[i,k]] = gp_c[[i,k]] #row swap
                M_c_l = stabilizer_literal(M_c, gp_c) #update M_c_l
                
            elif len(k) == 0:
                    
                k = uf.find_literal('Z',M_c_l[i:,j]) #indice of the rows where the j'th entry is Z, if exists
                
                if len(k) != 0: #there's at least a row with Z in j'th column
                
                    k = k[-1]+i #get index of last row containing Z
                    M_c[[i,k]] = M_c[[k,i]] #row swap
                    gp_c[[i,k]] = gp_c[[k,i]] #row swap
                    M_c_l = stabilizer_literal(M_c, gp_c) #update M_c_l
                    
                    if (len(uf.find_literal('X',M_c_l[i,j+1:])) + len(uf.find_literal('Y',M_c_l[i,j+1:])) + len(uf.find_literal('Z',M_c_l[i,j+1:]))) != 0: #if True apply Hadamard to qubit j
        
                        for n in range(N):#running lines
                    
                            if  M_c_l[n,j] == 'Y':
                                gp_c[n] = np.mod(gp_c[n]+2,4)
                            
                            temp_x = M_c[n,2*j]
                            M_c[n,2*j] = M_c[n,2*j+1]
                            M_c[n,2*j+1] = temp_x
                            
                        hadamard_block_1.append([1,j,0]) #append Hadamard gate to circuit C
                        M_c_l = stabilizer_literal(M_c, gp_c) #update M_c_l
                        
        i += 1
        
        if i == N:
            break
            
    #CNOT block
    
    
    for j in range(N): #run through rows
            
    
        for k in range(j+1,N): #run through columns
        
            if M_c_l[j,k] == 'X' or M_c_l[j,k] == 'Y': #if True apply CNOT between qubit j and qubit k
                
                for n in range(N): #running rows
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (M_c_l[n,j] == 'Y' and M_c_l[n,k] == 'Y') or (M_c_l[n,j] == 'X' and M_c_l[n,k] == 'Z'):
                        gp_c[n] = np.mod(gp_c[n]+2,4)
                        
        
                    M_c[n,2*k] = np.mod(M_c[n,2*k] + M_c[n,2*j] ,2)
                    M_c[n,2*j+1] = np.mod(M_c[n,2*j+1] + M_c[n,2*k+1],2)
            
                
                
                CNOT_block_2.append([3,j,k]) #append CNOT gate to circuit C
                M_c_l = stabilizer_literal(M_c, gp_c) #update M_c literal
                
    #CZ block
    
    
    for j in range(N): #run through rows
            
    
        for k in range(j+1,N): #run through columns
        
            if M_c_l[j,k] == 'Z': #if True apply CZ between qubit j and qubit k
            
            #the CZ gate between j and k is equivalent of applying a Hadamard to k, a CNOT to j and k and a Hadamard to k in this order
                
                #applying Hadamard to qubit k
                for n in range(N):#running lines
                    
                    if  M_c_l[n,k] == 'Y':
                        gp_c[n] = np.mod(gp_c[n]+2,4)
                    
                    temp_x = M_c[n,2*k]
                    M_c[n,2*k] = M_c[n,2*k+1]
                    M_c[n,2*k+1] = temp_x
                    
                M_c_l = stabilizer_literal(M_c, gp_c) #uptade M_c_l
            
    
                #applying CNOT to j and k
                for n in range(N): #running rows
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (M_c_l[n,j] == 'Y' and M_c_l[n,k] == 'Y') or (M_c_l[n,j] == 'X' and M_c_l[n,k] == 'Z'):
                        gp_c[n] = np.mod(gp_c[n]+2,4)
                        
        
                    M_c[n,2*k] = np.mod(M_c[n,2*k] + M_c[n,2*j] ,2)
                    M_c[n,2*j+1] = np.mod(M_c[n,2*j+1] + M_c[n,2*k+1],2)
                  
                    
                M_c_l = stabilizer_literal(M_c, gp_c) #uptade M_c_l
                    
                #applying Hadamard to qubit k
                for n in range(N):#running lines
                    
                    if  M_c_l[n,k] == 'Y':
                        gp_c[n] = np.mod(gp_c[n]+2,4)
                    
                    temp_x = M_c[n,2*k]
                    M_c[n,2*k] = M_c[n,2*k+1]
                    M_c[n,2*k+1] = temp_x
            
                
                CZ_block_3.append([1,k,0]) #append Hadamard gate to qubit k to circuit C
                CZ_block_3.append([3,j,k]) #append CNOT gate to circuit C
                CZ_block_3.append([1,k,0]) #append Hadamard gate to qubit k to circuit C
                M_c_l = stabilizer_literal(M_c, gp_c) #update M_c literal
    
    #Phase block
    
    
    for j in range(N): #run through rows
            
        if M_c_l[j,j] == 'Y': #if True apply Phase to qubit j
            
            for n in range(N):
                
                if  M_c_l[n,j] == 'Y':
                    gp_c[n] = np.mod(gp_c[n]+2,4)
            
                M_c[n,2*j+1] = np.mod(M_c[n,2*j]+M_c[n,2*j+1],2)
                
            phase_block_4.append([2,j,0]) #append Phase gate to circuit C
            M_c_l = stabilizer_literal(M_c, gp_c) #update M_c literal
            
    #Hadamard block
    
    for j in range(N): #run through rows
            
        if M_c_l[j,j] == 'X': #if True apply Hadamard to qubit j
            
            for n in range(N):#running lines
                    
                if  M_c_l[n,j] == 'Y':
                    gp_c[n] = np.mod(gp_c[n]+2,4)
                
                temp_x = M_c[n,2*j]
                M_c[n,2*j] = M_c[n,2*j+1]
                M_c[n,2*j+1] = temp_x
                
            hadamard_block_5.append([1,j,0]) #append Hadamard gate to qubit j to circuit C
            M_c_l = stabilizer_literal(M_c, gp_c) #update M_c literal
                    
                        
    #Eliminate trailing Z literals to ensure basis form
    
    M_c_l = stabilizer_literal(M_c, gp_c)
    
    for j in range(N): 
            
        for k in range(j+1,N):
        
            if M_c_l[k,j] == 'Z': #if True multiply row k and j
                
                gp_c[k] = np.mod(gp_c[k]+gp_c[j],4)
                        
                for n in range(N):
                        ans, signal = uf.pauli_mult(M_c[j,2*n:2*n+2],M_c[k,2*n:2*n+2])
                        M_c[k,2*n:2*n+2] = ans
                        gp_c[k] = np.mod(gp_c[k]+signal,4)
            
                M_c_l = stabilizer_literal(M_c, gp_c) #update M_c literal
    
    
    
    circuit_array = []
    if len(hadamard_block_1) != 0:
        for k in range(len(hadamard_block_1)):
            circuit_array.append(hadamard_block_1[k])
    if len(CNOT_block_2) != 0:
        for k in range(len(CNOT_block_2)):
            circuit_array.append(CNOT_block_2[k])
    if len(CZ_block_3) != 0:
        for k in range(len(CZ_block_3)):
            circuit_array.append(CZ_block_3[k])
    if len(phase_block_4) != 0:
        for k in range(len(phase_block_4)):
            circuit_array.append(phase_block_4[k])
    if len(hadamard_block_5) != 0:
        for k in range(len(hadamard_block_5)):
            circuit_array.append(hadamard_block_5[k])
    
    circuit_array = np.array(circuit_array)
    
    C = [hadamard_block_1 , CNOT_block_2 , CZ_block_3 , phase_block_4 , hadamard_block_5 , circuit_array]
    
    return C, M_c_l

def check_orthogonality(circ_1,circ_2,N):
    """
    Checks if circ_1 and circ_2 are orthogonal. The output is bool, True if
    circuits are orthogonal and False if not
    """
    
    M_1, gp_1 = stabilizer(N,circ_1)
    stab_group_1, gp_group_1, stab_lit_1 = stabilizer_group(N,circ_1)
    
    stab_group_1_list = [tuple(stab_lit_1[x,:N]) for x in range(len(stab_lit_1))]
    
    M_2, gp_2 = stabilizer(N,circ_2)
    M_l_2 = stabilizer_literal(M_2, gp_2)
    
    orthogonality = False
    
    for i in range(N):
        
        if tuple(M_l_2[i,:N]) in stab_group_1_list:
            
            idx = stab_group_1_list.index(tuple(M_l_2[i,:N]))
            
            if gp_group_1[idx] != gp_2[i]:
    
                orthogonality = True
                break
                
    return orthogonality

def inner_product(circ_1,circ_2,N):
    
    #testing if the states are orthogonal
    
    stop = False
    
    orthogonality = check_orthogonality(circ_1,circ_2,N)
    if orthogonality == True:
        inner_product = 0
        stop = True
    
    
    ######################################
    
    if stop == False:
        
        psi =  np.copy(circ_1)
        phi =  np.copy(circ_2)
        M_psi, gp_psi = stabilizer(N,psi)
        C_psi = basisnormcirc(M_psi, gp_psi)[0]
        
        
        if len(C_psi[-1]) == 0:
            phi_conj_c = phi
        elif len(phi) == 0:
            phi_conj_c = C_psi[-1]
        else:
            phi_conj_c = np.vstack([phi,C_psi[-1]])
    
        M_phi, gp_phi = stabilizer(N,phi_conj_c)
        M_phi_c, gp_phi_c = stabilizer_canonical(M_phi, gp_phi)
        M_phi_l = stabilizer_literal(M_phi_c, gp_phi_c)
        
        k = 0

        for i in range(N):
            
            X_Y_number = len(uf.find_literal('X',M_phi_l[i,:])) + len(uf.find_literal('Y',M_phi_l[i,:]))
            
            if X_Y_number != 0 : 
                k += 1
            
        inner_product = 2**(-1*k/2.)

    return inner_product

def orthogonal_states(N,circuit, global_phase_initial = None):
    
    #Since the tableau for all states is the same, only the
    #global phases change (and are all combinations of "+1" and "-1"),
    #it suffices to calculate 
    
    if str(type(global_phase_initial)) == "<class 'NoneType'>":
        
        #each element is a different combination of "+1" and "-1" global phase
        global_phase_initial = np.array(list(itertools.product([0, 2], repeat=N))[1:])
        
    stabilizer = np.zeros([N,2*N]) #each element share the same tableau
    for i in range(N):
        stabilizer[i,2*i+1] = 1
        
    #each line of global_phase is the final global phase of each orthogonal state
    global_phase = np.copy(global_phase_initial)
    
    if len(circuit) != 0:
        
        for i in range(len(circuit)): #running through circuit
        
            if circuit[i,0] == 1: #Hadamard
            
                for j in range(N):
                    
                    if  stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1:
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)
                        
                    temp_x = stabilizer[j,2*circuit[i,1]]
                    stabilizer[j,2*circuit[i,1]] = stabilizer[j,2*circuit[i,1]+1]
                    stabilizer[j,2*circuit[i,1]+1] = temp_x
                    
            elif circuit[i,0] == 2: #Phase
                
                for j in range(N):
                    
                    if  stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1:
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)
                
                    stabilizer[j,2*circuit[i,1]+1] = np.mod(stabilizer[j,2*circuit[i,1]]+stabilizer[j,2*circuit[i,1]+1],2)
            
            elif circuit[i,0] == 3: #CNOT (left)
            
                 for j in range(N):
                     
                    #this if tests if the generator pre-conjugation is Y_1*Y_2 or X_1*Z_2
                    if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1 and stabilizer[j,2*circuit[i,2]] == 1 and stabilizer[j,2*circuit[i,2]+1] == 1) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0 and stabilizer[j,2*circuit[i,2]] == 0 and stabilizer[j,2*circuit[i,2]+1] == 1):
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)
                        
        
                    stabilizer[j,2*circuit[i,2]] = np.mod(stabilizer[j,2*circuit[i,2]] + stabilizer[j,2*circuit[i,1]] ,2)
                    stabilizer[j,2*circuit[i,1]+1] = np.mod(stabilizer[j,2*circuit[i,1]+1] + stabilizer[j,2*circuit[i,2]+1],2)
            
            elif circuit[i,0] == 4: #X
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 0 and stabilizer[j,2*circuit[i,1]+1] == 1) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)
                        
            elif circuit[i,0] == 5: #Y
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0) or (stabilizer[j,2*circuit[i,1]] == 0 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)
                        
            elif circuit[i,0] == 6: #Z
            
                for j in range(N):
                    
                     if (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 0) or (stabilizer[j,2*circuit[i,1]] == 1 and stabilizer[j,2*circuit[i,1]+1] == 1):
                        global_phase[:,j] = np.mod(global_phase[:,j]+2,4)

    return stabilizer, global_phase


def circuit_reverser(C,M_c_l):
    
    N = len(M_c_l)

    x_block = []
    
    #forming x_block
    
    for i in range(N):
        if M_c_l[i,N] == -1: #if true, apply X to change |1> to |0>
            x_block.append([4,i,0])
            
    #notice that the Hadamards blocks don't need any operation to be reversed        
    
    #reversing order of CNOT block contained in C[1]
    
    CNOT_block_reversed = copy.deepcopy( C[1] )
    CNOT_block_reversed.reverse()
    
    #reversing order of CZ block contained in C[2]
    
    CZ_block_reversed = copy.deepcopy( C[2] )
    CZ_block_reversed.reverse()
    
    #reversing order of phase block contained in C[3]
    
    phase_block_reversed = []
    
    for i in range(len(C[3])):
        for j in range(3):
            phase_block_reversed.append(C[3][i])
            
    #mounting the reversed circuit
    
    C_array_reversed = []
    if len(x_block) != 0:
        for k in range(len(x_block)):
            C_array_reversed.append(x_block[k])
    if len(C[4]) != 0: #appending Hadamard block 5
        for k in range(len(C[4])):
            C_array_reversed.append(C[4][k])
    if len(phase_block_reversed) != 0:
        for k in range(len(phase_block_reversed)):
            C_array_reversed.append(phase_block_reversed[k])
    if len(CZ_block_reversed) != 0:
        for k in range(len(CZ_block_reversed)):
            C_array_reversed.append(CZ_block_reversed[k])
    if len(CNOT_block_reversed) != 0:
        for k in range(len(CNOT_block_reversed)):
            C_array_reversed.append(CNOT_block_reversed[k])
    if len(C[0]) != 0: #appending Hadamard block 1
        for k in range(len(C[0])):
            C_array_reversed.append(C[0][k])
    
    C_array_reversed = np.array(C_array_reversed)
    
    C_reversed = [x_block, C[4], phase_block_reversed, CZ_block_reversed, CNOT_block_reversed, C[0], C_array_reversed]
        
    return C_reversed
    
def generators2circuit(M,gp):
    
    C, M_c_l = basisnormcirc(M, gp)
    
    C_reversed = circuit_reverser(C,M_c_l)
    
    return C_reversed

def entropy_region(N,circuit,region):
    
    ''' 
    evaluates the mean von Neumann entropy of a given circuit
    '''
    
    M, gp = stabilizer(N,circuit)
    
    list_del_col = []
    
    for i in range(len(M)):
        if (i in region) == False:
            list_del_col.append(2*i)
            list_del_col.append(2*i+1)
            
    reduced_M = np.copy(M)
    reduced_M = np.delete(reduced_M,list_del_col,1)

    a,b = np.shape(reduced_M)
    I = uf.rank_mod2(reduced_M)
    entropy = I - b/2.
    
    return entropy

def mutual_information(N,circuit,region_A,region_B):
    
    entropy_A = entropy_region(N,circuit,region_A)
    entropy_B = entropy_region(N,circuit,region_B)
    entropy_AB = entropy_region(N,circuit,region_A+region_B)
    
    m_info = entropy_A + entropy_B - entropy_AB
    
    return m_info
    
def topological_entropy(stab, N, M, draws):
    
    rep = np.zeros([N,M]).astype(int)
    
    for i in range(N):
        for j in range(M):
            rep[i,j] = i*(2*M+1) + j
            
    p_row = ceil(N/2.)-1
    p_col = ceil(M/2.)-1
    
    #Mounting A region
    
    A = [rep[p_row,p_col]]
    
    for i in range(p_row):
        if i == 0:
            for j in range(p_col,M):
                if j == p_col:
                    A.append(A[-1] - (M + 1))
                    A.append(A[-1] - M)
                    A.append(A[-1] + M+1)
                else:
                    A.append(A[-1] - M)
                    A.append(A[-2] + M+1)
                    A.append(A[-1] - M)
            for j in range(0,p_col):
                if j == 0:
                    A.append(A[0]-1)
                    A.append(A[-1] - (2*M+1))
                    A.append(A[-1] + M)
                else:
                    A.append(A[-1] - (M+1))
                    A.append(A[-2] + M)
                    A.append(A[-1] - (M+1))
        else:
            
            A.append(rep[p_row-i,0] + M)
            
            for j in range(M):
    
                A.append(A[-1] - M)
                A.append(A[-1]+ M+1)

    #Mounting B region
    
    for i in range(p_col+1,M):
        if i == p_col+1:
            B = [rep[p_row,p_col] + M + 1]
            B.append(B[-1] + M +1)
            B.append(B[-1] - M)
            
            for j in range(p_row+1,N):
                B.append(B[-2]+M)
                B.append(B[-1] + M+1)
                B.append(B[-1] - M)
                
        else:
            for j in range(p_row,N):
                if j == p_row:
                    B.append(rep[p_row,p_col] + 2*M + +2 + 1*(i - (p_col+1)))
                    B.append(B[-1] - M)
                else:
                    B.append(B[-2]+2*M+1)
                    B.append(B[-1] - M)
                
    #Mounting C region
    
    for i in range(p_col+1):
        if i == 0:
            C = [rep[p_row,p_col] + M , rep[p_row,p_col] + 2*M+1]
            
            for j in range(p_row+1,N):
                C.append(C[-1] + M)
                C.append(C[-1] + M + 1)
                
        else:
            for j in range(p_row,N):
                if j == p_row:
                    C.append(rep[p_row,p_col] + M - 1*i)
                    C.append(C[-1] + M+1)
                else:
                    C.append(C[-1] + M)
                    C.append(C[-1] + M+1)
         
    S = []           
         
    for i in range(draws):
    
        A_size = np.random.randint(4,len(A)+1)
        B_size = np.random.randint(5,len(B)+1)
        C_size = np.random.randint(5,len(C)+1) #region C must have at least 2 qubits!
        
        
        Region_A = A[0:A_size]
        Region_B = B[0:B_size]
        Region_C = C[0:C_size]
        
        Region_AB = Region_A + Region_B
        Region_AC = Region_A + Region_C
        Region_BC = Region_B + Region_C
        Region_ABC = Region_A + Region_B + Region_C
        
        ent_A = entropy_region(stab, Region_A)
        ent_B = entropy_region(stab, Region_B)
        ent_C = entropy_region(stab, Region_C)
        ent_AB = entropy_region(stab, Region_AB)
        ent_BC = entropy_region(stab, Region_BC)
        ent_AC = entropy_region(stab, Region_AC)
        ent_ABC = entropy_region(stab, Region_ABC)
        
        S_temp = ent_A + ent_B + ent_C - ent_AB - ent_BC - ent_AC + ent_ABC 
        
        if S_temp != 0:
            S.append(S_temp)
            
    if len(S) != 0:
        S = set(S)
    else:
        S = 0
    
    return S

def error_operator(adj_mat,size):
    
    #returns an error circuit

    candidates = range(len(adj_mat))
    qubit_seed = random.sample(candidates,  1)[0]
    error = random.sample([4,6],1)[0]
    error_string = np.array([error,qubit_seed,0])
    
    last_qubit = error_string[1]
    b_last_qubit = -1
    
    for i in range(size-1):
        
        search = 0
        
        while search == 0:
            
            next_qubit = random.sample(candidates,  1)[0]
            
            if adj_mat[last_qubit,next_qubit] == 1 and next_qubit != b_last_qubit:
                
                error_string = np.vstack([error_string,[error,next_qubit,0]])
                search = 1
                
        last_qubit = error_string[-1,1]
        b_last_qubit = error_string[-2,1]
    
    
        
    return error_string
     
def mean_entropy(N,circ):
    
    ''' 
    evaluates the mean von Neumann entropy of a given circuit
    '''
    
    stab = stabilizer(N,circ)[0]

    mean_entropy = 0
    
    for i in range(N):
        list_del_col = []
        for j in range(i):
            list_del_col.append(2*j)
            list_del_col.append(2*j+1)
        stab_temp = np.delete(stab,list_del_col,1)
        a,b = np.shape(stab_temp)
        I = uf.rank_mod2(stab_temp)
        mean_entropy += I - b/2.
        
    mean_entropy = float(mean_entropy)/N
    
    return mean_entropy

def depth(N,circuit):

    init = np.zeros(N)
    
    for i in range(len(circuit)):
        
        if circuit[i,0] != 3 and circuit[i,0] != 0:
            init[circuit[i,1]] += 1
        if circuit[i,0] == 3:
            init[circuit[i,1]] = max(init[circuit[i,1]],init[circuit[i,2]]) + 1
            init[circuit[i,2]] = init[circuit[i,1]]
            
    D = max(init)
    
    return D

# def depth(N,circ):
    
#     ''' 
#     evaluates the depth of a given circuit
#     '''
    
#     c_qiskit = gene_express(N,circ)
#     D = c_qiskit.depth()
#     return D

def remove_id(circ):
    
    pos = range(len(circ))
    pos_with_id = [x for x in pos if circ[x,0] == 0]
    new_circ = np.copy(circ)
    new_circ = np.delete(new_circ,pos_with_id,0)
    
    return new_circ.astype(int)

def simplify_circuit(N,circ,loops):
    
    circuit = np.copy(circ[-1])
    circuit = remove_id(circuit)
    simplified_circuit = np.copy(circuit)
    simplified_circuit = simplified_circuit.astype(int)
    
    for i in range(loops):
        
        pos = range(len(simplified_circuit))
        p_1 = random.sample(pos,  1)[0]
        p_2 = random.sample(pos,  1)[0]
        original_gate_1 = simplified_circuit[p_1,0]
        original_gate_2 = simplified_circuit[p_2,0]
        simplified_circuit[p_1,0] = 0
        simplified_circuit[p_2,0] = 0
        if inner_product(circuit,simplified_circuit,N) == 1:
            simplified_circuit = remove_id(simplified_circuit)
            simplified_circuit = simplified_circuit.astype(int)
        else:
            simplified_circuit[p_1,0] = int(original_gate_1)
            simplified_circuit[p_2,0] = int(original_gate_2)
         
        # pos = range(len(simplified_circuit))
        # p_1 = random.sample(pos,  1)[0]
        # original_gate_1 = simplified_circuit[p_1,0]
        # simplified_circuit[p_1,0] = 0
        # if inner_product(circuit,simplified_circuit,N) == 1:
        #     simplified_circuit = remove_id(simplified_circuit)
        #     simplified_circuit = simplified_circuit.astype(int)
        # else:
        #     simplified_circuit[p_1,0] = int(original_gate_1)
            
    return simplified_circuit.astype(int)

def draw_circuit(N,circ, file_name, scale, fold):
    q_circ = gene_express(N,circ)
    circuit_drawer(q_circ, scale=scale, filename=file_name, style={'backgroundcolor': '#EEEEEE'}, output='mpl', interactive=False, plot_barriers=True, reverse_bits=False, justify=None, vertical_compression='high', idle_wires=True, with_layout=False, fold=fold, ax=None, initial_state=True, cregbundle=True)
    return 

def random_stabilizer(N,T):
    
    circuit = []
    
    for i in range(T):
        
        for j in range(N):
            
            operation = np.random.randint(0,4)
            
            if operation == 3:
                targets = list(range(N))
                targets.remove(j)
                target = random.sample(targets,  1)[0]
                
                circuit.append([operation,j,target])
                
            else:
                
                circuit.append([operation,j,0])
                
    circuit = np.array(circuit)
    
    return circuit