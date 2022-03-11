# -*- coding: utf-8 -*-

import numpy as np
import utility_functions as uf
import itertools


### LIST OF FUNTIONS ###
# syndrome_table(N,P) = s_table
# syndrome_operators(N,codewords) = P, P_bit
# logical_X(N,codewords) = logical_X_table
# phase_distance(N,codewords) = phase_distance
# flip_distance(N,codewords) = distance
# code_distance(N,codewords) = c_distance
# code_classification(N,codewords) = c_classification
# errors_list(N,t) = errors_literal, affected_qubits
# correctability_degree(N,codewords) = degree, undetectable_faults, undecidable_faults

def syndrome_table(N,P):
    
    s_table = []
    
    for i in range(3*N):
        
        s_table.append([])
    
    for i in range(N):
        
        X_error = ['I']*N
        Z_error = ['I']*N
        Y_error = ['I']*N
        
        X_error[i] = 'X'
        Z_error[i] = 'Z'
        Y_error[i] = 'Y'
        
        X_syndrome = ''
        Z_syndrome = ''
        Y_syndrome = ''
        
        for j in range(len(P)):
            
            if uf.pauli_commutation(X_error,P[j]) == -1:
                X_syndrome += '1'
            else:
                X_syndrome += '0'
                
            if uf.pauli_commutation(Z_error,P[j]) == -1:
                Z_syndrome += '1'
            else:
                Z_syndrome += '0'
                
            if uf.pauli_commutation(Y_error,P[j]) == -1:
                Y_syndrome += '1'
            else:
                Y_syndrome += '0'
                
        s_table[2*i] = [X_error,i,X_syndrome]
        s_table[2*i+1] = [Z_error,i,Z_syndrome]
        s_table[3*i+2] = ['Y',i,Y_syndrome]
        
    return s_table

def syndrome_table_t_errors(N,P,errors_literal,affected_qubits):
    
    
    s_table = []

    for i in range(len(errors_literal)):
        
        #if (tuple(errors_literal[i]) in common_stabilizers_set) == False: #error is not a common stabilizer
        if 1 == 1:
            syndrome = ''
            
            for j in range(len(P)):
                
                if uf.pauli_commutation(errors_literal[i],P[j]) == -1:
                    syndrome += '1'
                else:
                    syndrome += '0'
                    
                    
            s_table.append([errors_literal[i],affected_qubits[i],syndrome])

        
    return s_table

def common_stabilizers(N,stab_set):

    #INPUT: a list of the codewords circuits
    #OUTPUT: a list of common stabilizers of the codewords (the stabilizers must have the same signal to be considered common)

    #computing the group of stabilizers for each codeword
            
    common_stabilizers = stab_set[0].intersection(stab_set[1])
    
    common_stabilizers = list(common_stabilizers)

    common_stabilizers_list = [list(common_stabilizers[i][:N]) for i in range(len(common_stabilizers))]
    

    return common_stabilizers_list

def syndrome_operators(N,stab_set):
        
    common_stabilizers = stab_set[0].intersection(stab_set[1])
    common_stabilizers = list(common_stabilizers)
    common_stabilizers = [common_stabilizers[x] for x in range(len(common_stabilizers)) if common_stabilizers[x][-1] == '1']
          
    #computing the generators of common_commuting_stabilizers
    mat = []
    
    for i in range(len(common_stabilizers)):
        
        row = []
        
        for j in range(N):
            
            if common_stabilizers[i][j] == 'X':
                row.append(1)
                row.append(0)
            elif common_stabilizers[i][j] == 'Z':
                row.append(0)
                row.append(1)
            elif common_stabilizers[i][j] == 'Y':
                row.append(1)
                row.append(1)
            elif common_stabilizers[i][j] == 'I':
                row.append(0)
                row.append(0)
                
        mat.append(row)
        
    I = uf.rank_mod2(np.array(mat))
    
    generators = [mat[0]]
    rank_generators = 1
    
    for i in range(1,len(mat)):
        generators.append(mat[i])
        new_rank = uf.rank_mod2(np.array(generators))
        if new_rank == rank_generators: #rank didn't raise
            del(generators[-1])
        else:
            rank_generators = uf.rank_mod2(np.array(generators))
        
        if rank_generators == I:
            break
            
    #translating back to literals and forming the output
    P = []
    
    for i in range(len(generators)):
        P.append([])
    
    for i in range(len(generators)):
        syndrome_operator = []
        
        for j in range(N):
            syndrome_operator.append([])
        
        for j in range(N):
            if generators[i][2*j] == 0 and generators[i][2*j+1] == 0:
                syndrome_operator[j] = 'I'
            elif generators[i][2*j] == 0 and generators[i][2*j+1] == 1:
                syndrome_operator[j] = 'Z'
            elif generators[i][2*j] == 1 and generators[i][2*j+1] == 0:
                syndrome_operator[j] = 'X'
            elif generators[i][2*j] == 1 and generators[i][2*j+1] == 1:
                syndrome_operator[j] = 'Y'
        P[i] = syndrome_operator

    P_bit = generators    

    return P, P_bit

def syndrome_operators_topo(N,stab_set):

    #INPUT: a list of the codewords circuits
    #OUTPUT: a list of syndrome operators which 1) stabilizes all codewords
    #2) they commute with one another; 3) they are LI

    
    
    #computing a set of common stabilizers between the groups
    
    common_stabilizers = stab_set[0].intersection(stab_set[1])
    
    common_stabilizers = list(common_stabilizers)
    
    common_stabilizers = [common_stabilizers[x] for x in range(len(common_stabilizers)) if common_stabilizers[x][-1] == '1']
    
    notas = np.zeros(len(common_stabilizers))
    
    for i in range(len(common_stabilizers)):
        
        len_x = len(uf.find_literal("X",np.array(common_stabilizers[i])))
        len_y = len(uf.find_literal("Y",np.array(common_stabilizers[i])))
        len_z = len(uf.find_literal("Z",np.array(common_stabilizers[i])))
        
        if (len_y == 0 and len_z == 0) or (len_y == 0 and len_x == 0):
            
            notas[i] = 1
        
        elif len_y == max(len_x,len_y,len_z):
            
            notas[i] = (N - len_y)/N
        
        else:
            
            notas[i] = (N-min(len_x,len_z))/N
            
    css_group_degree = np.mean(notas)

    
    stab_only_X = []
    stab_only_Y = []
    stab_only_Z = []
    stab_XYZ = []
    set_X = []
    set_Y = []
    set_Z = []
    
    
    
    for i in range(len(common_stabilizers)):
        
        #only X
        if len(uf.find_literal("Y",np.array(common_stabilizers[i]))) == 0 and len(uf.find_literal("Z",np.array(common_stabilizers[i]))) == 0:
            
            stab_only_X.append([common_stabilizers[i],len(uf.find_literal("X",np.array(common_stabilizers[i])))])
            set_X.append(tuple(uf.find_literal("X",np.array(common_stabilizers[i]))))
            
        #only Y 
        elif len(uf.find_literal("X",np.array(common_stabilizers[i]))) == 0 and len(uf.find_literal("Z",np.array(common_stabilizers[i]))) == 0:
            
            stab_only_Y.append([common_stabilizers[i],len(uf.find_literal("Y",np.array(common_stabilizers[i])))])
            set_Y.append(tuple(uf.find_literal("Y",np.array(common_stabilizers[i]))))
            
            
        #only Z
        elif len(uf.find_literal("Y",np.array(common_stabilizers[i]))) == 0 and len(uf.find_literal("X",np.array(common_stabilizers[i]))) == 0:
            
            stab_only_Z.append([common_stabilizers[i],len(uf.find_literal("Z",np.array(common_stabilizers[i])))])
            set_Z.append(tuple(uf.find_literal("Z",np.array(common_stabilizers[i]))))
            
        else:
            
            stab_XYZ.append([common_stabilizers[i],len(uf.find_literal("X",np.array(common_stabilizers[i])))+len(uf.find_literal("Y",np.array(common_stabilizers[i])))+len(uf.find_literal("Z",np.array(common_stabilizers[i])))])
            
   
    
    stab_only_X = sorted(stab_only_X, key=lambda x: x[1])
    
    stab_only_Y = sorted(stab_only_Y, key=lambda x: x[1])
   
    stab_only_Z = sorted(stab_only_Z, key=lambda x: x[1])
   
    stab_XYZ = sorted(stab_XYZ, key=lambda x: x[1])
 
    
    stab_sorted = []
    
    for i in range(len(stab_only_X)):
        stab_sorted.append(stab_only_X[i][0])
    for i in range(len(stab_only_Z)):
        stab_sorted.append(stab_only_Z[i][0])      
    for i in range(len(stab_only_Y)):
        stab_sorted.append(stab_only_Y[i][0])
    for i in range(len(stab_XYZ)):
        stab_sorted.append(stab_XYZ[i][0])
        
    
    if len(set_X) == 0 or len(set_Z) == 0:
        
        notas = 0
    
    elif len(set_X) <= len(set_Z):
        
        notas = np.zeros(len(set_X))
        
        for i in range(len(set_X)):
            
            notas_temp = np.zeros(len(set_Z))
            x_temp = set(set_X[i])
            
            for j in range(len(set_Z)):
                
                z_temp = set(set_Z[j])
                
                int_degree = len(x_temp.intersection(z_temp))
                notas_temp[j] = int_degree
                
            notas[i] = np.max(notas_temp)
            
    else:
        
        notas = np.zeros(len(set_Z))
        
        for i in range(len(set_Z)):
            
            notas_temp = np.zeros(len(set_X))
            z_temp = set(set_Z[i])
            
            for j in range(len(set_X)):
                
                x_temp = set(set_X[j])
                
                int_degree = len(z_temp.intersection(x_temp))
                notas_temp[j] = int_degree
                
            notas[i] = np.max(notas_temp)
    
    topo_degree = 1+np.mean(notas)


    #computing the generators of common_commuting_stabilizers
    mat = []
    
    for i in range(len(stab_sorted)):
        
        row = []
        
        for j in range(N):
            
            if stab_sorted[i][j] == 'X':
                row.append(1)
                row.append(0)
            elif stab_sorted[i][j] == 'Z':
                row.append(0)
                row.append(1)
            elif stab_sorted[i][j] == 'Y':
                row.append(1)
                row.append(1)
            elif stab_sorted[i][j] == 'I':
                row.append(0)
                row.append(0)
                
        mat.append(row)
        
    I = uf.rank_mod2(np.array(mat))
    
    generators = [mat[0]]
    rank_generators = 1
    
    for i in range(1,len(mat)):
        generators.append(mat[i])
        new_rank = uf.rank_mod2(np.array(generators))
        if new_rank == rank_generators: #rank didn't raise
            del(generators[-1])
        else:
            rank_generators = uf.rank_mod2(np.array(generators))
        
        if rank_generators == I:
            break
        
    
    #translating back to literals and forming the output
    P = []
    
    for i in range(len(generators)):
        P.append([])
    
    for i in range(len(generators)):
        syndrome_operator = []
        
        for j in range(N):
            syndrome_operator.append([])
        
        for j in range(N):
            if generators[i][2*j] == 0 and generators[i][2*j+1] == 0:
                syndrome_operator[j] = 'I'
            elif generators[i][2*j] == 0 and generators[i][2*j+1] == 1:
                syndrome_operator[j] = 'Z'
            elif generators[i][2*j] == 1 and generators[i][2*j+1] == 0:
                syndrome_operator[j] = 'X'
            elif generators[i][2*j] == 1 and generators[i][2*j+1] == 1:
                syndrome_operator[j] = 'Y'
        P[i] = syndrome_operator

    P_bit = generators    
    
    css_degre = css_degree(P,N)

    return P, P_bit, topo_degree, css_group_degree, css_degre

# def css_degree(P):
    
#     faulty_stabilizers = 0
    
#     for i in range(len(P)):
        
#         if len(uf.find_literal("X",np.array(P[i]))) > 0 and len(uf.find_literal("Z",np.array(P[i]))) > 0:
            
#             faulty_stabilizers += 1
            
#         elif len(uf.find_literal("Y",np.array(P[i]))) > 0:
            
#             faulty_stabilizers += 1
            
#     degree = (len(P)-faulty_stabilizers)/len(P)
    
#     return degree

def css_degree(P,N):
    
    notas = np.zeros(len(P))
    
    for i in range(len(P)):
        
        len_x = len(uf.find_literal("X",np.array(P[i])))
        len_y = len(uf.find_literal("Y",np.array(P[i])))
        len_z = len(uf.find_literal("Z",np.array(P[i])))
        
        if (len_y == 0 and len_z == 0) or (len_y == 0 and len_x == 0):
            
            notas[i] = 1
        
        # elif len_y == max(len_x,len_y,len_z):
            
        #     notas[i] = (N - len_y)/N
        
        # else:
            
        #     notas[i] = (N-min(len_x,len_z))/N
            
    degree = np.mean(notas)
    
    return degree



# def locality_degree(P,lattice_graph):
    
#     local_stabilizers = 0
    
#     for i in range(len(P)):
        
#         nodes = list(uf.find_literal("X",np.array(P[i]))+1)+list(uf.find_literal("Y",np.array(P[i]))+1)+list(uf.find_literal("Z",np.array(P[i]))+1)
#         nodes_permut = list(itertools.permutations(nodes))
        
#         for j in range(len(nodes_permut)):
        
#             path_test = nx.is_simple_path(lattice_graph, list(nodes_permut[j]))
            
#             if path_test == True:
                
#                 local_stabilizers += 1
#                 break
            
#     degree = local_stabilizers/len(P)
    
#     return degree

def phase_distance(N,codewords,stab_set):
    
    #evaluating phase distance
    
        
    common_stabilizers = stab_set[0]
    dif_stabilizers = stab_set[0]
    
    for i in range(1,len(stab_set)):
        
        common_stabilizers = common_stabilizers.intersection(stab_set[i])
        dif_stabilizers = dif_stabilizers.union(stab_set[i])
     
    dif_stabilizers = list(dif_stabilizers - common_stabilizers)
    
    error_weight = [0]*len(dif_stabilizers)
    
    for i in range(len(dif_stabilizers)):
        for j in range(N):
            if dif_stabilizers[i][j] != 'I':
                error_weight[i] += 1
                
    phase_distance = min(error_weight)
    
    return phase_distance


def errors_list(N,t):
    
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

def correctability_degree(N,t,errors_literal,affected_qubits,stab_set):
    
    #P, P_bit,topo_degree,css_group_degree, css_degre = syndrome_operators_topo(N,stab_set)
    P, P_bit = syndrome_operators(N,stab_set)
    common_stab = common_stabilizers(N,stab_set)
    
    s_table = syndrome_table_t_errors(N,P,errors_literal,affected_qubits)
    
    P_set = set()
    for i in range(len(P)):
        P_set.add(tuple(P[i]))

    
    
    code_set = set()
        
    for i in range(len(s_table)):
        code_set.add(s_table[i][2])
            
    code_list = list(code_set)
    
    undetectable_faults = 0
    undecidable_faults = 0
    
    #counting lines that make the error undetectable
    for i in range(len(s_table)):
        
        if s_table[i][2] == '0'*len(P):
            
            undetectable_faults += 1
        
    #counting the syndromes that make the correction undecidable
    errors_per_code = []
    
    for i in range(len(code_list)):
        
        errors_per_code.append([])
    
    for i in range(len(code_list)):
        errors_per_code[i] = [code_list[i] , [s_table[x][0] for x in range(len(s_table)) if s_table[x][2] == code_list[i] ]]
    
    for i in range(len(errors_per_code)):
        
        if len(errors_per_code[i][1]) != 1 and errors_per_code[i][0] != '0'*len(P):
        
            errors = list(itertools.combinations(list(range(len(errors_per_code[i][1]))), 2))
            
            for j in range(len(errors)):
                
                ans,signal = uf.pauli_multi_list(errors_per_code[i][1][errors[j][0]],errors_per_code[i][1][errors[j][1]])
                
                if (ans in common_stab) == False:
                    undecidable_faults += 1
                    break

    c_degree = float(len(s_table) + len(code_list) - undetectable_faults - undecidable_faults)/float(len(s_table) + len(code_list))
        
    #return c_degree, topo_degree, css_group_degree, css_degre
    return c_degree