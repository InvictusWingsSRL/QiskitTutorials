from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit.primitives import BackendSampler
#from qiskit.extensions.simulator import snapshot
#from qiskit.tools.visualization import circuit_drawer
import numpy as np
import math as m
import scipy as sci
import random
import time
import matplotlib
import matplotlib.pyplot as plt
import ipywidgets as wd
S_simulator = Aer.backends(name='statevector_simulator')[0]
M_simulator = Aer.backends(name='qasm_simulator')[0]

S_simulator = Aer.backends(name='statevector_simulator')[0]
M_simulator = Aer.backends(name='qasm_simulator')[0] 

def execute(circuit, backend, **kwargs):
    s = 1024
    if 'shots' in kwargs:
        s = int( kwargs['shots'] )    
    new_circuit = transpile(circuit, backend)
    return backend.run(new_circuit, shots = s)

#Displaying Results 
def Wavefunction( obj , *args, **kwargs):
#Displays the wavefunction of the quantum system 
	if(type(obj) == QuantumCircuit ):
		statevec = execute( obj, S_simulator, shots=1 ).result().get_statevector()
	if(type(obj) == np.ndarray):
		statevec = obj
	sys = False
	NL = False
	dec = 5
	if 'precision' in kwargs:
		dec = int( kwargs['precision'] )
	if 'column' in kwargs:
		NL = kwargs['column']
	if 'systems' in kwargs:
		systems = kwargs['systems']
		sys = True
		last_sys = int(len(systems)-1)
		show_systems = []
		for s_chk in range(len(systems)):
			if( type(systems[s_chk]) != int ):
				raise Exception('systems must be an array of all integers')
		if 'show_systems' in kwargs:
			show_systems = kwargs['show_systems']
			if( len(systems)!= len(show_systems) ):
				raise Exception('systems and show_systems need to be arrays of equal length')
			for ls in range(len(show_systems)):
				if((show_systems[ls] != True) and (show_systems[ls] != False)):
					raise Exception('show_systems must be an array of Truth Values')
				if(show_systems[ls] == True):
					last_sys = int(ls) 
		else:
			for ss in range(len(systems)):
				show_systems.append(True)
	wavefunction = ''
	qubits = int(m.log(len(statevec),2))
	for i in range(int(len(statevec))):
		#print(wavefunction)
		value = round(statevec[i].real, dec) + round(statevec[i].imag, dec) * 1j
		if( (value.real != 0) or (value.imag != 0)):
			state = list(Binary(int(i),int(2**qubits)))
			state.reverse()
			state_str = ''
			#print(state)
			if( sys == True ): #Systems and SharSystems 
				k = 0 
				for s in range(len(systems)):
					if(show_systems[s] == True):
						if(int(s) != last_sys):
							state.insert(int(k + systems[s]), '>|' )
							k = int(k + systems[s] + 1)
						else:
							k = int(k + systems[s]) 
					else:
						for s2 in range(systems[s]):
							del state[int(k)]
			for j in range(len(state)):
				if(type(state[j])!= str):
					state_str = state_str + str(int(state[j])) 
				else:
					state_str = state_str + state[j]
			#print(state_str)
			#print(value)
			if( (value.real != 0) and (value.imag != 0) ):
				if( value.imag > 0):
					wavefunction = wavefunction + str(value.real) + '+' + str(value.imag) + 'j |' + state_str + '>   '
				else:
					wavefunction = wavefunction + str(value.real) + '' + str(value.imag) + 'j |' + state_str +  '>   '
			if( (value.real !=0 ) and (value.imag ==0) ):
				wavefunction = wavefunction  + str(value.real) + '  |' + state_str + '>   '
			if( (value.real == 0) and (value.imag != 0) ):
				wavefunction = wavefunction + str(value.imag)  + 'j |' + state_str + '>   '
			if(NL):
				wavefunction = wavefunction + '\n'
		#print(NL)
	
	print(wavefunction)
	return wavefunction


def Measurement(quantumcircuit, *args, **kwargs): 
	#Displays the measurement results of a quantum circuit 
	p_M = True
	S = 1
	ret = False
	NL = False
	if 'shots' in kwargs:
		S = int(kwargs['shots'])
	if 'return_M' in kwargs:
		ret = kwargs['return_M']
	if 'print_M' in kwargs:
		p_M = kwargs['print_M']
	if 'column' in kwargs:
		NL = kwargs['column']
	M1 = execute(quantumcircuit, M_simulator, shots=S).result().get_counts(quantumcircuit)
	M2 = {}
	k1 = list(M1.keys())
	v1 = list(M1.values())
	for k in range(len(k1)):
		key_list = list(k1[k])
		new_key = ''
		for j in range(len(key_list)):
			new_key = new_key+key_list[len(key_list)-(j+1)]
		M2[new_key] = v1[k]
	if(p_M):
		k2 = list(M2.keys())
		v2 = list(M2.values())
		measurements = ''
		for i in range(len(k2)):
			m_str = str(v2[i])+'|'
			for j in range(len(k2[i])):
				if(k2[i][j] == '0'):
					m_str = m_str + '0' 
				if(k2[i][j] == '1'):
					m_str = m_str + '1'
				if( k2[i][j] == ' ' ):
					m_str = m_str +'>|'
			m_str = m_str + '>   '
			if(NL):
				m_str = m_str + '\n'
			measurements = measurements + m_str
		print(measurements)
	if(ret):
		return M2

    
def Most_Probable(M,N):
    ''' 
    Input: M (Dictionary) N (integer)
    Returns the N most probable states accoding to the measurement counts stored in M 
    '''
    count = []
    state = []
    if( len(M) < N ):
        N = len(M)
    for k in range(N):
        count.append(0)
        state.append(0)
    for m in range(len(M)):
        new = True
        for n in range(N):
            if( (list(M.values())[m] > count[n]) and new ):
                for i in range( N-(n+1)):
                    count[-int(1+i)] = count[-int(1+i+1)]
                    state[-int(1+i)] = state[-int(1+i+1)]
                count[int(n)] = list(M.values())[m]
                state[int(n)] = list(M.keys())[m]
                new = False
    return count,state 

#Math Operations
def Oplus(bit1,bit2): 
	'''Adds too bits of O's and 1's (modulo 2)'''
	bit = np.zeros(len(bit1))
	for i in range( len(bit) ):
		if( (bit1[i]+bit2[i])%2 == 0 ):
			bit[i] = 0
		else: 
			bit[i] = 1
	return bit 


def Binary(number,total): 
#Converts a number to binary, right to left LSB 152 153 o
	qubits = int(m.log(total,2))
	N = number
	b_num = np.zeros(qubits)
	for i in range(qubits):
		if( N/((2)**(qubits-i-1)) >= 1 ):
			b_num[i] = 1
			N = N - 2 ** (qubits-i-1)
	B = [] 
	for j in range(len(b_num)):
		B.append(int(b_num[j]))
	return B

def BinaryL(number,total): 
    B = Binary(number, total)
    B.reverse()
    return B


def From_Binary(s):
    num = 0
    for i in range(len(s)):
        num = num + int(s[-(i+1)]) * 2**i
    return num



def From_BinaryLSB(S, LSB):
    num = 0
    for i in range(len(S)):
        if(LSB=='R'):
            num = num + int(S[-(i+1)]) * 2**i
        elif(LSB=='L'):
            num = num + int(S[i]) * 2**i
    return num

def B2D(in_bi):
    len_in = len(in_bi)
    in_bi = in_bi[::-1]
    dec = 0
    for i in range(0,len_in):
        if in_bi[i] != '0':
            dec += 2**i
    return dec

#  Custom Gates
def X_Transformation(qc, qreg, state): 
#Tranforms the state of the system, applying X gates according to as in the vector 'state' 
	for j in range(len(state)):
		if( int(state[j]) == 0 ):
			qc.x( qreg[int(j)] ) 


def n_NOT(qc, control, target, anc): 
#performs an n-NOT gate
	n = len(control)
	instructions = []
	active_ancilla = []
	q_unused = []
	q = 0
	a = 0
    
	while(n > 0):
		if(n >= 2):
			instructions.append( [control[q], control[q+1], anc[a]] )
			active_ancilla.append(a)
			a += 1
			q += 2
			n -= 2
		if(n == 1):
			q_unused.append(q)
			n -= 1
            
	while (len(q_unused) != 0):
		if(len(active_ancilla)!=1):
			instructions.append( [control[q], anc[active_ancilla[0]], anc[a]] )
			del active_ancilla[0]
			del q_unused[0]
			active_ancilla.append(a)
			a += 1
		else:
			instructions.append( [control[q], anc[active_ancilla[0]], target] )
			del active_ancilla[0]
			del q_unused[0]
                
	while(len(active_ancilla) != 0):
		if( len(active_ancilla) > 2 ):
			instructions.append( [anc[active_ancilla[0]], anc[active_ancilla[1]], anc[a]] )
			active_ancilla.append(a)
			del active_ancilla[0]
			del active_ancilla[0]
			a += 1
		elif( len(active_ancilla) == 2):
			instructions.append([anc[active_ancilla[0]], anc[active_ancilla[1]], target])
			del active_ancilla[0]
			del active_ancilla[0]
		elif( len(active_ancilla) == 1):
			instructions.append([anc[active_ancilla[0]], target])
			del active_ancilla[0]
        
                
	for i in range( len(instructions) ):
		if len(instructions[i]) == 2:
			qc.cx( instructions[i][0], instructions[i][1])
		else:
			qc.ccx( instructions[i][0], instructions[i][1], instructions[i][2] )
        
	del instructions[-1]
	
	for i in range( len(instructions) ):
		qc.ccx( instructions[0-(i+1)][0], instructions[0-(i+1)][1], instructions[0-(i+1)][2] )



def Control_Instruction( qc, vec ): 
#Ammends the proper quantum circuit instruction based on the input 'vec'
#Used for the function 'n_Control_U
	if( vec[0] == 'X' ):
		qc.cx( vec[1], vec[2] )
	elif( vec[0] == 'Z' ):
		qc.cz( vec[1], vec[2] )
	elif( vec[0] == 'CPHASE' ):
		qc.cu1( vec[2], vec[1], vec[3] )
	elif( vec[0] == 'SWAP' ):
		qc.cswap( vec[1], vec[2], vec[3] ) 


def sinmons_solver(E,N):
	'''Returns an array of s_prime candidates
	'''
	s_primes = []
	for s in np.ararge(1,2**N):
		sp = Binary( int(s), 2**N )
		candidate = True
		for e in range( len(E) ):
			value = 0
			for i in range( N ):
				value = value + sp[i]*E[e][i]
			if(value%2==1):
				candidate=False
		if(candidate):
			s_primes.append(sp)
	return s_primes


def Grover_Oracle(mark, qc, q, an1, an2): 
	'''
	picks out the marked state and applies a negative phase 
	'''
	qc.h( an1[0] )
	X_Transformation(qc, q, mark)
	if( len(mark) > 2 ):
		n_NOT( qc, q, an1[0], an2 )
	elif( len(mark) == 2 ):
		qc.ccx( q[0], q[1], an1[0] )
	X_Transformation(qc, q, mark)
	qc.h( an1[0] )

def Grover_Diffusion(mark, qc, q, an1, an2): 
	'''
	ammends the instructions for a Grover Diffusion Operation to the Quantum Circuit
	'''
	zeros_state = []
	for i in range( len(mark) ):
		zeros_state.append( 0 )
		qc.h( q[int(i)] )
	Grover_Oracle(zeros_state, qc, q, an1, an2)
	for j in range( len(mark) ):
		qc.h( q[int(j)] )



def Grover(Q, marked): 
	'''
	Amends all the instructions for a Grover Search 
	'''
	q = QuantumRegister(Q,name='q')
	an1 = QuantumRegister(1,name='anc')
	an2 = QuantumRegister(Q-2,name='nanc')
	c = ClassicalRegister(Q,name='c')
	qc = QuantumCircuit(q,an1,an2,c,name='qc')
	for j in range(Q):
		qc.h( q[int(j)] )
	qc.x( an1[0] )

	iterations = round( m.pi/4 * 2**(Q/2.0) )
	for i in range( iterations ):
		Grover_Oracle(marked, qc, q, an1, an2)
		Grover_Diffusion(marked, qc, q, an1, an2)

	return qc, q, an1, an2, c

def Multi_Grover(q, a1, a2, qc, marked, iters):
    '''
    Input: q (QuantumRegister) a1 (QuantumRegister) a2 (QuantumRegister) qc (QuantumCircuit)
    marked (array) iters (integer)
    Appends all of the gate operations for a multi-marked state Grover Search
    '''
    Q = int(len(marked))
    for i in np.arange( iters ):
        for j in np.arange(len(marked)):
            M = list(marked[j])
            for k in np.arange(len(M)):
                if(M[k]=='1'):
                    M[k] = 1
                else:
                    M[k] = 0
            Grover_Oracle(M, qc, q, a1, a2)
        Grover_Diffusion(M, qc, q, a1, a2)
    return qc, q, a1, a2


def n_Control_U(qc, control, anc, gates):
#Performs a list of single control gates, as an n-control operation

	instructions = []
	active_ancilla = []
	q_unused = []
	n = len(control)
	q = 0
	a = 0
        
	while(n > 0):
		if(n >= 2) :
			instructions.append([control[q], control[q+1], anc[a]])
			active_ancilla.append(a)
			a += 1
			q += 2
			n -= 2
		if(n == 1):
			q_unused.append( q )
			n -= 1
            
	while( len(q_unused) != 0 ) :
		if(len(active_ancilla)>1):
			instructions.append( [control[q] , anc[active_ancilla[0]], anc[a]])
			del active_ancilla[0]
			del q_unused[0]
			active_ancilla.append(a)
			a += 1
		else:
			instructions.append( [control[q] , anc[active_ancilla[0]], anc[a]]) 
			del active_ancilla[0]
			del q_unused[0]
			c_a = anc[a]
            
	while( len(active_ancilla) != 0 ) :
		if( len(active_ancilla) > 2 ) :
			instructions.append([anc[active_ancilla[0]], anc[active_ancilla[1]], anc[a]])
			active_ancilla.append(a)
			del active_ancilla[0]
			del active_ancilla[0]
			a += 1
		elif( len(active_ancilla)==2):
			instructions.append([anc[active_ancilla[0]], anc[active_ancilla[1]], anc[a]])
			del active_ancilla[0]
			del active_ancilla[0]
			c_a = anc[a]
		elif( len(active_ancilla)==1):
			c_a = anc[active_ancilla[0]]
			del active_ancilla[0]
                
	for i in range( len(instructions) ) :
		qc.ccx(instructions[i][0], instructions[i][1], instructions[i][2])
        
	for j in range(len(gates)):
		control_vec = [gates[j][0], c_a]
		for k in range( 1, len(gates[j])):
			control_vec.append( gates[j][k] )
		Control_Instruction( qc, control_vec )

	for i in range( len(instructions) ) :
		qc.ccx(instructions[0-(i+1)][0],instructions[0-(i+1)][1], instructions[0-(i+1)][2])


        
def Control_Instructions(qc, vec):
	if (vec[0] == 'X'):
		qc.cx(vec[1], vec[2])
	elif (vec[0] == 'Z'):
		qc.cz(vec[1], vec[2])


        
def Blackbox_g_D(qc, qreg):
	f_type=['f(0,1) -> (0,1)', 'f(0,1) -> (1,0)', 'f(0,1) -> 0', 'f(0,1) -> 1']
	r = int(m.floor(4*np.random.rand()))
	if (r == 0):
		qc.cx(qreg[0],qreg[1])
	if (r == 1):
		qc.x(qreg[0])
		qc.cx(qreg[0],qreg[1])
		qc.x(qreg[0])
	if (r == 2):
		qc.id(qreg[0])
		qc.id(qreg[1])
	if (r == 3):
		qc.x(qreg[1])
        
	return f_type[r]


def Deutsch(qc,qreg):
    qc.h(qreg[0])
    qc.h(qreg[1])
    f = Blackbox_g_D(qc, qreg)
    qc.h(qreg[0])
    qc.h(qreg[1])
    
    return f

def Blackbox_g_DJ(Q, qc, qreg, an1):
    f_type=['constant','balanced']
    f=[]
    r=int(m.floor(2**Q*np.random.rand()))
    
    if r==0:
        f.append(f_type[0])
    elif r == 1:
        qc.x(qreg[Q-1])
        f.append(f_type[0])
    else:
        control = []
        for i in range(Q):
            control.append(qreg[i])
        
        an2 = QuantumRegister(int(Q-2), name='nn_anc')
        qc.add_register(an2)
        f.append(f_type[1])
        
        S=[]
        for s in range(2**Q):
            S.append(s)
        
        for k in range(2**(Q-1)):
            S_num = S[int(m.floor(len(S)*np.random.rand()))]
            state = Binary(S_num,2**Q)
            S.remove(S_num)
            
            f_string = '|'
            for j in range(len(state)):
                f_string += str(int(state[j]))
                if (state[j] == 0):
                    qc.x(qreg[j])
            f.append(f_string + '>')
            
            n_NOT(qc, control, an1[0], an2)
            
            for j in range(len(state)):
                if (state[j] == 0):
                    qc.x(qreg[j])
            
    return f

def Deutsch_Josza(Q,qc,qreg,an1):
    for i in range(Q):
        qc.h(qreg[i])
    qc.h(an1[0])
    f = Blackbox_g_DJ(Q, qc, qreg, an1)
    for i in range(Q):
        qc.h(qreg[i])
    qc.h(an1[0])
    return f

def Blackbox_g_BV(Q,qc,qreg,an1):
    a = Binary(int(m.floor(2**Q*np.random.rand())),2**Q)
    control=[]
    for i in range(Q):
        control.append(qreg[i])
        
    an2 = QuantumRegister(Q-2,name='nn_anc')
    qc.add_register(an2)
    
    for s in range(2**Q):
        state = Binary(s,2**Q)
        dp = np.vdot(a, state)
        
        if (dp % 2 == 1):
            for j in range(len(state)):
                if int(state[j]) == 0:
                    qc.x(qreg[j])

            n_NOT(qc, control, an1[0], an2)
            
            for j in range(len(state)):
                if int(state[j]) == 0:
                    qc.x(qreg[j])
    return a

def Bernstein_Vazirani(Q,qc,qreg,an1):
    for i in range(Q):
        qc.h(qreg[i])
    qc.h(an1[0])
    a = Blackbox_g_BV(Q,qc,qreg,an1)
    for i in range(Q):
        qc.h(qreg[i])
    qc.h(an1[0])
    return a

def Blackbox_g_S(Q, qc, q, anc1):
    anc2 = QuantumRegister(Q-1,name='nU_anc')
    qc.add_register(anc2)
    s = np.zeros(Q)
    for i in range(Q):
        s[i] = int(m.floor(2*np.random.rand()))
        
    outputs=[]
    for o in range(2**Q):
        outputs.append(o)
    
    f = np.zeros(2**Q)
    for j in range(2**Q):
        out = outputs[int(m.floor(len(outputs)*np.random.rand()))]
        f[j] = out
        f[int(From_Binary(Oplus(Binary(j, 2**Q),s)))] = out
        outputs.remove(out)
        
    output_states=[]
    for k in range(2**Q):
        output_states.append(Binary(f[k],2**Q))
        
    for a in range(2**Q):
        c_ops=[]
        for b in range(Q):
            if output_states[a][b] == 1:
                c_ops.append(['X', anc1[b]])
                
        X_Transformation(qc, q, Binary(a, 2**Q))
        n_Control_U(qc, q, anc2, c_ops)
        # instead of n_Control_U it would work witn n_NOT as well, but the overhead would be much higher:
        #for b in range(Q):
        #    if output_states[a][b] == 1:
        #        n_NOT(qc, q, anc1[b], anc2)

        X_Transformation(qc, q, Binary(a, 2**Q))


    return qc, s, f

def Simons_Quantum(Q, qc, q, c, anc1):
    for i in range(Q):
        qc.h(q[i])   
    qc,s,f = Blackbox_g_S(Q,qc,q, anc1)
    for i in range(Q):
        qc.h(q[i])
    qc.measure(q,c)
    return qc, s

def Simons_Solver(E,N):
    s_primes = []
    for s in range(1, 2**N):
        sp = Binary(s, 2**N)
        candidate = True
        for e in range(len(E)):
            value = 0
            for i in range(N):
                value += sp[i] * E[e][i]
            if value%2 == 1:
                candidate = False
        if candidate:
            s_primes.append(sp)
    return s_primes


def Simons_Classical(Q, qc):
    run_quantum = True
    Equations = []
    Results = []
    quantum_runs = 0
    
    while(run_quantum):
        quantum_runs += 1
        M = Measurement(qc, shots = 20, return_M = True, print_M = False)
        new_result = True
        
        for r in range(len(Results)):
            if list(M.keys())[0] == Results[r]:
                new_result = False
                break
        
        if new_result:
            Results.append(list(M.keys())[0])
            eq = []
            for e in range(Q):
                eq.append(int(list(M.keys())[0][e]))
            Equations.append(eq)
            s_primes = Simons_Solver(Equations, Q)
            if len(s_primes) == 1:
                run_quantum = False
    return s_primes, Results, quantum_runs

def DFT(x, **kwargs):
    p = -1.0
    if 'inverse' in kwargs:
        P = kwargs['inverse']
        if P == True:
            p = 1.0
    L = len(x)
    X = []
    for i in range(L):
        value = 0
        for j in range(L):
            value += x[j] * np.exp(p * 2 * m.pi * 1.0j * i * j / L)
        X.append(value)
    for k in range(len(X)):
        re = round(X[k].real,5)
        im = round(X[k].imag,5)
        if abs(im) == 0 and abs(re) != 0:
            X[k] = re
        elif abs(re) == 0 and abs(im) != 0:
            X[k] = im * 1.0j
        elif abs(re) == 0 and abs(im) == 0:
            X[k] = 0
        else:
            X[k] = re + im * 1.0j
    return X

def QFT(qc, q, qubits, **kwargs):
    R_phis = [0]
    for i in range(2, qubits+1):
        R_phis.append( 2/(2**i) * m.pi )
    for j in range(qubits):
        qc.h( q[j] )
        for k in range(qubits-j-1):
            qc.cp( R_phis[k+1], q[j+k+1], q[j] )
    if 'swap' in kwargs:
        if(kwargs['swap'] == True):
            for s in range(m.floor(qubits/2.0)):
                qc.swap( q[s],q[qubits-1-s] )
                
def QFT_dgr(qc, q, qubits, **kwargs):
    if 'swap' in kwargs:
        if(kwargs['swap'] == True):
            for s in range(m.floor(qubits/2.0)):
                qc.swap( q[s],q[qubits-1-s] )
                
    R_phis = [0]
    for i in range(2,qubits+1):
        R_phis.append( -2/(2**i) * m.pi )
    for j in range(qubits):
        for k in range(j):
            qc.cp(R_phis[j-k], q[qubits-k-1], q[qubits-j-1] )
        qc.h( q[qubits-j-1] )
        
def Quantum_Adder(qc, Qa, Qb, A, B):
    Q = len(B)
    for n in range(Q):
        if( A[n] == 1 ):
            qc.x( Qa[n+1] )
        if( B[n] == 1 ):
            qc.x( Qb[n] )
    QFT(qc,Qa,Q+1)
    p = 1
    for j in range( Q ):
        qc.cp( m.pi/(2**p), Qb[j], Qa[0] )
        p = p + 1
    for i in range(1,Q+1):
        p = 0
        for jj in np.arange( i-1, Q ):
            qc.cp( m.pi/(2**p), Qb[jj], Qa[i] )
            p = p + 1
    QFT_dgr(qc,Qa,Q+1)
    

def QPE_phi(MP):
    ms = [[],[]]
    for i in range(2):
        for j in range(len(MP[1][i])):
            ms[i].append(int(MP[1][i][j]))
    n = int(len(ms[0]))
    MS1 = From_Binary(ms[0])
    MS2 = From_Binary(ms[1])
    
    estimatedProb = MP[0][0]
    aproxPhi = 0
    aproxProb = 1
    for k in np.arange(1,5000):
        phi = k/5000
        prob = 1/(2**(2*n)) * abs((-1 + np.exp(2.0j*m.pi*phi) )/(-1 + np.exp(2.0j*m.pi*phi/(2**n))))**2
        if abs(prob - estimatedProb) < abs(aproxProb - estimatedProb):
            aproxProb = prob
            aproxPhi = phi
            
    if( (MS1 < MS2) and ( (MS1!=0) and (MS2!=(2**n-1)) ) ):
        theta = (MS1+aproxPhi)/(2**n)
    elif( (MS1 > MS2) and (MS1!=0) ):
        theta = (MS1-aproxPhi)/(2**n)
    else:
        theta = 1+(MS1-aproxPhi)/(2**n)
        
    return aproxPhi,theta


def C_Oracle(qc, c, q, a1, a2, state):
    #qc.barrier()
    N = len(q)
    for i in np.arange(N):
        if( state[i]==0 ):
            qc.cx( c, q[int(i)] )

#---------------------------------
    qc.ccx( q[0], q[1], a1[0] )
    for j1 in np.arange(N-2):
        qc.ccx( q[int(2+j1)], a1[int(j1)], a1[int(1+j1)] )
    
    qc.ccx( c, a1[N-2], a2[0] )
    
    for j2 in np.arange(N-2):
        qc.ccx( q[int(N-1-j2)], a1[int(N-3-j2)], a1[int(N-2-j2)] )
    qc.ccx( q[0], q[1], a1[0] )

#---------------------------------
    for i2 in np.arange(N):
        if( state[i2]==0 ):
            qc.cx( c, q[int(i2)] )
    #qc.barrier()
            

def C_Diffusion(qc, c, q, a1, a2, ref):
    #qc.barrier()
    Q = len(q)
    N = 2**( Q )
    for j in np.arange(Q):
        qc.ch( c, q[int(j)] )
        
    if( ref ):
        for k in np.arange(1,N):
            C_Oracle(qc,c,q,a1,a2,Binary(int(k),N))
    else:
        C_Oracle(qc,c,q,a1,a2,Binary(0,N))
        
    for j2 in np.arange(Q):
        qc.ch( c, q[int(j2)] )
    #qc.barrier()

def C_Grover(qc, c, q, a1, a2, marked, **kwargs):
    #qc.barrier()
    Reflection=False
    if 'proper' in kwargs:
        Reflection = kwargs['proper']
        
    M = []
    for m1 in np.arange( len(marked) ):
        M.append( list(marked[m1]) )
        for m2 in np.arange( len(M[m1]) ):
            M[m1][m2] = int( M[m1][m2] )
            
    for i in np.arange(len(M)):        
        C_Oracle( qc,c,q,a1,a2,M[i] )
    
    C_Diffusion( qc,c,q,a1,a2,Reflection )
    #qc.barrier()

    
def GCD(a, b):
    gcd = 0
    if(a > b):
        num1 = a
        num2 = b
    elif(b > a):
        num1 = b
        num2 = a
    elif(a == b):
        gcd = a     
    while( gcd == 0 ):
        i = 1
        while( num1 >= num2*i ):
            i = i + 1
        if( num1 == num2*(i-1) ):
            gcd = num2
        else:
            r = num1 - num2*(i-1)
            num1 = num2
            num2 = r
    return gcd


def Euclids_Alg(a, b):
    if(a>=b):
        num1 = a
        num2 = b
    else:
        num1 = b
        num2 = a
    
    r_new = int( num1%num2 )
    r_old = int( num2 )
    while(r_new!=0):
        r_old = r_new
        r_new = int( num1%num2 )
        num1 = num2
        num2 = r_new
    
    gcd = r_old
    return gcd



def Modulo_f(Q, a, N):
    mods = []
    num = a%N
    for i in np.arange(1,2**Q):
        mods.append(num)
        num = (num*a)%N
            
    return mods

def Mod_Op(Q, qc, q1, q2, anc, a, N):
    #mods = Modulo_f(Q,a,N)
    num = a%N
    
    for j in np.arange( 2**Q ):
        q1_state = BinaryL( j, 2**Q )
        
        #q2_state = BinaryL( mods[j-1], 2**Q )
        q2_state = BinaryL(num, 2**Q )
        num = (num*a)%N
        
        X_Transformation(qc,q1,q1_state)
        
        gates = []
        for k in np.arange(Q):
            if(q2_state[k]==1):
                gates.append(['X',q2[int(k)]])
                
        n_Control_U(qc, q1, anc, gates)
        
        X_Transformation(qc,q1,q1_state)
        
def ConFrac(N, **kwargs):
    imax = 20
    r_a = False
    if 'a_max' in kwargs:
        imax = kwargs['a_max']
    if 'return_a' in kwargs:
        r_a = kwargs['return_a']
    a = []
    a.append( m.floor(N) )
    b = N - a[0]
    i = 1
    while( (round(b,10) != 0) and (i < imax) ):
        n = 1.0/b
        a.append( m.floor(n) )
        b = n - a[-1]
        i = i + 1
    #------------------------------
    a_copy = []
    for ia in np.arange(len(a)):
        a_copy.append(a[ia])
    for j in np.arange( len(a)-1 ):
        if( j == 0 ):
            p = a[-1] * a[-2] + 1
            q = a[-1]
            del a[-1]
            del a[-1]
        else:
            p_new = a[-1] * p + q
            q_new = p
            p = p_new
            q = q_new
            del a[-1]
    if(r_a == True):
        return q,p,a_copy
    return q,p


def r_Finder(a, N):
    value1 = a**1 % N
    r = 1
    value2 = 0
    while value1 != value2 or r > 1000:
        value2 = a**(int(1+r)) % N
        if( value1 != value2 ):
            r = r + 1
    return r

def Primality(N):
    is_prime = True
    if( (N==1) or (N==2) or (N==3) ):
        is_prime = True
    elif( (N%2==0) or (N%3==0) ):
        is_prime = False
    elif( is_prime==True ):
        p = 5
        while( (p**2 <= N) and (is_prime==True) ):
            if( (N%p==0) or (N%(p+2)==0) ):
                is_prime = False
            p = p + 6
    return is_prime

def Mod_r_Check(a, N, r):
    v1 = a**(int(2)) % N
    v2 = a**(int(2+r)) % N
    if( (v1 == v2) and (r<N) and (r!=0) ):
        return True
    return False



def Evaluate_S(S, L, a, N):
    Pairs = [[S,L]]
    for s in np.arange(3):
        S_new = int( S - 1 + s)
        for l in np.arange(3):
            L_new = int( L - 1 + l)
            if( ((S_new!=S) or (L_new!=L)) and (S_new!=L_new) ):
                Pairs.append( [S_new,L_new] )
#--------------------------- Try 9 combinations of S and L, plus or minus 1 from S & L
    period = 0
    r_attempts = []
    found_r = False
    while( (found_r==False) and (len(Pairs)!=0) ):
        order = 1
        S_o = Pairs[0][0]
        L_o = Pairs[0][1]
        q_old = -1
        q = 999
        while( q_old != q ):
            q_old = int( q )
            q,p = ConFrac(S_o/L_o,a_max=(order+1))
            new_r = True
            for i in np.arange(len(r_attempts)):
                if( q == r_attempts[i] ):
                    new_r = False
            if(new_r):
                r_attempts.append( int(q) )
                r_bool = Mod_r_Check(a,N,q)
                if( r_bool ):
                    found_r = True
                    q_old = q
                    period = int(q)
            order = order + 1
        del Pairs[0]
#--------------------------- Try higher multiples of already attempted r values
    r_o = 0
    while( (found_r == False) and (r_o < len(r_attempts)) ):
        k = 2
        r2 = r_attempts[r_o]
        while( k*r2 < N ):
            r_try = int(k*r2)
            new_r = True
            for i2 in np.arange(len(r_attempts)):
                if( r_try == r_attempts[i2] ):
                    new_r = False
            if(new_r):
                r_attempts.append( int(r_try) )
                r_bool = Mod_r_Check(a,N,r_try)
                if( r_bool ):
                    found_r = True
                    k = N
                    period = int(r_try)
            k = k + 1
        r_o = r_o + 1
#--------------------------- If a period is found, try factors of r for smaller periods
    if( found_r == True ):
        Primes = []
        for i in np.arange(2,period):
            if( Primality(int(i)) ):
                Primes.append(int(i))
        if( len(Primes) > 0 ):
            try_smaller = True
            while( try_smaller==True ):
                found_smaller = False
                p2 = 0
                while( (found_smaller==False) and (p2 < len(Primes)) ):
                    #print('p2: ',p2)
                    #print( 'period: ',period,' ',Primes[p2] )
                    try_smaller = False
                    if( period/Primes[p2] == m.floor( period/Primes[p2] ) ):
                        r_bool_2 = Mod_r_Check(a,N,int(period/Primes[p2]))
                        if( r_bool_2 ):
                            period = int(period/Primes[p2])
                            found_smaller = True
                            try_smaller = True
                    p2 = p2 + 1
    return period



def k_Data(k,n):
    Centers = []
    for i in np.arange(k):
        Centers.append( [1.5+np.random.rand()*5,1.5*np.random.random()*5] )
    count = int(round((0.7*n)/k))
    Data = []
    for j in range(len(Centers)):
        for j2 in range(count):
            r = np.random.random()*1.5
            x = Centers[j][0]+r*np.cos(np.random.random()*2*m.pi)
            y = Centers[j][1]+r*np.sin(np.random.random()*2*m.pi)
            Data.append([x, y])
                          
    for j2 in range(n - k * count):
        Data.append( [np.random.random()*8, np.random.random()*8] )
                         
    return Data
                          
def Initial_Centroids(k, D):
    D_copy = []
    for i in np.arange( len(D) ):
        D_copy.append( D[i] )
    Centroids = []
    for j in np.arange(k):
        p = np.random.randint(0,int(len(D_copy)-1))
        Centroids.append( [ D_copy[p][0] , D_copy[p][1] ] )
        D_copy.remove( D_copy[p] )
    return Centroids

def Update_Centroids(CT, CL):
    old_Centroids = []
    for c0 in np.arange(len(CT)):
        old_Centroids.append(CT[c0])
    Centroids = []
    for c1 in np.arange(len(CL)):
        mean_x = 0.
        mean_y = 0.
        for c2 in np.arange(len(CL[c1])):
            mean_x += CL[c1][c2][0]
            mean_y += CL[c1][c2][1]
        l = len(CL[c1])
        mean_x /= l
        mean_y /= l
        Centroids.append( [ mean_x,mean_y ] )
    return Centroids, old_Centroids

def Update_Clusters(D, CT, CL):
    old_Clusters = []
    for c0 in np.arange(len(CL)):
        old_Clusters.append(CL[c0])
    Clusters = []
    for c1 in np.arange( len(CT) ):
        Clusters.append( [] )
    for d in np.arange( len(D) ):
        closest = 'c'
        distance = 100000
        for c2 in np.arange( len(Clusters) ):
            Dist = m.sqrt( ( CT[c2][0] - D[d][0] )**2 + ( CT[c2][1] - D[d][1] )**2 )
            if( Dist < distance ):
                distance = Dist
                closest = int(c2)
        Clusters[closest].append( D[d] )
    return Clusters,old_Clusters


def Check_Termination(CL, oCL ):
    terminate = True
    for c1 in np.arange( len(oCL) ):
        for c2 in np.arange( len(oCL[c1]) ):
            P_found = False
            for c3 in np.arange( len(CL[c1]) ):
                if( CL[c1][c3] == oCL[c1][c2] ):
                    P_found = True
                    break
            if( P_found == False ):
                terminate = False
                break
    return terminate

def Draw_Data(CL, CT, oCT, fig, ax, colors, colors2 ):
    for j1 in np.arange( len(CL) ):
        ax.scatter( oCT[j1][0],oCT[j1][1], color='white', marker='s',s=80 )
    for cc in np.arange(len(CL)):
        for ccc in np.arange( len( CL[cc] ) ):
            ax.scatter( CL[cc][ccc][0],CL[cc][ccc][1], color=colors[cc],s=10 )
    for j2 in np.arange( len(CL) ):
        ax.scatter( CT[j2][0],CT[j2][1], color=colors2[j2], marker='x',s=50 )
    fig.canvas.draw()
    time.sleep(1)
    
    
def SWAP_Test( qc, control, q1, q2, classical, S ):
    qc.h( control )
    qc.cswap( control, q1, q2 )
    qc.h( control )
    qc.measure( control, classical )
    D = {'0':0}
    D.update( Measurement(qc,shots=S,return_M=True,print_M=False) )
    return D['0']

def Bloch_State( p,P ):
    x_min = P[0]
    x_max = P[1]
    y_min = P[2]
    y_max = P[3]
    theta = np.pi/2*( (p[0]-x_min)/(1.0*x_max-x_min) + (p[1]-y_min)/(1.0*y_max-y_min) )
    phi = np.pi/2*( (p[0]-x_min)/(1.0*x_max-x_min) - (p[1]-y_min)/(1.0*y_max-y_min) + 1 )
    return theta,phi



def Heatmap(data, show_text, show_ticks, ax, cmap, cbarlabel, **kwargs):
    valfmt="{x:.1f}"
    textcolors=["black", "white"]
    threshold=None
    cbar_kw={}
#----------------------------
    if not ax:
        ax = plt.gca()
    im = ax.imshow(data, cmap=cmap, **kwargs)
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    ax.grid(which="minor", color="black", linestyle='-', linewidth=1)
    if( show_ticks == True ):
        ax.set_xticks(np.arange(data.shape[1]))
        ax.set_yticks(np.arange(data.shape[0]))
        ax.tick_params(which="minor", bottom=False, left=False)
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
    if( show_text == True ):
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                kw.update(color=textcolors[int(im.norm(data[i, j]) < threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                
def Q_Update_Clusters(D, CT, CL, DS, shots):
    old_Clusters = []
    for c0 in np.arange(len(CL)):
        old_Clusters.append(CL[c0])
    Clusters = []
    for c1 in np.arange( len(CT) ):
        Clusters.append( [] )
#------------------------------------------------
    for d in np.arange( len(D) ):
        closest = 'c'
        distance = 0
        t,p = Bloch_State( D[d], DS )
        for c2 in np.arange( len(Clusters) ):
            t2,p2 = Bloch_State( CT[c2], DS )
            q = QuantumRegister( 3, name='q' )
            c = ClassicalRegister( 1, name='c' )
            qc= QuantumCircuit( q,c, name='qc' )
            qc.u( t, p, 0, q[1] )
            qc.u( t2, p2, 0, q[2] )
            IP = SWAP_Test( qc, q[0], q[1], q[2], c[0], shots )
            if( IP > distance ):
                distance = IP
                closest = int(c2)
        Clusters[closest].append(D[d] )
    return Clusters,old_Clusters


def E_Expectation_Value( qc, Energies ):
    SV = execute( qc, S_simulator, shots=1 ).result().get_statevector()
    EV = 0
    for i in range( len(SV) ):
        EV += Energies[i] *abs( SV[i] * np.conj(SV[i]) )
    EV = round(EV,4)
    return EV

def Top_States(States, Energies, SV, top):
    P = []
    S = []
    E = []
    for a in np.arange( top ):
        P.append(-1)
        S.append('no state')
        E.append('no energy')
    for i in range(len(States)):
        new_top = False
        probs = abs(SV[i]*np.conj(SV[i]))*100
        state = States[i]
        energ = Energies[i]
        j = 0
        while( (new_top == False) and (j < top) ):
            if( probs > P[j] ):
                for k in np.arange( int( len(P) - (j+1) ) ):
                    P[int( -1-k )] = P[int( -1-(k+1) )]
                    S[int( -1-k )] = S[int( -1-(k+1) )]
                    E[int( -1-k )] = E[int( -1-(k+1) )]
                P[j] = probs
                S[j] = state
                E[j] = energ
                new_top = True
            j = int(j+1)
    for s in range( top ):
        print('State ',S[s],' Probability: ',round(P[s],2),'%',' Energy: ',round(E[s],2))
        

def Ising_Energy(V, E, **kwargs):
    Trans = False
    if 'Transverse' in kwargs:
        if( kwargs['Transverse'] == True ):
            Trans = True
    Energies = []
    States = []
    for s in range( 2**len(V) ):
        B = BinaryL(int(s),2**len(V))
        B2 = []
        for i in range(len(B)):
            if( B[i] == 0 ):
                B2.append(1)
            else:
                B2.append(-1)
        state = ''
        
        energy = 0
        for s2 in range(len(B)):
            state = state+str(B[s2])
            energy -= V[s2][1]*B2[s2]
        States.append(state)
        
        for j in range( len(E) ):
            if( Trans == False ):
                energy -= B2[int(E[j][0])] * B2[int(E[j][1])]
            else:
                energy -= B2[int(E[j][0])] * B2[int(E[j][1])] * E[j][2]
                
        Energies.append(energy)
    return Energies,States

def Ising_Circuit(qc, q, V, E, beta, gamma, **kwargs):
    Trans = False
    if 'Transverse' in kwargs:
        if( kwargs['Transverse'] == True ):
            Trans = True
    Mixer = 1
    if 'Mixing' in kwargs:
        Mixer = int( kwargs['Mixing'] )

    p = 1
    if 'p' in kwargs:
        p = int(kwargs['p'])

    for c in range(p):
        Uc_Ising(qc,q,gamma,V,E,Trans)
        if( Mixer == 2 ):
            Ub_Mixer2(qc,q,beta,V)
        else:
            Ub_Mixer1(qc,q,beta,V)
        
def Uc_Ising(qc, q, gamma, Vert, Edge, T):
    for e in range( len(Edge) ): # ZZ
        if( T == False ):
            G = gamma
        else:
            G = gamma * Edge[e][2]
        qc.cx( q[int(Edge[e][0])], q[int(Edge[e][1])] )
        qc.rz( 2*G, q[int(Edge[e][1])] )
        qc.cx( q[int(Edge[e][0])], q[int(Edge[e][1])] )
    for v in range( len(Vert) ): # Z_gamma
        # WARNING: This was with error in the article, the 'on site' magnetic field was not taken into account
        qc.rz( gamma * Vert[v][1], q[int(Vert[v][0])] )
        

        
def Ub_Mixer1(qc, q, beta, Vert):
    for v in np.arange( len(Vert) ):
        qc.rx( beta, q[int(v)] )
        
        
def Ub_Mixer2(qc, q, beta, Vert):
    for v in range( len(Vert) ):
        qc.rx( beta, q[int(Vert[v][0])] )
    qc.cx( q[0], q[1] )
    qc.cx( q[1], q[2] )
    qc.cx( q[2], q[0] )
    for v2 in range( len(Vert) ):
        qc.ry( beta, q[int(Vert[v2][0])] )
        

def Ising_Gradient_Descent(qc, q, Circ, V, E, beta, gamma, epsilon, En, step, **kwargs):
    Trans = False
    if 'Transverse' in kwargs:
        if( kwargs['Transverse'] == True ):
            Trans = True
    Mixer = 1
    if 'Mixing' in kwargs:
        Mixer = int(kwargs['Mixing'])
    params = [ [beta+epsilon,gamma],[beta-epsilon,gamma],[beta,gamma+epsilon],[beta,gamma-epsilon] ]
    ev = []
    for i in np.arange( 4 ):
        q = QuantumRegister(len(V))
        qc= QuantumCircuit(q)
        for hh in np.arange(len(V)):
            qc.h( q[int(hh)] )
        Circ( qc, q, V, E, params[i][0], params[i][1], Transverse=Trans, Mixing=Mixer )
        ev.append( E_Expectation_Value( qc, En ) )
    beta_next = beta - ( ev[0] - ev[1] )/( 2.0*epsilon ) * step
    gamma_next = gamma - ( ev[2] - ev[3] )/( 2.0*epsilon ) * step
    return beta_next, gamma_next

def MaxCut_Energy(V, E):
    Energies = []
    States = []
    for s in np.arange( 2**len(V) ):
        B = BinaryL(int(s),2**len(V))
        B2 = []
        for i in np.arange(len(B)):
            if( B[i] == 0 ):
                B2.append(1)
            else:
                B2.append(-1)
        state = ''
        for s2 in np.arange(len(B)):
            state = state+str(B[s2])
        States.append(state)
        energy = 0
        for j in np.arange( len(E) ):
            energy = energy + 0.5* ( 1.0 - B2[int(E[j][0])]*B2[int(E[j][1])] )
        Energies.append(energy)
    return Energies,States

def MaxCut_Circuit(qc, q, V, E, beta, gamma):
    Uc_MaxCut( qc, q, gamma,E)
    Ub_Mixer1(qc,q,beta,V)
    
def Uc_MaxCut(qc, q, gamma, edge):
    for e in np.arange( len(edge) ):
        qc.cx( q[int(edge[e][0])], q[int(edge[e][1])] )
        qc.rz( gamma, q[int(edge[e][1])] )
        qc.cx( q[int(edge[e][0])], q[int(edge[e][1])] )
        
def p_Gradient_Ascent(qc, q, Circ, V, E, p, Beta, Gamma, epsilon, En, step):
    params = []
    for i in np.arange(2):
        for p1 in np.arange(p):
            if( i == 0 ):
                params.append( Beta[p1] )
            elif( i == 1 ):
                params.append( Gamma[p1] )
    ep_params = []
    for p2 in np.arange( len( params ) ):
        for i2 in np.arange( 2 ):
            ep = []
            for p3 in np.arange( len(params) ):
                ep.append( params[p3] )
            ep[p2] = ep[p2] + (-1.0)**(i2+1)*epsilon
            ep_params.append( ep )
    ev = []
    for p4 in np.arange( len( ep_params ) ):
        run_params = ep_params[p4]
        q = QuantumRegister(len(V))
        qc= QuantumCircuit(q)
        for hh in np.arange(len(V)):
            qc.h( q[int(hh)] )
        for p5 in np.arange(p):
            Circ( qc, q, V, E, run_params[int(p5)], run_params[int(p5+p)] )
        ev.append( E_Expectation_Value( qc, En ) )
    Beta_next = []
    Gamma_next = []
    for k in np.arange( len( params ) ):
        if( k < len( params )/2 ):
            Beta_next.append( params[k] - (ev[int(2*k)] - ev[int(2*k+1)])/( 2.0*epsilon ) * step )
        else:
            Gamma_next.append( params[k] - (ev[int(2*k)] - ev[int(2*k+1)])/( 2.0*epsilon ) * step )
    return Beta_next, Gamma_next







def Single_Qubit_Ansatz( qc, qubit, params ):
    qc.ry( params[0], qubit )
    qc.rz( params[1], qubit )
    
    
def VQE_Gradient_Descent(qc, q, H, Ansatz, theta, phi, epsilon, step, **kwargs):
    EV_type = 'measure'
    if 'measure' in kwargs:
        M_bool = kwargs['measure']
        if( M_bool == True ):
            EV_type = 'measure'
        else:
            EV_type = 'wavefunction'
    Shots = 1000
    if 'shots' in kwargs:
        Shots = kwargs['shots']
    params = [theta,phi]
    ep_params = [[theta+epsilon,phi],[theta-epsilon,phi],[theta,phi+epsilon],[theta,phi-epsilon]]
    Hk = list( H.keys() )
    EV = []
    for p4 in np.arange( len( ep_params ) ):
        H_EV = 0
        qc_params = ep_params[p4]
        for h in np.arange( len(Hk) ):
            qc_params = ep_params[p4]
            q = QuantumRegister(1)
            c = ClassicalRegister(1)
            qc= QuantumCircuit(q,c)
            Ansatz( qc, q[0], [qc_params[0], qc_params[1]] )
            if( Hk[h] == 'X' ):
                qc.ry(-m.pi/2,q[0])
            elif( Hk[h] == 'Y' ):
                qc.rx(m.pi/2,q[0])
            if( EV_type == 'wavefunction' ):
                sv = execute( qc, S_simulator, shots=1 ).result().get_statevector()
                H_EV = H_EV + H[Hk[h]]*( (np.conj(sv[0])*sv[0]).real - (np.conj(sv[1])*sv[1]).real )
            elif( EV_type == 'measure' ):
                qc.measure( q,c )
                M = {'0':0,'1':0}
                M.update( Measurement( qc, shots=Shots, print_M=False, return_M=True ) )
                H_EV = H_EV + H[Hk[h]]*(M['0']-M['1'])/Shots
        EV.append(H_EV)
    theta_slope = ( EV[0]-EV[1] )/(2.0*epsilon)
    phi_slope = ( EV[2]-EV[3] )/(2.0*epsilon)
    next_theta = theta - theta_slope*step
    next_phi = phi - phi_slope*step
    return next_theta,next_phi


def Two_Qubit_Ansatz(qc, q, params):
    Single_Qubit_Ansatz( qc, q[0], [params[0], params[1]] )
    Single_Qubit_Ansatz( qc, q[1], [params[2], params[3]] )
    qc.cx( q[0], q[1] )
    Single_Qubit_Ansatz( qc, q[0], [params[4], params[5]] )
    Single_Qubit_Ansatz( qc, q[1], [params[6], params[7]] )
    
    
def Calculate_MinMax(V, C_type):
    if( C_type == 'min' ):
        lowest = [V[0],0]
        for i in np.arange(1,len(V)):
            if( V[i] < lowest[0] ):
                lowest[0] = V[i]
                lowest[1] = int(i)
        return lowest
    elif( C_type == 'max' ):
        highest = [V[0],0]
        for i in np.arange(1,len(V)):
            if( V[i] > highest[0] ):
                highest[0] = V[i]
                highest[1] = int(i)
        return highest
    
    
def Compute_Centroid(V):
    points = len( V )
    dim = len( V[0] )
    Cent = []
    for d in np.arange( dim ):
        avg = 0
        for a in np.arange( points ):
            avg = avg + V[a][d]/points
        Cent.append( avg )
    return Cent

def Reflection_Point(P1, P2, alpha):
    P = []
    for p in np.arange( len(P1) ):
        D = P2[p] - P1[p]
        P.append( P1[p]+alpha*D )
    return P



def VQE_EV(params, Ansatz, H, EV_type, **kwargs):
    Shots = 10000
    if 'shots' in kwargs:
        Shots = int( kwargs['shots'] )
    Hk = list( H.keys() )
    H_EV = 0
    for k in range( len(Hk) ): # for each term in Hamiltonian
        L = list( Hk[k] )
        q = QuantumRegister(len(L))
        c = ClassicalRegister(len(L))
        qc= QuantumCircuit(q,c)
        Ansatz( qc, q, params )
        sv0 = execute( qc, S_simulator, shots=1 ).result().get_statevector()
        if( EV_type == 'wavefunction' ):
            for l in range( len(L) ):
                if( L[l] == 'X' ):
                    qc.x( q[int(l)] )
                elif( L[l] == 'Y' ):
                    qc.y( q[int(l)] )
                elif( L[l] == 'Z' ):
                    qc.z( q[int(l)] )
            sv = execute( qc, S_simulator, shots=1 ).result().get_statevector()
            H_ev = 0
            for l2 in range(len(sv)):
                H_ev = H_ev + (np.conj(sv[l2])*sv0[l2]).real
            H_EV = H_EV + H[Hk[k]] * H_ev
        elif( EV_type == 'measure' ):
            # apply Hamiltonian term
            for l in range( len(L) ):
                if( L[l] == 'X' ):
                    qc.ry(-m.pi/2,q[int(l)])
                elif( L[l] == 'Y' ):
                    qc.rx( m.pi/2,q[int(l)])
                    
            # measure it
            qc.measure( q,c )
            M = Measurement( qc, shots=Shots, print_M=False, return_M=True )
            Mk = list( M.keys() )
            
            # compute energy estimate
            H_ev = 0
            for m1 in range(len(Mk)):
                MS = list( Mk[m1] )
                e = 1
                for m2 in range(len(MS)):
                    if( MS[m2] == '1' ):
                        e = e*(-1)
                        
                H_ev = H_ev + e * M[Mk[m1]]
                
            H_EV = H_EV + H[Hk[k]]*H_ev/Shots
    return H_EV


def Nelder_Mead(H, Ansatz, Vert, Val, EV_type):
    alpha = 2.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5
    add_reflect = False
    add_expand = False
    add_contract = False
    shrink = False
    add_bool = False
#----------------------------------------
    hi = Calculate_MinMax( Val,'max' )
    Vert2 = []
    Val2 = []
    for i in np.arange(len(Val)):
        if( int(i) != hi[1] ):
            Vert2.append( Vert[i] )
            Val2.append( Val[i] )
    Center_P = Compute_Centroid( Vert2 )
    Reflect_P = Reflection_Point(Vert[hi[1]],Center_P,alpha)
    Reflect_V = VQE_EV(Reflect_P,Ansatz,H,EV_type)
#------------------------------------------------- # Determine if: Reflect / Expand / Contract / Shrink
    hi2 = Calculate_MinMax( Val2,'max' )
    lo2 = Calculate_MinMax( Val2,'min' )
    if( hi2[0] > Reflect_V >= lo2[0] ):
        add_reflect = True
    elif( Reflect_V < lo2[0] ):
        Expand_P = Reflection_Point(Center_P,Reflect_P,gamma)
        Expand_V = VQE_EV(Expand_P,Ansatz,H,EV_type)
        if( Expand_V < Reflect_V ):
            add_expand = True
        else:
            add_reflect = True
    elif( Reflect_V > hi2[0] ):
        if( Reflect_V < hi[0] ):
            Contract_P = Reflection_Point(Center_P,Reflect_P,rho)
            Contract_V = VQE_EV(Contract_P,Ansatz,H,EV_type)
            if( Contract_V < Reflect_V ):
                add_contract = True
            else:
                shrink = True
        else:
            Contract_P = Reflection_Point(Center_P,Vert[hi[1]],rho)
            Contract_V = VQE_EV(Contract_P,Ansatz,H,EV_type)
            if( Contract_V < Val[hi[1]] ):
                add_contract = True
            else:
                shrink = True
#------------------------------------------------- # Apply: Reflect / Expand / Contract / Shrink
    if( add_reflect == True ):
        new_P = Reflect_P
        new_V = Reflect_V
        add_bool = True
    elif( add_expand == True ):
        new_P = Expand_P
        new_V = Expand_V
        add_bool = True
    elif( add_contract == True ):
        new_P = Contract_P
        new_V = Contract_V
        add_bool = True
        
    if( add_bool ):
        del Vert[hi[1]]
        del Val[hi[1]]
        Vert.append( new_P )
        Val.append( new_V )
    
    if( shrink ):
        Vert3 = []
        Val3 = []
        lo = Calculate_MinMax( Val,'min' )
        Vert3.append( Vert[lo[1]] )
        Val3.append( Val[lo[1]] )
        for j in np.arange( len(Val) ):
            if( int(j) != lo[1] ):
                Shrink_P = Reflection_Point(Vert[lo[1]],Vert[j],sigma)
                Vert3.append( Shrink_P )
                Val3.append( VQE_EV(Shrink_P,Ansatz,H,EV_type) )
        for j2 in np.arange( len(Val) ):
            del Vert[0]
            del Val[0]
            Vert.append( Vert3[j2] )
            Val.append( Val3[j2] )
            
