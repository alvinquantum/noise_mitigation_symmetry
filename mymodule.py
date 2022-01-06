from copy import deepcopy
import numpy as np
from qiskit import execute, transpile, QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.circuit.library.standard_gates import (IGate, U1Gate, U2Gate, U3Gate, XGate,
                                                   YGate, ZGate, HGate, SGate, SdgGate, TGate,
                                                   TdgGate, RXGate, RYGate, RZGate, CXGate,
                                                   CYGate, CZGate, CHGate, CRZGate, CU1Gate,
                                                   CU3Gate, SwapGate, RZZGate,
                                                   CCXGate, CSwapGate)
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
from qiskit.visualization import circuit_drawer
from qiskit.opflow import X,Y,Z
import os, math, pickle, cirq, json
from cirq.contrib.qasm_import import circuit_from_qasm
from cirq.qis import fidelity
from os import listdir
from os.path import isfile

def random_circuit_cnot(num_qubits, num_cnots_required, seed=None):
    """Generate random circ: Generates a random circuit with num_qubits, 5 rz gates and num_cnots
    The circuit uses h, s, cz, cx, swap gates and no measurements. 

    Returns:
        QuantumCircuit: constructed circuit
    """
    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    qc=QuantumCircuit(QuantumRegister(num_qubits))
    cnot_count=0
    #Randomly generate cliffords until we have more than enough cnots
    #then truncate
    while cnot_count < num_cnots_required:
        qc_temp=decompose_clifford(random_clifford(num_qubits))
        cnot_count+=count_gate(qc_temp, "cx")
        qc.compose(qc_temp, inplace=True)
    
    # Too many cnots so trim.
    if cnot_count> num_cnots_required:
        qc=trim(qc, cnot_count, num_cnots_required)
        assert dict(qc.count_ops())["cx"] == num_cnots_required, "the circuit doesn't have the required cnots"
    # qc=add_rz_gates_prob(num_qubits, qc, rng) #Non deterministic # of rz.
    qc=add_rz_gates_det(num_qubits, qc, rng) #Deterministic # of rz
    return qc

def trim(qc, cnot_count, num_cnots_required):
    '''Generate random circ: Helper function. Remove the end of the circuit until we have enough cnots. Returns: QuantumCircuit'''
    assert cnot_count > num_cnots_required, "Number of cnots required is greater than actual. Cannot trim."
    qc_dag = circuit_to_dag(qc)
    layers= list(qc_dag.multigraph_layers())
    # Remove the necessary cnots.
    for layer in layers[::-1]:
        for node in layer:
            # Check if the node is cnot and reduce the cnot count
            if node.name=="cx":
                cnot_count-=1
            # Remove the node as long as it's an operation.
            if node.type=="op":
                qc_dag.remove_op_node(node)
            # We've reduced the circuit to required cnot so return.
            if cnot_count==num_cnots_required:
                return dag_to_circuit(qc_dag)

def insert_rz_gate(circ, rng, prob, operand):
    '''Generate random circ: Helper function. Add rz with probability.'''
    if rng.uniform(0,1)<=prob:
        #Add a random rz to the operand
        angle=rng.uniform(0, 2 * np.pi)
        circ.rz(angle, operand)

def add_rz_gates_det(num_qubits, qc, rng):
    '''Generate random circ: Helper function. Add a set number of rz gates randomly to the circuit. Returns: QuantumCircuit'''
    qc_dag=circuit_to_dag(qc)
    new_qc=QuantumCircuit(QuantumRegister(num_qubits))
    NUM_RZ=5

    # Each dag has layers and each layer has nodes.
    layers=list(qc_dag.multigraph_layers())
    # Get all the nodes. We will sample from this.
    all_nodes=[]
    for layer in layers:
        for node in layer:
            # We insert rz gates infront of the selected node. We distribute the required number of rzs
            # throught the circuit. We ignore the "in" nodes since they
            # would double the probability of getting an rz in the beginning. 
            # Also we need to add multiqubit nodes as many times as 
            # their number of arguments so that we have a uniform distribution.
            if node.type!="in":
                # print(node.type)
                # print(node.qargs)
                #Need to check that it's not an out node because out node returns an empty list for qargs.
                if node.type!="out" and len(node.qargs)>1:
                    # print(node.qargs)
                    for _ in range(len(node.qargs)):
                        all_nodes.append(node)
                else:
                    all_nodes.append(node)
    # Now we sample.
    chosen_nodes=rng.choice(all_nodes, replace=False, size=NUM_RZ).tolist()

    for layer in layers:
        for node in layer:
            # The node is part of the chosen set so add an rz.
            if node in chosen_nodes:
                if node.type == "out":
                    insert_rz_gate(new_qc, rng, 1, node.wire.index)

                elif node.type == "op":
                    add_rz_to_node_det(new_qc, chosen_nodes, node, rng)
            # Only copy op nodes
            if node.type == "op":
                copy_node(new_qc, node)    
    assert dict(new_qc.count_ops())["rz"] == NUM_RZ, "Number of rz's wasn't met."
    return new_qc

def add_rz_to_node_det(new_qc, chosen_nodes, node, rng):
    '''Generate random circ: Helper function for deterministic adding of rz gates. Process the given node and insert rz gates.'''
    # The number of qubits that the node/gate is operating on can be greater than 1.
    if len(node.qargs)==1:
        insert_rz_gate(new_qc, rng, 1, node.qargs[0].index)
    else:
        #Since it's a multi qubit gate there's a chance it was chosen multiple times.
        number_of_rzs_to_add=chosen_nodes.count(node)
        # Even though node.qargs is a list we can't iterate through it so it must be done
        # this way.
        remaining_qubits=[node.qargs[idx].index for idx in range(len(node.qargs))]
        # Add the gates.
        while number_of_rzs_to_add:
            rng.shuffle(remaining_qubits)
            operand=remaining_qubits[0]
            remaining_qubits = [q for q in remaining_qubits if q != operand]
            insert_rz_gate(new_qc, rng, 1, operand)
            number_of_rzs_to_add-=1

def add_rz_gates_prob(num_qubits, qc, rng):
    '''Generate random circ: Helper function. Add rz gates randomly to the circuit. Returns: QuantumCircuit'''
    const= 5 #Vary this parameter.
    cnot_count=count_gate(qc, "cx")
    swap_count=count_gate(qc, "swap")

    # We add the two qubit gates twice to the length and the number of qubits to account for
    # output nodes. This puts the expectation at const for any circuit.
    length=len(qc)+cnot_count+swap_count+ num_qubits
    prob=const/length

    qc_dag=circuit_to_dag(qc)
    new_qc=QuantumCircuit(QuantumRegister(num_qubits))

    # Each dag has layers and each layer has nodes.
    layers=list(qc_dag.multigraph_layers())
    for layer in layers:
        for node in layer:
            # We're adding rz gates to the front so we can ignore in nodes.
            if node.type == "out":
                insert_rz_gate(new_qc, rng, prob, node.wire.index)

            elif node.type == "op":
                # The number of qubits that the node/gate is operating on can be greater than 1.
                if len(node.qargs)==1:
                    insert_rz_gate(new_qc, rng, prob, node.qargs[0].index)
                else:
                    #Apply an rz gate to each qubit with some probability.
                    # Even though node.qargs is a list we can't iterate through it so it must be done
                    # this way.
                    remaining_qubits=[node.qargs[idx].index for idx in range(len(node.qargs))]
                    while remaining_qubits:
                        rng.shuffle(remaining_qubits)
                        operand=remaining_qubits[0]
                        remaining_qubits = [q for q in remaining_qubits if q != operand]
                        insert_rz_gate(new_qc, rng, prob, operand)
                #only copy op nodes.
                copy_node(new_qc, node)
    assert len(new_qc)>= len(qc), "Inserting RZ gates wasn't done properly."
    return new_qc

def count_gate(qc, gate):
    '''Generate random circ: Helper function: counts the number of occurrences of the gate in the qc.'''
    gates_counts=qc.count_ops()
    count=0
    if gate in dict(gates_counts):
        count+=dict(gates_counts)[gate]
    return count

def copy_node(new_qc, node):
    '''Generate random circ: Helper function: Copy the node into new_qc'''
    # Copy the node.
    if node.name=="x":
        new_qc.x(node.qargs[0].index)
    elif node.name=="y":
        new_qc.y(node.qargs[0].index)
    elif node.name=="z":
        new_qc.z(node.qargs[0].index)
    elif node.name=="h":
        new_qc.h(node.qargs[0].index)
    elif node.name=="s":
        new_qc.s(node.qargs[0].index)
    elif node.name=="sdg":
        new_qc.sdg(node.qargs[0].index)
    elif node.name=="cx":
        new_qc.cx(node.qargs[0].index, node.qargs[1].index)
    elif node.name=="swap":
        new_qc.swap(node.qargs[0].index, node.qargs[1].index)
    #Test
    elif node.name=="rz":
        # print(node.op.params)
        # print(type(node.op))
        new_qc.rz(node.op.params[0], node.qargs[0].index)
    elif node.name=="sx":
        new_qc.sx(node.qargs[0].index)   
        
    else:
        # We have overlooked a gate type.
        assert False, node.name + " gate wasn't matched in the DAG."

def random_circuit_depth(num_qubits, depth, seed=None):
    """Generate random circ: Generates a random circuit with num_qubits and depth (each wire has the specified depth). 
    The circuit uses h, s, t, cx gates and no measurements.
    Derived from https://qiskit.org/documentation/_modules/qiskit/circuit/random/utils.html#random_circuit.

    Returns:
        QuantumCircuit: constructed circuit
    """

    # This code is part of Qiskit.
    #
    # (C) Copyright IBM 2017.
    #
    # This code is licensed under the Apache License, Version 2.0. You may
    # obtain a copy of this license in the LICENSE.txt file in the root directory
    # of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
    #
    # Any modifications or derivative works of this code must retain this
    # copyright notice, and modified files need to carry a notice indicating
    # that they have been altered from the originals.

    max_operands=2 
    one_q_ops = [HGate, SGate, TGate]
    two_q_ops = [CXGate]

    qr = QuantumRegister(num_qubits, 'q')
    qc = QuantumCircuit(num_qubits)


    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    # apply arbitrary random operations at every depth
    for _ in range(depth):
        # choose either 1, 2, or 3 qubits for the operation
        remaining_qubits = list(range(num_qubits))
        while remaining_qubits:
            max_possible_operands = min(len(remaining_qubits), max_operands)
            num_operands = rng.choice(range(max_possible_operands)) + 1
            rng.shuffle(remaining_qubits) # We shuffle and then apply the operations to the necessary number of qubits
            #  at the beginning.
            operands = remaining_qubits[:num_operands]
            remaining_qubits = [q for q in remaining_qubits if q not in operands]
            if num_operands == 1:
                operation = rng.choice(one_q_ops)
            else:
                operation = rng.choice(two_q_ops)

            register_operands = [qr[i] for i in operands]
            op = operation()

            qc.append(op, register_operands)

    assert len(list(circuit_to_dag(qc).multigraph_layers()))-2 == depth, "The depth is wrong in the generated circuit."

    return qc

class PushOperator:
    '''Finding checks: Symbolic: push operations.'''
    # track of the direction at the start of the layer.
    @staticmethod
    def x(op2, temp_check_reversed):
        '''Finding checks: Symbolic: Push x through op2. 
        Forward: Circuit: [x]--[Op2]=[Op2]--[F] Math: (Op2)(x)=(F)(Op2)-->
        (F)=(Op2)(x)(Op2)^\dagger.

        Backwards: Circuit: [Op2]--[x]=[F]--[Op2] Math: (x)(Op2)=(Op2)(F)-->
        (F)=(Op2)^\dagger (x)(Op2)
        
        In return: first digit is phase.'''
        if op2=="X":
            return [1, "X"]
        elif op2=="Y":
            return [-1, "X"]
        elif op2=="Z":
            return [-1, "X"]
        elif op2=="H":
            return [1, "Z"]
        elif op2=="S":
            if temp_check_reversed.forward:
                return [1, "Y"]
            else:
                return [-1, "Y"]
        elif op2=="SDG":
            if temp_check_reversed.forward:
                return [-1, "Y"]
            else:
                return [1, "Y"]
        #Test
        elif op2=="RZ":
            temp_check_reversed.change_to_backwards=True
            # Increment the layer_idx so that we redo this layer when going
            # backwards. You can see this in the decrementing logic of find_checks_sym.
            # temp_check_reversed.layer_idx+=1
            return [1, "I"]
        else:
            # We have overlooked a gate type.
            assert False, op2 + " gate wasn't matched in the DAG."

    @staticmethod
    def y(op2, temp_check_reversed):
        '''Finding checks: Symbolic: Push y through op2.  
        Forward: Circuit: [y]--[Op2]=[Op2]--[F] Math: (Op2)(y)=(F)(Op2)-->
        (F)=(Op2)(y)(Op2)^\dagger.

        Backwards: Circuit: [Op2]--[y]=[F]--[Op2] Math: (y)(Op2)=(Op2)(F)-->
        (F)=(Op2)^\dagger (y)(Op2)
        
        In return: first digit is phase.'''
        if op2=="X":
            return [-1, "Y"]
        elif op2=="Y":
            return [1, "Y"]
        elif op2=="Z":
            return [-1, "Y"]
        elif op2=="H":
            return [-1, "Y"]
        elif op2=="S":
            if temp_check_reversed.forward:
                return [-1, "X"]
            else:
                return [1, "X"]
        elif op2=="SDG":
            if temp_check_reversed.forward:
                return [1, "X"]
            else:
                return [-1, "X"]
        elif op2=="RZ":
            temp_check_reversed.change_to_backwards=True
            # Increment the layer_idx so that we redo this layer when going
            # backwards. You can see this in the decrementing logic of find_checks_sym.
            # temp_check_reversed.layer_idx+=1
            return [1, "I"]
        else:
            # We have overlooked a gate type.
            assert False, op2 + " gate wasn't matched in the DAG." 

    @staticmethod        
    def z(op2):
        '''Finding checks: Symbolic: Push z through op2. 
        Forward: Circuit: [z]--[Op2]=[Op2]--[F] Math: (Op2)(z)=(F)(Op2)-->
        (F)=(Op2)(z)(Op2)^\dagger.

        Backwards: Circuit: [Op2]--[z]=[F]--[Op2] Math: (z)(Op2)=(Op2)(F)-->
        (F)=(Op2)^\dagger (z)(Op2)
        
        In return: first digit is phase.'''
        if op2=="X":
            return [-1, "Z"]
        elif op2=="Y":
            return [-1, "Z"]
        elif op2=="Z":
            return [1, "Z"]
        elif op2=="H":
            return [1, "X"]
        elif op2=="S":
            return [1, "Z"]
        elif op2=="SDG":
            return [1, "Z"]
        elif op2=="RZ":
            # Both operators are diagonal so they commute.
            return [1, "Z"]            
        else:
            # We have overlooked a gate type.
            assert False, op2 + " gate wasn't matched in the DAG." 

    @staticmethod
    def cx(op1):
        '''Finding checks: Symbolic: Push op1 through cx.  
        Forward: Circuit: [op1]--[cx]=[cx]--[F] Math: (cx)(op1)=(F)(cx)
        Backwards: Circuit: [cx]--[op1]=[F]--[cx] Math: (op1)(cx)=(cx)(F)
        -->
        (F)=(cx)(op1)(cx).
        
        In return: first digit is phase.'''

        if op1==["I", "I"]:
            return [1, "I", "I"]
        elif op1==["I", "X"]:
            return [1, "I", "X"]
        elif op1==["I", "Y"]:
            return [1, "Z", "Y"]
        elif op1==["I", "Z"]:
            return [1, "Z", "Z"]

        elif op1==["X", "I"]:
            return [1, "X", "X"]
        elif op1==["X", "X"]:
            return [1, "X", "I"]
        elif op1==["X", "Y"]:
            return [1, "Y", "Z"]
        elif op1==["X", "Z"]:
            return [-1, "Y", "Y"]

        elif op1==["Y", "I"]:
            return [1, "Y", "X"]
        elif op1==["Y", "X"]:
            return [1, "Y", "I"]
        elif op1==["Y", "Y"]:
            return [-1, "X", "Z"]
        elif op1==["Y", "Z"]:
            return [1, "X", "Y"]

        elif op1==["Z", "I"]:
            return [1, "Z", "I"]
        elif op1==["Z", "X"]:
            return [1, "Z", "X"]
        elif op1==["Z", "Y"]:
            return [1, "I", "Y"]
        elif op1==["Z", "Z"]:
            return [1, "I", "Z"]

        else:
            # We have overlooked a gate type.
            assert False, op1[0] + ", " + op1[1] + " wasn't a pauli element." 

    @staticmethod
    def swap(op1):
        '''Finding checks: Symbolic: pass op1 through.
        In return: first digit is phase.'''
        result_ops=list(reversed(op1))
        result=[1]
        return result+result_ops

class CheckOperator:
    '''Finding checks: Symbolic: Stores the check operation along with the phase. operations is a list of strings.'''
    __slots__=["phase", "operations"]
    def __init__(self, phase, operations):
        self.phase=phase
        self.operations=operations

class TempCheckOperator:
    '''Finding checks: Symbolic: Stores the check operation along with the phase. operations is a list of strings.'''
    __slots__=["phase", "operations", "change_to_backwards", "forward", "layer_idx"]
    def __init__(self, phase, operations):
        self.phase=phase
        self.operations=operations
        self.change_to_backwards=False
        self.forward=True
        self.layer_idx=1

class CircuitProperties:
    '''Circuit properties holder.'''
    __slots__=["number_of_qubits", "cnot_count", "number_of_circuits", "circ", "circ_operations"]
    def __init__(self, number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations,):
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count
        self.number_of_circuits=number_of_circuits
        self.circ=circ
        self.circ_operations=circ_operations

class TestCircuits:
    '''Testing circuits: For running the simulations.'''
    def __init__(self, cirq_circ_with_checks, cirq_circ_no_checks, number_of_qubits, rho_correct, sanity_check_fidelity, keep_qubits):
        self.cirq_circ_with_checks=cirq_circ_with_checks
        self.cirq_circ_no_checks=cirq_circ_no_checks
        self.number_of_qubits=number_of_qubits
        self.rho_correct=rho_correct
        self.keep_qubits=keep_qubits
        self.sanity_check_fidelity=sanity_check_fidelity

    def run_test(self, error_info):
        '''Testing circuits: runs the test.'''
        cirq_circ_with_checks=self.cirq_circ_with_checks
        cirq_circ_no_checks=self.cirq_circ_no_checks
        rho_correct=self.rho_correct
        number_of_qubits=self.number_of_qubits
        keep_qubits=self.keep_qubits
        sanity_check_fidelity=self.sanity_check_fidelity
        error_idx=error_info[0]
        single_qubit_error=error_info[1]

        '''Testing circuits: get testing results'''
        print("adding noise...")
        noisy_cirq_circ_with_checks=add_noise(cirq_circ_with_checks, single_qubit_error) 
        # print(noisy_cirq_circ_with_checks)   
        print("getting result rho...")    
        noisy_rho_with_checks=get_result_rho(noisy_cirq_circ_with_checks, number_of_qubits+1, keep_qubits)
        print("taking trace...")
        ancilla_zero_outcome_probability=np.trace(noisy_rho_with_checks)
        fidelity_noisy_rho_with_check=fidelity(np.around(noisy_rho_with_checks* 1/ancilla_zero_outcome_probability, 5), 
            rho_correct, qid_shape=(2**(number_of_qubits),), validate=False)
        
        noisy_cirq_circ_no_checks=add_noise(cirq_circ_no_checks, single_qubit_error)
        # print(noisy_cirq_circ_with_checks)
        noisy_rho_no_checks=get_result_rho(noisy_cirq_circ_no_checks, number_of_qubits, keep_qubits)
        fidelity_noisy_rho_no_check=fidelity(noisy_rho_no_checks, rho_correct, qid_shape=(2**(number_of_qubits),), validate=False)
        print("single qubit error rate: ", single_qubit_error)
        print("ancilla 0 prob outcome: ", ancilla_zero_outcome_probability)
        print("fidelity no check: ", fidelity_noisy_rho_no_check)
        print("fidelity with check: ", fidelity_noisy_rho_with_check)
        print()

        return {"num_results_before_postselect": 1, "num_results_after_postselect": ancilla_zero_outcome_probability, "error_idx": error_idx, 
        "one_qubit_err": single_qubit_error, "two_qubit_err": 10*single_qubit_error, 
        "state_fidelity_no_checks_with_errors": fidelity_noisy_rho_no_check, 
        "state_fidelity_with_checks_with_errors": fidelity_noisy_rho_with_check, 
        "state_fidelity_with_checks_no_errors": sanity_check_fidelity}

class ChecksProperties:
    '''Checks properties holder.'''
    __slots__=["count", "p2_weights", "pauli_str_p1s", "pauli_str_p2s"]
    def __init__(self, count, p2_weights, pauli_str_p1s, pauli_str_p2s):
        self.count=count
        self.p2_weights=p2_weights
        self.pauli_str_p1s=pauli_str_p1s
        self.pauli_str_p2s=pauli_str_p2s

def update_current_ops(op1, op2, temp_check_reversed, current_qubits):
    '''Finding checks: Symbolic: Finds the intermediate check. Always push op1 through op2. '''
    if len(op1)==1:
        if op1[0]=="X":
            result= PushOperator.x(op2, temp_check_reversed)
        elif op1[0]=="Y":
            result= PushOperator.y(op2, temp_check_reversed)
        elif op1[0]=="Z":
            result= PushOperator.z(op2)
        elif op1[0]=="I":
            result= [1, "I"]
        else:
            # Can expand to accomodate non pauli + I in the future.
            assert False, op1[0] + " is not I, X, Y, or Z."
    else:
        # Two qubit operations
        if op2=="CX":
            result = PushOperator.cx(op1)
        elif op2=="SWAP":
            result= PushOperator.swap(op1)
        else:
            assert False, op2 + " is not cx or swap."

    result_phase=result[0]
    result_ops=result[1::1]
    temp_check_reversed.phase=temp_check_reversed.phase*result_phase
    # Coppy the current ops into temp_check_reversed.
    for idx, op in enumerate(result_ops):
        temp_check_reversed.operations[current_qubits[idx]]=op

def get_check_strs(p1, p2):
    '''Finding checks: Symbolic: turns p1 and p2 to strings results.'''
    #P1s
    p1_operations=p1.operations
    p1_phase=str(p1.phase)
    if len(p1_phase)==1:#Add if + if positive phase
        p1_phase="+"+p1_phase
    p1_operations.insert(0, p1_phase)

    #P2s
    p2_operations=p2.operations
    p2_weight=get_weight(p2_operations)
    p2_phase=str(p2.phase)
    if len(p2_phase)==1:#Add + if positive phase
        p2_phase="+"+p2_phase
    p2_operations.insert(0, p2_phase)

    p1_str="".join(p1_operations)
    p2_str="".join(p2_operations)
    
    print("p1: ", p1_str)
    print("p2: ", p2_str)
    print("Pauli weight P2: ", p2_weight)
    print()

    return (p2_weight, p1_str, p2_str)

def can_continue(forward, op1, op2):
    '''Finding checks: Symbolic: Determine if can continue. If we're going backwards and op1 is not I or Z and op2 is RZ then don't continue.'''
    if forward==False and op2=="RZ" and op1!="I" and op1!="Z":
        return False
    else:
        return True

def get_current_qubits(node):
    '''Finding checks: Symbolic: get the current qubits whose operations that will be passed through.'''
    # We have to check for single or two qubit gates.
    if node.name in ["x", "y", "z", "h", "s", "sdg", "rz"]:
        return [node.qargs[0].index]
    elif node.name in ["cx", "swap"]:
        return [node.qargs[0].index, node.qargs[1].index]
    else:
        assert False, "Overlooked a node operation." 

def get_weight(pauli_string):
    '''Gets the weight of a Pauli string. Returns: int'''
    count=0
    for character in pauli_string:
        if character!="I":
            count+=1
    return count

def create_controlU(npmat, number_of_qubits):
    '''Testing circuits:
    Finding checks: Helper of numpy method.
    Returns: np.array'''
    return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)    

def check_p2(control_p1, control_p2, unitary, number_of_qubits):
    '''Finding checks: Helper for numpy method.
    Sanity check for p2. U\otimes I- ControlP2^\dagger(U\otimes I)ControlP1==0'''
    assert np.allclose(np.kron(np.eye(2),unitary.data)-control_p2.dot(np.kron(np.eye(2),unitary.data)).dot(control_p1), np.zeros(2**(number_of_qubits+1))), "wrong p2"

def layer_to_circ(layer, temp_layer_circ):
    '''Finding checks: Helper for transpile method. 
    Convert layer to circiut.'''
    for node in layer:
        if node.type == "op":
            copy_node(temp_layer_circ, node)

def pauli_to_circuit(pauli_str):
    '''Finding checks: Helper for transpile method. 
    Converts the Pauli string to a circuit. Returns: QuantumCircuit'''
    circ=QuantumCircuit(len(pauli_str))
    qubit_pos=len(pauli_str)-1
    for char in pauli_str:
        if char=="X":
            # print(phase_added, phase, char)
            circ.x(qubit_pos)
        elif char=="Y":
            circ.y(qubit_pos)
        elif char=="Z":
            circ.z(qubit_pos)
        qubit_pos-=1
    return circ

def append_checks_to_circ(circ_properties, checks_properties):
    '''Finding checks: combines everything into one circuit.'''
    circ=circ_properties.circ
    pauli_str_p1=checks_properties.pauli_str_p1s[0]
    pauli_str_p2=checks_properties.pauli_str_p2s[0]
    number_of_qubits=circ_properties.number_of_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits)
    ancilla_reg=QuantumRegister(1)
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)

    add_controlU(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    add_controlU(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg)
    return temp_circ

def find_checks_sym(pauli_group_elem, circ):
    '''Finding checks: Symbolic: Finds p1 and p2 elements symbolically.'''
    print(pauli_group_elem)

    # We will just iterate over the +1 phase elements of the pauli group since the 
    # scenarios can be recovered by just multiplying by the phase constant.
    pauli_group_elem_ops=list(pauli_group_elem)
    p1=CheckOperator(1, pauli_group_elem_ops)
    p2=CheckOperator(1, ["I" for _ in range(len(pauli_group_elem))])
    temp_check_reversed=TempCheckOperator(1, list(reversed(pauli_group_elem_ops)))

    # Iterate through the circuit. We manually keep track of the idx since
    # we can either go forward or backwards. This is kept track of inside the temp_check_reversed
    # We also track layer_idx in temp_check_reversed.
    # forward=True
    circ_dag = circuit_to_dag(circ)
    layers = list(circ_dag.multigraph_layers())
    num_layers=len(layers)
    # # We start index 1 since the first layer are just in nodes.
    # layer_idx = 1

    while True:
        # Get current layer
        layer=layers[temp_check_reversed.layer_idx]
        for node in layer:
            # Iterate through layers and nodes.
            # if found.value:
            #     print("exiting.")
            #     return
            # elif node.type=="op":
            if node.type=="op":
                current_qubits=get_current_qubits(node)
                current_ops=[temp_check_reversed.operations[qubit] for qubit in current_qubits]
                node_op= node.name.upper()

                # Update temp_check_reversed if possible
                if can_continue(temp_check_reversed.forward, current_ops[0], node_op):
                    update_current_ops(current_ops, node_op, temp_check_reversed, current_qubits)
                else:
                    return
        
        # See if we should start going backwards.
        if temp_check_reversed.change_to_backwards:
            temp_check_reversed.forward=False
            temp_check_reversed.change_to_backwards=False
            # we don't increment or decrement the layer_idx since we processed this layer.
            # We have to process the same layer going backwards.

        # Since we're not changing to backwards, either move forward or backards
        elif temp_check_reversed.forward:
            if temp_check_reversed.layer_idx==num_layers-1:
                p2.phase=temp_check_reversed.phase
                p2.operations=list(reversed(temp_check_reversed.operations))
                # Append operations.
                # with count.get_lock():
                result =get_check_strs(p1, p2)
                return result
            else:
                temp_check_reversed.layer_idx+=1            
        else:
            if temp_check_reversed.layer_idx==1:
                # We reached the first layer of operation nodes so move forward. Note
                # the zero index are all input nodes so we can skip.
                p1.phase=temp_check_reversed.phase
                p1.operations=list(reversed(temp_check_reversed.operations))
                temp_check_reversed.forward=True
            else:
                temp_check_reversed.layer_idx-=1

def find_checks_with_numpy(pauli_group_tuple, unitary, number_of_qubits, ABS_TOL, 
    pauli_labels, pauli_group_positive, table_length, count, p2_weights, pauli_str_p1s, pauli_str_p2s):
    '''Finding checks: Uses numpy. Probably depricate.'''
    idx1=pauli_group_tuple[0]
    # to_matrix() returns a list.
    p1=pauli_group_tuple[1].to_matrix()[0]
    if idx1 >= table_length:
        p1 = p1*-1
    #U.p1=p2.U ---->U.p1.U^\dagger=p2. Operator class so we need .data to access numpy array.
    p2=unitary.dot(p1).dot(unitary.adjoint()).data
    #Sanity check. Can comment out.
    control_p1=create_controlU(p1, number_of_qubits)
    control_p2=create_controlU(p2, number_of_qubits)
    if not check_p2(control_p1, control_p2, unitary, number_of_qubits):
        return
    # Check if p2 is traceless. All elements of the pauli group are traceless except identity.
    if not math.isclose(0.0,np.trace(p2), abs_tol=ABS_TOL):
        return
    #Only need p2 with +1 phase since the global phase can be absorbed into p1. Faster this way.
    for idx2, label_element in enumerate(pauli_group_positive):
        # to_matrix() returns a list.
        element=label_element.to_matrix()[0]
        # allclose checks if the values are within tolerance of atol=10^-8.
        if np.allclose(p2, element):
            #Have to check which part of the table p1 belongs to so we can print the correct phase.
            if idx1-table_length<0:
                p1_str="+1"+pauli_labels[idx1 % table_length]
            else:
                p1_str="-1"+pauli_labels[idx1 % table_length]
            # elif idx1-2*table_length<0:
            #     p1_str="-1"+pauli_labels[idx1 % table_length]
            # elif idx1-3*table_length<0:
            #     p1_str="+j"+pauli_labels[idx1 % table_length]
            # else:
            #     p1_str="-j"+pauli_labels[idx1 % table_length]
            #P2 is always +1 phase.
            p2_str="+1"+pauli_labels[idx2]
            print("p1: ", p1_str)
            print("p2: ", p2_str)

            # Print the weight. We care about P2 weight since we commute p1 through U.
            p2_weight=get_weight(pauli_labels[idx2])
            print("Pauli weight P2: ", p2_weight)
            print()
            #Need to lock the value so it doesn't change while checking. Since we write only
            #after locking the weight we don't need to lock the other values.
            with count.get_lock():
                count.value+=1
                if p2_weight>p2_weights[0]:
                    # We store the max values in the beginning.
                    p2_weights[0]=p2_weight
                    pauli_str_p1s[0]=p1_str
                    pauli_str_p2s[0]=p2_str
                pauli_str_p1s.append(p1_str)
                pauli_str_p2s.append(p2_str)
                p2_weights.append(p2_weight)

            # #Need to lock the value so it doesn't change. Nonatomic operation.
            # with count.get_lock():
            #     count.value+=1
            #We found p2 so just return.
            return

def find_checks_with_transpile(pauli_group_elem, circ):
    '''Finding checks: Uses transpile. May deprecate.'''
    # Temporay storage
    temp_p2_circ=pauli_to_circuit(pauli_group_elem)

    #Go through the layers of the circ
    circ_dag=circuit_to_dag(circ)
    layers=list(circ_dag.multigraph_layers())
    for layer in layers:
        #Convert layer to circuit.
        temp_layer_circ=QuantumCircuit(circ.num_qubits)
        layer_to_circ(layer, temp_layer_circ)
        # for node in layer:
        #     if node.type == "op":
        #         mymodule.copy_node(temp_layer_circ, node)
        temp_p2_circ=find_intermediate_p2(temp_p2_circ, temp_layer_circ)
    print("p2:", pauli_group_elem)
    print("Temp_p2: ")
    print(temp_p2_circ)
    return(temp_p2_circ)

def find_intermediate_p2(p1_circ, circ):
    '''Finding checks: Helper function for transpile method.'''
    circ_inverse=QuantumCircuit.inverse(circ)
    # U.p_1=c_2.U---->Up_1U^dagger=c_2
    p2_circ=circ_inverse.compose(p1_circ).compose(circ)
    # print(c2_circ)
    #Optimize the circuit
    basis_gates=['id', 'rz', 'sx', 'cx']
    p2_circ=transpile(p2_circ, basis_gates=basis_gates, optimization_level=2)
    # print(c2_circ)
    return p2_circ

# def create_controlU(npmat, number_of_qubits):
#     '''Testing cricuits: Returns a controlled operation. Type is np.array'''
#     return np.kron(np.array([[1,0],[0,0]]),np.eye(2**number_of_qubits))+np.kron(np.array([[0,0],[0,1]]),npmat)

def write_outputs(circ_properties, checks_properties, file_number, file_info_path, file_qasm_path, output_file):
    '''Finding checks: write outputs.'''
    pauli_str_p1s=checks_properties.pauli_str_p1s
    pauli_str_p2s=checks_properties.pauli_str_p2s
    p2_weights=checks_properties.p2_weights
    count=checks_properties.count
    circ_operations=circ_properties.circ_operations
    number_of_qubits=circ_properties.number_of_qubits
    circ=circ_properties.circ

    assert len(pauli_str_p1s)==len(pauli_str_p2s), "number of p1's and p2's don't match."
    assert len(pauli_str_p2s)==len(p2_weights), "number of p2's and weights don't match."
    # Count doesn't count the extra storage of max vales in the beginning of the lists of solutions (e.g. pauli_str_p1s)
    assert len(pauli_str_p2s)==count+1, "number of p2's and counts of solutions don't match."
    #Outputs
    if count==0:
        output_file.write("\n")
        output_file.write("nothing found: trivial solution\n")
        print("nothing found: trivial solution")
    else:
        for index1, strp1 in enumerate(pauli_str_p1s):
            if index1!=0: #The max stuff are stored in the beginning, which we print at the end.
                output_file.write("\n")
                output_file.write("p1: "+ str(strp1)+ "\n")
                output_file.write("p2: "+ str(pauli_str_p2s[index1])+ "\n")
                output_file.write("Pauli weight P2: "+ str(p2_weights[index1])+ "\n")
        output_file.write("\n")
        output_file.write("Found Matches: "+ str(count)+ "\n")
        output_file.write("Max Weight: "+ str(p2_weights[0]) + "\n")
        output_file.write("P1 that creates max P2: "+ str(pauli_str_p1s[0])+ "\n")
        output_file.write("Max P2: " + str(pauli_str_p2s[0])+ "\n")
        print("Found Matches: ", count)   
        print("Max P2 Weight: ", p2_weights[0])
        print("P1 that creates max P2: ", pauli_str_p1s[0])
        print("Max P2: ", pauli_str_p2s[0])

    cnot_count=0
    if "cx" in circ_operations:
        cnot_count=circ_operations["cx"]

    rz_count=0
    if "rz" in circ_operations:
        rz_count=circ_operations["rz"]
    output_file.write("Circuit no: " + str(file_number)+ "\n")
    output_file.write("Qubits: " + str(number_of_qubits)+ "\n")
    output_file.write("CNOT count: "+ str(cnot_count) + "\n")
    output_file.write("RZ count: "+ str(rz_count) + "\n")

    print("Qubits: ", str(number_of_qubits))
    print("Circuit no: ", str(file_number))
    print("CNOT count: ", cnot_count)
    print("RZ count: ", rz_count)

    # # Dump all the info into a pickle
    circ_file=open(file_info_path, "wb")
    pickle.dump({"cx": cnot_count, "rz": rz_count, "qubits": number_of_qubits, "circuit_num" : file_number, "found_matches": count,
            "max_pauli_weight": p2_weights[0], "max_pauli_str_p1": pauli_str_p1s[0], "max_pauli_str_p2": pauli_str_p2s[0]}, circ_file)
    circ.qasm(filename=file_qasm_path)
    # # Close the files. 
    circ_file.close()
    output_file.close()

def add_controlU(circ, pauli_str, number_of_qubits, quantum_register, ancilla_register):
    '''Testing cricuits: Adds a controlled Pauli to circ.'''
    print(pauli_str)
    c_minus_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(1.0592554e-13,3.288277,2.8116442) q[0];
    u3(pi,5.8991791,2.7575864) q[1];
    cx q[0],q[1];
    u3(2*pi,5.4995268,4.1085153) q[0];
    u3(pi,0.38400625,3.5255989) q[1];
    """.strip())

    c_minus_j_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(6.9759301e-13,4.3371506,0.076820611) q[0];
    u3(6.9089241e-09,5.477472,3.947306) q[1];
    cx q[0],q[1];
    u3(2*pi,0.60772643,2.832284) q[0];
    u3(5.5644668e-09,5.4477338,3.9770441) q[1];
    """.strip())

    c_plus_j_x =QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(pi,4.6477413,6.258122) q[0];
    u3(1.7599424,3*pi/2,3*pi/2) q[1];
    cx q[0],q[1];
    u3(pi,3.1666559,3.2062404) q[0];
    u3(4.9015351,3*pi/2,3*pi/2) q[1];
    """.strip())

    #Check the length of the pauli string. It maybe empty when there were no matches found.
    if len(pauli_str)>number_of_qubits:
        phase=pauli_str[:2]
        pauli_str=pauli_str[2:]
        qubit_pos=number_of_qubits-1
        phase_added=False
        for pos, char in enumerate(pauli_str):
            # First case is most common
            if phase_added or phase=="+1":
                if char=="X":
                    # print(phase_added, phase, char)
                    circ.cx(ancilla_register, quantum_register[qubit_pos])
                elif char=="Y":
                    # print(phase_added, phase, char)
                    circ.sdg(quantum_register[qubit_pos])
                    circ.cx(ancilla_register, quantum_register[qubit_pos])
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    # print(phase_added, phase, char)
                    circ.h(quantum_register[qubit_pos])
                    circ.cx(ancilla_register, quantum_register[qubit_pos])
                    circ.h(quantum_register[qubit_pos])
            # -1 phase
            else:
                if char=="X":
                    # print(phase_added, phase, char)
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                elif char=="Y":
                    # print(phase_added, phase, char)
                    circ.sdg(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    # print(phase_added, phase, char)
                    circ.h(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_register[0], quantum_register[qubit_pos]], inplace=True)
                    circ.h(quantum_register[qubit_pos])
                # Note the phase added needs this check because of identity terms in the pauli strings.
                # With no check, if there is an identity it will change phase_added to true and we won't get the
                # added necessary phase on of the pauli terms.
                if char!="I":
                    phase_added=True
            qubit_pos-=1
    return circ

def result_exists(base_path, file_name, result_num):
    '''Testing circuits: Checks if the result file exists.'''
    name_split=file_name.split("_")
    result_name="_".join(name_split[:-1])+"_result_"+str(result_num)+"_.txt"
    # print(file_name)
    # print(result_name)
    output_file_txt_path=os.path.join(base_path, result_name)
    # print(output_file_txt_path)
    return os.path.isfile(output_file_txt_path)

def get_number_of_results(job):
    '''Testing circuits: Get's the number of measurement outcomes. Before post selecting this should be number of tomography circuits*shots.'''
    tot=0
    for idx, _ in enumerate(job.results):
        # Each idx represents one tomography circuit. The number of shots is how much we repeat the circuit.
        # The list given by get_counts contains all the outcomes from running that specific tomography circuit. get_counts
        # returns a list of outcomes for that specific tomography circuit. The number of outcomes for that specific circuit
        # is equal to the number of shots.
        curr_exp=job.get_counts(idx)
        # print(curr_exp)
        tot+=sum(curr_exp.values())
    # print(tot)
    return tot

def run_job(job, basis_gates, noise_model, shots_arg, simulator):
    '''Testing cricuits: Runs the given job. Returns result.'''
    print("running job...")
    # simulator=Aer.get_backend('qasm_simulator')
    # simulator.set_option("method", "density_matrix")
    # simulator.set_option("max_parallel_experiments", psutil.cpu_count(logical=False))
    # simulator.set_option("max_parallel_shots", psutil.cpu_count(logical=False))
    # simulator.set_option("max_parallel_threads", psutil.cpu_count(logical=False))
    # result= execute(job, Aer.get_backend('qasm_simulator'),
    #             basis_gates=basis_gates,
    #             noise_model=noise_model, 
    #             shots=shots_arg, optimization_level=0).result()
    result= execute(job, simulator,
                basis_gates=basis_gates,
                noise_model=noise_model, 
                shots=shots_arg, optimization_level=0).result()
    print("job done.")
    
    return result

def store_results(circ, circ_full_with_ancilla, base_path, circ_file_name, number_of_qubits, results, cnot_count, rz_count, circuit_num,
    count, p2_weight, pauli_str_p1, pauli_str_p2):
    '''Testing circuits: Stores all the results.'''
    #File naming stuff. 
    #Strip the extension.
    split_circ_file_name=circ_file_name.split("_")
    file_name_no_extension="_".join(split_circ_file_name[:-1])

    # print(output_file_name_obj)
    # print(os.path.join(base_path,output_file_name_obj))
    temp_file_number=0
    output_file_name_obj=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.obj"
    while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
        temp_file_number+=1
        output_file_name_obj=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.obj"
    output_file_name_txt=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.txt"
    output_file_name_qasm=file_name_no_extension+"_result_"+ str(temp_file_number) +"_.qasm"
    circ_full_with_ancilla.qasm(filename=os.path.join(base_path, output_file_name_qasm))

    # Dump all the results into a pickle
    with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
        pickle.dump({"cx": cnot_count, "rz": rz_count, "qubits": number_of_qubits, "circuit_num" : circuit_num, "found_matches: ": count,
             "max_pauli_weight": p2_weight, "max_pauli_str_p1": pauli_str_p1, "max_pauli_str_p2": pauli_str_p2 ,"results": results}, circ_file)

    #Print text results to file
    output_file_txt_path=os.path.join(base_path, output_file_name_txt)
    circuit_drawer(circ_full_with_ancilla, filename=output_file_txt_path)
    with open(output_file_txt_path, "a") as output_file_txt:
        output_file_txt.write("\n")
        output_file_txt.write(json.dumps(circ.count_ops()))
        for result in results:
            output_file_txt.write("\n")
            output_file_txt.write("Error idx: "+str(result["error_idx"])+"\n")
            output_file_txt.write("One qubit error: "+str(result["one_qubit_err"])+"\n")
            # print("One_qubit_err:", one_qubit_err)
            output_file_txt.write("Two qubit error: "+str(result["two_qubit_err"])+"\n")
            # print("Two qubit error:", two_qubit_err)
            output_file_txt.write("State fidelity no checks and with errors: "+str(result["state_fidelity_no_checks_with_errors"])+"\n")
            output_file_txt.write("State fidelity with checks and with errors: "+str(result["state_fidelity_with_checks_with_errors"])+"\n")
            output_file_txt.write("Sanity check fidelity with checks and no errors: "+str(result["state_fidelity_with_checks_no_errors"])+"\n")

def get_files(base_path, number_of_qubits, cnot_count, start_circ_number, end_circ_number):
    '''Testing circuits: Get the desired files for testing.'''

    # Gets the files that match the string. Files include the path
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    rand_circ_files=[]
    circ_properties_files=[]
    for file in all_files:
        name_split=file.split("_")
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if (".qasm" in name_split and "result" not in name_split and name_split_nums[1]==cnot_count and name_split_nums[0]==number_of_qubits 
            and start_circ_number<=name_split_nums[2]<=end_circ_number):
            rand_circ_files.append(file)
            circ_properties_files.append("_".join(name_split[:-1])+"_.obj")

    return (rand_circ_files, circ_properties_files)

def add_rand_input_state(number_of_qubits, quantum_register, circ_with_checks, circ_no_checks):
    '''Testing circuits: Create a random state. Need to send both circs at the same time so they have the same random initial state.'''
    #Insert random state generator
    random_params = np.random.uniform(size=(number_of_qubits, 3))
    for i,qreg in enumerate(quantum_register):
        for j, pauli in enumerate([X,Y,Z]):
            rand_tuple=(random_params[i][j])
            circ_with_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
            circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)

def get_checks_properties(base_path, file_name):
    '''Testing circuits: Gets the check properties. Return type: CheckProperties'''
    with open(os.path.join(base_path, file_name), "rb") as circ_file:
        circ_info=pickle.load(circ_file)
            # circ=circ_info["circ"]
    return ChecksProperties(circ_info["found_matches"], circ_info["max_pauli_weight"], 
        circ_info["max_pauli_str_p1"], circ_info["max_pauli_str_p2"])

class Circs:
    '''Testing circuits: Holder for circ variables.'''
    def __init__(self, quantum_register, ancilla_register, circ_with_checks, circ_no_checks, circ_initial_state):
        self.quantum_register=quantum_register
        self.ancilla_register=ancilla_register
        self.circ_with_checks=circ_with_checks
        self.circ_no_checks=circ_no_checks
        self.circ_initial_state=circ_initial_state

def create_circs(number_of_qubits, circ_pieces, ancilla_name):
    '''Testing circuits: Creates circs no checks and with checks. The circs have the same random initial state.
    Return type: list'''
    #Create the complete circuit without a measurement
    quantum_register=QuantumRegister(number_of_qubits)
    ancilla_register=QuantumRegister(1, name=ancilla_name)
    
    all_qubits=[]
    for reg in quantum_register:
        all_qubits.append(reg)
    all_qubits.append(ancilla_register[0])
    circ_with_checks=QuantumCircuit(quantum_register, ancilla_register)
    #Create the circuit with no checks
    circ_no_checks=QuantumCircuit(quantum_register, ancilla_register)
    add_rand_input_state(number_of_qubits, quantum_register, circ_with_checks, circ_no_checks)
    # Save the initial state circuit.
    circ_initial_state=deepcopy(circ_no_checks)
    # print("random initial state", rand_initial_state)
    # circ_with_checks.barrier()
    circ_with_checks.h(ancilla_register)
    for elem in circ_pieces:
        # circ_with_checks.barrier()
        circ_with_checks.compose(elem, inplace=True)
    # circ_with_checks.barrier()
    circ_no_checks.compose(circ_pieces[1], inplace=True)
    circ_with_checks.h(ancilla_register)
    #Test
    # print("original")
    # for circ in circ_pieces:
    #     print(circ)
    # print("circ full")
    # print(circ_full)
    # print("circ no checks")
    # print(circ_no_checks)
    # print(circ_with_checks)
    #Test
    basis_gates=['u1', 'u2', 'u3', 'cx']
    circ_with_checks=transpile(circ_with_checks, basis_gates=basis_gates, optimization_level=0)
    circ_no_checks=transpile(circ_no_checks, basis_gates=basis_gates, optimization_level=0)
    circ_initial_state=transpile(circ_initial_state, basis_gates=basis_gates, optimization_level=0)
    print(circ_with_checks)
    return Circs(quantum_register, ancilla_register, circ_with_checks, circ_no_checks, circ_initial_state)

def split_circuit_by_barrier(qasm_file_path):
    '''Testing circuits: Split circuits by barrier.'''
    with open(qasm_file_path, "r") as file:
        qasm = file.read()
        prelude = []
        circuits = [[]]
        for line in qasm.splitlines():
            if any([line.startswith(t) for t in ['OPENQASM', 'include', 'qreg', 'creg']]):
                prelude.append(line)
            elif line.startswith('barrier'):
                circuits.append([])
            else:
                circuits[-1].append(line)
        circuits_with_prelude = [prelude+circuit for circuit in circuits]
        return list(map(lambda x: QuantumCircuit.from_qasm_str("\n".join(x)), circuits_with_prelude))

# def get_fidelity(circ, number_of_qubits, rho_compare):
#     '''Testing circuits: Uses Google Cirq. Returns the fidelity between the output of circ and rho_compare.'''
#     # Add noise to the circuit
#     rho_circ=simulator.simulate(circ).final_density_matrix
#     rho_circ_reduced=partial_trace(np.reshape(rho_circ, [2,2]*(number_of_qubits+1)), keep_indices=keep_qubits)
#     rho_circ_reduced=np.reshape(rho_circ_reduced, (2**number_of_qubits, 2**number_of_qubits))
#     result=fidelity(rho_circ_reduced, rho_compare, qid_shape=(2**(number_of_qubits),), validate=False)
#     return result

def get_result_rho(circ, number_of_qubits, keep_qubits):
    '''Testing circuits: Uses Google Cirq. 
    number_of_qubits: total number of qubits in circ. 
    keep_qubits: is a list. 
    Returns resulting rho from simulation of circ.'''
    simulator=cirq.DensityMatrixSimulator()
    print("simulating....")
    trial_result=simulator.simulate(circ)
    rho= np.around(trial_result.final_density_matrix, 5)
    final_size=len(keep_qubits)
    if number_of_qubits!=final_size:
        # Have to expand the indices for each qubit, e.g., for 2 qubits a_ij|i><j|-->a_ijkl|ij><kl|. Then 
        # to get the first state sum over indices j and l which are indices 1 and 3. For the second state
        # sum over i and k which are indices 0 and 2. In the cirq package, we specify which indices to keep.
        # But we only need to specify the first index not both. Thus, we can think of keep_indices as a list of
        # qubits to keep.
        print("taking partial trace...")
        # print(rho)
        rho_reduced=cirq.partial_trace(np.reshape(rho, [2,2]*(number_of_qubits)), keep_indices=keep_qubits)
        # We have to reshape back to a square matrix. We reasign to rho so we can just return rho.
        rho=np.reshape(rho_reduced, (2**final_size, 2**final_size))
        # print("reduced....")
        # print(rho)
        # print(rho.shape)
    return rho

def add_noise(circ, single_qubit_error):
    '''Testing circuits: Uses Google Cirq. Adds noise to circ.'''
    two_qubit_error=10*single_qubit_error
    twoqubit_noise_model=cirq.DepolarizingChannel(p=two_qubit_error, n_qubits=2)
    singlequbit_noise_model=cirq.DepolarizingChannel(p=single_qubit_error)
    # all_qubits=circ.all_qubits()
    moments=[]
    noisy_cirq=cirq.Circuit()
    # Iterate through the moments. For two qubit gates we add a two qubit depolarization at a two qubit error rate.
    # print(all_qubits)
    for idx, moment in enumerate(circ):
        # Skip the initial state.
        if idx <3:
            moments+=[moment]
        else:
            error_ops=[]
            # TODO: Should we do this? 
            # For the moment find all the qubits with no operations and noise to them. 
            # for qubit in all_qubits:
            #     if qubit not in moment.qubits:
            #         # print(qubit)
            #         error_ops += singlequbit_noise_model.on_each(qubit)

            # Go through the operations in the moment. For 2 qubit operations use the
            # two qubit error. For everything else, i.e. single qubit gates, use single
            # qubit error. 
            for operation in moment.operations:
                # print(operations)
                # print(type(operation.gate))
                # print(type(moment.qubits))
                # print("all moment qubits: ", moment.qubits)
                if isinstance(operation.gate, cirq.ops.common_gates.CXPowGate) or isinstance(operation.gate, cirq.ops.SwapPowGate):
                    # print("2 qubit ", operation)
                #     # print(twoqubit_noise_model.on_each(operation.qubits))
                    error_ops += twoqubit_noise_model.on_each(operation.qubits)
                    # error_ops += singlequbit_noise_model.on_each(operation.qubits)

                else:
                    # print("1 qubit ", operation)
                    # print(singlequbit_noise_model.on_each(operation.qubits))
                    error_ops += singlequbit_noise_model.on_each(operation.qubits)
            # print()
            # print(error_ops)
            moments+=[moment, cirq.ops.Moment(error_ops)]
        # print(moments)
    noisy_cirq+=cirq.Circuit(moments)
    # noisy_cirq=circ.with_noise(cirq.depolarize(p=single_qubit_error))
    return noisy_cirq

def post_select_on_ancilla(res, ancilla_value, new_nqubits): 
    """Testing circuits:
    strip the results where ancilla was not equal to `ancilla_value`
    This is some voodoo copied from 
    https://qiskit.org/documentation/tutorials/noise/8_tomography.html#2-Qubit-Conditional-State-Tomography  
    """
    assert(isinstance(ancilla_value, str))
    res_new = deepcopy(res)
    for resultidx, _ in enumerate(res.results):
        old_counts = res.get_counts(resultidx)
        new_counts = {}
        # print(res_new.results[resultidx].header.clbit_labels)
        res_new.results[resultidx].header.creg_sizes = [res_new.results[resultidx].header.creg_sizes[0]]
        res_new.results[resultidx].header.clbit_labels = res_new.results[resultidx].header.clbit_labels[0:-1]
        # print(res_new.results[resultidx].header.clbit_labels)
        res_new.results[resultidx].header.memory_slots = new_nqubits 
        for reg_key in old_counts:
            # print(reg_key)
            reg_bits = reg_key.split(' ')
            assert(len(reg_bits) == 2)
            assert(len(reg_bits[0]) == 1)
            if reg_bits[0]==ancilla_value:
                new_counts[reg_bits[1]]=old_counts[reg_key]
            res_new.results[resultidx].data.counts = new_counts
    return res_new

if __name__ == "__main__":
    pass