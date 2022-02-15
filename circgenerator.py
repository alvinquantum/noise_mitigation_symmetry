# from copy import deepcopy
from sys import maxsize
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import QuantumRegister
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
import checksfinder
from multiprocessing import Pool
import psutil
import time
import os
import utilities
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit.quantum_info import Operator
from qiskit.visualization import circuit_drawer
from cirq.contrib.qasm_import import circuit_from_qasm
import json
# import pickle

# def get_initial_state_file_names(raw_file_names):
#     initial_states_file_names=[]
#     for file_name in raw_file_names:
#         name_split=file_name.split("_")
#         base_file_name=name_split[:-2]
#         initial_states_file_names.append(f"{'_'.join(base_file_name)}_result_0_.qasm")
#     return initial_states_file_names

# class CircuitSplitter:
#     def __init__(self, circ):
#         self.unsplit_circ=circ
#         self.unsplit_circ_dag=circ

#     @staticmethod
#     def split_circs_by_qubits(base_path, raw_qasm_file, minimum_qubits, max_qubits):
#         '''splits the given cirq circuit into two circuits at the depth. 
#         The depth is part of the second circuit.'''
#         # for file_name in raw_qasm_files:
#         # print(os.path.join(base_path, file_name))
#         circ=QuantumCircuit.from_qasm_file(os.path.join(base_path, raw_qasm_file))
#         circ_dag=circuit_to_dag(circ)
#         layers=circ_dag.multigraph_layers()
#         cirq_circ=circuit_from_qasm(circ.qasm())
#         # while True:
#         qubit_sets=[]
#         layer_idx=0
#         found=False
#         for idx in range(len(cirq_circ)):
#             temp_qubit_sets=[]
#             independent_sets=cirq_circ[:idx].get_independent_qubit_sets()
#             if found:
#                 temp_qubit_sets.append()
#                 break
#             for elem in list(map(lambda x: len(x), independent_sets)):
#                 if elem>max_qubits:
#                     found=True
#             layer_idx=idx
#             qubit_sets=independent_sets
#         print(qubit_sets)
#         print(layer_idx)
#             #     layer.

#     # def traverse_forward(self, temp_qubits, temp_layer):
#     #     '''Go forward and expand number of qubits when hitting a 2 qubit gate.'''
#     #     for moment in self.unsplit_circ:
#     #         for ops in moment:

#     # def traverse_backwards(self, temp_qubits, temp_layer):
#     #     for moment in self.unsplit_circ:
#     #         for ops in moment:

class CircuitProperties:
    '''Circuit properties holder.'''
    __slots__=["number_of_qubits", "cnot_count", "number_of_circuits", "circ", "circ_operations"]
    def __init__(self, number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations,):
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count
        self.number_of_circuits=number_of_circuits
        self.circ=circ
        self.circ_operations=circ_operations

def generate_a_random_circuit(num_qubits, num_cnots_required, seed=None):
    """Generate random circ: Generates a random circuit with num_qubits, 5 rz gates and num_cnots
    The circuit uses h, s, cz, cx, swap gates and no measurements. 

    Returns:
        QuantumCircuit: constructed circuit
    """
    assert num_cnots_required>0, f"{num_cnots_required} needs to be greater than 0."
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

    #Transpile to expand the swaps to CNOTS.
    # basis_gates=['x', 'y', 'z', 's', 'sdg', 'h', 'cx']
    # qc=transpile(qc, basis_gates=basis_gates, optimization_level=0)
    # cnot_count=count_gate(qc, "cx")
    
    qc=remove_swap(qc)
    # Too many cnots so trim.
    if cnot_count> num_cnots_required:
        qc=trim(qc, cnot_count, num_cnots_required)
        assert dict(qc.count_ops())["cx"] == num_cnots_required, "the circuit doesn't have the required cnots"
    # qc=add_rz_gates_prob(num_qubits, qc, rng) #Non deterministic # of rz.
    qc=add_rz_gates_det(num_qubits, qc, rng) #Deterministic # of rz
    return qc

def remove_swap(qc):
    '''Generate random circ: Helper function. Remove swaps. Returns: QuantumCircuit'''
    cnot_count_init=dict(qc.count_ops())["cx"]
    qc_dag = circuit_to_dag(qc)
    layers= list(qc_dag.multigraph_layers())
    # Remove the necessary cnots.
    for layer in layers:
        for node in layer:
            # Remove the node as long as it's an operation.
            if node.name=="swap":
                # print("removing swap.")
                qc_dag.remove_op_node(node)
    new_circ=dag_to_circuit(qc_dag)
    assert cnot_count_init==dict(new_circ.count_ops())["cx"], "removing swap was done incorrectly."
    return new_circ

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
    new_qc=QuantumCircuit(QuantumRegister(num_qubits, name="q"))
    NUM_RZ=5
    # This number here seems to be "ok" (some circuits had no paulie checks)
    # for circuits with 200 cnots.
    # NUM_RZ=int(qc.count_ops()["cx"]/20)


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
    new_qc=QuantumCircuit(QuantumRegister(num_qubits, name="q"))

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
    elif node.name=="rz":
        new_qc.rz(node.op.params[0], node.qargs[0].index)
    elif node.name=="sx":
        new_qc.sx(node.qargs[0].index)   
        
    else:
        # We have overlooked a gate type.
        assert False, f"{node.name} gate wasn't matched in the DAG."

def generate_circuits_with_checks_to_files(number_of_qubits, cnot_count, number_of_circuits):
    '''Generates required circuits and saves them to files. Note that this doesn't
    save raw circuit, i.e., the circuit with no checks, but it can be easily recovered
    by removing the checks from the qasm file.'''
    time0=time.time()
    #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    pool=Pool(psutil.cpu_count(logical=False))
    if number_of_qubits<=5:
        PARALLEL=False
    else:
        PARALLEL=True
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    SUBDIR="data/results"
    BASE_PATH=os.path.join(CODE_DIR, SUBDIR)
    
    #Create +1 phase pauli group
    #Note that we can restrict the search for p2 to the +1 phase. The other results
    #can all be recovered by multiplying by the corresponding phase, e.g., -i. 
    pauli_table=pauli_basis(number_of_qubits)
    # pauli_group_positive=deepcopy(pauli_table)
    #Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()

    #main loop
    file_number=0
    for _ in range(number_of_circuits):
        time1=time.time()
        # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        output_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{file_number}_.txt"
        while os.path.isfile(os.path.join(BASE_PATH, output_file_name)):
            file_number+=1
            # output_file_name=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            output_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{file_number}_.txt"
        info_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{file_number}_.obj"
        qasm_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{file_number}_.qasm"
        #Random circuit.
        circ=generate_a_random_circuit(number_of_qubits, cnot_count)            
        #Save raw circuit
        write_circs_no_checks_to_qasm_file([circ], number_of_qubits, cnot_count, subdir="data/raw_circuits")

        print(circ)
        circ_operations=circ.count_ops()
        print(circ_operations)

        # max_pauli_str_p2=Array(ctypes.c_char, "+1"+"I"*number_of_qubits)

        circ_properties=CircuitProperties(number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations)
        checks_finder=checksfinder.ChecksFinder(number_of_qubits, circ_properties.circ)
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_checks_parallel(pool, pauli_labels)
            # with Pool(psutil.cpu_count(logical=False)) as pool:
        else:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_checks_linear(pauli_labels)
                # print(pauli_to_circuit(elem))
                #Test
                # cProfile.run("find_p1s_p2s(elem)", filename="rand_circ_stats.txt")
        checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
        # Combine the circuit and checks.
        if p2_weights[0]>0:
            circ_properties.circ=checksfinder.append_checks_to_circ(circ_properties, checks_properties)

        #Draw to file
        circuit_drawer(circ_properties.circ, filename=os.path.join(BASE_PATH, output_file_name))        
        output_file=open(os.path.join(BASE_PATH, output_file_name), "a")
        output_file.write("\n")
        # Count the number of CNOT gates
        output_file.write(json.dumps(circ_operations))

        checksfinder.write_outputs(circ_properties, checks_properties, file_number, os.path.join(BASE_PATH, info_file_name),os.path.join(BASE_PATH, qasm_file_name), output_file)

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    # if PARALLEL:
    pool.close() #Stops passing jobs to processes.
    pool.join() #Waits for processes to finish.
    print("done")

def generate_random_circuits(number_of_qubits, cnot_count, number_of_circuits):
    '''Returns a list of random cricuits each with the proper cnot count.'''
    print("running...")
    circuits=[]
    for _ in range(number_of_circuits):
        circ=generate_a_random_circuit(number_of_qubits, cnot_count)
        circuits.append(circ)
        print(circ)
    print("done")
    return circuits

def generate_random_circuits_to_qasm_files(number_of_qubits, cnot_count, number_of_circuits, subdir="data"):
    '''Generates random cricuits each with the proper cnot count and saves them
    as qasm files. Generated circuits contain no checks.'''
    circuits=[]
    for _ in range(number_of_circuits):
        circuits.append(generate_a_random_circuit(number_of_qubits, cnot_count))
    write_circs_no_checks_to_qasm_file(circuits, number_of_qubits, cnot_count, subdir)

def write_circs_no_checks_to_qasm_file(circuits, number_of_qubits, cnot_count, subdir="data"):
    '''Saves the circuits into qasm files in the specified subdirectory. Note that
    that this is for saving circuits that have no checks. This constraint exists
    because of the file names used.
    circuits: iterable'''
    print("running...")
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # subdir="data"
    BASE_PATH=os.path.join(CODE_DIR, subdir)
    circ_number=0
    for circ in circuits:
        # Increment the circuit number until we find one that's not taken.
        qasm_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{circ_number}_raw_.qasm"
        while os.path.isfile(os.path.join(BASE_PATH, qasm_file_name)):
            circ_number+=1
            # output_file_name=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            qasm_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{circ_number}_raw_.qasm"
        print(qasm_file_name)
        circ.qasm(filename=os.path.join(BASE_PATH, qasm_file_name))
    print("done")
    