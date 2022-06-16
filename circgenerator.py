import os
import numpy as np
import utilities
from qiskit import QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
from copy import deepcopy
import time
from qiskit.opflow import X,Y,Z

class CircuitProperties:
    '''Circuit properties holder.'''
    __slots__=["number_of_qubits", "cnot_count", "number_of_circuits", "circ", "circ_operations"]
    def __init__(self, number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations,):
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count
        self.number_of_circuits=number_of_circuits
        self.circ=deepcopy(circ)
        self.circ_operations=circ_operations

def generate_a_random_clifford_circuit(num_qubits, num_cnots_required, seed=None):
    """Generate random Clifford circ: Generates a random circuit with num_qubits and num_cnots
    The circuit uses h, s, cz, cx, swap gates and no measurements. 

    Returns:
        QuantumCircuit: constructed circuit
    """
    assert num_cnots_required>0, f"{num_cnots_required} needs to be greater than 0."
    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    # rng = np.random.default_rng(seed)

    qc=QuantumCircuit(QuantumRegister(num_qubits))
    cnot_count=0
    #Randomly generate cliffords until we have more than enough cnots
    #then truncate
    while cnot_count < num_cnots_required:
        qc_temp=decompose_clifford(random_clifford(num_qubits))
        # qc_temp=random_clifford(num_qubits).to_instruction()
        cnot_count+=count_gate(qc_temp, "cx")
        qc.compose(qc_temp, inplace=True)
    
    qc=remove_swap(qc)
    # Too many cnots so trim.
    if cnot_count> num_cnots_required:
        qc=trim(qc, cnot_count, num_cnots_required)
        assert dict(qc.count_ops())["cx"] == num_cnots_required, "the circuit doesn't have the required cnots"
    # qc=add_rz_gates_prob(num_qubits, qc, rng) #Non deterministic # of rz.
    # qc=add_rz_gates_det(num_qubits, qc, rng) #Deterministic # of rz
    return qc

def generate_a_random_circuit(num_qubits, num_cnots_required, rz_count, seed=None):
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
    
    qc=remove_swap(qc)
    # Too many cnots so trim.
    if cnot_count> num_cnots_required:
        qc=trim(qc, cnot_count, num_cnots_required)
        assert dict(qc.count_ops())["cx"] == num_cnots_required, "the circuit doesn't have the required cnots"
    # qc=add_rz_gates_prob(num_qubits, qc, rng) #Non deterministic # of rz.
    qc=add_rz_gates_det(num_qubits, qc, rng, rz_count) #Deterministic # of rz
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

def add_rz_gates_det(num_qubits, qc, rng, rz_count):
    '''Generate random circ: Helper function. Add a set number of rz gates randomly to the circuit. Returns: QuantumCircuit'''
    qc_dag=circuit_to_dag(qc)
    new_qc=QuantumCircuit(QuantumRegister(num_qubits, name="q"))
    # NUM_RZ=rz_count

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
    if rz_count>len(all_nodes):
        rz_count=len(all_nodes)
        chosen_nodes=all_nodes
    else:
        chosen_nodes=rng.choice(all_nodes, replace=False, size=rz_count).tolist()

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
    assert dict(new_qc.count_ops())["rz"] == rz_count, "Number of rz's wasn't met."
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

def add_rz_gates_prob(num_qubits, qc, rng, prob):
    '''Generate random circ: Helper function. Add rz gates randomly to the circuit. Returns: QuantumCircuit'''
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

def write_circs_no_checks_to_qasm_file(base_path, circuits, number_of_qubits, cnot_count):
    '''Saves the circuits into qasm files in the specified subdirectory. Note that
    that this is for saving circuits that have no checks.
    circuits: iterable'''
    circ_number=0
    for circ in circuits:
        # Increment the circuit number until we find one that's not taken.
        qasm_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{circ_number}_raw_.qasm"
        while os.path.isfile(os.path.join(base_path, qasm_file_name)):
            circ_number+=1
            qasm_file_name=f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_circuit_{circ_number}_raw_.qasm"
        print(qasm_file_name)
        circ.qasm(filename=os.path.join(base_path, qasm_file_name))
    print("done")

def filter_by_cnot_qubit(files, cnot_count, number_of_qubits):
    '''Only return files that have the specified cnot count and number of qubits.'''
    filtered_files=[]
    for file in files:
        file_nums_info=[int(elem) for elem in file.split("_") if elem.isdigit()]
        if file_nums_info[0]== number_of_qubits and file_nums_info[1]==cnot_count:
            filtered_files.append(file)
    return filtered_files
    
def gen_rand_circs(number_of_qubits, cnot_count, number_of_circuits, rz_count):
    '''Main function to call to create random input circs.'''
    time0=time.time()

    print("running...")
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{number_of_qubits}_rz_{rz_count}"
    # CHECKS_SUBDIR="checks"
    RAWS_SUBDIR="raws"
    # CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)

    # Generate random circuits
    # file_number=0
    for _ in range(number_of_circuits):
        time1=time.time()
        # circ=generate_a_random_circuit(number_of_qubits, cnot_count) 
        if rz_count>0:
            circ=generate_a_random_circuit(number_of_qubits, cnot_count, rz_count)             
        else:
            circ=generate_a_random_clifford_circuit(number_of_qubits, cnot_count)  
        #Save raw circuit
        write_circs_no_checks_to_qasm_file(RAWS_PATH, [circ], number_of_qubits, cnot_count)

        # print(circ)
        circ_operations=circ.count_ops()
        print(circ_operations)
        
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")

def gen_rand_input_state(number_of_qubits, cnot_count, rz_count):
    '''Main function to call to create random input states (input states
    are given by a random circuit).'''
    time0=time.time()
    print("running...")
    #Program parameters
    #Paths for outputs and inputs.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{number_of_qubits}_rz_{rz_count}"
    CHECKS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "checks")
    # RESULTS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "results")
    INITIAL_STATES_SUBDIR=os.path.join(MAIN_SUBDIR, "initial_states")
    RAWS_SUBDIR=os.path.join(MAIN_SUBDIR, "raws")
    CHECKS_PATH=os.path.join(CODE_DIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, RAWS_SUBDIR)
    # RESULTS_PATH=os.path.join(CODE_DIR, RESULTS_SUBDIR)
    INITIAL_STATES_PATH=os.path.join(CODE_DIR, INITIAL_STATES_SUBDIR)

    #We iterate over all the raw qasm files. For each file, we generate an input state.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:
        time1=time.time()

        # We strip the ending of the file name and add the appropriate ending.
        end_string="_raw_.qasm"
        output_file_name=f"{input_qasm_file_name[0:-len(end_string)]}_inputstate_0_.qasm"
        print(output_file_name)
        # Skip the raw qasm if we already have an output.
        if init_circ_exists(CHECKS_PATH, input_qasm_file_name):
            continue
        circ=make_rand_input_state_multilayer(number_of_qubits)
        circ.qasm(filename= os.path.join(INITIAL_STATES_PATH, output_file_name))

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    print("done")

def gen_rand_input_state_example(number_of_qubits, cnot_count, rz_count):
    '''Main function to call to create random input states (input states
    are given by a random circuit).'''
    time0=time.time()
    print("running...")
    #Program parameters
    #Paths for outputs and inputs.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=os.path.join(f"qubits_{number_of_qubits}_rz_{rz_count}","noiseless_checks_ex")
    CHECKS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "checks")
    # RESULTS_SUBDIR=os.path.join(MAIN_SUBDIR, "min_weight_checks", "results")
    INITIAL_STATES_SUBDIR=os.path.join(MAIN_SUBDIR, "initial_states")
    RAWS_SUBDIR=os.path.join(MAIN_SUBDIR, "raws")
    CHECKS_PATH=os.path.join(CODE_DIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, RAWS_SUBDIR)
    # RESULTS_PATH=os.path.join(CODE_DIR, RESULTS_SUBDIR)
    INITIAL_STATES_PATH=os.path.join(CODE_DIR, INITIAL_STATES_SUBDIR)

    #We iterate over all the raw qasm files. For each file, we generate an input state.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:
        time1=time.time()

        # We strip the ending of the file name and add the appropriate ending.
        end_string="_raw_.qasm"
        output_file_name=f"{input_qasm_file_name[0:-len(end_string)]}_inputstate_0_.qasm"
        print(output_file_name)
        # Skip the raw qasm if we already have an output.
        if init_circ_exists(CHECKS_PATH, input_qasm_file_name):
            continue
        circ=make_rand_input_state_multilayer(number_of_qubits)
        circ.qasm(filename= os.path.join(INITIAL_STATES_PATH, output_file_name))

        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    print("done")

def gen_rand_circs_example(number_of_qubits, cnot_count, number_of_circuits, rz_count):
    '''Main function to call to create random input circs.'''
    time0=time.time()

    print("running...")
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=os.path.join(f"qubits_{number_of_qubits}_rz_{rz_count}", "noiseless_checks_ex")
    # CHECKS_SUBDIR="checks"
    RAWS_SUBDIR="raws"
    # CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)

    # Generate random circuits
    # file_number=0
    for _ in range(number_of_circuits):
        time1=time.time()
        # circ=generate_a_random_circuit(number_of_qubits, cnot_count) 
        if rz_count>0:
            circ=generate_a_random_circuit(number_of_qubits, cnot_count, rz_count)             
        else:
            circ=generate_a_random_clifford_circuit(number_of_qubits, cnot_count)  
        #Save raw circuit
        write_circs_no_checks_to_qasm_file(RAWS_PATH, [circ], number_of_qubits, cnot_count)

        # print(circ)
        circ_operations=circ.count_ops()
        print(circ_operations)
        
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")

def make_rand_input_state_multilayer(number_of_compute_qubits):
    '''Testing circuits: Create a random state. Need to send both circs at the same time so they have the same random initial state.'''
    #Insert random state generator
    # temp_circ=deepcopy(circ)
    # total_number_of_qubits=circ.num_qubits
    quantum_register=QuantumRegister(number_of_compute_qubits, "a")
    # ancilla_register=AncillaRegister(total_number_of_qubits-number_of_compute_qubits)
    # idenity_circ=QuantumCircuit(quantum_register, ancilla_register)
    # for qubit_idx in range(total_number_of_qubits):
        # idenity_circ.i(qubit_idx)
    # print(idenity_circ)
    circ=QuantumCircuit(quantum_register)
    random_params = np.random.uniform(size=(number_of_compute_qubits, 3))
    for i,qreg in enumerate(quantum_register):
        for j, pauli in enumerate([X,Y,Z]):
            rand_tuple=(random_params[i][j])
            circ.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
            # circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
    return circ

def init_circ_exists(base_path, init_file_name):
    '''Checks if an initial state circuit exists.'''
    return os.path.isfile(os.path.join(base_path, init_file_name))