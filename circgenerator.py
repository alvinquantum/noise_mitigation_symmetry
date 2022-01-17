# from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
# import utilities
# import pickle

class CircuitProperties:
    '''Circuit properties holder.'''
    __slots__=["number_of_qubits", "cnot_count", "number_of_circuits", "circ", "circ_operations"]
    def __init__(self, number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations,):
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count
        self.number_of_circuits=number_of_circuits
        self.circ=circ
        self.circ_operations=circ_operations

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
    new_qc=QuantumCircuit(QuantumRegister(num_qubits, name="q"))
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