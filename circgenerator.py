# from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.converters.dag_to_circuit import dag_to_circuit
from qiskit.quantum_info import random_clifford, decompose_clifford
from qiskit.converters import circuit_to_dag
import utilities
import pickle

class CircuitProperties:
    '''Circuit properties holder.'''
    __slots__=["number_of_qubits", "cnot_count", "number_of_circuits", "circ", "circ_operations"]
    def __init__(self, number_of_qubits, cnot_count, number_of_circuits, circ, circ_operations,):
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count
        self.number_of_circuits=number_of_circuits
        self.circ=circ
        self.circ_operations=circ_operations

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
            assert False, f"{op2} gate wasn't matched in the DAG."

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
            assert False, f"{op2} gate wasn't matched in the DAG." 

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
            assert False, f"{op2} gate wasn't matched in the DAG." 

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
            assert False, f"{op1[0]} , {op1[1]} wasn't a pauli element." 

    @staticmethod
    def swap(op1):
        '''Finding checks: Symbolic: pass op1 through.
        In return: first digit is phase.'''
        result_ops=list(reversed(op1))
        result=[1]
        return result+result_ops

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
    elif node.name=="rz":
        new_qc.rz(node.op.params[0], node.qargs[0].index)
    elif node.name=="sx":
        new_qc.sx(node.qargs[0].index)   
        
    else:
        # We have overlooked a gate type.
        assert False, f"{node.name} gate wasn't matched in the DAG."

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
    
    print(f"p1: {p1_str}")
    print(f"p2: {p2_str}")
    print(f"Pauli weight P2: {p2_weight}")
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

def get_weight(pauli_string):
    '''Gets the weight of a Pauli string. Returns: int'''
    count=0
    for character in pauli_string:
        if character!="I":
            count+=1
    return count

def append_checks_to_circ(circ_properties, checks_properties):
    '''Finding checks: combines everything into one circuit.'''
    circ=circ_properties.circ
    pauli_str_p1=checks_properties.pauli_str_p1s[0]
    pauli_str_p2=checks_properties.pauli_str_p2s[0]
    number_of_qubits=circ_properties.number_of_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits)
    ancilla_reg=QuantumRegister(1)
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)

    temp_circ.h(ancilla_reg)
    utilities.add_controlU(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    utilities.add_controlU(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.h(ancilla_reg)
    
    return temp_circ

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
                output_file.write(f"p1: {strp1}\n")
                output_file.write(f"p2: {pauli_str_p2s[index1]}\n")
                output_file.write(f"Pauli weight P2: {p2_weights[index1]}\n")
        output_file.write("\n")
        output_file.write(f"Found Matches: {count}\n")
        output_file.write(f"Max Weight: {p2_weights[0]}\n")
        output_file.write(f"P1 that creates max P2: {pauli_str_p1s[0]}\n")
        output_file.write(f"Max P2: {str(pauli_str_p2s[0])}\n")
        print(f"Found Matches: {count}")   
        print(f"Max P2 Weight: {p2_weights[0]}")
        print(f"P1 that creates max P2: {pauli_str_p1s[0]}")
        print(f"Max P2: {pauli_str_p2s[0]}")

    cnot_count=0
    if "cx" in circ_operations:
        cnot_count=circ_operations["cx"]

    rz_count=0
    if "rz" in circ_operations:
        rz_count=circ_operations["rz"]
    output_file.write(f"Circuit no: {file_number}\n")
    output_file.write(f"Qubits: {number_of_qubits}\n")
    output_file.write(f"CNOT count: {cnot_count}\n")
    output_file.write(f"RZ count: {rz_count}\n")

    print(f"Qubits: {number_of_qubits}")
    print(f"Circuit no: {file_number}")
    print(f"CNOT count: {cnot_count}")
    print(f"RZ count: {rz_count}")

    # # Dump all the info into a pickle
    circ_file=open(file_info_path, "wb")
    pickle.dump({"cx": cnot_count, "rz": rz_count, "qubits": number_of_qubits, "circuit_num" : file_number, "found_matches": count,
            "max_pauli_weight": p2_weights[0], "max_pauli_str_p1": pauli_str_p1s[0], "max_pauli_str_p2": pauli_str_p2s[0]}, circ_file)
    circ.qasm(filename=file_qasm_path)
    # # Close the files. 
    circ_file.close()
    output_file.close()