import utilities
from qiskit.converters import circuit_to_dag
from qiskit.circuit import QuantumRegister
from qiskit import QuantumCircuit
import pickle
import itertools
import time
import sys

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
        self.layer_idx=0

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

    check_result=ChecksResult(p2_weight, p1_str, p2_str)
    return check_result

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

class ChecksResult:
    def __init__(self, p2_weight, p1_str, p2_str):
        self.p2_weight= p2_weight
        self.p1_str=p1_str
        self.p2_str=p2_str

class ChecksFinder:
    '''Finds checks symbolically.'''
    def __init__(self, number_of_qubits, circ):
        self.circ_reversed = circ.inverse()
        self.number_of_qubits= number_of_qubits

    def get_checks_linear(self, pauli_labels):
        count=0
        p2_weights=[0]
        pauli_str_p1s=[""]
        pauli_str_p2s=[""]
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=True)
        for elem in paulis_by_weight:
            # temp_p2_circ=find_p2s(elem, circ)
            checks_result=self.find_checks_sym(elem)
            if checks_result:
                count+=1
                append_result(checks_result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                # Terminate the pool when found is true, i.e., found=1.
                
                if checks_result.p2_weight==self.number_of_qubits:
                    break
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def _feed(self, pauli_labels, chunk_size):
        '''Generator of p2s. TODO: Instead of passing the pauli_labels
        we should generate the pauli_labels here and yield the ones needed.'''
        test_paulis=[]
        for pauli_elem in pauli_labels:
            test_paulis.append(pauli_elem)
            if len(test_paulis)==chunk_size:
                yield test_paulis
                #Reset the list
                test_paulis=[]

    def get_checks_parallel(self, pool, pauli_labels):
        '''Wrapper function. Get the checks in parallel.'''
        chunk_size=5
        count=0
        p2_weights=[0]
        pauli_str_p1s=[""]
        pauli_str_p2s=[""]
        optimal_weight=self.number_of_qubits
        # Highest weight first.
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=True)
        pauli_feed=self._feed(paulis_by_weight, chunk_size)
        found_solution=False
        #In some cases pool.imap_unordered needs to be wrapped in list in order to return properly. 
        #see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
        for chunk in pauli_feed:
            if found_solution:
                break
            # Set the possible optimal_weight. Since we traverse the possible p2s backwards
            # The best possible p2 weight is given by the weight of the first element in 
            # chunk.
            test_paulis_max_weight=get_weight(chunk[0])
            if test_paulis_max_weight<optimal_weight:
                optimal_weight=test_paulis_max_weight

            for result in pool.imap_unordered(self.find_checks_sym, chunk):
                if found_solution:
                    continue
                #Store the results.
                elif result:
                    count+=1
                    append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Terminate the pool when found is true, i.e., found=1.
                    if result.p2_weight>=optimal_weight:
                        print("terminating...")
                        found_solution=True
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def find_checks_sym(self, pauli_group_elem):
        '''Finding checks: Symbolic: Finds p1 and p2 elements symbolically.'''
        circ_reversed=self.circ_reversed
        print(pauli_group_elem)

        # We will just iterate over the +1 phase elements of the pauli group since the 
        # scenarios can be recovered by just multiplying by the phase constant.
        pauli_group_elem_ops=list(pauli_group_elem)
        p2=CheckOperator(1, pauli_group_elem_ops)
        p1=CheckOperator(1, ["I" for _ in range(len(pauli_group_elem))])
        temp_check_reversed=TempCheckOperator(1, list(reversed(pauli_group_elem_ops)))

        # Iterate through the circuit. We manually keep track of the idx since
        # we can either go forward or backwards. This is kept track of inside the temp_check_reversed
        # We also track layer_idx in temp_check_reversed.
        # forward=True
        circ_dag = circuit_to_dag(circ_reversed)
        layers = list(circ_dag.multigraph_layers())
        num_layers=len(layers)
        # # We start index 1 since the first layer are just out nodes.
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
                    p1.phase=temp_check_reversed.phase
                    p1.operations=list(reversed(temp_check_reversed.operations))
                    # Append operations.
                    # with count.get_lock():
                    checks_result =get_check_strs(p1, p2)
                    return checks_result
                else:
                    temp_check_reversed.layer_idx+=1            
            else:
                if temp_check_reversed.layer_idx==1:
                    # We reached the first layer of operation nodes so move forward. Note
                    # the zero index are all input nodes so we can skip.
                    p2.phase=temp_check_reversed.phase
                    p2.operations=list(reversed(temp_check_reversed.operations))
                    temp_check_reversed.forward=True
                else:
                    temp_check_reversed.layer_idx-=1

def append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s):
    p2_weight=result.p2_weight
    p1_str=result.p1_str
    p2_str=result.p2_str
    if p2_weight>p2_weights[0]:
        p2_weights[0]=p2_weight
        pauli_str_p1s[0]=p1_str
        pauli_str_p2s[0]=p2_str
    p2_weights.append(p2_weight)
    pauli_str_p1s.append(p1_str)
    pauli_str_p2s.append(p2_str)

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
    quant_comp_reg=QuantumRegister(number_of_qubits, name="q")
    ancilla_reg=QuantumRegister(1, name="a")
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)

    # temp_circ.h(ancilla_reg)
    utilities.add_controlU(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    utilities.add_controlU(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg)
    # temp_circ.h(ancilla_reg)
    
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