import utilities
from qiskit.converters import circuit_to_dag
from qiskit.circuit import QuantumRegister, AncillaRegister, ClassicalRegister
from qiskit import QuantumCircuit, transpile
import pickle
from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from multiprocessing import Pool
import circgenerator
import psutil
import os
import time
import math
from qiskit.visualization import circuit_drawer


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
        # self.change_to_backwards=False
        # self.forward=True
        self.layer_idx=1

class PushOperator:
    '''Finding checks: Symbolic: push operations.'''
    # track of the direction at the start of the layer.
    @staticmethod
    def x(op2, temp_check_reversed):
        '''Finding checks: Symbolic: Push x through op2. 
        Forward: Circuit: [x]--[Op2]=[Op2]--[F] Math: (Op2)(x)=(F)(Op2)-->
        (F)=(Op2)(x)(Op2)^\dagger.
        
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
            # if temp_check_reversed.forward:
            return [1, "Y"]
            # else:
                # return [-1, "Y"]
        elif op2=="SDG":
            # if temp_check_reversed.forward:
            return [-1, "Y"]
            # else:
                # return [1, "Y"]
        #Test
        # elif op2=="RZ":
        #     temp_check_reversed.change_to_backwards=True
        #     # Increment the layer_idx so that we redo this layer when going
        #     # backwards. You can see this in the decrementing logic of find_checks_sym.
        #     # temp_check_reversed.layer_idx+=1
        #     return [1, "I"]
        else:
            # We have overlooked a gate type.
            assert False, f"{op2} gate wasn't matched in the DAG."

    @staticmethod
    def y(op2, temp_check_reversed):
        '''Finding checks: Symbolic: Push y through op2.  
        Forward: Circuit: [y]--[Op2]=[Op2]--[F] Math: (Op2)(y)=(F)(Op2)-->
        (F)=(Op2)(y)(Op2)^\dagger.
        
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
            # if temp_check_reversed.forward:
            return [-1, "X"]
            # else:
            #     return [1, "X"]
        elif op2=="SDG":
            # if temp_check_reversed.forward:
            return [1, "X"]
            # else:
            #     return [-1, "X"]
        # elif op2=="RZ":
        #     temp_check_reversed.change_to_backwards=True
        #     # Increment the layer_idx so that we redo this layer when going
        #     # backwards. You can see this in the decrementing logic of find_checks_sym.
        #     # temp_check_reversed.layer_idx+=1
        #     return [1, "I"]
        else:
            # We have overlooked a gate type.
            assert False, f"{op2} gate wasn't matched in the DAG." 

    @staticmethod        
    def z(op2):
        '''Finding checks: Symbolic: Push z through op2. 
        Forward: Circuit: [z]--[Op2]=[Op2]--[F] Math: (Op2)(z)=(F)(Op2)-->
        (F)=(Op2)(z)(Op2)^\dagger.
        
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

class ChecksResult:
    def __init__(self, p2_weight, p1_str, p2_str):
        self.p2_weight= p2_weight
        self.p1_str=p1_str
        self.p2_str=p2_str

class ChecksFinder:
    '''Finds checks symbolically.'''
    def __init__(self, number_of_qubits, circ):
        # We are pushing p2 through U to get p1. Thus,
        # Circuit: [p1] [U] [p2] := [U].
        # Math : (p2)(U)(p1) = (U) -> (p1) = (U)^dagger(p2)(U).
        #
        # However, in the code we are pushing p2 from left to right through U.
        # We need to match the math given above. To do this, let
        # Circuit: [p2] [W] [p1] = [W], so that we are pushing left to right
        # and [W] is some circuit that we solve for to match the math above.
        # Math: (p1)(W)(p2) = (W) -> (p1) = (W)(p2)(W)^dagger.
        # Then, notice that (W)=(U)^\dagger:=[U].inverse().
        self.circ_reversed = circ.inverse()
        self.number_of_qubits= number_of_qubits

    @staticmethod
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
                # May be possible to expand to accomodate non pauli+I in the future.
                assert False, op1[0] + " is not I, X, Y, or Z."
        else:
            # Two qubit operations
            if op2=="CX":
                # print(f"op1 before pushing: {op1}")
                result = PushOperator.cx(op1)
                # print(f"result of pushing: {result}")
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

    @staticmethod
    def update_current_ops_approx(op1, op2, temp_check_reversed, current_qubits, rz_angle):
        '''Finding checks: Symbolic: Finds the intermediate check. Always push op1 through op2. '''
        # Rz = cos(lambda/2)I-isin(lambda/2)Z. We approximate by 
        # taking the operator with the larger magnitude coefficient. 
        # Note that the phases do not matter in determining the c_1 check
        # because they cancel out.
        # op2= Identity case:
        rz_angle=rz_angle[0]
        if abs(math.cos(rz_angle/2))>=abs(math.sin(rz_angle/2)):
            result=[1]
            result=result+op1
        # op2 = Z case:
        else:
            op2="Z"
            if op1[0]=="X":
                result= PushOperator.x(op2, temp_check_reversed)
            elif op1[0]=="Y":
                result= PushOperator.y(op2, temp_check_reversed)
            elif op1[0]=="Z":
                result= PushOperator.z(op2)
            elif op1[0]=="I":
                result= [1, "I"]
            else:
                # May be possible to expand to accomodate non pauli+I in the future.
                assert False, op1[0] + " is not I, X, Y, or Z."
        # else:
        #     # Two qubit operations
        #     if op2=="CX":
        #         # print(f"op1 before pushing: {op1}")
        #         result = PushOperator.cx(op1)
        #         # print(f"result of pushing: {result}")
        #     elif op2=="SWAP":
        #         result= PushOperator.swap(op1)
        #     else:
        #         assert False, op2 + " is not cx or swap."

        result_phase=result[0]
        result_ops=result[1::1]
        temp_check_reversed.phase=temp_check_reversed.phase*result_phase
        # Coppy the current ops into temp_check_reversed.
        for idx, op in enumerate(result_ops):
            temp_check_reversed.operations[current_qubits[idx]]=op

    @staticmethod
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
        
        # print(f"p1: {p1_str}")
        # print(f"p2: {p2_str}")
        # print(f"Pauli weight P2: {p2_weight}")
        # print()

        check_result=ChecksResult(p2_weight, p1_str, p2_str)
        return check_result

    @staticmethod
    def can_continue(op1, op2):
        '''Finding checks: Symbolic: Determine if can continue. If we're going backwards and op1 is not I or Z and op2 is RZ then don't continue.'''
        if op2=="RZ" and op1!="I" and op1!="Z":
            return False
        else:
            return True
        # if forward==False and op2=="RZ" and op1!="I" and op1!="Z":
        #     return False
        # else:
        #     return True


    @staticmethod
    def get_current_qubits(node):
        '''Finding checks: Symbolic: get the current qubits whose operations that will be passed through.'''
        # We have to check for single or two qubit gates.
        if node.name in ["x", "y", "z", "h", "s", "sdg", "rz"]:
            return [node.qargs[0].index]
        elif node.name in ["cx", "swap"]:
            return [node.qargs[0].index, node.qargs[1].index]
        else:
            assert False, "Overlooked a node operation."

    def get_min_checks_linear_approx(self, pauli_labels):
        count=0
        p2_weights=[]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=False)
        # print(paulis_by_weight)
        for elem in paulis_by_weight:
            # temp_p2_circ=find_p2s(elem, circ)
            if get_weight(elem)>0:
                checks_result=self.find_checks_sym_approx(elem)
                if checks_result:
                    count+=1
                    self._append_result(checks_result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    if count==3*self.number_of_qubits:
                        break
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def get_max_checks_linear(self, pauli_labels, max_layers):
        print('searching...')
        count=0
        p2_weights=[]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=True)
        # print(paulis_by_weight)
        for elem in paulis_by_weight:
            # temp_p2_circ=find_p2s(elem, circ)
            if get_weight(elem)>0:
                checks_result=self.find_checks_sym(elem)
                if checks_result:
                    count+=1
                    self._append_result(checks_result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    if count==max_layers:
                        break
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def get_min_checks_linear(self, pauli_labels, max_layers):
        print('searching...')
        count=0
        p2_weights=[]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=False)
        # print(paulis_by_weight)
        for elem in paulis_by_weight:
            # temp_p2_circ=find_p2s(elem, circ)
            if get_weight(elem)>0:
                checks_result=self.find_checks_sym(elem)
                if checks_result:
                    count+=1
                    self._append_result(checks_result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    if count==max_layers:
                        break
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def get_max_checks_parallel(self, pool, pauli_labels, max_layers):
        '''Wrapper function. Get the checks in parallel.'''
        print('searching...')
        chunk_size=100
        count=0
        p2_weights=[0]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        # Lowest weight first.
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
            # test_paulis_max_weight=get_weight(chunk[0])
            # if test_paulis_max_weight<optimal_weight:
            #     optimal_weight=test_paulis_max_weight

            for result in pool.imap(self.find_checks_sym, chunk):
                if found_solution:
                    continue
                #Store the results.
                elif result:
                    count+=1
                    self._append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    # if result.p2_weight>=optimal_weight:
                    if count>=max_layers:
                        print("terminating...")
                        found_solution=True
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def get_min_checks_parallel(self, pool, pauli_labels, max_layers):
        '''Wrapper function. Get the checks in parallel.'''
        print('searching...')
        chunk_size=100
        count=0
        p2_weights=[0]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        # Lowest weight first.
        paulis_by_weight=sorted(pauli_labels, key=get_weight, reverse=False)
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
            # test_paulis_max_weight=get_weight(chunk[0])
            # if test_paulis_max_weight<optimal_weight:
            #     optimal_weight=test_paulis_max_weight

            for result in pool.imap(self.find_checks_sym, chunk):
                if found_solution:
                    continue
                #Store the results.
                elif result:
                    count+=1
                    self._append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    # if result.p2_weight>=optimal_weight:
                    if count>=max_layers:
                        print("terminating...")
                        found_solution=True
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def get_min_checks_parallel_clifford(self, pool, pauli_labels):
        '''Wrapper function. Get the checks in parallel.'''
        chunk_size=1
        count=0
        p2_weights=[0]
        pauli_str_p1s=[]
        pauli_str_p2s=[]
        # optimal_weight=self.number_of_qubits
        paulis_by_weight=pauli_labels#sorted(pauli_labels, key=get_weight, reverse=False)
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
            # test_paulis_max_weight=get_weight(chunk[0])
            # if test_paulis_max_weight<optimal_weight:
            #     optimal_weight=test_paulis_max_weight

            for result in pool.imap(self.find_checks_sym, chunk):
                if found_solution:
                    continue
                #Store the results.
                elif result:
                    count+=1
                    self._append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s)
                    # Only get this number of layers. This limit (which works for weight one checks)
                    # comes from the fact that we only need an X and Y C_2 for each compute qubit.
                    # All the Pauli terms anticommute with X or Y. However, the search may
                    # pickup an X, Y, and Z for one qubit so we need more.
                    # if result.p2_weight>=optimal_weight:
                    if count>=3*self.number_of_qubits:
                        print("terminating...")
                        found_solution=True
        return count, p2_weights, pauli_str_p1s, pauli_str_p2s

    def _feed(self, pauli_labels, chunk_size):
        '''Generator of p2s. TODO: Instead of passing the pauli_labels
        we should generate the pauli_labels here and yield the ones needed.'''
        test_paulis=[]
        for pauli_elem in pauli_labels:
            if get_weight(pauli_elem)>=1:
                test_paulis.append(pauli_elem)
            if len(test_paulis)==chunk_size:
                yield test_paulis
                #Reset the list
                test_paulis=[]
        # We went through all the items but there might
        # not be enough to reach the chunk_size. yield
        # what we have.
        if test_paulis:
            yield test_paulis

    def find_checks_sym_approx(self, pauli_group_elem):
        '''Finding checks: Symbolic: Finds p1 and p2 elements symbolically.'''
        circ_reversed=self.circ_reversed
        # print(pauli_group_elem)

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
                    current_qubits=self.get_current_qubits(node)
                    # print(f"current qubits: {current_qubits}")
                    current_ops=[temp_check_reversed.operations[qubit] for qubit in current_qubits]
                    # print(f"current ops: {current_ops}")
                    node_op= node.name.upper()

                    # Update temp_check_reversed if possible
                    if self.can_continue(current_ops[0], node_op):
                        self.update_current_ops(current_ops, node_op, temp_check_reversed, current_qubits)
                    else:
                        rz_angle=node.op.params
                        self.update_current_ops_approx(current_ops, node_op, temp_check_reversed, current_qubits, rz_angle)

            # Reached the last layer so return result
            if temp_check_reversed.layer_idx==num_layers-1:
                p1.phase=temp_check_reversed.phase
                p1.operations=list(reversed(temp_check_reversed.operations))
                # Append operations.
                # with count.get_lock():
                checks_result =self.get_check_strs(p1, p2)
                # Exit point with checks.
                return checks_result
            else:
                temp_check_reversed.layer_idx+=1

            # # See if we should start going backwards.
            # if temp_check_reversed.change_to_backwards:
            #     temp_check_reversed.forward=False
            #     temp_check_reversed.change_to_backwards=False
            #     # we don't increment or decrement the layer_idx since we processed this layer.
            #     # We have to process the same layer going backwards.

            # # Since we're not changing to backwards, either move forward or backards
            # elif temp_check_reversed.forward:
            #     if temp_check_reversed.layer_idx==num_layers-1:
            #         p1.phase=temp_check_reversed.phase
            #         p1.operations=list(reversed(temp_check_reversed.operations))
            #         # Append operations.
            #         # with count.get_lock():
            #         checks_result =self.get_check_strs(p1, p2)
            #         # Exit point with checks.
            #         return checks_result
            #     else:
            #         temp_check_reversed.layer_idx+=1            
            # else:
            #     if temp_check_reversed.layer_idx==1:
            #         # We reached the first layer of operation nodes so move forward. Note
            #         # the zero index are all input nodes so we can skip.
            #         p2.phase=temp_check_reversed.phase
            #         p2.operations=list(reversed(temp_check_reversed.operations))
            #         temp_check_reversed.forward=True
            #     else:
            #         temp_check_reversed.layer_idx-=1  

    def find_checks_sym(self, pauli_group_elem):
        '''Finding checks: Symbolic: Finds p1 and p2 elements symbolically.'''
        circ_reversed=self.circ_reversed
        # print(pauli_group_elem)

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
                    current_qubits=self.get_current_qubits(node)
                    # print(f"current qubits: {current_qubits}")
                    current_ops=[temp_check_reversed.operations[qubit] for qubit in current_qubits]
                    # print(f"current ops: {current_ops}")
                    node_op= node.name.upper()

                    # Update temp_check_reversed if possible
                    if self.can_continue(current_ops[0], node_op):
                        self.update_current_ops(current_ops, node_op, temp_check_reversed, current_qubits)
                    else:
                        #Exit point with no checks.
                        return None

            # Reached the last layer so return result
            if temp_check_reversed.layer_idx==num_layers-1:
                p1.phase=temp_check_reversed.phase
                p1.operations=list(reversed(temp_check_reversed.operations))
                # Append operations.
                # with count.get_lock():
                checks_result =self.get_check_strs(p1, p2)
                # Exit point with checks.
                return checks_result
            else:
                temp_check_reversed.layer_idx+=1

            # # See if we should start going backwards.
            # if temp_check_reversed.change_to_backwards:
            #     temp_check_reversed.forward=False
            #     temp_check_reversed.change_to_backwards=False
            #     # we don't increment or decrement the layer_idx since we processed this layer.
            #     # We have to process the same layer going backwards.

            # # Since we're not changing to backwards, either move forward or backards
            # elif temp_check_reversed.forward:
            #     if temp_check_reversed.layer_idx==num_layers-1:
            #         p1.phase=temp_check_reversed.phase
            #         p1.operations=list(reversed(temp_check_reversed.operations))
            #         # Append operations.
            #         # with count.get_lock():
            #         checks_result =self.get_check_strs(p1, p2)
            #         # Exit point with checks.
            #         return checks_result
            #     else:
            #         temp_check_reversed.layer_idx+=1            
            # else:
            #     if temp_check_reversed.layer_idx==1:
            #         # We reached the first layer of operation nodes so move forward. Note
            #         # the zero index are all input nodes so we can skip.
            #         p2.phase=temp_check_reversed.phase
            #         p2.operations=list(reversed(temp_check_reversed.operations))
            #         temp_check_reversed.forward=True
            #     else:
            #         temp_check_reversed.layer_idx-=1            

    @staticmethod
    def _append_result(result, p2_weights, pauli_str_p1s, pauli_str_p2s):
        '''Helper for appending results for found checks.'''
        p2_weight=result.p2_weight
        p1_str=result.p1_str
        p2_str=result.p2_str
        # if p2_weight>p2_weights[0]:
        #     p2_weights[0]=p2_weight
        #     pauli_str_p1s[0]=p1_str
        #     pauli_str_p2s[0]=p2_str
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
    '''Finding checks: combines everything into one circuit.
    circ_properties: circgenerator.CircuitProperties
    checks_properties: circgenerator.ChecksProperties'''
    circ=circ_properties.circ
    pauli_str_p1=checks_properties.pauli_str_p1s[0]
    pauli_str_p2=checks_properties.pauli_str_p2s[0]
    number_of_qubits=circ_properties.number_of_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits, name="q")
    ancilla_reg=QuantumRegister(1, name="a")
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)

    temp_circ.h(ancilla_reg)
    temp_circ.barrier()
    utilities.add_controlU(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    utilities.add_controlU(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
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

def write_outputs_multilayer(circ_properties, file_qasm_path):
    '''Finding checks: write outputs.'''
    circ=circ_properties.circ
    circ.qasm(filename=file_qasm_path)

def append_multilayer_checks_to_circ_clifford(circ, layers_count, check_layers):
    # print(layers_count)
    # print(check_layers)
    '''Finding checks: combines everything into one circuit.
    circ_properties: circgenerator.CircuitProperties
    checks_properties: circgenerator.ChecksProperties'''
    # layers_count=5
    # check_layers=check_layers[:layers_count]
    number_of_qubits=circ.num_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits, name="a")
    ancilla_reg=AncillaRegister(layers_count, "b")
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)

    #Add the p1s. Note that we need to reverse the checks layers when either
    # adding the p1s or p2s. Ex: the for the innermost pair (p11, p21) p11 has to be added
    # last for the p1s, but the p21 should be added first for the p2s .
    for ancilla_idx, pair in reversed(list(enumerate(check_layers))):
        pauli_str_p1=pair[0]
        # pauli_str_p2=pair[1]
        temp_circ.h(ancilla_reg[ancilla_idx])
        utilities.add_controlU_multilayer_clifford(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])

    # Add the main circuit.
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)


    #Add the p2s.
    for ancilla_idx, pair in enumerate(check_layers):
        # pauli_str_p1=pair[0]
        pauli_str_p2=pair[1]
        # print(pauli_str_p2)
        utilities.add_controlU_multilayer_clifford(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])
        # temp_circ.barrier()
        temp_circ.h(ancilla_reg[ancilla_idx])
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)
    # print(temp_circ)
    
    return temp_circ

def append_multilayer_checks_to_circ(circ, layers_count, check_layers):
    # print(layers_count)
    # print(check_layers)
    '''Finding checks: combines everything into one circuit.
    circ_properties: circgenerator.CircuitProperties
    checks_properties: circgenerator.ChecksProperties'''
    # layers_count=5
    # check_layers=check_layers[:layers_count]
    number_of_qubits=circ.num_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits, name="a")
    ancilla_reg=AncillaRegister(layers_count, "b")
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)

    #Add the p1s. Note that we need to reverse the checks layers when either
    # adding the p1s or p2s. Ex: the for the innermost pair (p11, p21) p11 has to be added
    # last for the p1s, but the p21 should be added first for the p2s .
    for ancilla_idx, pair in reversed(list(enumerate(check_layers))):
        pauli_str_p1=pair[0]
        # pauli_str_p2=pair[1]
        temp_circ.h(ancilla_reg[ancilla_idx])
        utilities.add_controlU_multilayer(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])

    # Add the main circuit.
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)


    #Add the p2s.
    for ancilla_idx, pair in enumerate(check_layers):
        # pauli_str_p1=pair[0]
        pauli_str_p2=pair[1]
        # print(pauli_str_p2)
        utilities.add_controlU_multilayer(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])
        # temp_circ.barrier()
        temp_circ.h(ancilla_reg[ancilla_idx])
    temp_circ.barrier()
    utilities.insert_identity_layer(temp_circ)
    # print(temp_circ)
    
    return temp_circ

def append_multilayer_checks_to_circ_no_identities(circ, layers_count, check_layers):
    # print(layers_count)
    # print(check_layers)
    '''Finding checks: combines everything into one circuit.
    circ_properties: circgenerator.CircuitProperties
    checks_properties: circgenerator.ChecksProperties'''
    # layers_count=5
    # check_layers=check_layers[:layers_count]
    number_of_qubits=circ.num_qubits
    quant_comp_reg=QuantumRegister(number_of_qubits, name="a")
    ancilla_reg=AncillaRegister(layers_count, "b")
    temp_circ=QuantumCircuit(quant_comp_reg, ancilla_reg)
    temp_circ.barrier()
    # utilities.insert_identity_layer(temp_circ)

    #Add the p1s. Note that we need to reverse the checks layers when either
    # adding the p1s or p2s. Ex: the for the innermost pair (p11, p21) p11 has to be added
    # last for the p1s, but the p21 should be added first for the p2s .
    for ancilla_idx, pair in reversed(list(enumerate(check_layers))):
        pauli_str_p1=pair[0]
        # pauli_str_p2=pair[1]
        temp_circ.h(ancilla_reg[ancilla_idx])
        utilities.add_controlU_multilayer(temp_circ, pauli_str_p1, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])

    # Add the main circuit.
    temp_circ.barrier()
    # utilities.insert_identity_layer(temp_circ)
    temp_circ.compose(circ, quant_comp_reg, inplace=True)
    temp_circ.barrier()
    # utilities.insert_identity_layer(temp_circ)


    #Add the p2s.
    for ancilla_idx, pair in enumerate(check_layers):
        # pauli_str_p1=pair[0]
        pauli_str_p2=pair[1]
        # print(pauli_str_p2)
        utilities.add_controlU_multilayer(temp_circ, pauli_str_p2, number_of_qubits, quant_comp_reg, ancilla_reg[ancilla_idx])
        # temp_circ.barrier()
        temp_circ.h(ancilla_reg[ancilla_idx])
    temp_circ.barrier()
    # utilities.insert_identity_layer(temp_circ)
    # print(temp_circ)
    
    return temp_circ


def get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits):
    number_of_startingIs=qubit_idx
    number_of_finalIs=number_of_qubits-number_of_startingIs-1
    return ["I"*number_of_startingIs+"X"+"I"*number_of_finalIs,
            # "I"*number_of_startingIs+"Y"+"I"*number_of_finalIs,
            "I"*number_of_startingIs+"Z"+"I"*number_of_finalIs]

def checks_file_exists(checks_path, qasm_file_name, layers_count):
    '''Returns: True or False'''
    return os.path.isfile(os.path.join(checks_path, generate_checks_file_name_from_raw_qasm(qasm_file_name, layers_count)))

def generate_checks_file_name_from_raw_qasm(qasm_file_name, number_of_layers):
    '''Returns: Checks qasm file name or raises an error.'''
    str_ending="_raw_.qasm"
    if qasm_file_name.endswith(str_ending):
        return f"{qasm_file_name[:-len(str_ending)]}_layers_{number_of_layers}_.qasm"
    raise AssertionError(f"{qasm_file_name} does not end with '_raw_.qasm'")

def generate_obj_file_name_from_checks_file_name(output_file_name):
    '''Returns: Checks qasm file name or raises an error.'''
    str_ending=".qasm"
    if output_file_name.endswith(str_ending):
        return f"{output_file_name[:-len(str_ending)]}.obj"
    raise AssertionError(f"{output_file_name} does not end with '.qasm'")

def generate_txt_file_name_from_checks_file_name(output_file_name):
    '''Returns: Checks qasm file name or raises an error.'''
    str_ending=".qasm"
    if output_file_name.endswith(str_ending):
        return f"{output_file_name[:-len(str_ending)]}.txt"
    raise AssertionError(f"{output_file_name} does not end with '.qasm'")

def get_circuit_number(file_name):
    '''Gets the circuit number from the file name or raises an error.
    Return: type int'''
    split_name=file_name.split("_")
    for idx, elem in enumerate(split_name):
        if elem == "circuit":
            if idx+1<len(file_name):
                return split_name[idx+1]
    raise AssertionError(f"{file_name} does not contain 'circuit_*'")

def get_p2_qubit(p2):
    '''Helper function for removing redundant checks.'''
    # print(p2)
    for idx, elem in enumerate(reversed(p2[2:])):
        if elem !="I":
            # print(idx)
            return idx

def post_process(count, p2_weights, pauli_str_p1s, pauli_str_p2s):
    '''If there are three weight one p2 checks acting on the same qubit,
    remove one since it is redundant.'''
    weight_one_checks=[]
    weight_one_qubits=[]
    for p2_idx, p2 in enumerate(pauli_str_p2s):
        p2_weight=get_weight(p2[2:])
        if p2_weight==1:
            p2_qubit=get_p2_qubit(p2)
            # p2_idx is check idx in case we want to remove it.
            # p2 qubit is the qubit being acted on by the check p2.
            weight_one_checks.append((p2_idx, p2_qubit, p2))
            weight_one_qubits.append(p2_qubit)

    remove_check_for_qubit=[]
    for p2_qubit in weight_one_qubits:
        if weight_one_qubits.count(p2_qubit)==3:
            remove_check_for_qubit.append(p2_qubit)
            weight_one_qubits.remove(p2_qubit)


    remove_check_idxs=[]
    # print(remove_check_for_qubit)
    for p2_info in weight_one_checks:
        if p2_info[1] in remove_check_for_qubit:
            remove_check_idxs.append(p2_info[0])
            remove_check_for_qubit.remove(p2_info[1])
            # print(remove_check_for_qubit)

    for idx in sorted(remove_check_idxs, reverse=True):
        # print("REMOVED...")
        del p2_weights[idx] 
        del pauli_str_p1s[idx] 
        del pauli_str_p2s[idx]
        count=count-1
    # print(f"final {pauli_str_p1s} {pauli_str_p2s}")

    return count, p2_weights, pauli_str_p1s, pauli_str_p2s

def create_low_weight_checks_circs(number_of_qubits, cnot_count, target_layers, max_layers, rz_count):
    '''Main function to call to create circuits with checks. Follows a rigid folder structure. Lowe
    weight checks.
    '''
    time0=time.time()
    print("running...")
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))

    #Program parameters
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    MAIN_SUBDIR=f"qubits_{number_of_qubits}_rz_{rz_count}"
    CHECKS_SUBDIR=os.path.join(MAIN_SUBDIR,"low_weight_checks","checks")
    RAWS_SUBDIR=os.path.join(MAIN_SUBDIR, "raws")
    CHECKS_PATH=os.path.join(CODE_DIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, RAWS_SUBDIR)
    
    #Create +1 phase pauli group
    #Note that we can restrict the search for p2 to the +1 phase. The other results
    #can all be recovered by multiplying by the corresponding phase, e.g., -i.
    # This can probably be made more efficient by defining this once somewhere
    # instead of generating it everytime. 
    pauli_table=pauli_basis(number_of_qubits, pauli_list=True)
    # pauli_group_positive=deepcopy(pauli_table)
    # Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()

    if number_of_qubits<=5:
        PARALLEL=False
    else:
        PARALLEL=True
    # PARALLEL=True
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:
        # print(input_qasm_file_name)
        # print(pauli_labels)

        circ=QuantumCircuit.from_qasm_file(os.path.join(RAWS_PATH, input_qasm_file_name))
        circ_operations=circ.count_ops()
        circ_properties=circgenerator.CircuitProperties(number_of_qubits, cnot_count, NUMBER_OF_CIRCUITS, circ, circ_operations)
        checks_finder=ChecksFinder(number_of_qubits, circ_properties.circ)

        layers_count=0
        check_layers=[]
        # for qubit_idx in range(NUMBER_OF_QUBITS):
        time1=time.time()
        # pauli_labels=checksfinder.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, NUMBER_OF_QUBITS)
        # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_parallel(pool, pauli_labels, max_layers)
            # with Pool(psutil.cpu_count(logical=False)) as pool:
        else:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels, max_layers)

        if count>0:
            # Remove one of the single weight checks if there are three acting on the same qubits.
            count, p2_weights, pauli_str_p1s, pauli_str_p2s= post_process(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            count=0 #Reset the count. Since we only max of 2 checks.
            for idx, elem in enumerate(pauli_str_p1s):
                if idx<=target_layers-1:
                    check_layers.append((elem, pauli_str_p2s[idx]))
                    count+=1
            layers_count+=count

        if layers_count>0:
            # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            # loop through the layers.
            for final_layers_idx in range(layers_count):
                temp_layers_count=final_layers_idx+1
                # Skip the raw qasm if we already have an output.
                if checks_file_exists(CHECKS_PATH, input_qasm_file_name, temp_layers_count):
                    continue

                circ_properties.circ=append_multilayer_checks_to_circ(circ, temp_layers_count, check_layers[0:temp_layers_count])
                
                # if final_layers_idx==layers_count-1:
                print(circ_properties.circ)
                #Draw to file
                output_qasm_file_name=generate_checks_file_name_from_raw_qasm(input_qasm_file_name, temp_layers_count)
                # print(output_qasm_file_name)
                # info_file_name=generate_obj_file_name_from_checks_file_name(output_qasm_file_name)
                txt_file_name=generate_txt_file_name_from_checks_file_name(output_qasm_file_name)
                circuit_drawer(circ_properties.circ, filename=os.path.join(CHECKS_PATH, txt_file_name))        
                output_txt_file=open(os.path.join(CHECKS_PATH, txt_file_name), "a")
                # circuit_number=get_circuit_number(output_qasm_file_name)
                write_outputs_multilayer(circ_properties, os.path.join(CHECKS_PATH, output_qasm_file_name))
                output_txt_file.close()
                    
                print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print("done")

def create_high_weight_checks_circs(number_of_qubits, cnot_count, target_layers, max_layers, rz_count):
    '''Main function to call to create circuits with checks. Follows a rigid folder structure. High
    weight checks.
    '''
    time0=time.time()
    print("running...")
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))

    #Program parameters
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    MAIN_SUBDIR=f"qubits_{number_of_qubits}_rz_{rz_count}"
    CHECKS_SUBDIR=os.path.join(MAIN_SUBDIR,"high_weight_checks","checks")
    RAWS_SUBDIR=os.path.join(MAIN_SUBDIR, "raws")
    CHECKS_PATH=os.path.join(CODE_DIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, RAWS_SUBDIR)
    
    #Create +1 phase pauli group
    #Note that we can restrict the search for p2 to the +1 phase. The other results
    #can all be recovered by multiplying by the corresponding phase, e.g., -i.
    # This can probably be made more efficient by defining this once somewhere
    # instead of generating it everytime. 
    pauli_table=pauli_basis(number_of_qubits, pauli_list=True)
    # pauli_group_positive=deepcopy(pauli_table)
    # Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()

    if number_of_qubits<=5:
        PARALLEL=False
    else:
        PARALLEL=True
    # PARALLEL=True
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:
        # print(input_qasm_file_name)
        # print(pauli_labels)

        circ=QuantumCircuit.from_qasm_file(os.path.join(RAWS_PATH, input_qasm_file_name))
        circ_operations=circ.count_ops()
        circ_properties=circgenerator.CircuitProperties(number_of_qubits, cnot_count, NUMBER_OF_CIRCUITS, circ, circ_operations)
        checks_finder=ChecksFinder(number_of_qubits, circ_properties.circ)

        layers_count=0
        check_layers=[]
        # for qubit_idx in range(NUMBER_OF_QUBITS):
        time1=time.time()
        # pauli_labels=checksfinder.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, NUMBER_OF_QUBITS)
        # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
        # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
        # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
        if PARALLEL:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_max_checks_parallel(pool, pauli_labels, max_layers)
            # with Pool(psutil.cpu_count(logical=False)) as pool:
        else:
            count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_max_checks_linear(pauli_labels, max_layers)

        if count>0:
            # Remove one of the single weight checks if there are three acting on the same qubits.
            count, p2_weights, pauli_str_p1s, pauli_str_p2s= post_process(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            count=0 #Reset the count. Since we only max of 2 checks.
            for idx, elem in enumerate(pauli_str_p1s):
                if idx<=target_layers-1:
                    check_layers.append((elem, pauli_str_p2s[idx]))
                    count+=1
            layers_count+=count

        if layers_count>0:
            # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            # loop through the layers.
            for final_layers_idx in range(layers_count):
                temp_layers_count=final_layers_idx+1
                # Skip the raw qasm if we already have an output.
                if checks_file_exists(CHECKS_PATH, input_qasm_file_name, temp_layers_count):
                    continue

                circ_properties.circ=append_multilayer_checks_to_circ(circ, temp_layers_count, check_layers[0:temp_layers_count])
                
                # if final_layers_idx==layers_count-1:
                print(circ_properties.circ)
                #Draw to file
                output_qasm_file_name=generate_checks_file_name_from_raw_qasm(input_qasm_file_name, temp_layers_count)
                # print(output_qasm_file_name)
                # info_file_name=generate_obj_file_name_from_checks_file_name(output_qasm_file_name)
                txt_file_name=generate_txt_file_name_from_checks_file_name(output_qasm_file_name)
                circuit_drawer(circ_properties.circ, filename=os.path.join(CHECKS_PATH, txt_file_name))        
                output_txt_file=open(os.path.join(CHECKS_PATH, txt_file_name), "a")
                # circuit_number=get_circuit_number(output_qasm_file_name)
                write_outputs_multilayer(circ_properties, os.path.join(CHECKS_PATH, output_qasm_file_name))
                output_txt_file.close()
                    
                print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print("done")

def create_weight_one_check_circs(number_of_qubits, cnot_count, target_layers, rz_count):
    '''Main function to call to create circuits with checks. Follows a rigid folder structure. 
    Weight one checks.
    '''
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))

    time0=time.time()
    #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=f"qubits_{number_of_qubits}_rz_{rz_count}"
    CHECKS_SUBDIR=os.path.join("low_weight_checks","checks")
    RAWS_SUBDIR="raws"
    CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)
    
    # if number_of_qubits<=5:
    #     PARALLEL=False
    # else:
    #     PARALLEL=True
    PARALLEL=True
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:

        # Skip the raw qasm if we already have an output.
        # if checksfinder.checks_file_exists(CHECKS_PATH, input_qasm_file_name):
        #     continue

        circ=QuantumCircuit.from_qasm_file(os.path.join(RAWS_PATH, input_qasm_file_name))
        circ_operations=circ.count_ops()
        circ_properties=circgenerator.CircuitProperties(number_of_qubits, cnot_count, NUMBER_OF_CIRCUITS, circ, circ_operations)
        checks_finder=ChecksFinder(number_of_qubits, circ_properties.circ)
        # print(circ)
        # # Iterate over the qubits and try to get two p2 checks for each qubit
        layers_count=0
        check_layers=[]
        for qubit_idx in range(number_of_qubits):
            time1=time.time()
            pauli_labels=get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
            # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
            # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
            if PARALLEL:
            # if False:
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_parallel_clifford(pool, pauli_labels)
                # with Pool(psutil.cpu_count(logical=False)) as pool:
            else:
                # count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)


            if count>0:
                # count=0 #Reset the count. Since we only max of 2 checks.
                for idx, elem in enumerate(pauli_str_p1s):
                    # if idx<=1:
                    check_layers.append((elem, pauli_str_p2s[idx]))
                        # count+=1
                layers_count+=count
        if layers_count>0:
            # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            # loop through the layers.
            increment_size=1
            start_size=0
            # for final_layers_idx in list(range(layers_count))[4:MAX_LAYERS:5]:
            #     temp_layers_count=final_layers_idx+1
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ(circ, temp_layers_count, check_layers[0:temp_layers_count])
            sets_of_check_layers=[check_layers[0:final_layers_idx+1] for final_layers_idx in list(range(layers_count))[start_size:target_layers:increment_size]]
            circs=[circ]*len(sets_of_check_layers)
            temp_layers_counts=[len(elem) for elem in sets_of_check_layers]
            result_circs =pool.starmap(append_multilayer_checks_to_circ_clifford, zip(circs, temp_layers_counts, sets_of_check_layers))
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ_clifford(circ, temp_layers_count, check_layers[0:temp_layers_count])

            # print(circ_properties.circ)
            for idx, circ in enumerate(result_circs):
                circ_properties.circ=circ
                temp_layers_count=temp_layers_counts[idx]
                #Draw to file
                output_qasm_file_name=generate_checks_file_name_from_raw_qasm(input_qasm_file_name, temp_layers_count)
                print(output_qasm_file_name)
                info_file_name=generate_obj_file_name_from_checks_file_name(output_qasm_file_name)
                txt_file_name=generate_txt_file_name_from_checks_file_name(output_qasm_file_name)
                circuit_drawer(circ_properties.circ, filename=os.path.join(CHECKS_PATH, txt_file_name))        
                output_txt_file=open(os.path.join(CHECKS_PATH, txt_file_name), "a")
                circuit_number=get_circuit_number(output_qasm_file_name)
                write_outputs_multilayer(circ_properties, os.path.join(CHECKS_PATH, output_qasm_file_name))
                
                output_txt_file.close()
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print("done")

def create_weight_one_check_circs_example(number_of_qubits, cnot_count, target_layers, rz_count):
    '''Main function to call to create circuits with checks. Follows a rigid folder structure. 
    Weight one checks.
    '''
    NUMBER_OF_CIRCUITS=1#int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))

    time0=time.time()
    #Disable qiskit parallel.
    os.environ['QISKIT_IN_PARALLEL'] = 'TRUE'
    print("running...")
    #Program parameters
    # func_params=FunctionParams()
    #Paths for outputs and pickle file of circuit. sys.path[0] on laptop and the other on hpc.
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    MAIN_SUBDIR=os.path.join(f"qubits_{number_of_qubits}_rz_{rz_count}", "noiseless_checks_ex")
    CHECKS_SUBDIR=os.path.join("low_weight_checks","checks")
    RAWS_SUBDIR="raws"
    CHECKS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, CHECKS_SUBDIR)
    RAWS_PATH=os.path.join(CODE_DIR, MAIN_SUBDIR, RAWS_SUBDIR)
    
    # if number_of_qubits<=5:
    #     PARALLEL=False
    # else:
    #     PARALLEL=True
    PARALLEL=True
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))
    #Load qasm files and find checks.
    input_qasm_file_names=utilities.get_files_from_dir_by_extension(RAWS_PATH, ".qasm")
    input_qasm_file_names=circgenerator.filter_by_cnot_qubit(input_qasm_file_names, cnot_count, number_of_qubits)
    for input_qasm_file_name in input_qasm_file_names:

        # Skip the raw qasm if we already have an output.
        # if checksfinder.checks_file_exists(CHECKS_PATH, input_qasm_file_name):
        #     continue

        circ=QuantumCircuit.from_qasm_file(os.path.join(RAWS_PATH, input_qasm_file_name))
        circ_operations=circ.count_ops()
        circ_properties=circgenerator.CircuitProperties(number_of_qubits, cnot_count, NUMBER_OF_CIRCUITS, circ, circ_operations)
        checks_finder=ChecksFinder(number_of_qubits, circ_properties.circ)
        # print(circ)
        # # Iterate over the qubits and try to get two p2 checks for each qubit
        layers_count=0
        check_layers=[]
        for qubit_idx in range(number_of_qubits):
            time1=time.time()
            pauli_labels=get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
            # output_file_path=BASE_FILE_PATH+ "circuit_" + str(file_number) +"_.txt"
            # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
            # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
            if PARALLEL:
            # if False:
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_parallel_clifford(pool, pauli_labels)
                # with Pool(psutil.cpu_count(logical=False)) as pool:
            else:
                # count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)
                count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels)


            if count>0:
                # count=0 #Reset the count. Since we only max of 2 checks.
                for idx, elem in enumerate(pauli_str_p1s):
                    # if idx<=1:
                    check_layers.append((elem, pauli_str_p2s[idx]))
                        # count+=1
                layers_count+=count
        if layers_count>0:
            # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
            # loop through the layers.
            increment_size=1
            start_size=0
            # for final_layers_idx in list(range(layers_count))[4:MAX_LAYERS:5]:
            #     temp_layers_count=final_layers_idx+1
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ(circ, temp_layers_count, check_layers[0:temp_layers_count])
            sets_of_check_layers=[check_layers[0:final_layers_idx+1] for final_layers_idx in list(range(layers_count))[start_size:target_layers:increment_size]]
            circs=[circ]*len(sets_of_check_layers)
            temp_layers_counts=[len(elem) for elem in sets_of_check_layers]
            result_circs =pool.starmap(append_multilayer_checks_to_circ_clifford, zip(circs, temp_layers_counts, sets_of_check_layers))
            # circ_properties.circ=checksfinder.append_multilayer_checks_to_circ_clifford(circ, temp_layers_count, check_layers[0:temp_layers_count])

            # print(circ_properties.circ)
            for idx, circ in enumerate(result_circs):
                circ_properties.circ=circ
                temp_layers_count=temp_layers_counts[idx]
                #Draw to file
                output_qasm_file_name=generate_checks_file_name_from_raw_qasm(input_qasm_file_name, temp_layers_count)
                print(output_qasm_file_name)
                info_file_name=generate_obj_file_name_from_checks_file_name(output_qasm_file_name)
                txt_file_name=generate_txt_file_name_from_checks_file_name(output_qasm_file_name)
                circuit_drawer(circ_properties.circ, filename=os.path.join(CHECKS_PATH, txt_file_name))        
                output_txt_file=open(os.path.join(CHECKS_PATH, txt_file_name), "a")
                circuit_number=get_circuit_number(output_qasm_file_name)
                write_outputs_multilayer(circ_properties, os.path.join(CHECKS_PATH, output_qasm_file_name))
                
                output_txt_file.close()
        print(f"file execution time {time.time()-time1}")
    print(f"total execution time {time.time()-time0}")
    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print("done")
    
############################
def get_error_mitig_circ(qasmstr, target_layers):
    '''Main function to call to create circuits with low weight checks.
    qasmstr is the input circuit (the input circuit can't have measurements or barriers 
    and the qubits must be on a single quantum register) and target_layers is the 
    desired number of Pauli check layers. Returns an error mitigated qasm str 
    with the number of Pauli check layers <= target_layers. Note that to implement 
    the protocol you still have to postselect on the ancilla 
    measurements (Only keep the state when you measure all zeros on the 
    ancillas. Discard otherwise.)
    
    The number of input qubits (tested with up to 10) can't be too many since we generate 
    a list of all the possible +1Pauli checks to iterate through. 
    This memory limitation can easily be removed by using a generator instead of a list, 
    but is not yet implemented. 
    
    qasmstr: str
    target_layers: int
    returns: str or None if nothing found
    '''
    time0=time.time()
    print("running...")
    circ=QuantumCircuit.from_qasm_str(qasmstr)
    number_of_qubits=circ.num_qubits
    max_layers=target_layers+target_layers/3
    NUMBER_OF_CIRCUITS=1

    # Transpile to the base set we are using.
    basis_gates=["x", "y", "z", "h", "s", "sdg", "rz", "cx", "id"]
    circ=transpile(circ, basis_gates=basis_gates, optimization_level=0)
    cnot_count=dict(circ.count_ops())["cx"]
    
    #Create +1 phase pauli group
    #Note that we can restrict the search for p2 to the +1 phase. The other results
    #can all be recovered by multiplying by the corresponding phase, e.g., -i.
    # This can probably be made more efficient by defining this once somewhere
    # instead of generating it everytime. 
    pauli_table=pauli_basis(number_of_qubits, pauli_list=True)
    # pauli_group_positive=deepcopy(pauli_table)
    # Labels will be used to print pauli strings in the loop
    pauli_labels=pauli_table.to_labels()

    if number_of_qubits<=5:
        PARALLEL=True
    else:
        PARALLEL=True
    # PARALLEL=True
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    circ_operations=circ.count_ops()
    circ_properties=circgenerator.CircuitProperties(number_of_qubits, cnot_count, NUMBER_OF_CIRCUITS, circ, circ_operations)
    checks_finder=ChecksFinder(number_of_qubits, circ_properties.circ)

    layers_count=0
    check_layers=[]

    # Doing pool this way is faster when the circuits become large since the cpus will be fully utilized
    # each time. If we parallelize across individual circuits, each generation of circuit will be slow.
    if PARALLEL:
        count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_parallel(pool, pauli_labels, max_layers)
        # with Pool(psutil.cpu_count(logical=False)) as pool:
    else:
        count, p2_weights, pauli_str_p1s, pauli_str_p2s=checks_finder.get_min_checks_linear(pauli_labels, max_layers)

    if count>0:
        # Remove one of the single weight checks if there are three acting on the same qubits.
        count, p2_weights, pauli_str_p1s, pauli_str_p2s= post_process(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
        count=0 #Reset the count. Since we only max of 2 checks.
        for idx, elem in enumerate(pauli_str_p1s):
            if idx<=target_layers-1:
                check_layers.append((elem, pauli_str_p2s[idx]))
                count+=1
        layers_count+=count

    if layers_count>0:
        # checks_properties=checksfinder.ChecksProperties(count, p2_weights, pauli_str_p1s, pauli_str_p2s)
        # loop through the layers.
        for final_layers_idx in range(layers_count):
            temp_layers_count=final_layers_idx+1

            #If you want identity gates after the barriers use this func instead: append_multilayer_checks_to_circ
            circ_properties.circ=append_multilayer_checks_to_circ_no_identities(circ, temp_layers_count, check_layers[0:temp_layers_count])
            
            # print(circ_properties.circ)
        circ=circ_properties.circ
        ancillas=circ.ancillas
        classical_reg=ClassicalRegister(len(ancillas))
        circ.add_register(classical_reg)
        for idx, qancilla in enumerate(ancillas):
            circ.measure(qancilla,classical_reg[idx])
        print(circ)
        circ_qasm=circ.qasm()

    if PARALLEL:
        pool.close() #Stops passing jobs to processes.
        pool.join() #Waits for processes to finish.
    print(f"total execution time {time.time()-time0}")
    print("done")
    if layers_count:
        return circ_qasm
    else:
        return None