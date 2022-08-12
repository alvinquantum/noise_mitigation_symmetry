# Quantum Error Mitigation by Pauli Check Sandwiching
arXiv paper: [link](https://doi.org/10.48550/arXiv.2206.00215)

Important files:
- circgenerator.py: module for generating random circuits.  
- checksfinder.py: module for finding checks for a given circuit.  
- circtester.py: module for testing circuits.
- test_units.py: unit tests.
- tested with qiskit 0.36.1 and cirq 0.14.1.

To produce an error mitigated circuit use the function 
checksfinder.get_error_mitig_circ(qasmstr, target_layers). qasmstr is the input 
circuit (the input circuit can't have measurements or barriers and the qubits 
must be on a single quantum register) and target_layers is the desired number of 
Pauli check layers. Returns an error mitigated qasm str with the number of Pauli 
check layers <= target_layers or if no checks are found returns None.
Note that to implement the protocol you still have 
to postselect on the ancilla measurements (Only keep the state when you 
measure all zeros on the ancillas. Discard otherwise.)

The number of input qubits (tested with up to 10) can't be too many since we generate 
a list of all the possible +1Pauli checks to iterate through. This memory limitation 
can easily be removed by using a generator instead of a list, but is not yet implemented. 

To replicate the graphs in the paper extract the contents of data zip folder to 
the base directory of the project. Then execute exec_gen_test.py.

To redo all the numerical experiments from scratch including the random input circuits
delete all the file contents (leave the folder structure intact) from the unzipped data. 
Uncomment out line 113 in exec_gen_test.py and run. Note that this will take an extremely long time
(multiple days).
