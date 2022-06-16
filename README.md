# Quantum Error Mitigation by Pauli Check Sandwiching
arXiv paper: [link](https://doi.org/10.48550/arXiv.2206.00215)

Important files:
- circgenerator.py: module for generating random circuits.  
- checksfinder.py: module for finding checks for a given circuit.  
- circtester.py: module for testing circuits.
- test_units.py: unit tests.
- tested with qiskit 0.36.1 and cirq 0.14.1.

To produce an error mitigated circuit use the function 
checksfinder.get_error_mitig_circ(qasmstr, target_layers). qasmstr
is the input qasm string and target_layers is the number of check layers desired. 
Returns an error mitigated qasm string.

To replicate the graphs in the paper extract the contents of data zip folder to 
the base directory of the project. Then execute exec_gen_test.py.

To redo all the numerical experiments from scratch including the random input circuits
delete all the file contents (leave the folder structure intact) from the unzipped data. 
Uncomment out line 113 in exec_gen_test.py and run. Note that this will take an extremely long time
(multiple days).
