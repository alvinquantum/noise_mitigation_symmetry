#!/usr/bin/env python3

import numpy as np
from qiskit import QuantumCircuit
from matplotlib import pyplot as plt

if __name__ == "__main__":
    print("hello")
    circ = QuantumCircuit(5)
    circ.h(0)
    circ.h(3)
    circ.cx(0,1)
    circ.cx(0,2)
    circ.cx(3,4)
    circ.draw('mpl')
    plt.show()