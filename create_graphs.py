#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt
import operator

def get_lines(file_path):
    '''Return the lines of file but in reverse order.'''
    output=[]
    with open(file_path, "r", errors="ignore") as file:
        lines=file.readlines()
        for line in reversed(lines):
            if isinstance(line, str):
                # print(line)
                # print(type(line))
                output.append(line)
    return output

# def get_cnot_counts(file_lines):
#     '''Get CNOT counts'''
#     for line in file_lines:
#         if line.startswith("CNOT count"):
#             print(line)
#             return

def get_count(file_lines, str):
    '''Get count associated with str'''
    for line in file_lines:
        if line.startswith(str):
            return line[len(str):]
    return "0"

# File path stuff
CODE_DIR=os.path.dirname(os.path.realpath('__file__'))
SUBDIR="/data/"
PATH=CODE_DIR+SUBDIR
os.chdir(PATH)

def create_cnot_plot():
    cnot_counts=[]
    max_weight=[]
    for file in os.listdir():
        if file.endswith(".txt") and not file.endswith("result.txt"):
            print(file)
            file_lines_reverse=get_lines(file)
            cnot_counts.append(get_count(file_lines_reverse, "CNOT count:  "))
            max_weight.append(get_count(file_lines_reverse, "Max Weight:  "))
    cnot_counts=list(map(int, cnot_counts))
    max_weight=list(map(int, max_weight))
    # X axis is cnot_counts and y is max_weight of Pauli P2.
    plt.scatter(cnot_counts,max_weight)
    plt.xlabel("CNOT count")
    plt.ylabel("Pauli Weight")
    plt.savefig(PATH+"weight_vs_cnot_count.png")
    plt.show()

def create_fidelity_plot():
    fidelity_no_checks=[]
    fidelity_with_checks=[]
    for file in os.listdir():
        if file.endswith("result.txt"):
            print(file)
            file_lines_reverse=get_lines(file)
            fidelity_no_checks.append(get_count(file_lines_reverse, "State fidelity no checks with errors:  "))
            fidelity_with_checks.append(get_count(file_lines_reverse, "State fidelity with checks and errors:  "))
    fidelity_no_checks=list(map(float, fidelity_no_checks))
    fidelity_with_checks=list(map(float, fidelity_with_checks))
    # X axis is cnot_counts and y is max_weight of Pauli P2.
    plt.scatter(fidelity_no_checks, list(map(operator.sub,fidelity_with_checks,fidelity_no_checks)))
    plt.xlabel("Fidelity no checks")
    plt.ylabel("Fidelity gain")
    plt.savefig(PATH+"fidelity.png")
    plt.show()

if __name__ == "__main__":
    create_cnot_plot()
    create_fidelity_plot()
