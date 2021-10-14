#!/usr/bin/env python3

import os, sys
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
CODE_DIR=sys.path[0]
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
            for item in file_lines_reverse:
                if item.startswith("CNOT count: "):
                    cnot_counts.append(item.split(" ")[2].strip("\n"))
            # for item in file_lines_reverse:
                elif item.startswith("Max Weight: "):
                    max_weight.append(item.split(" ")[2].strip("\n"))
    # print(len(cnot_counts))
    # print(len(max_weight))
            # print(file_lines_reverse)
            # print(next(obj for obj in file_lines_reverse if obj.startswith("CNOT count: ")))
    #         cnot_counts.append(next(obj for obj in file_lines_reverse if obj.startswith("CNOT count: ")).split(" ")[1])
    #         max_weight.append(next(obj for obj in file_lines_reverse if obj.startswith("Max weight: ")).split(" ")[1])
    cnot_counts=list(map(int, cnot_counts))
    max_weight=list(map(int, max_weight))
    # print(cnot_counts)
    # print(max_weight)
    # X axis is cnot_counts and y is max_weight of Pauli P2.
    plt.boxplot(cnot_counts,max_weight)
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
    plt.boxplot(fidelity_no_checks, list(map(operator.sub,fidelity_with_checks,fidelity_no_checks)))
    plt.xlabel("Fidelity no checks")
    plt.ylabel("Fidelity gain")
    plt.savefig(PATH+"fidelity.png")
    plt.show()

if __name__ == "__main__":
    create_cnot_plot()
    # create_fidelity_plot()
