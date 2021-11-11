#!/usr/bin/env python3
import os, sys, shutil

def rename_file(file):
    '''adds _ between numbers and strings except before "."'''
    fn_old=file.split("_")
    fn_new=""
    for elem_old in fn_old:
        if elem_old.startswith("qubits"):
            new_elem=elem_old[:6]+"_"+elem_old[6:]+"_"
        elif elem_old.startswith("depth"):
            new_elem=elem_old[:5]+"_"+elem_old[5:]+"_"
        elif elem_old.startswith("erroridx"):
            new_elem=elem_old[:8]+"_"+elem_old[8:]+"_"
        elif elem_old.startswith("result"):
            new_elem=elem_old[:6]+"_"+elem_old[6:]
        elif elem_old!=fn_old[-1]:
            new_elem=elem_old+"_"
        else:
            new_elem=elem_old
        fn_new+=new_elem
    # Add _ before period.
    temp = fn_new.split(".")
    return temp[0]+"_."+temp[1]

if __name__ == "__main__":
    print("running...")
    CODE_DIR=sys.path[0]
    OLD_DIR=CODE_DIR+"\qubits5_depth_generated_old\\"
    NEW_DIR=CODE_DIR+"\qubits5_depth_generated_renamed\\"
    files=os.listdir(OLD_DIR)
    for fn_old in files:
        if fn_old.endswith(".png"):
            shutil.copy(os.path.join(OLD_DIR+fn_old), os.path.join(NEW_DIR+fn_old))
        else:
            fn_new=rename_file(fn_old)
            # print(fn_new)
            shutil.copy(os.path.join(OLD_DIR+fn_old), os.path.join(NEW_DIR+fn_new))
    print("done.")