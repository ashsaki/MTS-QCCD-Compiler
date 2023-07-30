import numpy as np
import subprocess as sp
import os
import time


# PROG=["programs/sup64_new.qasm"]
PROG=["programs/maxcut.qasm"] # qaoa
# PROG=["programs/square_root_clean_new.qasm"] # sqrt
# PROG=["programs/qft64_trimmed.qasm"] #qft
# PROG=["programs/quadratic_form_clean.qasm"] # quadratic_form


output_file = open('output.log','w')

MACHINE=["L6"]

IONS = ["15"]

lookahead = "1"
mapper = "Greedy"
reorder = "Naive"
p = PROG[0]
num_iter = 1
moves_arr = np.zeros(num_iter, dtype=int)
f = "num_moves.txt"
if os.path.exists(f):
    os.remove(f)

for idx in range(num_iter):
    for m in MACHINE:
        for i in IONS:
            ts = time.time()
            sp.call(["python", "run.py", p, m, i, mapper, reorder, "1", "0", "0", "FM", "GateSwap", lookahead], stdout=output_file)
            te = time.time()
            print(f"Benchmark {p} time taken {te-ts} sec")

            with open(f, "r") as rd_file:
                moves = rd_file.read()[:-1]
                moves_arr[idx] = int(moves)
            print(f"Iter {idx} Moves {moves}")

print(f"MIN {moves_arr.min()} MAX {moves_arr.max()} AVG {moves_arr.mean()}")
