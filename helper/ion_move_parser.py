import numpy as np
import re

filename = "../output.log"

def find(s, ch="T"):
    return [i for i, ltr in enumerate(s) if ltr == ch]

nqubits = 78
moves = np.zeros(nqubits, dtype=int)
with open(filename, "r") as rd_file:
    for line in rd_file:
        if line.startswith("MOV"):
            result = re.search(r"\[([A-Za-z0-9_]+)\]", line)
            ion = int(result.group(1))
            indices_of_T = find(line, "T")
            traps = [int(line[i+1]) for i in indices_of_T]
            moves[ion] += abs(traps[1] - traps[0])

mov_tuples = []
for ion, nmoves in enumerate(moves):
    # print(f"Ion {ion}, moved {nmoves} times")
    mov_tuples.append((ion, nmoves))

mov_tuples.sort(key=lambda x:x[1], reverse=True)

for ion, nmoves in mov_tuples:
    print(f"Ion {ion}, moved {nmoves} times")

print(f"total moves {sum(moves)}")
print(f"max moves {max(moves)}")