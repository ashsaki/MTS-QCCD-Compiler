import sys
# m = 1454
# m = 970
m = sys.argv[1]
rd_file = f"output_{m}moves.log"
wr_file = f"output_{m}moves_parsed.log"

with open(wr_file, "w") as fw:
    with open(rd_file, "r") as fr:
        moves = 0
        total = 0
        flag = 0
        for line in fr:
            if line.startswith("MOV ["):
                moves += 1
                t = [i for i, ltr in enumerate(line) if ltr == "T"]
                trap_a, trap_b = int(line[t[0]+1]), int(line[t[1]+1])
                total += abs(trap_a - trap_b)

            if line.startswith("#") and flag == 0:
                fw.write(f"Moves: {moves} Total {total}\n")
                moves = 0
                total = 0
                flag = 1
                continue

            if line.startswith("#") and flag == 1:
                flag = 0
                continue

            # fw.write(line)
