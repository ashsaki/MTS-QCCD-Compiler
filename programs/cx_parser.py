nq = 64
# filename = f"qv{nq}.qasm"
filename = "quadratic_form"
with open(f"{filename}.qasm", "r") as f_rd:
    with open(f"{filename}_clean.qasm", "w") as f_wr:
        for line in f_rd:
            if line.startswith(("cx", "qreg", "OPENQASM", "include")):
                f_wr.write(line)