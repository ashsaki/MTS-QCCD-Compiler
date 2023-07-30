# MTS-QCCD-Compiler
Muzzle-the-Shuttle (MTS) Compiler for linearly connected multi-trap QCCD Trapped Ion Quantum Computers.

The repository contains the compiler code for the paper "Muzzle the Shuttle: Efficient Compilation for Multi-Trap Trapped-Ion Quantum Computers" (Accepted in 2022 Design, Automation and Test in Europe (DATE) Conference and Exhibition)  (Link: https://ieeexplore.ieee.org/abstract/document/9774619).

The MTS compiler is based on the excellent work of Murali et al. "Architecting noisy intermediate-scale trapped ion quantum computers" (ISCA'2020)
  - Paper Link: https://dl.acm.org/doi/10.1109/ISCA45697.2020.00051
  - GitHub Link: https://github.com/prakashmurali/QCCDSim

The MTS compiler improves upon the previous QCCDSim compiler from Murali et al. by including three heuristic architectural policies:
  - Future operations-based shuttle direction policy
  - Opportunistic gate re-ordering, and
  - Nearest neighbor first trap re-balancing

The codebase was last tested with
```
python           3.10.12
networkx         2.5.0 (important)
numpy            1.25.0
decorator        5.1.1
joblib           1.3.1
scikit-learn     1.3.0
scipy            1.11.1
sklearn          0.0
threadpoolctl    3.2.0
qiskit           0.44.0
```

To run, go to the `mts` directory and run the following command inside the directory:
```
python run_batch.py
```

Please cite the work using the following:
```
@INPROCEEDINGS{saki2022mts,
  author={Saki, Abdullah Ash and Topaloglu, Rasit Onur and Ghosh, Swaroop},
  booktitle={2022 Design, Automation & Test in Europe Conference & Exhibition (DATE)}, 
  title={Muzzle the Shuttle: Efficient Compilation for Multi-Trap Trapped-Ion Quantum Computers}, 
  year={2022},
  volume={},
  number={},
  pages={322-327},
  doi={10.23919/DATE54114.2022.9774619}
}
```

```
A. A. Saki, R. O. Topaloglu and S. Ghosh, "Muzzle the Shuttle: Efficient Compilation for Multi-Trap Trapped-Ion Quantum Computers," 2022 Design, Automation & Test in Europe Conference & Exhibition (DATE), Antwerp, Belgium, 2022, pp. 322-327, doi: 10.23919/DATE54114.2022.9774619.
```
