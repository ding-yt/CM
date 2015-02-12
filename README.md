# CM

Tumor growth simulation basic model 2D

main.cpp: main
ParFile.cpp: parameter file processing
randgen.cpp: random number generater (by Jeet)
CellType.cpp: define properties of different type of cells
GenealogyNode.cpp: a single cell/agent/genealogy node, with a pointer to its parent
Lattices.cpp: 2D space, container for genealogyNode, cell fate and interaction control
Topology.cpp: tree structure, pointers to both parents and children, sequence simulation based on lineage

######
Tree.cpp: duplicate functions with topology.cpp
SeqGen.cpp: duplicate functions with topology.cpp



