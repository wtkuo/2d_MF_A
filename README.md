# Mutual Information between conduction electron and impurity electron in 2d Mean Field Single Impurity Anderson Model   
This code is to calculate the MI between c-electron and f-electron in 2d Mean Field single impurity Anderson model. 

The lattice is square lattice. The impurity is beside the center lattice point and the transition between conduction site and impurity is local and only happens at the center lattice point. The chemical potential for conduction electrons is half-filled case.

My EE calculation is based on the local correlation functions to construct the reduced density matrix. Therefore, there are three main codes.

eecal is to calculate entanglement entropy from a given correlation function matrix (eecal(Correlation matrix)). 

A2d_model is to construct the Hamiltonian and correlation function matrix between every sites and impurity.

Third one is to choose the region for us to calculate the mutual information.



