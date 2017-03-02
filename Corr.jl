function Corr(n,C,M)
# This subroutine is to construct the two-point correlation function matrix <c*_{i}c_{j}> for ground state
# From Pauli exclusion principle, we know only one electron can live in a certain energy state. In our case,
# since there is no spin interaction term, the spectrum is 2-fold degenerate. If d_{i} and d*_{i} are the 
# annihlation and creation operator for energy eigenstate i. <d*{i}d{j}> = \delta_{ij} for i,j <= Fermi energy
# Fermi energy will be determined by the total number of electrons which is related to chemical potential.
# Here, I focus on the half-filled case. For each site, it only has one electron and the maximum number for 
# each site is two. To connect d operator and c operator, we need the matrix(M) which can diagonalize our 
# Hamiltonian. c{i} = \sum^{n}_{1}M_{ij}d{j} => <c*_{i}c_{j}> = M*_{ik}M_{jk}<d*_{k}d_{k}> = M*_{ik}M_{jk}
# 
# Parameter 
# n : the number of the conduction sites
# C : the full two-point correlation function matrix
# M : the matrix which can diagonalize M^{T}*H*M = E => d = M^{T}c => c = M d

for i = 1:2n+2
    for j = 1:n
       C[i,i] += M[i,j]* M[i,j]; 
    end
end

for i = 1:2n+1
    for j = i+1:2n+2
            
            for k = 1:n
                C[i,j] += M[i,k]* M[j,k];
            end
        
        C[j,i] = C[i,j];
    end
end

end
