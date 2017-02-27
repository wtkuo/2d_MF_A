# This script is to generate a 2D single impurity Anderson model
include("eecal.jl");

n = 9;          # the number of conduction electon sites
N = 2*(n+1);    # the size of the full Hamiltonian 2x(n c.sites + 1 f.site)
l = Int(sqrt(n));
t = 1;          # transition amplitude
H = zeros(N,N); # Hamiltonian
Eu = 0;
Ed = 0;
V = 0;
C = zeros(N,N); # Correlation function matrix

;# Set up the c-electron

for i = 1:n-1

        if  i % l == 0
            continue 
        end

        H[i,i+1] = t;
        H[i+1,i] = t;
        H[i+n,i+n+1] = t;
        H[i+n+1,i+n] = t;

end

for i = 1:l-1
        
        H[i,i+l] = t;
        H[i+l,i] = t;
        H[i+n,i+n+l] = t;
        H[i+n+l,i+n] = t;

end

;# Set up the f-electron

H[2n+1,2n+1] = Eu;
H[2n+2,2n+2] = Ed;

;# Set up the interaction between c-electron and f-electron only at the center of the lattice

c = round(Int,n/2)+1;
H[2n+1,c] = V;
H[c,2n+1] = V;
H[2n+2,c+n] = V;
H[c+n,2n+2] = V;


E,M = eig(H);

for i = 1:N
    for j = 1:n
       C[i,i] += M[i,j]* M[i,j]; 
    end
end

for i = 1:N-1
    for j = i+1:N
        for k = 1:n
            C[i,j] += M[i,k]* M[j,k];
        end
        C[j,i] = C[i,j];
    end
end




