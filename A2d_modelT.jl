# This script is to generate a 2D single impurity Anderson model on a square lattice
# To simplify, the possible transition between c-electron and f-electron happens at the center of the lattice (only one site) and the impurity site with transition amplitude V.
using Plots;
include("eecal.jl");
include("findcorT.jl");


n = 169;      # the number of conduction electon sites (WARNING: this number should be the square of ODD number.)
N = 2*(n+1);     # the size of the full Hamiltonian (n c.sites + 1 f.site) x 2 (spin up and spin down)
l = Int(sqrt(n));# the width of the side
iter = Int((l+1)/2);  # the radius range
t = 1;           # transition amplitude between the neighboring sites for conduction electron
H = zeros(N,N);  # Hamiltonian matrix; 1~n: spin up, n+1~2n: spin down, 2n+1: spin-up f, 2n+2: spin-down f
U = 10000;

V = -10;           # Hybrization energy between conduction electron and impurity electron
Eu = -V*V;          # On-site energy for spin up electron at the impurity
Ed = U-V*V;          # On-site energy for spin down electron at the impurity
C = zeros(N,N);  # Correlation function matrix
C_im = zeros(2,2); # Correlation function for impurity electron
S_im = 0;          # EE for impurity 
MI = zeros(iter); # Mutual information 



# Set up the c-electron

# This step means the transition between the horizontally neighrboring sites.

for i = 1:n-1

        # For the site at boundary, there is no transtion. Thus, leave the loop.
        
        if  i % l == 0
            continue 
        end

        # Spin up conduction electrons
        H[i,i+1] = t;
        H[i+1,i] = t; # symmetric construction for Hamiltonian
        
        # Spin down conduction electrons 
        H[i+n,i+n+1] = t;
        H[i+n+1,i+n] = t;
end

# This step is to set up the transition between the vertically neighboring sites

for i = 1:n-l-1
        
        # Spin up conduction electron
        H[i,i+l] = t;
        H[i+l,i] = t;
        
        # Spin down conduction electron
        H[i+n,i+n+l] = t;
        H[i+n+l,i+n] = t;

end

# Set up the f-electron

# On-site energy for the impurity site

    H[2n+1,2n+1] = Eu;
    H[2n+2,2n+2] = Ed;

# Set up the interaction between c-electron and f-electron only at the center of the lattice

# find the central site 
c = Int((1+n)/2);

H[2n+1,c] = V;
H[c,2n+1] = V;
H[2n+2,c+n] = V;
H[c+n,2n+2] = V;

# Solve the eigenstate 
E,M = eig(H);



# Construct the correlation function for each site including c-electron site and impurity site 
# Here, i set up the half-filled case. Thus, the second loop is only summed to n (roughly equal to N/2)

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

# Calculate the impurity EE

C_im[1,1] = C[N-1,N-1];
C_im[1,2] = C[N-1,N];
C_im[2,1] = C_im[1,2];
C_im[2,2] = C[N,N];
S_im = eecal(C_im);


# Calculate the Mutual Information 


for i = 1:iter
    
    size = Int((2*i-1)*(2*i-1));
    C_red = zeros(N,N);    
    C_con = zeros(N,N);
    S_total = 0;
    S_con = 0;
    findcorT(C,C_red,C_con,l,n,i,size);
    S_con = eecal(C_con);
    S_total = eecal(C_red);    
    MI[i] = S_con + S_im - S_total;

end


plot(MI)
