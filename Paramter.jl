# This script is to generate a 2D single impurity Anderson model on a square lattice
# To simplify, the possible transition between c-electron and f-electron happens at the center of the lattice (only one site) and the impurity site with transition amplitude V.

include("eecal.jl");
include("findcor.jl");
include("Hamiltonian.jl");
include("Corr.jl");
n = 121;      # the number of conduction electon sites (WARNING: this number should be the square of ODD number.)
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
MI1 = zeros(iter); # Mutual information 

# This subroutine can help us construct the Hamiltonian
Hamiltonian(H,n,t,V,Eu,Ed);
# Solve the eigenstate 
E,M = eig(H);
# Construct the two-point correlation function
Corr(n,C,M);

# Calculate the impurity EE
C_im[1,1] = C[N-1,N-1];
C_im[1,2] = C[N-1,N];
C_im[2,1] = C_im[1,2];
C_im[2,2] = C[N,N];
S_im = eecal(C_im);

# Calculate the Mutual Information 

for i = 2:iter
    
    size = Int((2*i-3)*(2*i-3));
    Dim = Int(2*size+2);
    dim = Int(2*size+2);
    C_red = zeros(Dim,Dim);    
    C_con = zeros(dim,dim);
    S_total = 0;
    S_con = 0;
    findcor(C,C_red,C_con,l,n,i,size);
    S_con = eecal(C_con);
    S_total = eecal(C_red);    
    MI1[i] = S_con + S_im - S_total;

end

MI1[1] = 0;


