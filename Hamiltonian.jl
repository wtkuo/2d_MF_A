function Hamiltonian(H,n,t,V,Eu,Ed)
# This subroutine is to setup the 2D single impurity Mean Field Approximation Anderson model. Lattice type is 
# square lattice. The single impurity is at the center. The interaction between conduction electron sites and 
# impurity only shows up at the center.
# 
# Parameters: 
# H : the Hamiltonian
# n : the number of conduction sites
# t : the hopping amplitude for conduction electrons between the neighboring sites
# V : the hopping amplitude between the central conduction site and the impurity
# Eu: the on-site energy at impurity for spin-up f electron 
# Ed: the on-site energy at impurity for spin-down f electron 
#

# Define the size of the side first 

len_side = Int(sqrt(n));

# Set up the c-electron
# This step is to set up the transition between the horizontally neighrboring sites.

for i = 1:n-1

        # For the site at boundary, there is no transtion. Thus, leave the loop.
        
        if  i % len_side == 0
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

for i = 1:n-len_side-1
        
        # Spin up conduction electron
        H[i,i+len_side] = t;
        H[i+len_side,i] = t;
        
        # Spin down conduction electron
        H[i+n,i+n+len_side] = t;
        H[i+n+len_side,i+n] = t;

end

# Set up the f-electron

# On-site energy for the impurity site

    H[2n+1,2n+1] = Eu;
    H[2n+2,2n+2] = Ed;

# Set up the interaction between c-electron and f-electron only at the center of the lattice

# find the center site 
center = Int((1+n)/2);

H[2n+1,center] = V;
H[center,2n+1] = V;
H[2n+2,center+n] = V;
H[center+n,2n+2] = V;

end
