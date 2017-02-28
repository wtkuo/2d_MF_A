# C is the original correlation function
# C_red is the reduced correlation function which would be created in this subroutine
# D is the dimension of full hamiltonian (also the dimension of full correlation function)
# l is the length of the side in the original lattice
# N is the number of conduction electron site
# k is the index number for the range we are calculating


function findcor(C,C_red,l,N,k)

    size = Int((2k-1)*(2k-1));    # the number of sites in this calculation
    len = sqrt(size);               # the length of each side
    center = (N+1)/2;             # center coordinate 
    v = zeros(size+2);            # index for us to use 
    

    v[size+1] = N+1;
    v[size+2] = N+2;
    v[1] = c-(k-1)*l-(k-1);
    
    

    for j = 0:len-1
        for i = 1:len
            v[i+j*len] = v[1] + (i-1) + j*l;
        end
    end


    for i = 1:size+2
        for j = 1:size+2 
            C_red[i,j] = C[v[i],v[j]];
        end
    end
           
end
