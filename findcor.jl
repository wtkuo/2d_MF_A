# C is the original correlation function
# C_red is the reduced correlation function which would be created in this subroutine
# l is the length of the side in the original lattice
# n is the number of conduction electron site
# k is the index number for the range we are calculating


function findcor(C,C_red,C_con,l,n,k,size)

    len = sqrt(size);               # the length of each side in kth iteration
    center = (n+1)/2;       # center coordinate 
    v = zeros(2*size+2);    # this v vector stores the index of the sites which we need for the kth iteration
                            
    v[size+1] = 2*n+1;      # This is spin up impurity site 
    v[size+2] = 2*n+2;      # This is spin down impurity site 
    
    v[1] = c-(k-1)*l-(k-1); # This is the left-bottom site for this kth iteration
   
# I set up for spin up condution sites first

    for j = 0:len-1
        for i = 1:len
            v[i+j*len] = v[1] + (i-1) + j*l;
        end
    end

# This is for spin down conduction sites

    for i = 1:size
        v[i+size] = v[i] + n;
    end

# Happily construct the reduced density matrix

    for i = 1:2*size+2
        for j = 1:2*size+2 

            if i < 2*size+1 && j < 2*size+1
                C_con[i,j] = C[v[i],v[j]];
            end
            
            C_red[i,j] = C[v[i],v[j]];
        end
    end
           
end
