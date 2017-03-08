function findcor(C,C_red,C_con,n,k,size)
# This subroutine is to construct the reduced matrix for conduction electron only and the total reduced matrix
# for conduction electrons and impurity electron. This counting method is from the center to the outside. Basically, 
# the mutual information will keep increasing. However, the increament will become less and less as we go to 
# the outer region. 
# 
# Parameter :
# C is the original correlation function
# C_red is the correlation function for reduced region
# C_con is the correlation function for conduction electrons
# n is the number of conduction electron site
# k is the index number for the range we are calculating
# size is the number of sites in the kth iteration



    len_side = Int(sqrt(n)); # the size of the side in the original space 
    len = Int(sqrt(size));   # the size of the side in kth iteration space
    center = (n+1)/2;        # center coordinate 
    N = 2*size+2;            # the dimension for correlation function for reduced region
    v = zeros(N);            # this v vector stores the index of the sites which we need for the kth iteration
                            
    v[N-1] = 2*n+1;             # This is spin up impurity site 
    v[N] = 2*n+2;               # This is spin down impurity site 
    
    v[1] = center-(k-1)*len_side-(k-1);     # This is the left-bottom site for this kth iteration


# I set up for spin up condution sites first

    for j = 0:len-1
        for i = 1:len
            v[i+j*len] = v[1] + (i-1) + j*len_side;
        end
    end

# This is for spin down conduction sites

    for i = 1:size
        v[i+size] = v[i] + n;
    end

# Happily construct the reduced density matrix


    for i = 1:N
        for j = 1:N

            if i < N-1 && j < N-1
                C_con[i,j] = C[Int(v[i]),Int(v[j])];
            end
        
            C_red[i,j] = C[Int(v[i]),Int(v[j])];
        
        end
    end    
end

