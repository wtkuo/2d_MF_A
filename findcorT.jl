# C is the original correlation function
# C_red is the correlation function for reduced region
# C_con is the correlation function for conduction electrons
# l is the length of the side in the original lattice
# n is the number of conduction electron site
# k is the index number for the range we are calculating
# size is the abandoned number of sites for conduction electrons

function findcorT(C,C_red,C_con,l,n,k,size)

    len = Int(sqrt(size));  # the length of each side in kth iteration
    center = (n+1)/2;       # center coordinate 
    N = 2*size;             # the total elimanated DOF 
    v = zeros(N);        # this v vector stores the index of the sites which we need for the kth iteration
    RN = 2*n+2;        # the dimension for correlation function for reduced region                           
    v[1] = c-(k-1)*l-(k-1);     # This is the left-bottom site for this kth iteration

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

    
# The situation is a little bit different here. v[i] lists the points which I don't want to include 
# Happily construct the reduced density matrix

    for i = 1:RN
        for j = 1:RN
    
            if i in v || j in v
                C_red[i,j] = 0;
                continue
            end

            C_red[i,j] = C[i,j];
            
            if i < RN-1 && j < RN-1
                C_con[i,j] = C_red[i,j];
            end

        end
    end
    

    
    
    
end

