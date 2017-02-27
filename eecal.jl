
function eecal(C)

    d = size(C,1);
    eva_C = zeros(d);
    epi = zeros(d);
    entropy = 0.0;
    eva_C = eigvals(C);
    
    for i = 1:d
            if eva_C[i] > 0.999 || eva_C[i]  < 0.001
                eva_C[i] = 0.000000000000000000000001
            end
        
        epi[i] = log(abs(eva_C[i])/(1 - abs(eva_C[i])));

        entropy += log( 1 + e^(-epi[i])) + (epi[i])*e^(-epi[i])/(1+e^(-epi[i]));
    end

    
return entropy
end
