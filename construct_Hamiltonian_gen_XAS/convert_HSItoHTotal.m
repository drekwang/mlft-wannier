function [HTotal,binaryindices] = convert_HSItoHTotal(...
    HSI,numorbitals,numelectrons)
    
    TM3dindices = 1:10;
    L2pindices = 11:46;
    TM2pindices = 47:52;
    binaryindices = gen_binaryindices(numorbitals,numelectrons);
    %Number of many body states
    numstates = size(binaryindices,1);
    HTotal = zeros(numstates,numstates);
    for bra = 1:numstates
        %Only calculate the diagonal and upper triangle to halve calc time
        for ket = bra:numstates
            a = binaryindices(bra,:);
            b = binaryindices(ket,:);
            %Find indices of all occupied orbitals
            differentorbitals = sum(mod(a+b,2))/2;
            %If bra and ket are identical, 
            if differentorbitals == 0
                %Sum across all occupied orbitals
                indices_a = find(a);
                for alpha = indices_a
                    beta = alpha;
                    c = a;
                    d = b;
                    %Annhilate e- in 'alpha' & 'beta' orbitals
                    c(alpha) = c(alpha) - 1;
                    d(beta) = d(beta) - 1;
                    %Count number of times orbital creation operators must be
                    %switched to move adjacent to annihilation operators of
                    %the Hamiltonian
                    numswitches = sum(c(1,1:alpha)) + sum(d(1,1:beta));
                    %'~any(mod(c+d,2))' is 0 if the remaining orbitals are not
                    %identically occupied and 1 if they are.
                    HTotal(bra,ket) = HTotal(bra,ket) + ~any(mod(c+d,2)) * ...
                        (-1)^numswitches*HSI(alpha,beta);
                end
            elseif differentorbitals == 1
                %Sum across all occupied orbitals
                indices_a = find(a);
                indices_b = find(b);
                alpha = setdiff(indices_a,indices_b);
                beta = setdiff(indices_b,indices_a);
                c = a;
                d = b;
                %Annhilate e- in 'alpha' & 'beta' orbitals
                c(alpha) = c(alpha) - 1;
                d(beta) = d(beta) - 1;
                %Count number of times orbital creation operators must be
                %switched to move adjacent to annihilation operators of
                %the Hamiltonian
                numswitches = sum(c(1,1:alpha)) + sum(d(1,1:beta));
                %'~any(mod(c+d,2))' is 0 if the remaining orbitals are not
                %identically occupied and 1 if they are.
                HTotal(bra,ket) = HTotal(bra,ket) + ~any(mod(c+d,2)) * ...
                    (-1)^numswitches*HSI(alpha,beta);
            end
        end
    end
    %Delete elements on the diagonal
    HTotal_copy = HTotal;
    for state = 1:numstates
        HTotal_copy(state,state) = 0;
    end
    %Mirror 'HTotal' over the diagonal
    HTotal = HTotal_copy + rot90(fliplr(HTotal),1);
end