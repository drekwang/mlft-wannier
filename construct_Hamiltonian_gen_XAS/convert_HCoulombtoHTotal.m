function [HTotal,binaryindices] = convert_HCoulombtoHTotal(...
    V, numorbitals,numelectrons)
    L2pindices = 11:46;
    TM2pindices = 47:52;
    
    binaryindices = gen_binaryindices(numorbitals,numelectrons);
    %Since all Coulomb matrix elements involving L 2p e- (indices 11 to 36
    %in one-particle basis) are zero completely omit calculations 
    %involving these elements
    occupied = ones(1,36);
    [~,idx] = ismember(binaryindices(:,L2pindices), occupied, 'rows');
    %Number of many body states
    numstates = size(binaryindices,1);
    HTotal = zeros(numstates,numstates);
    for bra = 1:numstates %binaryindices_occupiedL2p
        %Only calculate the diagonal and upper triangle to halve calc time
        for ket = bra:numstates %setdiff(binaryindices_occupiedL2p,1:bra-1)
            state_bra = binaryindices(bra,:);
            state_ket = binaryindices(ket,:);
            %Find indices of all occupied orbitals
            indices_bra = find(state_bra);
            indices_ket = find(state_ket);
            %Calculate the differential number of orbital occupancy
            differentorbitals = sum(mod(state_bra+state_ket,2))/2;
            if differentorbitals == 0
                %Sum across all occupied orbitals, except L 2p
                indices_bra_noL2p = setdiff(indices_bra,L2pindices);
                for alpha = indices_bra_noL2p
                    %Spin-orbitals are single occupancy, so there's no point
                    %in considering double annihilation of an orbital
                    indices_bra_noL2p_noalpha = setdiff...
                        (indices_bra_noL2p,alpha);
                    %There are no Coulomb repulsions for L 2p-L 2p, so no
                    %need to consider 'beta' as L 2p if 'alpha' already is
                    if alpha >= TM2pindices(1)
                        indices_bra_noL2p_noalpha_noTM2p = ...
                            setdiff(indices_bra_noL2p_noalpha,TM2pindices);
                    else
                        indices_bra_noL2p_noalpha_noTM2p = ...
                            indices_bra_noL2p_noalpha;
                    end
                    for beta = indices_bra_noL2p_noalpha_noTM2p
                        for gamma = [alpha beta]
                            for delta = setdiff([alpha beta], gamma)
                                %[alpha beta gamma delta]
                                state_bra_ann = state_bra;
                                state_ket_ann = state_ket;
                                numswitches = 0;
                                %Annhilate e- in 'alpha' & 'beta' orbitals of bra
                                state_bra_ann(alpha) = state_bra_ann(alpha) - 1;
                                %Count number of times orbital creation operators must be
                                %switched to move adjacent to annihilation operators of
                                %the Hamiltonian
                                numswitches = numswitches + sum(state_bra_ann(1,1:alpha));
                                state_bra_ann(beta) = state_bra_ann(beta) - 1;
                                numswitches = numswitches + ...
                                    sum(state_bra_ann(1,1:beta));
                                %Annhilate e- in 'gamma' and 'delta' orbitals of bra
                                state_ket_ann(delta) = state_ket_ann(delta) - 1;
                                numswitches = numswitches+sum(state_ket_ann(1,1:delta));
                                state_ket_ann(gamma) = state_ket_ann(gamma) - 1;
                                numswitches = numswitches+sum(state_ket_ann(1,1:gamma));
                                %'~any(mod(c+d,2))' is 0 if the remaining orbitals are not
                                %identically occupied and 1 if they are.
                                if ~any(mod(state_bra_ann+state_ket_ann,2))
                                    HTotal(bra,ket) = HTotal(bra,ket)+...
                                        (-1)^numswitches*1/2*...
                                        get_V(alpha,beta,gamma,delta,V);
                                end
                            end
                        end
                    end
                end
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