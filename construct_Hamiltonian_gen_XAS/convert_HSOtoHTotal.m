function [HTotal_SO,binaryindices] = convert_HSOtoHTotal(...
    H_SO,numorbitals,numelectrons)
    % numorbitals = 5;
    % numelectrons = 3;
    TM3dindices = 1:10;
    L2pindices = 11:46;
    TM2pindices = 47:52;
    binaryindices = gen_binaryindices(numorbitals,numelectrons);
    %Number of many body states
    numstates = size(binaryindices,1);
    HTotal_SO = zeros(numstates,numstates);
    for bra = 1:numstates
        %Only calculate the diagonal and upper triangle to halve calc time
        for ket = bra:numstates
            state_bra = binaryindices(bra,:);
            state_ket = binaryindices(ket,:);
            differentorbitals = sum(mod(state_bra+state_ket,2))/2;
            if differentorbitals <= 1
                %Find indices of all occupied orbitals except L 2p orbitals
                %because spin-orbit coupling is omitted for these terms
                indices_bra = setdiff(find(state_bra),L2pindices);
                indices_bra;
                %Sum across all occupied orbitals
                for alpha = indices_bra
                    indices_ket = setdiff(find(state_ket),L2pindices);
                    %If 'alpha' annihilates a 3d / 2p orbital, 'beta' must also
                    %annihilate a 3d / 2p orbital, respectively because
                    %there are not spin-orbit coupling terms between a 3d
                    %electron and a 2p electron
                    if alpha <= max(TM3dindices)
                        indices_ket = setdiff(indices_ket, TM2pindices);
                    else
                        indices_ket = setdiff(indices_ket,TM3dindices);
                    end
                    for beta = indices_ket
                        %[alpha beta]
                        state_bra_ann = state_bra;
                        state_ket_ann = state_ket;
                        %Annhilate e- in 'alpha' & 'beta' orbitals
                        state_bra_ann(alpha) = state_bra_ann(alpha) - 1;
                        state_ket_ann(beta) = state_ket_ann(beta) - 1;
                        %Count number of times orbital creation operators must be
                        %switched to move adjacent to annihilation operators of
                        %the Hamiltonian
                        numswitches = sum(state_bra_ann(1,1:alpha)) + ...
                            sum(state_ket_ann(1,1:beta));
                        %'~any(mod(c+d,2))' is 0 if the remaining orbitals are not
                        %identically occupied and 1 if they are.
                        if ~any(mod(state_bra_ann+state_ket_ann,2))
                            HTotal_SO(bra,ket) = HTotal_SO(bra,ket) + ...
                             (-1)^numswitches*H_SO(alpha,beta);
                        end
                    end
                end
            end
        end
    end
    %Delete elements on the diagonal
    HTotal_copy = HTotal_SO;
    for state = 1:numstates
        HTotal_copy(state,state) = 0;
    end
    %Mirror 'HTotal' over the diagonal
    HTotal_SO = HTotal_copy + rot90(fliplr(HTotal_SO),1);
end