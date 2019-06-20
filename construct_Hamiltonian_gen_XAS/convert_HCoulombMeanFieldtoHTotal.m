function [HTotal_CoulombMeanField] = convert_HCoulombMeanFieldtoHTotal(...
    numorbitals, numelectrons, V, rho)
    %Use Hartree mean-field approximation on HCoulomb and convert into 
    %HTotal in many-body basis
    L2pindices = 11:46;
    TM2pindices = 47:52;
    TM3dindices = 1:10;
    %Generate indexed list of many-body basis states
    binaryindices = gen_binaryindices(numorbitals,numelectrons);
    numstates = size(binaryindices,1);
    
    %Since all Coulomb matrix elements involving L 2p e- (indices 11 to 36
    %in one-particle basis) are zero completely omit calculations 
    %involving these elements
    occupied = ones(1,36);
    [~,idx] = ismember(binaryindices(:,L2pindices), occupied, 'rows');
    
    %Number of many body states
    numstates = size(binaryindices,1);
    
    %Initialize HTotal, which will later be set to HTotal_CoulombMeanField
    HTotal = zeros(numstates,numstates);
    
    for bra = 1:numstates
        bra
        %Only calculate the diagonal and upper triangle to halve calc time
        for ket = bra:numstates
            state_bra = binaryindices(bra,:);
            state_ket = binaryindices(ket,:);
            %Find indices of all occupied orbitals
            indices_bra = find(state_bra);
            indices_ket = find(state_ket);
            %Calculate the differential number of orbital occupancy
            differentorbitals = sum(mod(state_bra+state_ket,2))/2;
            %V_ijkl is non-zero only identical orbitals
            if differentorbitals == 0
                %Sum across only TM 3d orbitals
                %Original DFT calculation includes mean-field interactions of 
                %TM 3d, L 2p electrons, and TM 2p
                %MLFT corrected for TM 3d and TM 2p
                %Thus, the only double counted interactions are the TM 3d
                %and TM 2p orbitals
                %Still unsure whether to subtract out F0, F2, and F4, or
                %just F2 and F4. Technically it should be all if we had
                %calculated all F(k) parameters from first-principles
                %but we fit F0 to experiment, so it's possible that doesn't
                %need to be subtracted out? I think it does, but Haverkort
                %says it doesn't
                indices_bra_TM3d2p = setdiff(indices_bra,...
                    L2pindices);
                for alpha = indices_bra_TM3d2p
                    %Spin-orbitals are single occupancy, so there's no point
                    %in considering double annihilation of an orbital
                    indices_bra_TM3d2p_noalpha = setdiff...
                        (indices_bra_TM3d2p,alpha);
                    for beta = indices_bra_TM3d2p_noalpha
                        for gamma = [alpha beta]
                            for delta = setdiff([alpha beta], gamma)
                                %[alpha beta gamma delta]
                                % Two cases: 
                                % Case 1: alpha=gamma and beta=delta
                                if alpha == gamma
                                    %-a_i_dag*a_k*<a_j_dag*a_l>
                                    %Number of switches n_ik starts at 1 
                                    %because sign in front of operator pair 
                                    %is -1
                                    n_ik = 1;
                                    %Initialize a bra state to annihilate
                                    state_bra_ann = state_bra;
                                    %Annhilate e- in 'alpha' orbital of bra
                                    state_bra_ann(alpha) = state_bra_ann(alpha) - 1;
                                    %Count up the number of switches
                                    %necessary to move the alpha creation
                                    %operator next to the annihilation
                                    n_ik = n_ik + sum(state_bra_ann(1,1:alpha));
                                    %Initialize a ket state to annihilate
                                    state_ket_ann = state_ket;
                                    %Annhilate e- in 'gamma' orbital of bra
                                    state_ket_ann(gamma) = state_ket_ann(gamma) - 1;
                                    %Count up the number of switches
                                    %necessary to move the gamma creation
                                    %operator next to the annihilation
                                    n_ik = n_ik + sum(state_ket_ann(1,1:gamma));
                                   
                                    %-a_j_dag*a_l*<a_i_dag*a_k>
                                    %Number of switches n_jl starts at 1 
                                    %because sign in front of operator pair 
                                    %is -1
                                    n_jl = 1;
                                    %Initialize a bra state to annihilate
                                    state_bra_ann = state_bra;
                                    %Annhilate e- in 'beta' orbital of bra
                                    state_bra_ann(beta) = state_bra_ann(beta) - 1;
                                    %Count up the number of switches
                                    %necessary to move the beta creation
                                    %operator next to the annihilation
                                    n_jl = n_jl + sum(state_bra_ann(1,1:beta));
                                    %Initialize a ket state to annihilate
                                    state_ket_ann = state_ket;
                                    %Annhilate e- in 'delta' orbital of bra
                                    state_ket_ann(delta) = state_ket_ann(delta) - 1;
                                    %Count up the number of switches
                                    %necessary to move the delta creation
                                    %operator next to the annihilation
                                    n_jl = n_jl + sum(state_ket_ann(1,1:delta));
                                    
                                    %-a_i_dag*a_k*<a_j_dag*a_l>-
                                    %a_j_dag*a_l*<a_i_dag*a_k>+
                                    %<a_j_dag*a_l><a_i_dag*a_k>
                                    HTotal(bra,ket) = HTotal(bra,ket)+...
                                    1/2*get_V(alpha,beta,gamma,delta,V)*...
                                    ((-1)^n_ik*rho(beta,delta)+(-1)^n_jl*rho(alpha,gamma)...
                                    +rho(beta,delta)*rho(alpha,gamma));
                                % Case 2: alpha=delta and beta=gamma
                                else
                                    %a_i_dag*a_l*<a_j_dag*a_k>
                                    %Number of switches n_il starts at 0 
                                    %because sign in front of operator pair 
                                    %is +1
                                    n_il = 0;
                                    %Initialize a bra state to annihilate
                                    state_bra_ann = state_bra;
                                    %Annhilate e- in 'alpha' orbital of bra
                                    state_bra_ann(alpha) = state_bra_ann(alpha) - 1;
                                    %Count up the number of switches
                                    %necessary to move the alpha creation
                                    %operator next to the annihilation
                                    n_il = n_il + sum(state_bra_ann(1,1:alpha));
                                    %Initialize a ket state to annihilate
                                    state_ket_ann = state_ket;
                                    %Annhilate e- in 'delta' orbital of bra
                                    state_ket_ann(delta) = state_ket_ann(delta) - 1;
                                    %Count up the number of switches
                                    %necessary to move the delta creation
                                    %operator next to the annihilation
                                    n_il = n_il + sum(state_ket_ann(1,1:delta));
                                   
                                    %+a_j_dag*a_k*<a_i_dag*a_l>
                                    %Number of switches n_jk starts at 1 
                                    %because sign in front of operator pair 
                                    %is -1
                                    n_jk = 0;
                                    %Initialize a bra state to annihilate
                                    state_bra_ann = state_bra;
                                    %Annhilate e- in 'beta' orbital of bra
                                    state_bra_ann(beta) = state_bra_ann(beta) - 1;
                                    %Count up the number of switches
                                    %necessary to move the beta creation
                                    %operator next to the annihilation
                                    n_jk = n_jk + sum(state_bra_ann(1,1:beta));
                                    %Initialize a ket state to annihilate
                                    state_ket_ann = state_ket;
                                    %Annhilate e- in 'gamma' orbital of bra
                                    state_ket_ann(gamma) = state_ket_ann(gamma) - 1;
                                    %Count up the number of switches
                                    %necessary to move the gamma creation
                                    %operator next to the annihilation
                                    n_jk = n_jk + sum(state_ket_ann(1,1:gamma));
                                    
                                    %a_i_dag*a_l*<a_j_dag*a_k>+
                                    %a_j_dag*a_k+1*<a_i_dag*a_l>-
                                    %<a_j_dag*a_k>*<a_i_dag*a_l>
                                    HTotal(bra,ket) = HTotal(bra,ket)+...
                                    1/2*get_V(alpha,beta,gamma,delta,V)*...
                                    ((-1)^n_il*rho(beta,gamma)+(-1)^n_jk*rho(alpha,delta)-...
                                    rho(beta,gamma)*rho(alpha,delta));
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