function [w_final, K_final] = gen_XASspectra(HTotal,numinitial,nvector,qn,...
    gauntcoefficients,fullwidthhalfmax)

    numstates = size(HTotal,1);
    nx = nvector(1,1); %Polarization of incoming x-ray in x direction
    ny = nvector(1,2); %Polarization of incoming x-ray in y direction
    nz = nvector(1,3); %Polarization of incoming x-ray in z direction
    n = [0 0 complex(nx,ny)/sqrt(2) nz complex(-1*nx,ny)/sqrt(2) 0 0];
    binaryindices = gen_binaryindices(52,50);
    %full width half maximum for the Lorentzian broadening
    FWHM = fullwidthhalfmax; 
    maxE = 2000; %maximum energy in eV of calculated spectra
    c = gauntcoefficients;

    % Use the below script to 
    % zero off diagonal terms that are less than some value
    %Find eigenstates 'V' and eigenvalues 'D'
    [V,D] = eig(HTotal);
    %Zero all elements that are less than 0.01 for convenience
    V_simp = V;
    V_simp(abs(V) < 0.01) = 0;
    %Removing low values retains nearly all of the electron density
    %sum(V_diag_simp.^2,1);
    %Make V_diag_simp new eigenstate matrix
    V = V_simp;

    %Generate an array of eigenenergies for convenience
    E = zeros(1,numstates);
    for i = 1:numstates
        E(i) = D(i,i);
    end
    %Find index of eigenstate with energy below maxE
    index_maxstate = max(find(E<maxE));

    %Here goes the giant loop to calculate the spectrum
    w = []; %frequency omega
    K = []; %intensity
    count = 1;
    %For all initial states
    for u = 1:numinitial
        %For all final states
        for v = numinitial+1:index_maxstate
            summand = 0;
            %Quantum number m for d orbitals
            for md = -2:2
                %Quantum number m for p orbitals
                for mp = -1:1
                    %For spin down (0) and up (1)
                    for spin = 0:1
                        %Add 4 because n_tilde_neg3 occurs at index 1
                        %Should have used a dictionary but couldn't figure it
                        %out lol
                        n_mm = n(md-mp+4);
                        c_mm = get_c(2,md,1,mp,1,c);
                        %Get index of spin orbital
                        qn_index_2p = get_qn_index(2,1,mp,spin,1,qn);
                        qn_index_3d = get_qn_index(3,2,md,spin,1,qn);
                        %Find indices of many-particle basis state in
                        %eigenstates u and v
                        bases_indices_v = find(V(:,v))';
                        bases_indices_u = find(V(:,u))';
                        %For each many-particle basis state in u
                        for index_u = bases_indices_u
                            %Find occupancy of each spin orbital in basis_u
                            basis_u = binaryindices(index_u,:);
                            %If (2,1,mp,spin) orbital is occupied
                            if basis_u(1,qn_index_2p) == 1
                                %Annihilate e- in (2,1,mp,spin) orbital;
                                basis_u(1,qn_index_2p) = basis_u(1,qn_index_2p)-1;
                                %For each many-particle basis state in v
                                for index_v = bases_indices_v
                                    %Find occupancy of each spin orbital in
                                    %basis_v
                                    basis_v = binaryindices(index_v,:);
                                    %If (3,2,md,spin) orbital is occupied
                                    if basis_v(1,qn_index_3d) == 1
                                        %Annihilate e- from (3,2,md,spin)
                                        %orbital
                                        basis_v(1,qn_index_3d)=basis_v(1,qn_index_3d)-1;
                                        sgn = sum(basis_u(1,1:qn_index_2p))+...
                                            sum(basis_v(1,1:qn_index_3d));
                                        if basis_u == basis_v
                                            summand = summand+V(index_u,u)*...
                                                V(index_v,v)*n_mm*c_mm*(-1)^sgn
    %                                         %Constants to print for debugging
    %                                         n_mm
    %                                         c_mm
    %                                         md
    %                                         mp
    %                                         spin
    %                                         index_u
    %                                         index_v
    %                                         summand
    %                                         qn_index_2p
    %                                         qn_index_3d
    %                                         V(index_u,u)
    %                                         V(index_v,v)
    %                                         V(index_u,u)*...
    %                                             V(index_v,v)*n_mm*c_mm*(-1)^sgn
    %                                         index_u
    %                                         index_v
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            w(count)  = E(v) - E(u);
            K(count) = summand*conj(summand);%*exp(-1*(E(u)-E(v))/kT);
            count = count + 1;
        end
    end

    %Sort by increasing w
    XAS = [w; K];
    XAS_sorted = sort(XAS,2);
    w = XAS(1,:);
    K = XAS(2,:);
    %w has duplicate values, so combine corresponding K into single value
    count = 1;
    w_merge(1,count) = w(1,1);
    K_merge(1,count) = K(1,1);
    for i = 2:length(w)
        if w(1,i) == w_merge(1,count)
            K_merge(1,count) = K_merge(1,count)+K(1,i);
        else
            count = count+1;
            K_merge(1,count) = K(1,i);
            w_merge(1,count) = w(1,i);
        end
    end

    %Apply Lorentzian broadening to each peak and sum individual peaks together
    %to get the total spectrum
    w_final = linspace(0,max(w_merge+15),10000);
    K_final = zeros(1,length(w_final));
    for peak = 1:length(w_merge)
        K_temp = K_merge(1,peak)*(0.5*FWHM)./((w_final-w_merge(1,peak)).^2+...
            (0.5*FWHM)^2);
        K_final = K_final + K_temp;
    end
end