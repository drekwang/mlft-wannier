function [H_SO] = gen_HSO(numorbitals,qn, lambda_3d, lambda_2p)
    indices_TM3d = 1:10;
    indices_L2p = 11:46;
    indices_TM2p = 47:52;
    coupling_3d = lambda_3d;%eV
    coupling_2p = lambda_2p;%eV
    n_3d = 3;
    l_3d = 2;
    n_2p = 2;
    l_2p = 1;
    
    H_SO = zeros(numorbitals,numorbitals);
    %Fill in matrix elements for TM 3d orbitals
    for bra = indices_TM3d
        for ket = indices_TM3d
            summand = 0;
            %Equation 13 in Eder
            for m = -1*l_3d:l_3d
                %'c1' annihilates bra with n = 3, l = 2, m = m, s = +
                c1 = (bra==find(ismember(qn(:,1:5),[n_3d l_3d m 1 1], 'rows'),1));
                %'c2' annihilates ket with n = , l = 2, m = m, s = + 
                c2 = (ket==find(ismember(qn(:,1:5),[n_3d l_3d m 1 1], 'rows'),1));
                %'c3' annihilates bra with n = 3, l = 2, m = m, s = - 
                c3 = (bra==find(ismember(qn(:,1:5),[n_3d l_3d m 0 1], 'rows'),1));
                %'c4' annihilates ket with n = 3, l = 2, m = m, s = -
                c4 = (ket==find(ismember(qn(:,1:5),[n_3d l_3d m 0 1], 'rows'),1));
                summand = summand + m*(c1*c2-c3*c4);
            end
            %Equation 14 in Eder
            for m = -1*l_3d:l_3d-1
                %'c5' annihilates bra with n = 3, l = 2, m = m+1, s = -
                c5 = (bra==find(ismember(qn(:,1:5),[n_3d l_3d m+1 0 1], 'rows'),1));
                %'c6' annihilates ket with n = 3, l = 2, m = m, s = +
                c6 = (ket==find(ismember(qn(:,1:5),[n_3d l_3d m 1 1], 'rows'),1));
                %'c7' annihilates bra with n = 3, l = 2, m = m, s = +
                c7 = (bra==find(ismember(qn(:,1:5),[n_3d l_3d m 1 1], 'rows'),1));
                %'c8' annihilates ket with n = 3, l = 2, m = m+1, s = -
                c8 = (ket==find(ismember(qn(:,1:5),[n_3d l_3d m+1 0 1], 'rows'),1));
                summand = summand + sqrt((l_3d-m)*(l_3d+m+1))*(c5*c6+c7*c8);
            end
            H_SO(bra,ket) = H_SO(bra,ket)+coupling_3d/2*summand;
        end     
    end
    
    %Fill in matrix elements for TM 2p orbitals
    for bra = indices_TM2p
        for ket = indices_TM2p
            summand = 0;
            for m = -1*l_2p:l_2p
                %'c1' annihilates bra with n = 2, l = 1, m = m, s = +
                c1 = (bra==find(ismember(qn(:,1:5),[n_2p l_2p m 1 1], 'rows'),1));
                %'c2' annihilates ket with n = 2, l = 1, m = m, s = + 
                c2 = (ket==find(ismember(qn(:,1:5),[n_2p l_2p m 1 1], 'rows'),1));
                %'c3' annihilates bra with n = 2, l = 1, m = m, s = - 
                c3 = (bra==find(ismember(qn(:,1:5),[n_2p l_2p m 0 1], 'rows'),1));
                %'c4' annihilates ket with n = 1, l = 1, m = m, s = -
                c4 = (ket==find(ismember(qn(:,1:5),[n_2p l_2p m 0 1], 'rows'),1));
                summand = summand + m*(c1*c2-c3*c4);
            end
            for m = -1*l_2p:l_2p-1
                %'c5' annihilates bra with n = 3, l = 2, m = m+1, s = -
                c5 = (bra==find(ismember(qn(:,1:5),[n_2p l_2p m+1 0 1], 'rows'),1));
                %'c6' annihilates ket with n = 3, l = 2, m = m, s = +
                c6 = (ket==find(ismember(qn(:,1:5),[n_2p l_2p m 1 1], 'rows'),1));
                %'c7' annihilates bra with n = 3, l = 2, m = m, s = +
                c7 = (bra==find(ismember(qn(:,1:5),[n_2p l_2p m 1 1], 'rows'),1));
                %'c8' annihilates ket with n = 3, l = 2, m = m+1, s = -
                c8 = (ket==find(ismember(qn(:,1:5),[n_2p l_2p m+1 0 1], 'rows'),1));
                summand = summand + sqrt((l_2p-m)*(l_2p+m+1))*(c5*c6+c7*c8);
            end
            H_SO(bra,ket) = H_SO(bra,ket)+coupling_2p/2*summand;
        end     
    end
end