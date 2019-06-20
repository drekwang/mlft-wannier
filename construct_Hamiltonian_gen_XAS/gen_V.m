function [V] = gen_V(numorbitals, gauntcoefficients, qn, Rk)
%Generate a list of values of V(alpha,beta,gamma,delta) (where they are
%listed in the first four columns and the value is in in the fifth column)
    %numorbitals = 52; 
    %5 TM 3d + 18 L 2p + 3 TM 3p = 26 or 52 spin orbitals
    %'n' is number of ijkl permutations, excluding i~=k or l and j~=k or l
    c = gauntcoefficients;
    n = numorbitals*(numorbitals-1)*2;
    %Column 1 is i, 2 is j, 3 is k, 4 is l, 5 is the value of V_ijkl
    V = zeros(n,5);
    count = 1;
    for i = 1:numorbitals
        for j = 1:numorbitals
            if i~=j
                %Exchange (G)
                V(count,1) = i;
                V(count,2) = j;
                V(count,3) = i;
                V(count,4) = j;
                count = count + 1;
                %Coulombic (F)
                V(count,1) = i;
                V(count,2) = j;
                V(count,3) = j;
                V(count,4) = i;
                count = count + 1;
            end
        end
    end

    for ijkl = 1:n
        i = V(ijkl,1);
        j = V(ijkl,2);
        k = V(ijkl,3);
        l = V(ijkl,4);
        %In 'qn', column 1 is n, 2 is l, 3 is m, 4 is s, 
        %5 is the atom (1 for metal, 0 for ligand)
        %Calc delta(si,sl)*delta(sj,sk)*delta(mi+mj,ml+mk)
        deltafactor = (qn(i,4)==qn(l,4))*(qn(j,4)==qn(k,4))*...
            (qn(i,3)+qn(j,3)==qn(l,3)+qn(k,3));
        %Calc c^(k)(li,mi,ll,ml)*c^(k)(lk,mk,lj,mj)*
        %R^(k)(ni,li,nj,lj,nk,lk,nl,ll)
        %Rk is constructed such that only TM 2p-3d and TM 3d-3d
        %interactions are considered, not TM 2p-2p, TM 2p-L 2p,
        %TM 3d-L 2p, or L 2p-L 2p.
        sum = 0;
        for order = 0:4
            sum = sum+...
                get_c(qn(i,2),qn(i,3),qn(l,2),qn(l,3),order,c)...
                *get_c(qn(k,2),qn(k,3),qn(j,2),qn(j,3),order,c)...
                *get_R(qn(i,1),qn(i,2),qn(j,1),qn(j,2),...
                qn(k,1),qn(k,2),qn(l,1),qn(l,2),qn(k,5),qn(l,5),order,Rk);
        end
        V_ijkl = deltafactor*sum;
        V(ijkl,5) = V_ijkl;
end