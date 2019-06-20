function [rho] = gen_rho(H_TB, numelectrons)
    %Generates density matrix 'rho' by determining occupation of each spin
    %orbital
    [V,~] = eig(H_TB);
    numorbitals = size(H_TB,1);
    V_2 = V.^2;
    rho_array = sum(V_2(:,1:numelectrons),2);
    rho = zeros(numorbitals, numorbitals);
    for i = 1:length(rho_array)
        rho(i,i) = rho_array(i);
    end    
end
    