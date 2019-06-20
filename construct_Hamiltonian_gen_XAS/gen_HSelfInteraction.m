function HSI = gen_HSelfInteraction(numspinorbitals, numelectrons, V, rho)
    %quanty.org/documentation/language_reference/functions/meanfieldoperator
    HSI = zeros(numspinorbitals,numspinorbitals);
    %We double count self-interaction only for d orbitals
    for i = 1:10
        %summand is U_ijji-U_ijij
        summand = 0;
        for j = 1:10
            if i ~= j
                summand = get_V(i,j,j,i,V) - get_V(i,j,i,j,V);
            end
        end
        HSI(i,i) = 2/(numelectrons*(numelectrons-1))*summand;
    end
end