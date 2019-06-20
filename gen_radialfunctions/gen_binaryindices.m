function [binaryindices] = gen_binaryindices(numorbitals,numelectrons)
    numstates = factorial(numorbitals)/(factorial(numelectrons)*...
        factorial(numorbitals-numelectrons));
    binaryindices = zeros(round(numstates),numorbitals);
    state = zeros(1,numorbitals);
    for i = numorbitals - numelectrons + 1: numorbitals
        state(1,i) = 1;
    end
    binaryindices(1,:) = state;
    for i = 2:numstates
        state = gen_nextstate(state);
        binaryindices(i,:) = state;
    end
    %Change the ordering so that the lowest indices are occupied first
    %binaryindices = fliplr(binaryindices')';
end