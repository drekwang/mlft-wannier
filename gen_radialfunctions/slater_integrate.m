function [X_k] = slater_integrate(r1_old, R1_old, r2_old, R2_old, co_or_ex, k)
    %written on 4/17/19 by Derek Wang
    %input: 'R1'/'R2' are radial functions (Nx2), 'l1'/'l2' are angular 
    %momentum quantum number (0 for s, 1 for l, 2 for d, etc.), 'co_or_ex' is 1
    %for Coulomb integral (F) and 0 for exchange (G), 'k' is k in Haverkort et
    %al. Phys. Rev. B 86, 2012. 
    %output: F_{l1,l2}^{(k)} for 'co_or_ex' = 1 (Coulomb coefficient)
    %or G_{l1,l2}^{(k)} for 'co_or_ex' = 0 (exchange coefficient)

    %Constants
    e = 1.6*10^-19; %elementary charge in C
    eps_0 = 8.854187*10^-12; %vacuum permittivity in F/m = s^2*C^2*m^-2*kg^-1
    %r1 and r2 likely don't have the same spacing or max value, so re-cast onto
    %uniform grid with 'griddata'
    r1 = linspace(max(min(r1_old),min(r2_old)),min(max(r1_old),max(r2_old)),100)';
    r2 = linspace(max(min(r1_old),min(r2_old)),min(max(r1_old),max(r2_old)),100)';
    R1 = interp1(r1_old,R1_old,r1,'pchip');
    R2 = interp1(r2_old,R2_old,r2,'pchip');
    %Create a 2D meshgrid using r1 and r2
    [r1_mesh, r2_mesh] = meshgrid(r1,r2);
    %List out all the coordinates
    r_points = [r1_mesh(:) r2_mesh(:)];
    values = zeros(size(r_points,1),1);
    if co_or_ex == 1
        for i = 1:size(r_points,1)
            %Equation 1 in Haverkort et al, Phys. Rev. B 85, 2012.
            values(i,1) = min(r_points(i,1),r_points(i,2))^k/...
                (max(r_points(i,1),r_points(i,2)))^(k+1)*r_points(i,1)^2*...
                r_points(i,2)^2*R1(find(r1==r_points(i,1)))^2*...
                R2(find(r2==r_points(i,2)))^2;
        end
        %Only for (r1,r2)=(0,0) should MATLAB output NaN. Physically, we know it
        %should be 0.
        values(isnan(values))=0;
        values_mesh = zeros(size(r1,1),size(r2,1));
        count = 1;
        %Because 'r_points' holds 'r1' constant and increments 'r2', the values 
        %correspond to going down columns in 'value_mesh', which is why the index
        %ordering is weird
        for j = 1:size(r1,1)
            for i = 1:size(r2,1)
                values_mesh(i,j) = values(count,1);
                count = count+1;
            end
        end
    else
        for i = 1:size(r_points,1)
            %Equation 1 in Haverkort et al, Phys. Rev. B 85, 2012.
            values(i,1) = min(r_points(i,1),r_points(i,2))^k/...
                (max(r_points(i,1),r_points(i,2)))^(k+1)*r_points(i,1)^2*...
                r_points(i,2)^2*R1(find(r1==r_points(i,1)))*...
                R2(find(r2==r_points(i,2)))*R1(find(r2==r_points(i,2)))*...
                R2(find(r1==r_points(i,1)));
        end
        %Only for (r1,r2)=(0,0) should MATLAB output NaN. Physically, we know it
        %should be 0.
        values(isnan(values))=0;
        values_mesh = zeros(size(r1,1),size(r2,1));
        count = 1;
        %Because 'r_points' holds 'r1' constant and increments 'r2', the values 
        %correspond to going down columns in 'value_mesh', which is why the index
        %ordering is weird
        for j = 1:size(r1,1)
            for i = 1:size(r2,1)
                values_mesh(i,j) = values(count,1);
                count = count+1;
            end
        end
    end
    %F_k in atomic units (Angstrom^-1)
    X_k_au = trapz(r2,trapz(r1,values_mesh,2));
    %Convert from atomic units to eV
    X_k = X_k_au*e^2/(4*pi()*eps_0)/e/10^-10;
end


            
