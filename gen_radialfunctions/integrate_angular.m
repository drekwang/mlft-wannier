function [r_plot, R] = integrate_angular(sphcoords_ordered,...
    values_ordered)

phi = sphcoords_ordered(:,1);
theta = sphcoords_ordered(:,2);
r = sphcoords_ordered(:,3);
original = [phi theta r values_ordered];

%'phi' ranges from 0 to 2pi and 'theta' ranges from 0 to pi
%If a point lies within the cutoff distance from either theta/phi limit,
%it will be mirrored to avoid NaN number errors during griddata 
%interpolation
phi_cutoff = pi();%max is pi()
theta_cutoff = pi()/2;%max is pi()/2
duplicate = zeros(size(sphcoords_ordered,1),4);
count = 1;
for i = 1:size(sphcoords_ordered,1)
    if phi(i)<phi_cutoff || 2*pi()-phi(i)<phi_cutoff || ...
            theta(i)<theta_cutoff || pi()-theta(i)<theta_cutoff 
        duplicate(count,:) = [phi(i) theta(i) r(i) values_ordered(i)];
        count = count+1;
    end
end
%Remove excess rows of zeroes at the end
duplicate = duplicate(any(duplicate,2),:);
%Counting shows that there will be 14 versions of 'duplicate'
mirrored = zeros(size(sphcoords_ordered,1)+14*size(duplicate,1),4);
count = 1;
for i = -1:1
    for j = -1:1
        %No need to duplicate into the original main cell
        if i~= 0 || j~=0
            %Rotate 'phi' by -2*pi, 2*pi, or 0
            new_phi = duplicate(:,1)+i*2*pi();
            %If 'theta' is not to be rotated
            if j==0
                translated = [new_phi duplicate(:,2) duplicate(:,3) duplicate(:,4)];
                mirrored((count-1)*size(duplicate,1)+1:count*size(duplicate,1),:)=translated;
                count = count+1;
            %If 'theta' is to be rotated
            else
                %Rotate 'phi' by +-pi() and adjust 'theta' to place the
                %point in the same spot in the unit sphere
                new_theta = (j+1)*pi()-duplicate(:,2);
                new_phi1 = new_phi + pi();
                new_phi2 = new_phi - pi();
                translated1 = [new_phi1 new_theta duplicate(:,3) duplicate(:,4)];
                translated2 = [new_phi2 new_theta duplicate(:,3) duplicate(:,4)];
                mirrored((count-1)*size(duplicate,1)+1:(count)*size(duplicate,1),:) = translated1;
                count = count+1;
                mirrored((count-1)*size(duplicate,1)+1:(count)*size(duplicate,1),:) = translated2;
                count = count+1;
            end
        end  
    end
end

mirrored((count-1)*size(duplicate,1)+1:end,:) = original;

mirrored_phi = mirrored(:,1);
mirrored_theta = mirrored(:,2);
mirrored_r = mirrored(:,3);
mirrored_values = mirrored(:,4);

numpoints = 100;
phi_range = linspace(0,2*pi(),numpoints);
theta_range = linspace(0,pi(),numpoints);
r_range = linspace(min(mirrored_r)+0.1,max(mirrored_r)-0.1,numpoints);
[phi_mesh, theta_mesh, r_mesh] = meshgrid(phi_range,theta_range,r_range);
values_mesh = griddata(mirrored_phi,mirrored_theta,mirrored_r,...
    mirrored_values.^2.*abs(sin(mirrored_theta)),phi_mesh,theta_mesh,r_mesh,'natural');

% %Count the number of NaN in 'values_mesh'. Was necessary when I was
% %testing methods to fill in the missing values. Taking advantange of
% %spherical continuity and re-casting onto a uniform grid of spherical 
% %coordinates using griddata fixed everything.
%nnz(isnan(values_mesh))

% %This was a cheat used before all the NaN were filled in. I just filled
% %in all NaN with zeroes, which only works with orbitals that are zero
% %near low/high phi/theta/r.
%values_mesh(isnan(values_mesh))=0;

%Integrate across phi and theta and remove extra dimensions
%'R_unnormalized' corresponds to R(r)^2*N^2 
% where N is a normalization constant
%Add physically intuitive point (0,0)
r_plot = [r_range'];
R2_unnormalized = [squeeze(trapz(theta_range,trapz(phi_range,values_mesh,2)))];
%total_int = N^2
total_int = trapz(r_plot,R2_unnormalized.*r_plot.^2)
%Calculate normalized R(r)
R = sqrt(R2_unnormalized/total_int);
r_plot = [0; r_plot];
R = [0; R];
%trapz(r_plot,R.^2.*r_plot.^2)
end
