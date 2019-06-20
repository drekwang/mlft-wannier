%written on 4/12/19 by Derek Wang
%input: 'cartcoords' from 'gen_cartcoords.m'
%output: 'sphcoords', a 3 x N array

%CONSIDER CUTTING OUT ALL POINTS BEYOND lengths/2 TO MAINTAIN SPHERICAL
%ARRAY OF POINTS

function [sphcoords] = gen_sphcoords(cartcoords)

%Transpose for cart2sph comptability
cartcoords = cartcoords';
[phi, theta, r] = cart2sph(cartcoords(1,:),cartcoords(2,:),cartcoords(3,:));
%In cart2sph, phi is -pi->pi, theta is -pi/2->pi/2, and z is 0->inf
%phi is the counterclockwise angle from positive x axis and theta is angle
%from xy plane to positive z axis.
sphcoords = [phi;theta;r]';
%Convert spherical coordinates to standard convention ('theta' is angle
%from positive z axis and ranges from 0 to pi, and 'phi' is counterclockwise angle
%from positive x axis and ranges from 0 to 2pi.
for i = 1:size(sphcoords,1)
    if sphcoords(i,1)<0
        sphcoords(i,1)=sphcoords(i,1)+2*pi();
    end
    sphcoords(i,2) = pi()/2 - sphcoords(i,2);
end
            
