%written on 4/12/19 by Derek Wang
%input: 'npts' and 'lengths' in each direction, 'origin'. All 1x3 arrays
%output: 'x/y/z_mesh' is necessary input for 'isosurface' visualization.
%cartcoords is a 3 x N array of all coordinates in column-major mode

function [ x_mesh, y_mesh, z_mesh, cartcoords ] = gen_cartcoords(npts,...
    lengths, origin,atomcenter)
%origin = [-4.013472 -4.013472 -4.013472]
%lengths = [  11.6065277   11.6065277   11.6065277];
%npts = [108 108 108];
%atomcenter = [1.95 1.95 1.95];

%Create 3 x npts matrix of possible x, y, and z values
cartpts = [linspace(origin(1)-atomcenter(1),origin(1)+lengths(1)-atomcenter(1),npts(1));...
    linspace(origin(2)-atomcenter(2),origin(2)+lengths(2)-atomcenter(2),npts(2));...
    linspace(origin(3)-atomcenter(3),origin(3)+lengths(3)-atomcenter(3),npts(3))];
[x_mesh,y_mesh,z_mesh] = meshgrid(cartpts(1,:),cartpts(2,:),cartpts(3,:));
%For some reason, to output in column-major mode (see xsf file format spec
%in xcrysden), y(:) has to go first. We could swap columns 1 and 2, but
%since in this situation, lattice vectors in x and y are the same, we do it
%the lazy way by leaving it unchanged. Tranpose for cart2sph compatibility
cartcoords = [y_mesh(:),x_mesh(:),z_mesh(:)];

end

