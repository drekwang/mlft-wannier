function [ values, values_3D ] = gen_values(xsffile, npts)
%gen_values : read Wannier function values from Wannier90 output
%input: 'xsffile' output by Wannier90 with 3D Cartesian scalar field
%output: 1 x N(umber of scalar field values) and cube of npts^3
%representing where each value is for 'isosurface' visualization

%'data' contains the scalar 'values' in 'xsffile', in which all lines up 
%to the scalar values have been deleted
data=importdata(xsffile,' ');
values = [];
%The scalar values are listed in columns (probably to minimize number of
%lines used), so this loop concatenates them all together
for i = 1:size(data,1)
    values = [values,data(i,:)];
end
%This converts the list of 'values' into a 3D array to plot with isosurface
values_3D = zeros(npts(1),npts(2),npts(3));
count = 1;
for k = 1:npts(1)
    for j = 1:npts(2)
        for i = 1:npts(3)
            values_3D(i,j,k) = values(count);
            count = count+1;
        end
    end
end
values = values';

end

