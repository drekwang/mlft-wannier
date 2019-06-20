%written on 3/30/19

%points is an N x 3 array of k points generated with k-path selection in
%xcrysden, e.g. [0.5 0 0; 0 0 0] is 2 x 3 array of real form of k-point
%coordinates M and Gamma

%This code is not very smart, so make sure steps is divisible by
%N-1
function kpath = gen_kpath(points,steps)

numpoints = size(points,1);

steps_per_k = steps/(numpoints-1);
for i = 1:numpoints-1
    for j = 1:steps_per_k
    disp([num2str(points(i,1)+j/steps_per_k*(points(i+1,1)-points(i,1)))...
    ' ' num2str(points(i,2)+j/steps_per_k*(points(i+1,2)-points(i,2)))... 
    ' ' num2str(points(i,3)+j/steps_per_k*(points(i+1,3)-points(i,3)))... 
    ' 1'])
    end
end

        
            

