function area = voronoidens(k);
%
% function area = voronoidens(k);
%
% input:  k = kx + i ky is the  k-space trajectory
% output: area of cells for each point 
%           (if point doesn't have neighbors the area is NaN)

kx = real(k);
ky = imag(k);

[row,column] = size(kx);

% uncomment these to plot voronoi diagram
%[vx, vy] = voronoi(kx,ky);
%plot(kx,ky,'r.',vx,vy,'b-'); axis equal

% scale k-space points to prevent "non-unique points excluded" operation


kxy = [kx(:),ky(:)]; 
% returns vertices and cells of voronoi diagram
[V,C] = voronoin(kxy); 
area = [];
for j = 1:length(kxy)
  x = V(C{j},1); y = V(C{j},2); lxy = length(x);
  A = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
  area = [area A];
end

area = reshape(area,row,column);
