function [dcf_out, ind_duplicates] = voronoi_volume(ktrajs)
    % Calculate voronoi density correction function for 3D k-space trajectory
    % Input: 
    %      ktrajs - c1 x c2 x 3
    % Output:
    %      dcf_out - c1 x c2
    %      ind_duplicates - indices in flattened array of values, discarded
    %      (their dcf values set to zero in dcf_out)
    % Gehua Tong, 02/20/2020
    
    d = length(size(ktrajs));
    if size(ktrajs,d) ~= 3
        disp("This is not 3D k-space!")
        return
    end
    
    M = size(ktrajs,1);
    N = size(ktrajs,2);
    % Ndim = 3 (ro, line)
    X = reshape(ktrajs, M*N, 3);
      
    % Find unique points and store
    dcf_out = zeros(size(X,1),1);
    % Calculate volume 
    % https://www.mathworks.com/matlabcentral/answers/354850-voronoi-tessellation-volume-calculation
    [C,IA,IC] = unique(X,'rows'); % All non unique points are removed
    B = 1:size(X,1);
    missingvalues = setdiff(B,IA); 
    [Y,~] = removerows(X,'ind', missingvalues);
   
    % Declare variable for storing voronoi volumes
    volumes = zeros(1,size(Y,1));
    disp(size(C))
    
%    [V,C] = voronoin(Y, {'Qbb'}); % Calculate Tessellation
    [V,CC] = voronoin(C, {'Qbb'});
    disp("Voronoin step complete")
    
    V(~isfinite(V)) = 0;
    
    % Iterate through all outputs and store 3D volume2
     for j=1:length(CC)
         disp("Calculating vol. " + j)
         X2 = [V(CC{j},1) V(CC{j},2) V(CC{j},3)];
         if (isempty(X2))
             vol = 0;
         else
             [~, vol] = convhulln(X2);
         end
         volumes(j) = vol;
     end
     
     % Recover into orignal shape
     % Make use of IA : ) 
     dcf_out(IA,1) = volumes;
     dcf_out = reshape(dcf_out, M,N);
     ind_duplicates = missingvalues;
end