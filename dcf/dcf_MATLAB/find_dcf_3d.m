% Parametrized version of find_dcf_voronoi to directly support ktraj
% 3D version
% Gehua Tong, Oct 2019 

function dcf_areas = find_dcf_3d(filepath,opt_save,savename)
    % Only works for 2D right now 
    % Default save option
    if nargin < 3
        savename = 'dcf.mat';
    end
    if nargin < 2
        opt_save = true;
    end
    
    % Load k-space trajectory
    data = load(filepath);    
    if ~isfield(data,'ktrajs')
        error("ktrajs should be a field in loaded mat file");
    end
    
    ktrajs = data.ktrajs;
    sizek = size(ktrajs);

    if sizek(end) ~= 3
        error("Last dimension of k-trajs should be 3; use find_dcf() for 2D")
    end
    
    % Convert into complex array for input into voronoi_area()
    if length(sizek) == 4
        kk = reshape(ktrajs, sizek(1)*sizek(2)*sizek(3),sizek(4));
    elseif length(sizek) == 3
        kk = reshape(ktrajs, sizek(1)*sizek(2),sizek(3));
    elseif length(sizek) == 2
        kk = ktrajs;
        
    end
    
    kxy = (kk(:,1)+1i*kk(:,2))';

    % Scale k-space points to prevent "non-unique points excluded" operation
    dcf_out = voronoi_area(kxy);

    % Normalize and reshape
    dcf_out = dcf_out/max(dcf_out(:));
    if length(sizek) == 3
        dcf_out = reshape(dcf_out', sizek(1),sizek(2));
    end
    
    if opt_save
        save(strcat(savename), 'dcf_out');
    end
    
    dcf_areas = dcf_out;
end