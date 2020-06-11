function dcf_areas = find_dcf(filepath,opt_save,savename)
% FIND_DCF    returns density correction function (dcf) for arbitrary 2D
% k-space trajectory based on voronoi areas. 
% Parametrized version of find_dcf_voronoi to directly support ktraj
% Gehua Tong, Sept 2019 
% Inputs:
%   - filepath : path to kspace trajectory file (.mat) 
%                The loaded data should have a field called ktrajs that
%                has at most 3 dimensions with the last dimension = 2
%   - opt_save : whether to save the dcf areas as a .mat file. Default is
%                true.
%   - savename : name by which to save the dcf file.

% Dependencies : VORONOI_AREA() by : 
    % Created Florian Wiesinger
    % Modified 7/2006 Rolf Schulte
    

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

    if sizek(end) ~= 2
        error("Last dimension of k-trajs should always be 2")
    end
    
    % Convert into complex array for input into voronoi_area()
    
    if length(sizek) == 3
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