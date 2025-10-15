function [OutputFields, pp] = RSoft_import2(path,filename)

%% Import output fields from RSoft into MATLAB

% file is: path/filename
% file must have .fld extension
% file must be in RSoft "OUTPUT_REAL_IMAG_3D" format
% pp is the pixel pitch in m/px

AA = dir(join([path,filename],""));


filename = AA(1).name; % read the name of the first field file to get size data
sz_f_px = dlmread([path,filename],'',[3 0 3 0]);               % size of one field matrix [px]
sz_f = 2*1e-6*abs(dlmread([path,filename],'',[3 1 3 1]));      % spatial extent of one field [m]
pp = sz_f/(sz_f_px-1);  % pixel pitch [m/px]

%% uniterleave output fields
 
        M = dlmread([path,filename],'',4,0); % save mode data points to M
    
        %% uninterleave
        for j = 1:sz_f_px
            for k = 1:sz_f_px
                OutputFields(k,j) = M(j,2*k-1) + 1i*M(j,2*k);    % k & j are swapped because in MATLAB first index is row index, i.e. y-coordinate
            end
        end
  
        OutputFields(:,:) = flipud(OutputFields(:,:));     % flip upside-down is necessary because in RSoft y-coordinate increases from bottom to top
end
