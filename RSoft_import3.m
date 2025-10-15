%function [OutputFields, pp] = RSoft_import2(path,filename_pre)

path = '/work3/INNOFSPEC/mdiab/Documents/MATLAB/Closure/';

filename_pre = 'Closure_EP.fld';

%% Load output fields
%path = '/home/INNOFSPEC/mdiab/Documents/my-RSoft/Designs/3M/1x7PLx27SMF-asy-Lopt-degen-Rev1//';
%AA = dir(join([path,'1x7PLx27SMF-asy-Lopt-Degen-Rev1-0*.*'],""));

path = path;
AA = dir(join([path,filename_pre],""));


filename = AA(1).name; % read the name of the first field file to get size data
sz_f_px = dlmread([path,filename],'',[3 0 3 0]);               % size of one field matrix [px]
sz_f = 2*1e-6*abs(dlmread([path,filename],'',[3 1 3 1]));      % spatial extent of one field [m]
pp = sz_f/(sz_f_px-1);  % pixel pitch [m/px]

%% uniterleave output fields

        filename = [filename_pre];
 
        M = dlmread([path,filename],'',4,0); % save mode data points to M
    
        %% uninterleave
        for j = 1:sz_f_px
            for k = 1:sz_f_px
                OutputFields(k,j) = M(j,2*k-1) + 1i*M(j,2*k);    % k & j are swapped because in MATLAB first index is row index, i.e. y-coordinate
            end
        end
  
        OutputFields(:,:) = flipud(OutputFields(:,:));     % flip upside-down is necessary because in RSoft y-coordinate increases from bottom to top
%end
