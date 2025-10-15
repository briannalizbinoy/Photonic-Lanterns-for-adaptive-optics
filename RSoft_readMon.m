function mon = RSoft_readMon(path,filename_pre)

%path = 'C:\Users\momen\Documents\MATLAB\GC\';

%filename_pre = 'x1.txt';

%% Load output fields
%path = '/home/INNOFSPEC/mdiab/Documents/my-RSoft/Designs/3M/1x7PLx27SMF-asy-Lopt-degen-Rev1//';
%AA = dir(join([path,'1x7PLx27SMF-asy-Lopt-Degen-Rev1-0*.*'],""));

%path = path;
AA = dir(join([path,filename_pre],""));


filename = AA(1).name; % read the name of the first field file to get size data
C = readcell([path,filename]);
mon = -C{2,2};
%sz_f_px = dlmread([path,filename],'',[1 0 1 0]);               % size of one field matrix [px]
% sz_f = 1*1e-6*abs(C{2,3} - C{2,2});      % spatial extent of one field [m]
% pp = sz_f/(sz_f_px-1);  % pixel pitch [m/px]
% 
% %% uniterleave output fields
% 
%         filename = [filename_pre];
%  
%         M = dlmread([path,filename],'',2,0); % save mode data points to M
%     
%         %% uninterleave
%         for j = 1:sz_f_px
%             %for k = 1:2
%                 OutputFields(j) = M(j,1)*cosd(M(j,2)) + 1i.*M(j,1)*sind(M(j,2));    % k & j are swapped because in MATLAB first index is row index, i.e. y-coordinate
%             %end
%         end
%   
%         %OutputFields(:) = flipud(OutputFields(:));     % flip upside-down is necessary because in RSoft y-coordinate increases from bottom to top
% %end
