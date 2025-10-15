function Output_to_Rsoft2(field,filename,part,size,lambda) %field=imaginary field distribution, filename, size must be the same in x and y, 
% field = A_temp;
% filename = 'field_dist5';
% part = 'REAL_IMAG';
% size = focal_field_m;
% lambda =lambda;
%filename = fullfile(savepath, 'myFile.txt');
    %field(:,:) = flipud(field(:,:));     % flip upside-down is necessary because in RSoft y-coordinate increases from bottom to top
    
    fid = fopen(filename, 'wt');
    % create header here
    part1='/rn,a,b/nx0/ls1'; %This is a comment for RSoft defining how the data shall be plotted
    part2='/r,qa,qb';
    %example of Rsoft header
    %441 -44 44 0 OUTPUT_AMPLITUDE_3D 1.498135354 -1.133978591e-012 Wavelength=1.5
    %441 -44 44
    part3=[num2str(length(field(:,1))) ' ' '-' num2str(size/2*1E6) ' ' num2str(size/2*1E6) ' 0 OUTPUT_' part '_3D Wavelength=' num2str(lambda*1E6)];
    part4=[num2str(length(field(1,:))) ' ' '-' num2str(size/2*1E6) ' ' num2str(size/2*1E6)];
            
    fprintf(fid, '%s\n%s\n%s\n%s\n', part1,part2,part3,part4);  % header
    for i=1:length(field(1,:))
        for j=1:length(field(:,1))
            fprintf(fid,'%E\t%E\t',real(field(j,i)),imag(field(j,i)));
%            fprintf(fid,'%E\t%E\t',real(field(j,i)),imag(field(j,i)),'delimiter','\t');
            %disp(['array:' num2str(real(field(j,i)))]); %just for debugging
            
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    %dlmwrite(filename,field,'delimiter','\t','-append');