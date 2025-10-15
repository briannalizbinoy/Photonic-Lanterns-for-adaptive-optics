%% 03-06-2023: 

% generate screens, slice into subapertures, propagate to focal plane,
% launch into grating coupler in RSoft, combine output fields (uses 2x1 bc), 
% calculate "SR"

dbstop if error
close all
%% Intialization
%close all
clearvars %-except WF_t
%clearvars %-except field_nT WF_n_ii
%set(0,'DefaultFigureWindowStyle','normal');
set(0,'DefaultFigureWindowStyle','docked');
t_start = tic;

%% Global inputs
lambda = 1.55e-6;       % [m] wavelength
k0 = 2*pi/lambda;       % [m^-1] wavenumber


%% Simulation settings
newWF = 1;              % generate new WF or use one from workspace 
scinti = 1;             % turn scintillation on/off 
FSM = 1;                % turn TT reomval on/off 
nZerm = 3;              % no. of Zernike modes to remove (mathematically)  
bc = 0;                 % switch RSoft beam combiner on/off 
gc = 0;                 % switch RSoft grating coupler on/off 
nP = 32;                % number of processors to use for FDTD calc 
hide = 1;               % Hide the Rsoft simulations 
gpuArrayUse = 0;        % turn on/off gpuArray for fourier transform 

%% RSoft path
%indFile_out = 'grating_coupler_FDTD3_upright_toy2.ind'; % .ind file name
% indFile_out = 'MSPL6_2.ind'; % .ind file name
indFile_out = 'free_space_test.ind'; % .ind file name
path = '';
% outfilename_pre = 'mspl6_bpm_ey.dat';                      % RSoft file for calculated output
outfilename_pre = 'bptmp.fld'; 
% outfilename_pre_mon = 'mspl6.txt';
outfilename_pre_mon = 'mspl6.mon';

% [a1,hostname] = system('hostname');                      % find hostname
% hostname = strip(hostname);

%RSoftCAD = ['fwmpirun -np',num2str(nP),' ', indFile_out,' wait=0 prefix=mspl6'];     % construct RSoft command
RSoftCAD = ['bsimw32 ', indFile_out,' wait=0']; % prefix=',outfilename_pre];     % construct RSoft command

% f = 2*0.5e-3;           % [m] focal length
f = 5;
    
% Telescope
D = 1;%(0.4/scFactor);                      % [m]  aperture diameter (def: 1.2)
Ds = 1*D;
    
% Atmosphere
r0 = 1*0.1559e0;          % [m] r0 at zenith
 
    
Dp = 4.*Ds;
% M = 2*round(Ds*subN*200);
M = 100;                    % [px] size of the wavefront
N = 4*M;                    % Number of points in the turbulence realization side = number of points of the focal plane side
PR = 1;                     % pupil ratio
N_ps = 1;               % number of phase screens used.
   
D_r0 = D./max(r0);      %[] turbulence strenght

    
% Fraunhofer diffraction
NN = (D./1)*4;       % [px] focal field size
zpf_d1 = 1.5; %3;           % desired zero padding factor for aperture function
zpf_d = 1; %3;             % desired zero padding factor for modes
    
%% Telescope aperture
%[x1,y1]   = meshgrid(linspace(-D_ML/2,D_ML/2,M));
[x1,y1]   = meshgrid(linspace(-D/2,D/2,M));
cx1 = 0;%MM/2;    % coordinate x of the circle centre [m]
cy1 = 0;%MM/2;    % coordinate y of the circle centre [m]
%r1 = 1*D_ML/2;     % radius of the circle [m]
r1 = 1*D/2;     % radius of the circle [m]
r_c = 0*r1;                  % obstruction radius
CA1 = (((x1-cx1).^2 + (y1-cy1).^2) <= r1^2) & (((x1-cx1).^2 + (y1-cy1).^2) >= r_c^2);
    
%% Reference plane wave
WF_0 = ones(M);
figure; imagesc(WF_0); colorbar
WF_0 = WF_0 .* CA1;
figure; imagesc(WF_0); colorbar
    
% construct complex field
field_0T = CA1.*exp(1i.*WF_0.*1);
total_power0T = sum(sum(abs(field_0T.^ 2))); % calculate total power
    
%% propagating the subfields
MM = size(WF_0,1);
% MM = size(subWF1,1);
zpf_m = (2^nextpow2(zpf_d1*MM))./MM;         % screens zero padding factor recalculated to get a matrix with a power of 2 dimensions
                                                              
% size_m = 4*2.44*lambda*f/Dsub;
size_m = 4*2.44*lambda*f/D;
size_m_x = NN;                      % size of focal field matrix [px]

U_zp_m = zeroPadMK(field_0T,round(size(field_0T,1)*zpf_m),0,'center');  % zeropadding: function out=zeroPadMK(in,padSize,padValue,type)
figure; imagesc(abs(U_zp_m)); colorbar
%%
dx = size_m/(size_m_x-1);
dxf = lambda*f/(D*zpf_m);

% mesh for zero-padded field at pupil plane 
xi = linspace(-(D*zpf_m)/2,(D*zpf_m)/2,MM*zpf_m);
yi = linspace(-(D*zpf_m)/2,(D*zpf_m)/2,MM*zpf_m);
[XI,YI] = meshgrid(xi,yi);

% mesh for focal plane field
xxi = linspace(-dx*size(U_zp_m,1)/2,dx*size(U_zp_m,1)/2,size(U_zp_m,1));
yyi = linspace(-dx*size(U_zp_m,2)/2,dx*size(U_zp_m,2)/2,size(U_zp_m,2));
[xx,yy] = meshgrid(xxi,yyi);

focus_cut = 30;
focal_field_size = size(U_zp_m,1);

range = floor(focal_field_size/2)+1 - floor(focus_cut/2) : floor(focal_field_size/2) + ceil(focus_cut/2);

% calculate focal fields
focal_field = ifftshift(fft2(ifftshift((1i*lambda*f)^-1.*(exp(-1i*(2*pi/lambda)*(f+((xx.^2+yy.^2)./(2*f))))).*U_zp_m)));
figure; imagesc(abs(focal_field)); colorbar

% Cutting the propagated fields
focal_field1 = focal_field(range,range);
figure; imagesc(abs(focal_field1)); colorbar

% Interpolation (higher resolution image)
xxxi = linspace(-(D*zpf_m)/2,(D*zpf_m)/2,size(focal_field1,1));
yyyi = linspace(-(D*zpf_m)/2,(D*zpf_m)/2,size(focal_field1,2));
[xxx,yyy] = meshgrid(xxxi,yyyi);

focal_field_cut_i_0 = interp2(xxx,yyy,focal_field1,XI,YI,'cubic',0);
figure; imagesc(abs(focal_field_cut_i_0).^2); colorbar

pp_i = (dxf*focus_cut)/(size(focal_field_cut_i_0,1)-1);   % [m/px] pixel pitch

psf_I = abs(focal_field_cut_i_0).^2;
% psf_I_0 = psf_I;


%%
xa = linspace(-(pp_i.*size(focal_field_cut_i_0(:,:,1,1),1))/2, (pp_i.*size(focal_field_cut_i_0(:,:,1,1),1))/2, size(focal_field_cut_i_0(:,:,1,1),1));

% filename_1D = ['PS_MLA_GC_PSF_1D_',num2str(1),'_',num2str(1),'.fld'];
filename_2D = 'focus_2_pl.fld';
        
[fieldMax, maxRow] = max(max(abs(focal_field_cut_i_0)));

% Output_to_Rsoft2_1D(focal_field_cut_i_0(maxRow,:) + 1e-16,filename_1D,'REAL_IMAG',pp_i.*size(focal_field_cut_i_0(:,:),1),lambda);
Output_to_Rsoft2(focal_field_cut_i_0, filename_2D, 'REAL_IMAG', pp_i.*size(focal_field_cut_i_0,1), lambda);

status = 1;
while status ~= 0
    [status,results] = system(RSoftCAD);

    if status ~= 0
        disp('RSoft error');
        disp(status);
        %disp(results);
        %pause(); 
    end
end 

% [OutputFields_0, pp_out] = RSoft_import3_1D(path,outfilename_pre);
%[OutputFields_0, pp_out] = RSoft_import2(path,[outfilename_pre,'.fld']);
[OutputFields_0, pp_out] = RSoft_import2(path,outfilename_pre);

% monPower_0 = RSoft_readMon2(path,outfilename_pre_mon);  % if using R2020 or later use RSoft_readMon2
figure; imagesc(abs(OutputFields_0)); colorbar

%Isolate each mode

%Mode1
range1 =  200 : 300;
range2 =  350 : 450;
mode1 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode1)); colorbar

%Mode2
range1 =  320 : 420;
range2 =  270 : 370;
mode2 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode2)); colorbar

%Mode3
range1 =  350 : 400;
range2 =  150 : 250;
mode3 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode3)); colorbar

%Mode4
range1 =  200 : 300;
range2 =  70 : 150;
mode4 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode4)); colorbar

%Mode5
range1 =  100 : 200;
range2 =  150 : 250;
mode5 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode5)); colorbar

%Mode6
range1 =  100 : 200;
range2 =  300 : 400;
mode6 = OutputFields_0(range1,range2);
figure; imagesc(abs(mode6)); colorbar

range1list = {200, 320, 350, 200, 100, 100};
range2list = {350, 270, 150, 70, 150, 300};

% Store all modes in a cell array
modes = {mode1, mode2, mode3, mode4, mode5, mode6};
rowIndices = [];
colIndices = [];
maxValues = [];

% Loop through each mode
for i = 1:length(modes)
    % Get the current mode matrix
    currentMode = modes{i};
    
    [maxValue, rowIdx] = max(max(abs(currentMode))); % Get maximum value and row index
    [~, colIdx] = max(abs(currentMode(:, rowIdx))); % Get column index for that row
    
    % Display results
    disp(['Mode ', num2str(i)]);
    disp(['  Maximum Amplitude: ', num2str(maxValue)]);
    disp(['  Row Index: ', num2str(rowIdx)]);
    disp(['  Column Index: ', num2str(colIdx)]);

    rowIndices = [rowIndices, rowIdx];
    colIndices = [colIndices, colIdx];
    maxValues = [maxValues, maxValue];
end


% parameters
M = size(OutputFields_0, 1);  
D = 500;                      

%meshgrid for coordinates
[x1,y1]   = meshgrid(linspace(0,D,M));
powers = zeros(1, length(modes));
for i = 1:length(modes)
    % center coordinates and radius for the bright spot
    cx1 = range2list{i}+rowIndices(i);  
    cy1 = range1list{i}+colIndices(i);  
    r1 = 55;            
    %circular mask centered at (cx1, cy1) with radius r1
    CA1 = ((x1 - cx1).^2 + (y1 - cy1).^2) <= r1^2;
    %circular mask applied
    maskedImage = abs(OutputFields_0) .* CA1;
    power = sum(sum(abs(maskedImage).^ 2));
    powers(i) = power;
    %figure
    figure;
    imagesc(maskedImage); axis image;
    colorbar;
end

figure; imagesc(abs(focal_field_cut_i_0)); colorbar
figure; imagesc(abs(OutputFields_0)); colorbar

[xx,yy] = meshgrid(xxi,yyi);
F1 = abs(focal_field_cut_i_0).^2;
input_power = trapz(yyi,trapz(xxi,F1,2));

[OutputFields_0, pp_out] = RSoft_import2(path,outfilename_pre);


%create the meshgrid

[rows, cols] = size(OutputFields_0);

xxo = linspace(-pp_out*size(OutputFields_0,1)/2,pp_out*size(OutputFields_0,1)/2,size(OutputFields_0,1));
yyo = linspace(-pp_out*size(OutputFields_0,2)/2,pp_out*size(OutputFields_0,2)/2,size(OutputFields_0,2));

% Create the meshgrid
[xx2, yy2] = meshgrid(xxo, yyo);

F2 = abs(OutputFields_0).^2;
output_power = trapz(yyo,trapz(xxo,F2,2));

disp(['Total Input Power: ', num2str(input_power)])
disp(['Total Output Power: ', num2str(output_power)])

disp(['Input Power: ',num2str(sum(sum(abs(focal_field_cut_i_0).^2)))])
disp(['Total Power for Each Mode: ', num2str(powers)])
disp(['Total Output Power: ', num2str(sum(powers))])



%{
outPhase_0 = angle(OutputFields_0(peakInd));
outAmp_0 = abs(OutputFields_0(peakInd));
zField_0 = OutputFields_0(peakInd);
%}        