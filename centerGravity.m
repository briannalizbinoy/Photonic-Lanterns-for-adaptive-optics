
% SH = Ir;
% Npixels = sideML;
%
% % circular mask
% cx1 = (M/2)+1;      % coordinate x of the circle centre [px]
% cy1 = (M/2)+1;      % coordinate y of the circle centre [px]
% r1 = M/2*PR;        % radius of the circle [px]
% [x1,y1] = meshgrid(1:M+1,1:M+1);
% pupil = (x1-cx1).^2 + (y1-cy1).^2 <= r1^2;

%M = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 3 4];

function [Cx,Cy] = centerGravity(M)

%% centroid.m
%% centroiding algorithm for the evaluation of Shack-Harmanngram
%% input:
%% M = matrix contianing the Shack-Hartmanngram
%% Npixels = number of pixel in the side of a ROI (should be even)

[Nrow,Ncol]=size(M);

gridx=1:Ncol; %coordinate x of the computational grid
gridy=1:Nrow; %coordinate y of the computational grid

[xm,ym]=meshgrid(gridx,gridy);

%%pixel (1,1) is on the top left of the input matrix
Tp=sum(sum(M));
Cx =sum(sum(M.*xm))/(Tp); %centroid x
Cy =sum(sum(M.*ym))/(Tp); %centroid y


