function [target, n, mm] = zernikeGenerator2(Zxx, M)

%Zxx = 1;

%load ('ALPAO DM97-15', 'wfsMask');

% load look table for OSA- Noll conversion
C = table2array(readtable('b176988.txt'));
j = C(C(:,2)==Zxx);

% find n and mm indicies from j
n = ceil((-3+sqrt(9+8*j))./2);
mm = 2*j-n*(n+2);

%% generate coordinate grid
[X,Y]   = meshgrid(linspace(-1,1,M));
r       = sqrt(X.^2+Y.^2);  
t       = angle(X+Y*1i);

m = abs(mm);
R = zeros(size(r));
for S = 0 : (n-m)/2
    R = R + ((-1).^S.*factorial(n-S).*r.^(n-2.*S))./(factorial(S).*factorial(((n+m)./2)-S).*factorial(((n-m)./2)-S)); 
end

if mm == 0
    target = sqrt(n+1).*R;
else
    switch mm > 0
        case 1
            target = sqrt(n+1).*R.*sqrt(2).*cos(m.*t);
        case 0
            target = sqrt(n+1).*R.*sqrt(2).*sin(m.*t);
    end
end
 

% switch Zxx
%     case 1
%         target = zeros(size(wfsMask,1));
%     case 2
%         target = X;
%     case 3
%         target = Y;
%     case 4
%         target = sqrt(3)*(2*r.^2-1);     
%     case 5
%         target = sqrt(6)*r.^2.*sin(2*t);
%     case 6
%         target = sqrt(6)*r.^2.*cos(2*t);
%     case 11
%         target = sqrt(5)*(6*r.^4-6*r.^2+1);
%     case 9
%         target = sqrt(8)*r.^3.*sin(3*t);        
%     case 7
%         target = sqrt(8)*(3*r.^3 - 2*r).*sin(t);
%     case 8
%         target = sqrt(8)*(3*r.^3 - 2*r).*cos(t);
%     case 10
%         target = sqrt(8)*r.^3.*cos(3*t);
%     case 15
%         target = sqrt(10)*r.^4.*sin(4*t);
%     case 13
%         target = sqrt(10)*(4*r.^4 - 3*r.^2).*sin(2*t);
%     case 12
%         target = sqrt(10)*(4*r.^4 - 3*r.^2).*cos(2*t);
%     case 14
%         target = sqrt(10)*r.^4.*cos(4*t);
% end
% generate target wavefront
%target = Zxx;

% remove data outside pupil
cx1 = 0;%round(M/2)+1;      % coordinate x of the circle centre [px]
cy1 = 0;%round(M/2)+1;      % coordinate y of the circle centre [px]
r1 = 1;%M/2;        % radius of the circle [px]
%[x1,y1] = meshgrid(1:M,1:M);
pupil = (X-cx1).^2 + (Y-cy1).^2 <= r1^2;

target = target.*pupil;
%target(~wfsMask) = 0; 

% scale target
%target = target ./ p2v(target)*40;
% remove piston from the target wavefront
%target = target - mean(target(wfsMask));
%target = target - mean(target(pupil));
target = target ./ max(max(target));

