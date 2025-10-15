% function out=zeroPadMK(in,padSize,padValue,type)
%
% Zero Padding of the input matrix
% Input
%  in       Input matrix
%  padSize  Size of the output matrix, either a scalar, if both dimensions
%           are equal, or a vector with 2 elements [x,y].
%  padValue Padding values (optional, default zero)
%  type     Padding type (optional, default 'center')
%           'center' - IN matrix in the center of the output
%           'lefttop' - IN matrix in the left top corner of the output
%           'righttop' - IN matrix in the right top corner of the output
%           'rightbottom' - IN matrix in the right bottom corner of the output
%           'leftbottom' - IN matrix in the left bottom corner of the output
% Return
%  out      Output matrix of size padSize
% MK 19.2.08, 11.2.09,15.5.09

function out=zeroPadMK(in,padSize,padValue,type)

if (nargin==2) type='center'; padValue=0; end;                                  %nargin returns the number of input arguments passed in the call (in this case in "zeroPadMK")
if (nargin==3) type='center'; end;                                              %Example: if only "in" and "padSize" nargin=2 then asign by default padValue=0 and type='center'

if (strcmp(type,'center')),
    if (length(padSize)==1), padX=padSize; padY=padSize; end;                   %padSize has been defined with equal [X,Y] dimensions
    if (length(padSize)==2), padX=padSize(1); padY=padSize(2); end;             %padSize has been defined with [X,Y] dimensions independently 
    s=size(in);
    if (s(1)>padY || s(2)>padX) out=in; fprintf('Warning: input image larger than pad size\n'); return; end;
    
    out=ones(padY,padX)*padValue;       %out=ones(padY,padX)*padValue; out=[257,257]*padValue                     %out=ones(padY,padX)*padValue; out=[257,257]*padValue 
    out(floor(padY/2-s(1)/2)+1:floor(padY/2-s(1)/2)+s(1),floor(padX/2-s(2)/2)+1:floor(padX/2-s(2)/2)+s(2))=in;    %out(1:256,1:256) because floor(padY/2-s(1)/2)+1 =floor(1.5) : floor(padY/2-s(1)/2)+s(1)=floor(256.5)
end;                                                                                                              %Idem with PadX 
    
if (strcmp(type,'lefttop')),
    if (length(padSize)==1), padX=padSize; padY=padSize; end;
    if (length(padSize)==2), padX=padSize(1); padY=padSize(2); end;
    s=size(in);
    if (s(1)>padY || s(2)>padX) out=in; fprintf('Warning: input image larger than pad size\n'); return; end;
    
    out=ones(padY,padX)*padValue;
    out(1:s(1),1:s(2))=in;                                                   %  out(1:s(1),1:s(2))=in  --> out(1:256,1:256)=in    if S=size(in) and in=phase (phase comes from the branchpoint_function.m )        
end;  

if (strcmp(type,'rightbottom')),
    if (length(padSize)==1), padX=padSize; padY=padSize; end;
    if (length(padSize)==2), padX=padSize(1); padY=padSize(2); end;
    s=size(in);
    if (s(1)>padY || s(2)>padX) out=in; fprintf('Warning: input image larger than pad size\n'); return; end;
    
    out=ones(padY,padX)*padValue;                                            %  out=ones(padY,padX)*padValue; out=[257,257]*padValue  
    out(padY-s(1)+1:padY,padX-s(2)+1:padX)=in;                               %  out(padY-s(1)+1:padY,padX-s(2)+1:padX)=in -->   out(256+1-256+1:256+1,256+1-256+1:256+1)=out(2:257,2:257)=in
end;  

if (strcmp(type,'righttop')),
    if (length(padSize)==1), padX=padSize; padY=padSize; end;
    if (length(padSize)==2), padX=padSize(1); padY=padSize(2); end;
    s=size(in);
    if (s(1)>padY || s(2)>padX) out=in; fprintf('Warning: input image larger than pad size\n'); return; end;
    
    out=ones(padY,padX)*padValue;                                            %  out=ones(padY,padX)*padValue; out=[257,257]*padValue    
    out(1:s(1),padX-s(2)+1:padX)=in;                                         %  out(1:s(1),padX-s(2)+1:padX)=in  -->  out(1:256,257-256+1:257)=out(1:256,2:257)=in
end;  

if (strcmp(type,'leftbottom')),
    if (length(padSize)==1), padX=padSize; padY=padSize; end;
    if (length(padSize)==2), padX=padSize(1); padY=padSize(2); end;
    s=size(in);
    if (s(1)>padY || s(2)>padX) out=in; fprintf('Warning: input image larger than pad size\n'); return; end;
    
    out=ones(padY,padX)*padValue;                                            %  out=ones(padY,padX)*padValue; out=[257,257]*padValue   
    out(padY-s(1)+1:padY,1:s(2))=in;                                         %  out(padY-s(1)+1:padY,1:s(2))=in  -->  out(2:257,1:256)=in
end;  