% Copyright (c) 2010, Felix Totir
% Copyright (c) 2009, Ian Howat
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function enviwrite(D,info,datafile,hdrfile)
%enviread: Read binary image files using ENVI header information
%[D,info,x,y]=enviread(datafile) where the header file is named
%filename.hdr and exists in the same directory. Otherwise use
%[D,info,x,y]=enviread(datafile,hdrfile) The output contains D, info, x and
%y containing the images data, header info, x-coordinate vector and
%y-coordinate vector, respectively. D will be in whatever number format
%(double, int, etc.) as in the ENVI file.

if nargin < 3
    error('You must specify at least three inputs');
elseif nargin < 4
    hdrfile=[deblank(datafile),'.hdr']; %implicit name
end

envihdrwrite(info,hdrfile);
envidatawrite(D,datafile,info);


function envidatawrite(D,datafile,info)
%enviread: Read binary image files using ENVI header information
%D=envidataread(datafile,info) The output contains D, info, x and
%y containing the images data, header info, x-coordinate vector and
%y-coordinate vector, respectively. D will be in whatever number format
%(double, int, etc.) as in the ENVI file.

%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Thanks to Yushin Ahn and Ray Jung
%Heavily modified by Felix Totir.

if nargin < 2
    error('You must specify the information about data (input of envihdrwrite)');
end

%-improved- assert(info.header_offset==0,'Non-zero header offset not supported');

%% Set binary format parameters
switch info.byte_order
    case {0}
        machine = 'ieee-le';
    case {1}
        machine = 'ieee-be';
    otherwise
        machine = 'n';
end

iscx=false; %if it is complex
switch info.data_type
    case {1}
        format = 'uint8';
    case {2}
        format= 'int16';
    case{3}
        format= 'int32';
    case {4}
        format= 'single';
    case {5}
        format= 'double';
    case {6}
        iscx=true;
        format= 'single';
    case {9}
        iscx=true;
        format= 'double';
    case {12}
        format= 'uint16';
    case {13}
        format= 'uint32';
    case {14}
        format= 'int64';
    case {15}
        format= 'uint64';
    otherwise
        error(['File type number: ',num2str(dtype),' not supported']);
end

% *** BSQ Format ***
% [Band 1]
% R1: C1, C2, C3, ...
% R2: C1, C2, C3, ...
%  ...
% RN: C1, C2, C3, ...
%
% [Band 2]
%  ...
% [Band N]

% *** BIL Format ***
% [Row 1]
% B1: C1, C2, C3, ...
% B2: C1, C2, C3, ...
%
%  ...
% [Row N]

% *** BIP Format ***
% Row 1
% C1: B1 B2 B3, ...
% C2: B1 B2 B3, ...
% ...
% Row N

switch lower(info.interleave)
    case {'bsq'}
        %D = reshape(D,[info.samples,info.lines,info.bands]);
        %D = permute(D,[2,1,3]);
        D = permute(D,[2,1,3]); %adica: D = ipermute(D,[2,1,3]);
    case {'bil'}
        %D = reshape(D,[info.samples,info.bands,info.lines]);
        %D = permute(D,[3,1,2]);
        D = permute(D,[2,3,1]); %adica: D = ipermute(D,[3,1,2]);
    case {'bip'}
        %D = reshape(D,[info.bands,info.samples,info.lines]);
        %D = permute(D,[3,2,1]);
        D = permute(D,[3,2,1]); %adica: D = ipermute(D,[3,2,1]);
end

D=D(:);
if iscx
    D = [real(D), imag(D)]; %2 columns
    D=D.';
    D=D(:);
end

fid=fopen(datafile,'w');
fwrite(fid,zeros(info.header_offset,1),'uint8',0,machine); %here we write a dummy header in order to respect the offset
fwrite(fid,D,format,0,machine); %alternatively, multibandwrite (matlab function) might be used
fclose(fid);

function envihdrwrite(info,hdrfile)
% ENVIHDRWRITE read and return ENVI image file header information.
%   INFO = ENVIHDRWRITE(info,'HDR_FILE') writes the ENVI header file
%   'HDR_FILE'. Parameter values are provided by the fields of structure
%   info.
%
%   Example:
%   >> info = envihdrread('my_envi_image.hdr')
%   info =
%          description: [1x101 char]
%              samples: 658
%                lines: 749
%                bands: 3
%        header_offset: 0
%            file_type: 'ENVI Standard'
%            data_type: 4
%           interleave: 'bsq'
%          sensor_type: 'Unknown'
%           byte_order: 0
%             map_info: [1x1 struct]
%      projection_info: [1x102 char]
%     wavelength_units: 'Unknown'
%           pixel_size: [1x1 struct]
%           band_names: [1x154 char]
%     info = enviwrite('my_envi_image2.hdr');
%
% Felix Totir

params=fieldnames(info);

fid = fopen(hdrfile,'w');
fprintf(fid,'%s\n','ENVI');

for idx=1:length(params)
    param = params{idx};
    value = getfield(info,param);
    param(findstr(param,'_')) = ' '; %automatic name
    
    line=[param,' = ',num2str(value)];
    
    fprintf(fid,'%s\n',line);
    
end

fclose(fid);

