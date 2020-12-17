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

function [D,info]=enviread(datafile,hdrfile)
%ENVIREAD Reads ENVI image file.
%[D,INFO]=ENVIREAD(DATAFILE,HDRFILE) provides images data (D) and a
%structure (INFO) whose fields correspond to header items.
%[D,INFO]=ENVIREAD(DATAFILE) assumes the header file is named
%"DATAFILE.hdr" and exists in the same directory.
%D will be in whatever number format (double, int, etc.) as in the ENVI
%file.

%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Thanks to Yushin Ahn and Ray Jung
%Heavily modified by Felix Totir.

if nargin < 1
    error('You must specify at least one input');
elseif nargin <2
    hdrfile=[deblank(datafile),'.hdr']; %implicit name
end

info=envihdrread(hdrfile);
D=envidataread(datafile,info);


function D=envidataread(datafile,info)
%ENVIDATAREAD Reads data of ENVI image file according to ENVI image file header information.
%D=ENVIDATAREAD(DATAFILE,INFO) provides image data (D).
%D will be in whatever number format (double, int, etc.) as in the ENVI
%file.

%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Thanks to Yushin Ahn and Ray Jung
%Heavily modified by Felix Totir.

if nargin < 2
    error('You must specify the information about data (output of envihdrread)');
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

n = info.lines*info.samples*info.bands;

fid=fopen(datafile,'r');
fread(fid,info.header_offset,'uint8',0,machine); %we skip the header
if ~iscx
    D = fread(fid,n,format,0,machine); %alternatively, multibandreader (matlab function) might be used
else
    D = fread(fid,2*n,format,0,machine); %alternatively, multibandreader (matlab function) might be used
    D = complex(D(1:2:end),D(2:2:end));
end
fclose(fid);
    
switch lower(info.interleave)
    case {'bsq'}
        D = reshape(D,[info.samples,info.lines,info.bands]);
        D = permute(D,[2,1,3]);
    case {'bil'}
        D = reshape(D,[info.samples,info.bands,info.lines]);
        D = permute(D,[3,1,2]);
    case {'bip'}
        D = reshape(D,[info.bands,info.samples,info.lines]);
        D=permute(D,[3,2,1]);
end


function info = envihdrread(hdrfile)
% ENVIHDRREAD Reads header of ENVI image.
%   INFO = ENVIHDRREAD('HDR_FILE') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters.
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
%
%   NOTE: This function is used by ENVIREAD to import data.

% Ian M. Howat, Applied Physics Lab, University of Washington
% ihowat@apl.washington.edu
% Version 1: 19-Jul-2007 00:50:57
% Modified by Felix Totir

fid = fopen(hdrfile);
while true
    line = fgetl(fid);
    if line == -1
        break
    else
        eqsn = findstr(line,'=');
        if ~isempty(eqsn)
            param = strtrim(line(1:eqsn-1));
            param(findstr(param,' ')) = '_';
            value = strtrim(line(eqsn+1:end));
            if isempty(str2num(value))
                if ~isempty(findstr(value,'{')) && isempty(findstr(value,'}'))
                    while isempty(findstr(value,'}'))
                        line = fgetl(fid);
                        value = [value,strtrim(line)];
                    end
                end
                eval(['info.',param,' = ''',value,''';'])
            else
                eval(['info.',param,' = ',value,';'])
            end
        end
    end
end
fclose(fid);

if isfield(info,'map_info')
    line = info.map_info;
    line(line == '{' | line == '}') = [];
    
    %originally: line = strtrim(split(line,','));
    %replaced by
    line=textscan(line,'%s',','); %behavior is not quite the same if "line" ends in ','
    line=line{:};
    line=strtrim(line);
    %
    
    info.map_info = [];
    info.map_info.projection = line{1};
    info.map_info.image_coords = [str2num(line{2}),str2num(line{3})];
    info.map_info.mapx = str2num(line{4});
    info.map_info.mapy = str2num(line{5});
    info.map_info.dx  = str2num(line{6});
    info.map_info.dy  = str2num(line{7});
    if length(line) == 9
        info.map_info.datum  = line{8};
        info.map_info.units  = line{9}(7:end);
    elseif length(line) == 11
        info.map_info.zone  = str2num(line{8});
        info.map_info.hemi  = line{9};
        info.map_info.datum  = line{10};
        info.map_info.units  = line{11}(7:end);
    end
    
    %part below comes form the original enviread
    %% Make geo-location vectors
    xi = info.map_info.image_coords(1);
    yi = info.map_info.image_coords(2);
    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    %adjust points to corner (1.5,1.5)
    if yi > 1.5
        ym =  ym + ((yi*info.map_info.dy)-info.map_info.dy);
    end
    if xi > 1.5
        xm = xm - ((xi*info.map_info.dy)-info.map_info.dx);
    end
    
    info.x= xm + ((0:info.samples-1).*info.map_info.dx);
    info.y = ym - ((0:info.lines-1).*info.map_info.dy);
end

if isfield(info,'pixel_size')
    line = info.pixel_size;
    line(line == '{' | line == '}') = [];
    
    %originally: line = strtrim(split(line,','));
    %replaced by:
    line=textscan(line,'%s',','); %behavior is not quite the same if "line" ends in ','
    line=line{:};
    line=strtrim(line);
    
    info.pixel_size = [];
    info.pixel_size.x = str2num(line{1});
    info.pixel_size.y = str2num(line{2});
    info.pixel_size.units = line{3}(7:end);
end

% function split is only used when replacements above do not work
% function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
% A = {};
% while (~isempty(s))
%     [t,s] = strtok(s,d);
%     A = {A{:}, t};
% end

