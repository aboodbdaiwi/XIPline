function [a, header] = read_MR_image(filename, width, height, depth)
%read_MR_image - read MR images extracted from Signa image database
%
% Images may be extracted by either the LISTSELECT tool, direct file copy,
% or the XIMG tool.
%
%   [A, HEADER] = read_MR_image(FILENAME);    returns the
%   MR image in FILENAME into A, and returns the headers into HEADER.
%
%   A = read_MR_image(FILENAME);    returns the MR image in FILENAME.
%   Only the pixel header, needed for image bitmap info, is read.  No
%   No header are returned, making exectuion much faster.
%
%   A = read_MR_image(FILENAME, WIDTH, HEIGHT, DEPTH );    reads the
%   MR image from FILENAME using the WIDTH, HEIGHT, and DEPTH parameters.
%   No headers are read, making for the fastest execution.

% Copyright (c) 2020 by General Electric Company. All rights reserved.

switch nargin
    case 1    % Read headers to get image info
        switch nargout
            case 2    % If header will be returned, read ALL the headers.
                header = read_MR_headers( filename, 'all','image' );
            case 1    % If not, read only the pixel header (needed for image bitmap info)
                header = read_MR_headers( filename, 'qimg','image' );
        end
        if ischar(header)
            if strcmp(header,'DICM')
                a = dicomread(filename);
                return
            end
        end
        width  = header.pixel.width;
        height = header.pixel.height;
        depth  = header.pixel.depth;
        offset = header.total_length;
        
    case 4      % Use input params for image info
        offset = 8432;		% Assume offset -- this works for either XIMG or LISTSELECT format
        
    otherwise
        fprintf('Incorrect number of input arguments -- type "help read_MR_image".\n');
        a = [];
        return;
end

% Read image data
fid = fopen(filename,'r',header.endian);
fseek(fid, offset, 'bof');
switch depth
    case 8
        a = fread(fid, [ width, height ], '256*int8');
    case 16
        a = fread(fid, [ width, height ], '256*int16');
    case 32
        a = fread(fid, [ width, height ], '256*int32');
    otherwise
        a = [];
end

a = transpose(a);		% Transpose since MATLAB reads data down columns, not accross rows.

fclose(fid);

%{
Modification History
Rev 1.7  2003-MAY-07  Matthew Eash
Rev 1.8  29 Sept 2003 B.Sivalingam
MRIhc37346 UJT Moved symbolic link location of gems-dicom-dict.txt to /usr/g/tools

2020 Nov 18 Jen Taylor
HCSDM00635322
Remove dependency on isHost for easier source code use by internal MR teams
%}
