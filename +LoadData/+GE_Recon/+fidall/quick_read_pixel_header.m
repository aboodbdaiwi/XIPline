function a = quick_read_pixel_header( fid )
%quick_read_pixel_header - read only the fields needed to locate raw data
%
% QUICK_READ_PIXEL_HEADER reads only the fields needed to locate
% the raw data in the rawfile, for decreased execution time.
%
%  function a = quick_read_pixel_header( fid )
%    fid - open file handle
%    a - structure with necessary fields

% Copyright (c) 2012 by General Electric Company. All rights reserved.

% Modification History
% MRIhc31374  11 Jan 2008 Puneet This function is taken from 12.0 build 
%   to read the images from ximg output

a.magic         = freadc( fid, 4 );
a.hdr_length    =  fread( fid, 1,  'int32' );
a.width         =  fread( fid, 1,  'int32' );
a.height        =  fread( fid, 1,  'int32' );
a.depth         =  fread( fid, 1,  'int32' );
