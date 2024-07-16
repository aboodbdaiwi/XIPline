function long_uid = uncompress_uid(short_uid,minlen,verb)
%UNCOMPRESS_UID Uncompress dicom uid
% compressed UID shares 8bits in unsigned char
% first 4bits->odd uncompressed UID
% last 4bits->even uncompressed UID
% unsigned char not existing in Matlab, hence using uint16 and doing bit
% conversion manually
%  long_uid = uncompress_uid(short_uid,verb,minlen)
% short_uid   Compressed UID (as storied in p-file header)
%    minlen   Minimum length of UID; if below, generate dicom UID   ('')
%      verb   Verbose mode                                          (false)
%  long_uid   Uncompressed UID (as required for dicom header)
%
% 7/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~exist('minlen','var'), minlen = []; end
if isempty(minlen),        minlen = 8; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 0; end


%% uncompress UID
long_uid = '';
translation_used = false;
for l=1:(2*length(short_uid))
    sval1 = uint16(short_uid(ceil(l/2)));
    if sval1>255
        translation_used = true;
        if verb>2 
            fprintf('l=%d %d (>255)\n',l,sval1);
        end
        switch sval1         % correct for unsigned char
            case 352,  sval1 = 138;
            case 353,  sval1 = 154;
            case 402,  sval1 = 131;
            case 710,  sval1 = 136;
            case 732,  sval1 = 152;
            case 8211, sval1 = 150;
            case 8212, sval1 = 151;    
            case 8216, sval1 = 145;
            case 8217, sval1 = 146;
            case 8218, sval1 = 130;
            case 8220, sval1 = 147;
            case 8221, sval1 = 148;
            case 8222, sval1 = 132;
            case 8224, sval1 = 134;
            case 8225, sval1 = 135;
            case 8226, sval1 = 149;
            case 8230, sval1 = 133;
            case 8240, sval1 = 137;
            case 8249, sval1 = 178;
            case 8250, sval1 = 155;
            case 8364, sval1 = 128;
            case 8482, sval1 = 153;
            otherwise
                warning('sval1=%d not in correction list; setting to 17 (->''00'')',sval1);
                sval1 = 17;
                % get correct number via bin2dec([dec2bin(n1,4),dec2bin(n2,4)])
                % with n1/n2 being the numbers according to the conversion table
        end
        %         in = [352,353,402,710,732,8211,8216,8217,8218,...
        %             8220,8224,8225,8226,8230,8240,8250,8364,8482];
        %         out = [138;154;131;136;152;150;145;146;...
        %             130;147;134;135;149;133;137;155;128;153].';
    end
    sval = dec2bin(sval1,8); % convert to 0-1 bits
    
    if verb>2, fprintf('l=%d  \t%d %s ',l,sval1,sval); end
    if isodd(l)
        val1 = sval(1:4);    % odd->using first four bits
    else
        val1 = sval(5:8);    % even->using last four bits
    end
    val = bin2dec(val1);     % convert back to decimal number
    
    switch val
        case 0
            if l~=2*length(short_uid)
                if verb>2, fprintf('\n'); end
                % warning('break prematurely: %d/%d',l,2*length(short_uid));
            end
            break;
        case {1,2,3,4,5,6,7,8,9,10}
            lval = num2str(val-1);
        case 11, lval = '.';
        otherwise
            warning('val=%d; setting lval=''0''',val); 
            lval = '0';
    end
    if verb>2, fprintf('  %s %g %s\n',val1,val,lval); end
    long_uid(l) = lval;
end
if translation_used
    fprintf('Attention: sval1>255; switch character set before reading header via\n');
    fprintf('\tfeature(''DefaultCharacterSet'',''ISO-8859-1'')\n');
end

%% check UIDs
if minlen>0
    if length(long_uid)<minlen
        warning('length(long_uid)(=%d)<minlen(=%d)',length(long_uid),minlen);
        fprintf('short_uid = %s\n',short_uid);
        fprintf('long_uid = %s\n',long_uid);
        fprintf('Generating dicom UID\n');
        rootuid = '1.2.840.113619.2.156.';
        tmp = dicomuid;
        long_uid = [rootuid tmp(24:length(tmp))];
    end
end
if length(long_uid)>7
    if isempty(strfind(long_uid,'1.2.840.'))
        warning('Uncompressed UID (=''%s'') not starting with ''1.2.840.''',...
            long_uid);
    end
else
    warning('length(long_uid(=''%s''))(=%d)<8',long_uid,length(long_uid));
end



end           % sub_uncompress_uid
