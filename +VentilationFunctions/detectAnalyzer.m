function analyzer = detectAnalyzer()

hostname = lower(strtrim(string(getenv('COMPUTERNAME'))));

switch hostname
    case "ew19-03703"
        analyzer = 'CBM';
    otherwise
        analyzer = 'MRM';
end

end