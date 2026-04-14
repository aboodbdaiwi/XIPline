function info = detectAnalyzer()

    user = lower(strtrim(string(getenv('USERNAME'))));
    host = lower(strtrim(string(getenv('COMPUTERNAME'))));

    info = struct();
    info.code = "";
    info.name = "";
    info.role = "";                 % 'primary' or 'backup'
    info.primaryPolarizer = "";
    info.otherPolarizer = "";

    % Carter
    if user == "mcm5bk" || host == "ew19-03703"
        info.code = "CBM";
        info.name = "McMaster, Carter";
        info.primaryPolarizer = "McMaster, Carter";
        info.otherPolarizer = "McNeill, Mary";
        return
    end

    % 
    if user == "mcn4vc" 
        info.code = "MRM";
        info.name = "McNeill, Mary";
        info.primaryPolarizer = "McNeill, Mary";
        info.otherPolarizer = "McMaster, Carter";
        return
    end


    error("detectAnalyzer:UnknownIdentity", ...
        "Unknown user '%s' on host '%s'. Update detectAnalyzer.m", user, host);
end