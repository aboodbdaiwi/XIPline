function matlab_scripter(varargin)
%MATLAB_SCRIPTER Matlab interface to run miscelaneous functions on scanner
% matlab_scripter(varargin)
%     1st input  function name
% further input  arguments for function
%                automatically converts double + characters
%  2/2021  Rolf Schulte

%% set startup defaults
feature('DefaultCharacterSet','ISO-8859-1');     % for short UID in header

%% interactive part: opens a prompt and runs compiled in matlab functions
interactive = false;
dummy_false = (10<0);               % for mcc V9.5
if nargin==0, interactive = true; 
else
    if strcmp(varargin{1},'interactive'), interactive = true; end
end


if interactive
    while (1)
        cmd = input('>> ','s');
        if (strcmp(cmd,'return') || strcmp(cmd,'break')), break; end
        if regexpi(cmd,'help ')==1, deployhelp(cmd(6:end)); 
        else
            try 
                eval(cmd);
                if regexp(lastwarn,'HELP','ONCE')
                    deployhelp(cmd);
                    lastwarn('');
                end
            catch
                err = lasterror;
                fprintf('Matlab scripter: cannot evaluate ''%s''\n',cmd);
                fprintf('Error: %s\n',err.message);
            end
        end
    end
    return
end

%% execute function passed as argument
cmdstr = varargin{1};
if nargin>1
    cmdstr = [cmdstr '('];
    for l=2:nargin
        tmp = varargin{l};
        if isempty(str2num(tmp))
            if (isempty(strfind(tmp,'[')) && isempty(strfind(tmp,']')))
                tmp = ['''' tmp ''''];
            end
        end
        cmdstr = [cmdstr tmp];
        if l<nargin, cmdstr = [cmdstr ',']; end
    end
    cmdstr = [cmdstr ')'];
end
cmdstr = [cmdstr ';'];
fprintf('Matlab command:\n>> %s\n',cmdstr);
try
    outp=evalc(cmdstr);
catch
    err = lasterror;
end

%% write function output to logfile
fid = 1;
logdir = getenv('MASCRILOGDIR');
if exist(logdir,'dir')
    fid = fopen([logdir '/matlab_scripter' ...
        datestr(now,'yyyymmdd_hhMMss') '.txt'],'w+');
else
    if exist([getenv('homedir') '/log/'],'dir')
        fid = fopen([getenv('homedir') '/log/matlab_scripter' ...
            datestr(now,'yyyymmdd_hhMMss') '.txt'],'w+');
    end
end
fprintf(fid,'************************************\n');
fprintf(fid,'Matlab command:\n>> %s\n',cmdstr);
fprintf(fid,'************************************\n');
if ~exist('err','var')
    outp = regexprep(outp,'%','%%');       % fix formatting of % sign
    fprintf('evalc output = \n');
    fprintf(outp);
    fprintf(fid,'evalc output = \n');
    fprintf(fid,outp);
else
    fprintf('Matlab scripter: cannot evaluate ''%s''\n',cmdstr);
    fprintf('Error:\n%s\n',err.message);
    fprintf(fid,'Matlab scripter: cannot evaluate ''%s''\n',cmdstr);
    fprintf(fid,'Error:\n%s\n',err.message);
end



%% name all desired functions here for inclusion into ctf archive
if dummy_false
    dicom2mat;eval_pol_mr;gyrogamma;hyper_flip;NSA;pfile2mat;
    fid2spec;read_p;mns_prescan;freqs13c;coil_characterise;
    recon_spiral;recon_spsp_profile;stack_plot;help;deployhelp;
    read_fidall_gating;design_spiral;recon_csi;image_overlay;
    memo;grep_hdr;dbm2mw;cart2b0map;SliceBrowser;
    Create5xxRfExtFile;grep_struct;
    read_ak_wav;plot_spec;pns;read_rho;Read5xxRfExtFile;
    design_radial;recon_grid3d;plot_csi;blosi_b1map;recon_zte;
    relative_b1map;segment_sphere;read_dicom;b1optimiser;
    circular_mask;dicom2png;fidall_b1shimming.m;zte_b1shimming;
    imagesc_ind3d;b1map_phase;GERecon;read_archive;zte_image_fid;
    read_archive_header;read_fdl;write_fdl;characterise_grad_noise;
    ernst_angle;recon_epi;recon_snr_map_2Dgr;design_cart;recon_cart;
    ssfp_phase_cycling;design_epi;read_plotter;design_rf4fidall;
    psf_mrsi;recon_mrs;recon_mrsi;design_mrsi;design_spiral3d;design_cones;
    design_hardpulse4fidall;mri_coil_combine;design_corings;recon_corings;
    design_spsp2d;plot_mrsi;plot_mrsi_file;design_spiral_pns;design_floret;
    ssfp_angle;mns_freq_converter;create_spsp_cv_pars;analyse_trajectory;
    golden_angle;vfa;ssfp_combine_phase_cycling;integrate_mrsi;
    recon_mrsi_xenon;mns_get_f0;design_qti.m;anonymise_mat;
end

if isdeployed      % when running on scanner    
    for l=1:90     % wait for 45 minutes before closing figures by force
        if ~any(ishandle(1:100)), break; end
        pause(30);
    end
    fprintf('Quitting matlab_scripter\n');
    fprintf(fid,'Quitting matlab_scripter\n');
    if fid~=1, fclose(fid); end
    quit;          % exit matlab to avoid matlab_crash_dump.XXX
else
    if fid~=1, fclose(fid); end
end

end   % main function matlab_scripter.m
