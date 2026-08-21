function wfn = get_wfname(hh,vrsn,dname,verb)
%GET_WFNAME Return ak-grad waveform names for preinstalled MNSRP trajectories
%   wfn = get_wfname(hh,vrsn,dname,verb)
%    hh   Header structure or directly cv3 value
%  vrsn   Waveform version                                              (3)
% dname   Directory name
%  verb   Verbose mode                                              (false)
%
% 11/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters + defaults
if ~exist('vrsn','var'), vrsn = []; end
if isempty(vrsn),        vrsn = 3; end
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = false; end
if ~exist('dname','var'), dname = []; end
if isempty(dname)
    if isdeployed
        dname = '/export/home/sdc/fidall/waveforms/';
    else
        switch vrsn
            case 1,     dname = 'C:\BoxSync\data\waveforms_old\';
            case {2,3}, dname = 'C:\BoxSync\data\waveforms\';
            otherwise,  dname = '';
        end
    end
end
if ~exist(dname,'dir'), warning('directory (=''%s'') not found',dname); end
if isstruct(hh)
    cv3 = hh.image.user3;
else
    cv3 = hh;
end
if verb, fprintf('get_wfname: cv3=%g; waveform version=%g\n',cv3,vrsn); end


%% waveform names
switch vrsn
    case 1
        ak{2}='spiral/spiral_13C_fov80_npix32_arms1_ksamp16_gmax22_smax73';
        ak{3}='spiral/spiral_13C_fov120_npix40_arms1_ksamp16_gmax23_smax77';
        ak{4}='spiral/spiral_13C_fov200_npix50_arms1_ksamp8_gmax23_smax76';
        ak{5}='spiral/spiral_13C_fov300_npix60_arms1_ksamp8_gmax23_smax77';
        ak{6}='spiral/spiral_13C_fov400_npix66_arms1_ksamp8_gmax23_smax74';
        ak{7}='spiral/spiral_13C_fov51_npix34_arms1_ksamp16_gmax39_smax135';
        ak{8}='spiral/spiral_13C_fov90_npix45_arms1_ksamp8_gmax40_smax140';
        ak{9}='spiral/spiral_13C_fov135_npix54_arms1_ksamp8_gmax40_smax139';
        ak{10}='spiral/spiral_13C_fov180_npix60_arms1_ksamp4_gmax38_smax136';
        ak{11}='spiral/spiral_13C_fov300_npix75_arms1_ksamp4_gmax40_smax145';
        ak{12}='spiral/spiral_13C_fov400_npix80_arms1_ksamp4_gmax37_smax134';
        ak{13}='spiral/spiral_1H_fov90_npix60_arms1_ksamp4_gmax23_smax77';
        ak{14}='spiral/spiral_1H_fov80_npix80_arms1_ksamp4_gmax40_smax150';
        ak{15}='spiral/spiral_1H_fov120_npix80_arms1_ksamp4_gmax40_smax150';
        ak{16}='spiral/spiral_1H_fov180_npix90_arms1_ksamp4_gmax40_smax150';
        ak{17}='spiral/spiral_1H_fov400_npix100_arms4_ksamp2_gmax33_smax120';
        ak{19}='radial/radial3D_23Na_fov130_mtx70_intlv15460_kdt18_gmax33_smax119_dur3_coda_rew2_1';
        ak{28}='mrsi/mrsi3D_2H_scheme4_fov100_mtx10_bw5000_npts512_nexc2504_gmax32_smax163';
    case {2,3}
        ak{2}='spiral/spiral_13C_fov80_mtx28_arms1_kdt24_gmax33_smax119_dur24p6';
        ak{3}='spiral/spiral_13C_fov120_mtx34_arms1_kdt16_gmax33_smax119_dur25p3';
        ak{4}='spiral/spiral_13C_fov200_mtx42_arms1_kdt8_gmax33_smax119_dur25p7';
        if vrsn==2
            ak{5}='spiral/old/spiral_13C_fov300_mtx48_arms1_kdt4_gmax30_smax119_dur25p5';
            ak{6}='spiral/old/spiral_13C_fov400_mtx52_arms1_kdt4_gmax27_smax119_dur24p9';            
        else
            ak{5}='spiral/spiral_13C_fov300_mtx48_arms1_kdt8_gmax30_smax119_dur25p5';
            ak{6}='spiral/spiral_13C_fov400_mtx52_arms1_kdt8_gmax27_smax119_dur24p9';
        end
        ak{7}='spiral/spiral_1H_fov200_mtx50_arms1_kdt2_gmax19_smax119_dur16p6';
        ak{8}='spiral/spiral_1H_fov200_mtx200_arms50_kdt2_gmax29_smax119_dur3p3';
        ak{9}='spiral/spiral_1H_fov400_mtx100_arms4_kdt2_gmax15_smax119_dur9';
        ak{10}='mrp/radial2D_1H_fov200_mtx300_intlv1000_kdt4_gmax25_smax100_dur4p5_fsca_rew1_GA7';
        ak{11}='mrp/radial4D_1H_fov200_53_mtx300_20_intlv40000_kdt4_gmax25_smax100_dur4p5_fsca_rew1_GA7_zPE3';
        ak{15}='cart/cart2D_129Xe_fov300x240_mtx100x80_nexc80_kdt20_gmax32_smax99_dur4p92';
        ak{16}='cart/cart3D_129Xe_fov300x240x144_mtx100x80x24_nexc1928_kdt20_gmax33_smax92_dur5p76_bal_co';
        ak{19}='radial/radial3D_23Na_fov159_mtx80_intlv20230_kdt24_gmax33_smax120_dur16p1_coda_rew2_1';
        if vrsn==2
            ak{17}='mrsi/old/mrsi3D_129Xe_scheme4_fov280_280_140_mtx14_14_7_bw20000_npts88_nexc1799_gmax18_smax114';
            ak{28}='mrsi/old/mrsi3D_2H_scheme4_fov200_mtx12_bw5000_npts380_nexc2405_gmax27_smax113';
        else
            ak{17}='mrsi/mrsi3D_129Xe_scheme2_3_fov270_240_150_mtx18_16_10_bw20000_npts88_nexc1586_gmax16_smax87_crush2';
            ak{28}='mrsi/mrsi3D_2H_scheme4_fov140_mtx10_bw5000_npts700_nexc1678_gmax31_smax116';
        end
    otherwise
        error('vrsn(=%g) unkown',vrsn);
end
ak{18}='radial/radial3D_1H_fov170_mtx200_intlv61712_kdt4_gmax32_smax120_dur1p3_coca_rew2_1';
ak{20}='spsp_profile/gz_fov200_mtx512_13c_gmax23_smax75';
ak{21}='spsp_profile/gz_fov200_mtx512_1h_gmax23_smax75';
ak{22}='epsi/flyback_13C_x_fov82_gmax50_smax200';
ak{23}='epsi/flyback_13C_y_fov82_gmax50_smax200';
ak{24}='epsi/flyback_13C_z_fov82_gmax50_smax200';
ak{25}='mrp/spiral225_225_MRP_16Gmax_80Slew_1H_377int';
ak{26}='mrp/spiral225_225_MRP_33Gmax_120Slew_1H_89int';
ak{27}='grad_noise/grad_noise_spiral_13C_fov200_ksamp2_gmax33_smax120';
ak{29}='mrp/3Dspiral_fov_256_npix_200_gmax_18_smax_60_nphi_880';


%% checks
if cv3>length(ak)
    warning('cv3(=%d) not in ak list (length(ak)=%d); return wfn=''''',cv3,length(ak));
    wfn = '';
    return;
end
if cv3<1
    if exist('fidspiral.grad','file')
        wfn = 'fidspiral.grad';
    else
        warning('cv3(=%d)<1; return wfn=''''',cv3);
        wfn = '';
    end
    return
end
if isempty(ak{cv3}), error('ak{%d} is empty',cv3); end


%% assemble waveform name
wfn = [dname ak{cv3} '.mat'];
if ispc, wfn = strrep(wfn,'/','\'); end


%% final check
if ~exist(wfn,'file'), warning('wfname (=''%s'') not found',wfn); end
if verb, fprintf('wfn = ''%s''\n',wfn); end


end      % get_wfname.m
