function y = fiex3d(se0,fi_order,fi_nr,fi_omega,...
    ex_order,ex_nr,ex_omega,ex_rounds,verb)
%FIEX3D  Voxel-wise polynomial FItting and EXtrapolation of 2D & 3D maps
%  1. fits sensitivity in a smaller area within the mask
%  2. extrapolates exceeding mask
%
%y = fiex3d(se0,fi_order,fi_nr,fi_omega,ex_order,ex_nr,ex_omega,ex_rounds)
%       se0   masked raw sensitivity [nx,ny,nz]
% optional parameters:                                          defaults
% fitting parameters (inside masked area)
%  fi_order   order of fit (1=linear; 2=quadratic; 3=cubic)     1
%     fi_nr   kernel width for fitting                          6
%  fi_omega   kernel width for Gaussian weighting               4
% extrapolation parameters (around masked area)
%  ex_order   order                                             1
%     ex_nr   kernel width                                      4
%  ex_omega   Gaussian factor                                   4
% ex_rounds   cycles around (0=off)                             2
%      verb   verbose mode                                      true
%
%         y   fitted and extrapolated maps
%
% Literature: Pruessmann, SENSE: MRM1999;42:952-962.
%
% Created (2D) 10.11,2 | Florian Wiesinger
% 10/2019 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% parameters + input checks
mtx = size(se0);            % matrix size
dim = length(mtx);          % dimensions of input image
if ((dim==2) && (mtx(2)==1)), dim = 1; end
if ((dim>1) && any(mtx==1)), warning('singleton dimension'); end
if dim>3, error('dimensions of se0 (=%g) > 3',dim); end
if any(isnan(se0(:)))
    warning('se0 contains NaN; setting to zero');
    se0(isnan(se0)) = 0;
end
if any(isinf(se0(:)))
    warning('se0 contains Inf; setting to zero');
    se0(isinf(se0)) = 0;
end


%% defaults for input parameters
% FITTING PARAMETERS (orig 1 8 6)
% order of fit (1=linear; 2=quadratic; 3=cubic)
if ~exist('fi_order','var'), fi_order = []; end
if isempty(fi_order),        fi_order = 1;  end
% kernel width for fitting
if ~exist('fi_nr','var'),    fi_nr = []; end
if isempty(fi_nr),           fi_nr = 6;  end
% kernel width for Gaussian weighting
if ~exist('fi_omega','var'), fi_omega = []; end
if isempty(fi_omega),        fi_omega = 4; end

% EXTRAPOLATION PARAMETERS (1 5 5 1)
if ~exist('ex_order','var'), ex_order = []; end  % order
if isempty(ex_order),        ex_order = 1;  end
if ~exist('ex_nr','var'),    ex_nr = []; end     % kernel width
if isempty(ex_nr),           ex_nr = 4;  end
if ~exist('ex_omega','var'), ex_omega = []; end  % Gaussian factor
if isempty(ex_omega),        ex_omega = 4; end
if ~exist('ex_rounds','var'),ex_rounds = []; end % cycles around
if isempty(ex_rounds),       ex_rounds = 2; end

if ~exist('verb','var'),     verb = []; end      % verbose mode
if isempty(verb),            verb = true; end


%% print info
if verb
    fprintf('fi_order=%g, fi_nr=%g, fi_omega=%g | ',fi_order,fi_nr,fi_omega);
    fprintf('ex_order=%g, ex_nr=%g, ex_omega=%g, ex_rounds=%g\n',...
        ex_order,ex_nr,ex_omega,ex_rounds);
end


%% extend matrix for fitting
n_ext = max([fi_nr ex_nr]);       % extend matrix by 2*n_ext
mtx_ext = mtx+2*n_ext;            % extended matrix size
if dim==1, mtx_ext(2) = 1; end
mtx_ext = [mtx_ext ones(1,3-dim)];
ie{2} = 1; ie{3} = 1;             % array for looping
for l=1:dim, ie{l} = n_ext+(1:mtx(l)); end
se1  = zeros(mtx_ext);            % extended sensitivity map
se1(ie{1},ie{2},ie{3}) = se0;


%% 1) fit sub-area
if verb, fprintf('Polynomial smoothing inside\n'); end
se2 = zeros(mtx_ext);
if ((fi_nr>0)&&(fi_order>0))
    [X,W] = sub_construct_polynomial_coefficients(dim,fi_nr,fi_order,fi_omega);
    
    fi_arr{2} = 0; fi_arr{3} = 0;     % index array
    for l=1:dim, fi_arr{l} = (-fi_nr:fi_nr); end
    
    warning('off');
    nwarn = 0;
    for l1=ie{1}
        for l2=ie{2}
            for l3=ie{3}
                if (se1(l1,l2,l3)~=0)
                    lastwarn('','');
                    % select and vectorise data
                    s = se1(l1+fi_arr{1},l2+fi_arr{2},l3+fi_arr{3});
                    s = s(:);         % vectorise
                    isi = (s~=0);     % exclude masked points
                    X1  = X(isi,:);   % polynomial coefficients
                    W1  = W(:,isi);   % Gaussian weighting
                    % if sum(isi)<100, fprintf('%g ',sum(isi)); end
                    c = (W1.*X1'*X1)\((W1.*X1')*s(isi));  % actual fit
                    % c = ((bsxfun(@times,X1',W1)*X1)\(bsxfun(@times,X1',W1)*s(isi)));
                    % c = pinv(W1.*X1'*X1)*(W1.*X1')*s(isi);
                    [msgstr,msgid]=lastwarn;
                    if ~isempty(msgstr)
                        nwarn = nwarn+1;
                        if ~strcmp(msgid,'MATLAB:singularMatrix')
                            fprintf('Warning: %s\n',msgstr);
                        end
                        se2(l1,l2,l3) = se1(l1,l2,l3);
                    else
                        se2(l1,l2,l3) = c(1);  % assign to smoothed map
                    end
                end
            end
        end
    end
    warning('on');
    if nwarn>1
        fprintf('%g warnings: %s\n',nwarn,msgstr);
        fprintf('\t-> voxels removed\n');
    end
else
    if verb, fprintf('Skipping polynomial fitting inside\n'); end
    se2 = se1;
end


%% 2) extrapolation of sensitivities
if verb, fprintf('Extrapolation of sensitivites\n'); end
se3 = se2;
if sum(abs(se1)>eps)==0
    if verb, fprintf('mask all zeroes -> skipping extrapolation\n'); end
    ex_rounds = 0;
end
if ((ex_nr>0)&&(ex_order>0)&&(ex_rounds>0))
    mask1 = (se1~=0);                 % mask
    [X,W] = sub_construct_polynomial_coefficients(dim,ex_nr,ex_order,ex_omega);
    
    ex_arr{2} = 0; ex_arr{3} = 0;     % index array
    for l=1:dim, ex_arr{l} = (-ex_nr:ex_nr); end
    
    % creat spherical structuring element
    r = 1;                            % radius
    switch dim
        case 1, [x1] = ndgrid(-r:r); x2 = zeros(size(x1)); x3 = x2;
        case 2, [x1,x2] = ndgrid(-r:r); x3 = zeros(size(x1));
        case 3, [x1,x2,x3] = ndgrid(-r:r);
    end
    str_ele_sphere1 = strel(sqrt(x1.^2 + x2.^2 + x3.^2) <=r);

    warning('off');
    nwarn = 0;
    for lr = 1:ex_rounds
        % extend mask by one in each direction (requires image_toolbox)
        mask2 = imdilate(mask1,str_ele_sphere1);
        
        diff_mask = mask2 - mask1;    % mask for extended area
        mask1 = mask2;
        for l1=1:mtx_ext(1)
            for l2=1:mtx_ext(2)
                for l3=1:mtx_ext(3)
                    if diff_mask(l1,l2,l3)~=0
                        lastwarn('','');
                        % data selection
                        s = se1(l1+ex_arr{1},l2+ex_arr{2},l3+ex_arr{3});
                        s = s(:);     % vectorise
                        isi = (s~=0); % exclude zero values
                        X1 = X(isi,:);
                        W1 = W(:,isi);
                        c = (W1.*X1'*X1)\((W1.*X1')*s(isi));   % actual LSQ fit
                        % c = pinv(W1.*X1'*X1)*(W1.*X1')*s(isi);
                        [msgstr,msgid]=lastwarn;
                        if ~isempty(msgstr)
                            nwarn = nwarn+1;
                            % if ~strcmp(msgid,'MATLAB:singularMatrix')
                            %     fprintf('Warning: %s\n',msgstr);
                            % end
                            se3(l1,l2,l3) = se2(l1,l2,l3);
                        else
                            se3(l1,l2,l3)=c(1);  % assign to extrapolatied map
                        end
                    end
                end
            end
        end
    end
    warning('on');
    if nwarn>1
        fprintf('%g warnings: %s\n',nwarn,msgstr);
        fprintf('\t-> voxels removed\n');
    end
else
    if verb, fprintf('Skipping polynomial extrapolation outside\n'); end
end


%% return map
y = se3(ie{1},ie{2},ie{3});       % crop smoothed map to original size

if any(isnan(y(:)))
    warning('y contains NaN (%g); setting to zero',sum(isnan(y(:))));
    y(isnan(y)) = 0;
end

if any(isinf(y(:)))
    warning('y contains Inf (%g); setting to zero',sum(isinf(y(:))));
    y(isinf(y)) = 0;
end


end      % main funciton fiex3d


%% subfunction
function [X,W] = sub_construct_polynomial_coefficients(dim,nr,order,omega)
%SUB_CONSTRUCT_POLYNOMIAL_COEFFICIENTS  
% Generate polynomial coefficients for smoothing sensitivity maps 
% and corresponding Gaussian weighting function
if order>3, error('polynomial order (=%g)>3',order); end
if order<1, error('polynomial order (=%g)<3',order); end
if dim>3, error('dim (=%g)>3',dim); end

% generate grid
switch dim
    case 1
        [x1] = ndgrid(-nr:nr);
        x1=x1(:); x2 = zeros(size(x1)); x3 = x2;
    case 2
        [x1,x2] = ndgrid(-nr:nr);
        x1=x1(:); x2=x2(:); x3 = zeros(size(x1));
    case 3
        [x1,x2,x3] = ndgrid(-nr:nr);
        x1=x1(:); x2=x2(:); x3=x3(:);
end

% polynomial fitting coefficients
X = [ones(size(x1)), x1];    % linear terms
if dim>1, X = [X, x2]; end
if dim>2, X = [X, x3]; end
if order>1                   % quadratic terms
    X = [X, x1.^2];
    if dim>1, X = [X, x2.^2, x1.*x2]; end
    if dim>2, X = [X, x3.^2, x1.*x3, x2.*x3]; end
end
if order>2                   % cubic terms
    X = [X, x1.^3];
    if dim>1, X = [X, x2.^3, x1.^2.*x2, x1.*x2.^2]; end
    if dim>2, X = [X, x3.^3, x1.^2.*x3, x1.*x3.^2,...
            x2.^2.*x3, x2.*x3.^2, x1.*x2.*x3]; end
end

% weighting function (Gaussian)
W = meshgrid(exp(-(x1.^2+x2.^2+x3.^2)/omega^2), 1:size(X,2) );

end      % sub_construct_polynomial_coefficients
