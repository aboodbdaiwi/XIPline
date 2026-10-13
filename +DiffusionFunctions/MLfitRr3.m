function [x,fval] = MLfitRr3(datacol,b,s2,A,mymodel,ws,lb,ub,D0,deltaT,gas)
% 129Xe: exact original MLfitRr2 code and optimization.
% 3He: additional separate model branch.
if nargin < 11 || isempty(gas), gas = '129Xe'; end
if any(strcmpi(strtrim(char(gas)),{'3He','he','helium'}))
    [x,fval] = fitHelium(datacol,b,s2,A,mymodel,ws,lb,ub,D0,deltaT);
    return;
end

% Fits model 'mymodel" to data 'datacol' using maximum likelihood method
% for Rician distributed noisy data
% Juan Parra-Robles, CPIR-CCHMC, 09/17/2023
% Edited by Abood Bdaiwi, CPIR-CCHMC, 011/10/2023

% ws are weights used in weighted ML or with repeated measures
% D0 free diff coeff
% deltaT diffusion time 
% s2 noise variance
% A initial values
% b b-values
% datacol data vector

if nargin<5
    ws=1;
    mymodel="expn";
elseif nargin<6
    ws=1;
    if mymodel=="cylRandr"
        D0=0.14;
        deltaT=3.5e-3;
        lb=[0.0200 0.0020]; ub=[0.0440 0.0250];
    elseif mymodel=="cylRandrAnimal"
        D0=0.14;
        deltaT=2.5e-3;
        lb=[0.0080 0.0010]; ub=[0.0160 0.0100];        
    elseif mymodel=="strexp"
        lb=[0.00001 0.005 0.5]; ub=[datacol(1)+3*sqrt(s2) 0.12 1.1];
        %lb=[]; ub=[];
    end
elseif nargin<10
    D0=0.14;
    deltaT=3.5e-3;
end

if mymodel=="expn"
    fun=@expnfun;
    [x,fval] = fminsearch(fun,A);

elseif mymodel=="strexp"
    fun=@strexpfun;
    options = optimoptions(@fmincon,'Display','off');
    [x,fval] = fmincon(fun,A,[],[],[],[],lb,ub,[],options);
    %options=optimset('display','off');
    %[x,fval] = fminsearch(fun,A,options);

elseif mymodel=="cylRandr"
    options = optimoptions(@fmincon,'Display','off');
    fun=@cylRandHfun;
    ld=sqrt(2*D0*deltaT);  % Diffusion length
    ld2=sqrt(4*D0*deltaT);
    b=max(b,1e-5);   % to avoid error when evaluating at b==0

    Aineq=[0.38,-1,0];                                       %  using r/R >= 0.38 instead of r/R >= 0.3   
    bineq=0;                                                % can change both for [] when calling fmincon to remove condition

    [x,fval] = fmincon(fun,A,Aineq,bineq,[],[],lb,ub,[],options);

elseif mymodel=="cylRandrAnimal"
    options = optimoptions(@fmincon,'Display','off');
    fun=@cylRandHfunAnimal;
    ld=sqrt(2*D0*deltaT);  % Diffusion length
    ld2=sqrt(4*D0*deltaT);
    b=max(b,1e-5);   % to avoid error when evaluating at b==0

    Aineq=[0.38,-1,0];                                       %  using r/R >= 0.38 instead of r/R >= 0.3   
    bineq=0;                                                % can change both for [] when calling fmincon to remove condition

    [x,fval] = fmincon(fun,A,Aineq,bineq,[],[],lb,ub,[],options);
end


% Nested functions that compute the objective function for minimization

function logml = expnfun(A)                   % monoexponential decay, fits S0 and ADC
    fun2fit=A(1)*exp(-(b.*A(2))); % s0*exp(-b*adc) 
    logbessel1=log(besseli(0,fun2fit.*datacol/s2));
    logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
end

function logml = strexpfun(A)                 % Stretched exponential funtion, fits S0, DDC and alpha parameters
    fun2fit=A(1)*exp(-(b.*A(2)).^A(3)); % s0*exp(-b*ddc)^alpha
    logbessel1=log(besseli(0,fun2fit.*datacol/s2));
    logbessel2=(fun2fit.*datacol/s2)-0.5*log(fun2fit.*datacol/s2);         %  Alternative calculation to use when besseli function returns infinite numbers
    logbessel1(isinf(logbessel1))=logbessel2(isinf(logbessel1));  
    logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
end

function logml = cylRandHfun(A)               % morphometry model, fits S0, R and H parameters
    R=A(1);
    r=A(2);
    So=A(3);
    
    DLo=D0*exp(-2.81*(1-r./R)^1.76);
    betaL=(21.2*(R./ld)^1.5)*exp(-3.65*(1-r/R)^(-0.5));

    c1=1.13*(R/ld2)-1.40*(R/ld2)^2;
    c2=3.55-11.27*(R/ld2)+7.44*(R/ld2)^2;
    u=c1*(1-r/R)+c2*(1-r/R)^2;
    DTo=D0*exp(-0.74*(ld2/R)^1.47)*(1+u);

    ct0=1.89-5.82*(R/ld2)+3.6*(R/ld2)^2;
    ct1=-4.21+14.1*(R/ld2)-9.52*(R/ld2)^2;
    ct2=2.14-7.35*(R/ld2)+5.03*(R/ld2)^2;
    betaT=ct0+ct1*(r/R)+ct2*(r/R)^2;

    DT=DTo*(1+betaT.*b*DTo);
    DL=DLo*(1-betaL.*b*DLo);
    Dan=DL-DT;
    if Dan>0 & isreal(Dan)
        % This may not be needed if additional conditions added to ensure Dan>0 and real
        fun2fit=So*exp(-b.*DT).*(pi./(4.*b.*Dan)).^(1/2).*erf((b.*Dan).^(1/2));

        logbessel1=log(besseli(0,fun2fit.*datacol/s2));
        logbessel2=(fun2fit.*datacol/s2)-0.5*log(fun2fit.*datacol/s2);
        logbessel1(isinf(logbessel1))=logbessel2(isinf(logbessel1));
        logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
    else
        %disp('Dan<=0')
        logml= 0;
    end        
end

function logml = cylRandHfunAnimal(A)               % morphometry model, fits S0, R and H parameters
    R=A(1);
    r=A(2);
    So=A(3);

    DLo= D0*exp(-2.89*((1-(r/R))^1.78));
    betaL=(35.6*(R./ld)^1.5)*exp(-4/sqrt(1-r/R));
    A = 1.3+0.25*exp(14*(R/ld)^2);
    u = exp(-A*(1-r/R)^2) *(exp(-5*(1-r/R)^2) +5*((1-r/R)^2)-1);
    DTo=D0*exp(-0.73*(ld2/R)^1.4)*(1+u);
    betaT = 0.06;

    DT=DTo*(1+betaT.*b*DTo);
    DL=DLo*(1-betaL.*b*DLo);
    Dan=DL-DT;

    if Dan>0 & isreal(Dan)                % This may not be needed if additional conditions added to ensure Dan>0 and real
        fun2fit=So*exp(-b.*DT).*(pi./(4.*b.*Dan)).^(1/2).*erf((b.*Dan).^(1/2));

        logbessel1=log(besseli(0,fun2fit.*datacol/s2));
        logbessel2=(fun2fit.*datacol/s2)-0.5*log(fun2fit.*datacol/s2);
        logbessel1(isinf(logbessel1))=logbessel2(isinf(logbessel1));
        logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
    else
        logml= 0;
    end
end

end

function [x,fval] = fitHelium(datacol,b,s2,A,mymodel,ws,lb,ub,D0,deltaT)
% 3He-only extension. Does not affect the original 129Xe path.
b=b(:); datacol=datacol(:); A=A(:).';
if numel(b)~=numel(datacol), error('b/data length mismatch'); end
if strcmp(char(mymodel),'cylRandrAnimal')
    mouse=true;
else
    mouse=false;
end
Aineq=[0.4,-1,0;-1,1,0]; bineq=[0;0];
opts=optimoptions('fmincon','Display','off');
fun=@(p) heliumObjective(p,mouse);
[x,fval]=fmincon(fun,A,Aineq,bineq,[],[],lb,ub,[],opts);
    function value = heliumObjective(p,mouse)
        R=p(1); r=p(2); S0=p(3);
        if R<=0 || r<=0 || r>=R || S0<=0
            value=1e100; return;
        end
        v=1-r/R;
        L1=sqrt(2*D0*deltaT);
        L2=sqrt(4*D0*deltaT);
        xR=R/L2;
        if mouse
            % Osmanagic et al., J Appl Physiol 2010, Appendix A1-A2.
            % Validated for mouse/rat geometry and Delta~0.3-0.6 ms.
            DL0=D0*exp(-2.99*v^1.88);
            betaL=113.5*(R/L1)^1.5*exp(-5.11/sqrt(v));
            F=(7/16)*xR^(4-3.63*xR^0.6);
            AA=1.1+3.14*xR-18.39*xR^2;
            u=v^1.75*exp(AA*v);
            DT0=D0*F*(1+u);
            betaT=0;
        else
            % Sukstanskii & Yablonskiy, J Magn Reson 2008, Eqs 3-4.
            L = 2*R*sin(pi/8);  % Standard 8-alveoli model
            DL0=D0*(1-sqrt(R/L)*(1-exp(-2.5*v^1.8)));
            betaL=43.5*(R/L1)^2*exp(-4/sqrt(v));
            DT0=D0*(7/16)*xR^(4*(1-xR^0.7));
            betaT=0.09*xR^(0.1*(L2/R)^1.5);
        end
        DL=DL0*(1-betaL*b*DL0);
        DT=DT0*(1+betaT*b*DT0);
        if any(~isfinite(DL)) || any(~isfinite(DT)) || ...
                any(DL<=0) || any(DT<0) || any(DL<DT) || ...
                any(DL>D0*1.01) || any(DT>D0*1.01)
            value=1e100; return;
        end
        value=objective(cylinderSignal(S0,DL,DT));
    end

    function signal=cylinderSignal(S0,DL,DT)
        % Orientation-averaged acinar-cylinder signal (both gases).
        dan=DL-DT;
        z=b.*dan;
        factor=ones(size(z));
        mask=z>1e-8;
        factor(mask)=sqrt(pi)./(2*sqrt(z(mask))).*erf(sqrt(z(mask)));
        factor(~mask)=1-z(~mask)/3+z(~mask).^2/10;
        signal=S0*exp(-b.*DT).*factor;
    end

    function value=objective(signal)
        if any(~isfinite(signal)) || any(signal<0)
            value=1e100; return;
        end
        z=signal.*datacol/s2;
        % Scaled modified Bessel function prevents overflow.
        logI0=log(besseli(0,z,1))+z;
        terms=logI0-signal.^2/(2*s2);
        value=-sum(ws.*terms);
        if ~isfinite(value), value=1e100; end
    end
end
