function [x,fval] =  MLricianfitcon(datacol,b,s2,A,mymodel,ws,D0,deltaT)
% Fits model 'mymodel" to data 'datacol' using maximum likelihood method
% for Rician distributed noisy data
% Juan Parra-Robles, CPIR-CCHMC, 03/17/2023

% ws are weights used in weighted ML or with repeated measures
% D0 free diff coeff
% deltaT diffusion time
% s2 noise variance
% A initial values
% b b-values
% datacol data vector

% datacol = datacol./max(datacol(:));
if nargin<5
    ws=1;
    mymodel="expn";
elseif nargin<6
    ws=1;
    if mymodel=="cylRandH"
        D0=0.14;
        deltaT=5e-3;
    end
end

if mymodel=="expn"
    fun=@expnfun;
    [x,fval] = fminsearch(fun,A);

elseif mymodel=="strexp"
    fun=@strexpfun;
    [x,fval] = fminsearch(fun,A);

%     modeloptions = fitoptions('Method', 'NonlinearLeastSquares',...
%                'Lower', [0, -Inf -Inf],...
%                'Upper', [2, Inf Inf],...
%                'StartPoint', [0 -1 0]);
%     fun = fittype('a*exp((x*b)^c)', 'coeff', {'a', 'b', 'c'}, 'options', modeloptions);
% %     fun = fittype('a*exp(-(x*b)^c)+d', 'coeff', {'a', 'b', 'c','d'}, 'options', modeloptions);
%     f= fit(b, datacol, fun);
%     x = coeffvalues(f);
%     fval = 0;



elseif mymodel=="cylRandH"
    options = optimoptions(@fmincon,'Display','off');
    fun=@cylRandHfun;
    ld=sqrt(2*D0*deltaT);  % Diffusion length
    ld2=sqrt(2)*ld;

    Aineq=[0,-0.7,1];                                       % this is equivalent to r/R <= 0.3   ----->    -0.7 R + H <= 0
    bineq=0;                                                % can change both for [] when calling fmincon to remove condition

    %lb=[0 0 0]; ub=[datacol(1).*2+3*sqrt(s2) 0.1 0.04];    % lower and upper bounds for estimated parameter
    lb=[0  0  0]; 
    ub=[2  0.0800 0.0300];    % lower and upper bounds for estimated parameter
    %A = [datacol(1), 0.030 0.015];
%     try
        [x,fval] = fmincon(fun,A,Aineq,bineq,[],[],lb,ub,[],options);
%     catch
%         x = [0 0 0];
%         disp('Objective function is undefined at initial point. Fmincon cannot continue.')
%     end

%     fun=@cylRandHfun;
%     fitoptions = optimoptions('lsqcurvefit','Display','off');
%     lb=[0  0.0100  0.0050]; 
%     ub=[2  0.0600 0.0300];    % lower and upper bounds for estimated paramete
%     x = lsqcurvefit(fun,A,b,datacol,lb,ub,fitoptions)


elseif mymodel=="cylRandHfunAnimal"
    options = optimoptions(@fmincon,'Display','off');
    fun=@cylRandHfun;
    ld=sqrt(2*D0*deltaT);  % Diffusion length
    ld2=sqrt(2)*ld;

    Aineq=[0,-0.7,1];                                       % this is equivalent to r/R <= 0.3   ----->    -0.7 R + H <= 0
    bineq=0;                                                % can change both for [] when calling fmincon to remove condition

    %lb=[0 0 0]; ub=[datacol(1).*2+3*sqrt(s2) 0.1 0.04];    % lower and upper bounds for estimated parameter
    lb=[0  0  0]; 
    ub=[2  0.0200 0.0100];    % lower and upper bounds for estimated parameter
    %A = [datacol(1), 0.030 0.015];
%     try
        [x,fval] = fmincon(fun,A,Aineq,bineq,[],[],lb,ub,[],options);
%     catch
%         x = [0 0 0];
%         disp('Objective function is undefined at initial point. Fmincon cannot continue.')
%     end

%     fun=@cylRandHfun;
%     fitoptions = optimoptions('lsqcurvefit','Display','off');
%     lb=[0  0.0100  0.0050]; 
%     ub=[2  0.0600 0.0300];    % lower and upper bounds for estimated paramete
%     x = lsqcurvefit(fun,A,b,datacol,lb,ub,fitoptions)


end


% Nested functions that compute the objective function

    function logml = expnfun(A)                   % monoexponential decay, fits S0 and ADC
        
        fun2fit=A(1)*exp(-(b.*A(2))); % s0*exp(-b*adc) 
        logbessel1=log(besseli(0,fun2fit.*datacol/s2));
        logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
    end
    
    function logml = strexpfun(A)                 % Stretched exponential funtion, fits S0, DDC and alpha parameters

        fun2fit=A(1)*exp(-(b.*A(2)).^A(3)); % s0*exp(-b*ddc)^alpha
        logbessel1=log(besseli(0,fun2fit.*datacol/s2));
        logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
    end

    function logml = cylRandHfun(A)               % morphometry model, fits S0, R and H parameters
        R=A(2);
        r=R-A(3);

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

        if Dan>0 & isreal(Dan)                % This may not be needed if additional conditions added to ensure Dan>0 and real
            fun2fit=A(1)*exp(-b.*DT).*(pi./(4.*b.*Dan)).^(1/2).*erf((b.*Dan).^(1/2));
            if b(1)==0
                fun2fit(b==0)=A(1);       % avoid division by zero for b=0;
            end
            logbessel1=log(besseli(0,fun2fit.*datacol/s2));
            logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
        else
            logml= inf;
        end
    end

    function logml = cylRandHfunAnimal(A)               % morphometry model, fits S0, R and H parameters
        R=A(2);
        r=R-A(3);

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
            fun2fit=A(1)*exp(-b.*DT).*(pi./(4.*b.*Dan)).^(1/2).*erf((b.*Dan).^(1/2));
            if b(1)==0
                fun2fit(b==0)=A(1);       % avoid division by zero for b=0;
            end
            logbessel1=log(besseli(0,fun2fit.*datacol/s2));
            logml=-(sum(ws.*(logbessel1-(fun2fit.^2)/(2*s2))));
        else
            logml= inf;
        end
    end







end


