function sse = cylM(x,b,ydata,Do,deltaT)
So = x(1);
R = x(2);
h = x(3);

r=R-h;
ld=sqrt(2*Do*deltaT);  % Diffusion length
ld2=sqrt(2)*ld;

DLo=Do*exp(-2.81*(1-r./R)^1.76);
betaL=(21.2*(R./ld)^1.5)*exp(-3.65*(1-r/R)^(-0.5));

c1=1.13*(R/ld2)-1.40*(R/ld2)^2;
c2=3.55-11.27*(R/ld2)+7.44*(R/ld2)^2;
u=c1*(1-r/R)+c2*(1-r/R)^2;
DTo=Do*exp(-0.74*(ld2/R)^1.47)*(1+u);

ct0=1.89-5.82*(R/ld2)+3.6*(R/ld2)^2;
ct1=-4.21+14.1*(R/ld2)-9.52*(R/ld2)^2;
ct2=2.14-7.35*(R/ld2)+5.03*(R/ld2)^2;
betaT=ct0+ct1*(r/R)+ct2*(r/R)^2;

DT=DTo*(1+betaT.*b*DTo);
DL=DLo*(1-betaL.*b*DLo);
Dan=DL-DT;

if Dan>0 & isreal(Dan)                % This may not be needed if additional conditions added to ensure Dan>0 and real
    A = So*exp(-b.*DT).*(pi./(4.*b.*Dan)).^(1/2).*erf((b.*Dan).^(1/2));
    sse = sum((ydata-A).^2);
else
    sse = 0;
end









