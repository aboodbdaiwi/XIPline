function sse = fmaxlikelihood_ADC(x,b,ydata,sigma)
So = x(1);
ADC = x(2);
A=So*exp(-ADC*b);
B=(sigma.^2)./(2*A);
sse = sum((ydata-A-B).^2);
