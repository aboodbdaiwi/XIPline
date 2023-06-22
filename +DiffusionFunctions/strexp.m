function sse = strexp(x,b,ydata)
So = x(1);
DDC = x(2);
alpha = x(3);
A=So.*exp(-(DDC.*b).^alpha);
sse = sum((ydata-A).^2);