function [mu,mask]=kmeans3D(image3D,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   kmeans image segmentation
%
%   Input:
%          image3D: grey color 3D image
%          k: Number of classes
%   Output:
%          mu: vector of class means 
%          mask: clasification image mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check image
image3D=double(image3D);
imageVector=image3D(:);       % vectorize image3D
mi=min(imageVector);      % deal with negative 
imageVector=imageVector-mi+1;     % and zero values
s=length(imageVector);
% create image histogram
m=max(imageVector)+1;
h=zeros(1,m);
hc=zeros(1,m);
for i=1:s
  if(imageVector(i)>0) h(imageVector(i))=h(imageVector(i))+1; end
end
ind=find(h);
hl=length(ind);
% initiate centroids
mu=(1:k)*m/(k+1);
% start process
iter = 1; 
while(true)  
      oldmu=mu;
      % current classification 
      for i=1:hl
          c=abs(ind(i)-mu);
          cc=find(c==min(c));
          hc(ind(i))=cc(1);
      end
      %recalculation of means  
      for i=1:k
          a=find(hc==i);
          mu(i)=sum(a.*h(a))/sum(h(a));
      end
      if(mu==oldmu) 
          break;
      end
    if iter > 1000
      break;
    end
    iter = iter + 1;
end
% calculate mask
s=size(image3D);
mask=zeros(s);
for i=1:s(1)
    for j=1:s(2)
        for slice=1:s(3)
          c=abs(image3D(i,j,slice)-mu);
          a=find(c==min(c));  
          mask(i,j,slice)=a(1);
        end
    end
end
mu=mu+mi-1;   % recover real range

