function [ventP_stats,VentByLobes] = classify_ventilation_levles(ventilation_mask,lungsmask,lobemask)
% Calculate percent ventilation distribution
% W. Zha @ August 2016
if nargin<3
    lobemask = [];
    VentByLobes = [];
end
ventLevels = unique(ventilation_mask);
ventLevels(ventLevels<1)=[];
if length(ventLevels)<max(ventLevels)
    ventLevels = 1:max(ventLevels);
else
    ventLevels = ventLevels.';
end

lungVent  = ventilation_mask.*lungsmask;

lungVol = sum(lungVent(:)>0);
ventPixels = zeros(1,length(ventLevels));
for nv = ventLevels%1:length(ventLevels)
    level = sprintf('v%d',nv);
    
    ventPixels(nv) = sum(lungVent(:)== nv);
end
if max(ventLevels)==6
    nVentLevels = 6;
    ventNames = {'VDP','vent1','vent2','vent3','vent4','vent5'};
else
nVentLevels = 4;
ventNames = {'VDP','LVP','MVP','HVP'};
end
% nVentLevles = 2;
% ventNames = {'LVP','HVP'};


ventP_stats = zeros(1,nVentLevels);
% if max(ventLevels) ==6
% ventP_stats(1) = (ventPixels(2)+ventPixels(3))/lungVol*100; % low ventilated
% ventP_stats(2) = ventPixels(4)/lungVol*100; % normal ventilated
% ventP_stats(3) = (ventPixels(5)+ventPixels(6))/lungVol*100; %highly ventilated
% else
    for nv = 1:nVentLevels
        ventP_stats(nv) = ventPixels(nv)/lungVol*100;
    end
    
% end
disp(ventNames);
disp(ventP_stats);
%%
if ~isempty(lobemask)
lobevalues = unique(lobemask);
lobevalues(lobevalues<2 | lobevalues>6) = [];
nLobes = length(lobevalues);
lobemask = double(lobemask).*lungsmask;

if nLobes == 2% left and right
    lobelabels = {'RUL','LUL'};
    lobevalues = [2,5];
elseif nLobes == 5
    lobelabels = {'RUL','RML','RLL','LUL','LLL'};      
        lobevalues = 2:6;

end
[nRows,nCols,nSlices] = size(lungVent);
for nv = 1: nVentLevels
    binary_vent = zeros(nRows,nCols,nSlices);
    ventname = ventNames{nv};
    binary_vent(lungVent == nv) = 1;
    
    lobeVent = binary_vent.*lobemask;
    VentByLobes.(ventname) = zeros(1,nLobes);
    for n_lobe = 1:nLobes
        this_lobe = lobemask;
        this_lobe( this_lobe ~= lobevalues(n_lobe)) = 0;
        per_lobe_vent = binary_vent .* this_lobe;
        VentByLobes.(ventname)(n_lobe) = sum(per_lobe_vent(:) == lobevalues(n_lobe))*100/sum(this_lobe(:)==lobevalues(n_lobe));
    end
    
end
end