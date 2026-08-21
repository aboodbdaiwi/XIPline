% *******************************************************
% set slider to allowed position
% 7/2006 Rolf Schulte
function set_slider(hdl,val)
hmax = get(hdl,'Max');
hmin = get(hdl,'Min');
if length(val)>1, 
    fprintf('Warning: setting slider to first value\n');
    val = val(1);
end
if val>hmax,
  val = hmax;
end
if val<hmin
  val = hmin;
end
set(hdl,'Value',val);

