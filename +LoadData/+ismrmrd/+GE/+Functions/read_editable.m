function val = read_editable(hdl)
% READ_EDITABLE  Read in double from editable text
% val = read_editable(hdl)
% 7/2006 Rolf Schulte
val = str2num(get(hdl,'String'));
if isnan(val)
  fprintf('Warning: Invalid input in: %s\n',get(hdl,'Tag'));
  name = get(get(hdl,'Parent'),'Name');
  switch name(1:3),
   case 'css', val = get_cssf_default(hdl);
   case 'jpr', val = get_jpress_default(hdl);
   otherwise,  val = 0;
  end
  
  fprintf('Setting value to default of %g\n',val);
  set(hdl,'String',num2str(val));
end
