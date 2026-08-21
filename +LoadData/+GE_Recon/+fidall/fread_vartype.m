function a = fread_vartype( fid )
  a.length =  fread( fid, 1, 'uint32');    % length of the data
  a.pointer = fread( fid, 1, 'uint32');    % pointer to the data -- is uint32 correct?