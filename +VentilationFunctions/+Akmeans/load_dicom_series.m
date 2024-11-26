function [image_stack,acquisition_order,image_header,IPP] = load_dicom_series(dicom_path,common_string)
%
% W. Zha@ 2015
if nargin<2
    common_string = [];
end
datasets = dicom_folder_info(dicom_path);
folder_content = datasets.Filenames;
if ~isempty(folder_content)
    nSlices = numel(folder_content);
    image_header = dicominfo(folder_content{end});
    
   % ADJUST THIS BACK FOR UW DATA
    %nRows=image_header.Width;
    %nCols = image_header.Height;
    nCols=image_header.Width;
    nRows = image_header.Height;
    image_stack = zeros(nRows,nCols,nSlices);
    IPP=zeros(nSlices,3);
    image_no = zeros(nSlices,1);
    instanceNumbers = zeros(nSlices,1);
    echoNumbers = zeros(nSlices,1);
    echoTimes = zeros(nSlices,1);
    for n_sl=1:numel(folder_content)
        image_stack(:,:,n_sl) = dicomread(folder_content{n_sl});
        image_no(n_sl) = str2double(strrep(folder_content{n_sl},common_string,''));
        image_header = dicominfo(folder_content{n_sl});
        IPP(n_sl,:) = double(image_header.ImagePositionPatient);
        instanceNumbers(n_sl) = image_header.InstanceNumber; % acquisition order
        if isfield(image_header,'EchoNumber')
            
        echoNumbers(n_sl) = image_header.EchoNumber;
        echoTimes(n_sl) = image_header.EchoTime;
        else
            echoNumbers=[];
            echoTimes = [];
        end
    end
%     if instanceNumbers(2)-instanceNumbers(1) >1
%         corrected_instance = instanceNumbers./(instanceNumbers(2)-instanceNumbers(1));
%     else
       [~, corrected_instance] = sort(instanceNumbers,'ascend');
%     end
    IPP = IPP(corrected_instance,:);
    
    
  col_index= IPP(1,:)-IPP(2,:)~=0;
  if IPP(1,col_index)<IPP(end,col_index)
      acquisition_order = 'AP';
  else
      acquisition_order = 'PA';
  end
  nEchos  = length(unique(echoNumbers));
  fprintf('%d echos detected: %6.2f.\n',nEchos,unique(echoTimes));
if nEchos>1 
    image_tmp = zeros(nRows,nCols,round(nSlices/nEchos),nEchos);
    for ne = 1:nEchos
        image_tmp(:,:,:,ne) = image_stack(:,:,echoNumbers==ne);
    end
    clear image_stack
    image_stack = image_tmp;
end

end

