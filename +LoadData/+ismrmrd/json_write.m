function json_write(dcm_file,filename)

myinfo = dicominfo(dcm_file);

myjson.Manufacturer = myinfo.Manufacturer;
myjson.ManufacturersModelName = myinfo.ManufacturerModelName;
myjson.SoftwareVersions = myinfo.SoftwareVersions;
myjson.MagneticFieldStrength = myinfo.MagneticFieldStrength;
try
    myjson.ReceiveCoilName = myinfo.ReceiveCoilName;
catch
    myjson.ReceiveCoilName = myinfo.TransmitCoilName;
end

myjson.ScanningSequence = myinfo.ScanningSequence;
myjson.SequenceVariant = myinfo.SequenceVariant;
myjson.SequenceName = myinfo.SequenceName;


%% Let's see if we can guess what type of sequence this is from the dicom:

prot = myinfo.ProtocolName;
st = myinfo.SliceThickness;

isdiff = 0;
if contains(prot,'diff')
    isdiff = 1;
elseif st > 15
    isdiff = 1;
end

prompt = {'Enter type of Pulse Sequence:','Enter specific details regarding pulse sequence'};
dlgtitle = 'Input';
fieldsize = [1 45; 1 45];
if isdiff
    definput = {'Diffusion Weighted Gradient Echo','John Mugler Diffusion Sequence - HPG_Diffusion_2301'};
else
    definput = {'Gradient Echo','John Mugler Ventilation Sequence - gre_hpg_2201'};
end
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
myjson.PulseSequenceType = answer{1};
myjson.PulseSequenceDetails = answer{2};

myjson.MRAcquisitionType = myinfo.MRAcquisitionType;
myjson.FlipAngle = myinfo.FlipAngle;
myjson.ImagingFrequency = myinfo.ImagingFrequency;
myjson.ImagedNucleus = myinfo.ImagedNucleus;

json = jsonencode(myjson,PrettyPrint=true);

%% Write
fid = fopen([filename '.json'], 'w');
fprintf(fid, '%s ', json); % one row
fclose(fid);

