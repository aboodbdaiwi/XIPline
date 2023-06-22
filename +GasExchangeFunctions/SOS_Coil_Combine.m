function Image = SOS_Coil_Combine(Image_Stack)

%Function to take a stack of images 
%(Size =ImageSize(1),ImageSize(2),ImageSize(3),NCoils)
%Return sum of squares combination of the images following GPI's
%implementation of this

nCoil = size(Image_Stack,4);

tempSOS = zeros(size(Image_Stack));
for i = 1:nCoil
    tempSOS(:,:,:,i) = abs(Image_Stack(:,:,:,i).^2);
end

Image = squeeze(sqrt(sum(tempSOS,4)));