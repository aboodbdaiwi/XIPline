% read the data from complex data pair with an extra1 parameter
%   files accompanying .data/.list
%   
%   
%   
%
%
function [ImgK,NoiK,kx_oversample_factor] = load_philips_extr1_2D(filename,ch_range)

%chIndexKluge ='yes';

FileNameData = [filename,'.data'];
FileNameList = [filename,'.list'];

fidList = fopen(FileNameList);
fidData = fopen(FileNameData);

DataType = 'float32';
DataBytesPerType = 4;
ch_size=ch_range(2)-ch_range(1)+1;
NoiK=[];
ImgK=[];

while ~feof(fidList)
    textline = fgetl(fidList);
    if ((length(textline)>54)&&(textline(1)=='.'))
        if (strcmp(textline(1:54),'.    0    0    0  number_of_locations                :'))
            sl_size=str2num(textline(55:end));
        end
    
        if (strcmp(textline(1:54),'.    0    0    0  number_of_extra_attribute_1_values :'))
            extr1_size=str2num(textline(55:end));
        end
    end

    if (strcmp(textline,'# mix  echo n.a.  k-space coordinate ranges            start  end'))
        textline=fgetl(fidList);
        textline=fgetl(fidList);
        kx_range=str2num(textline(55:end));
        textline=fgetl(fidList);
        ky_range=str2num(textline(55:end));
        
        kx_size=kx_range(2)-kx_range(1)+1;
        ky_size=ky_range(2)-ky_range(1)+1;    
        ImgK=zeros([ky_size kx_size sl_size ch_size extr1_size]);    
    end
    
    if (strcmp(textline,'# mix  echo n.a.  k-space oversample factors           value'))
        textline=fgetl(fidList);
        textline=fgetl(fidList);
        kx_oversample_factor=str2num(textline(58:end));
    end
    
    if (textline(1)==' ')
        if (textline(1:5)=='  NOI')
            dataIndex=str2num(textline(6:end));
            
%             if strcomp(chIndexKluge,'yes')==1  
%                 chIndex = 1;
%             else
%                 chIndex=dataIndex(6)+1;
%             end
            chIndex=dataIndex(6)+1;
            if ((chIndex<=ch_range(2))&&(chIndex>=ch_range(1)))
                fseek(fidData, dataIndex(end),'bof');
                datasizeBytes=dataIndex(end-1);
                temp=fread(fidData, datasizeBytes/DataBytesPerType, DataType);
                NoiK(:,chIndex-ch_range(1)+1)=complex(temp(1:2:length(temp)), temp(2:2:length(temp)));
            end
        end
            
        if (textline(1:5)=='  STD')
            dataIndex=str2num(textline(6:end));
            
            slIndex=dataIndex(5)+1;
            chIndex=dataIndex(6)+1;
            extr1Index=dataIndex(7)+1;
            
            if ((chIndex<=ch_range(2))&&(chIndex>=ch_range(1))&&(slIndex<=sl_size)&&(extr1Index<=extr1_size))
                kyIndex=dataIndex(9)-ky_range(1)+1;
                fseek(fidData, dataIndex(end),'bof');
                datasizeBytes=dataIndex(end-1);
                temp=fread(fidData, datasizeBytes/DataBytesPerType, DataType);
                ImgK(kyIndex,:,slIndex,chIndex-ch_range(1)+1,extr1Index)=squeeze(complex(temp(1:2:2*kx_size), temp(2:2:2*kx_size)));
            end
        end
    end
end
fclose(fidList);
fclose(fidData);


