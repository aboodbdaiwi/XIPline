function xemri_mrd()

Subj_ID = inputdlg('Participant ID','Input Participant ID',[1 50]); % Prompt User for Participant ID
Subj_ID = Subj_ID{1};

Contrasts = {'Ventilation','Diffusion','Gas Exchange'};

Files_Acquired = listdlg('PromptString',{'What Contrasts were Collected?',...
    'Use Shift-Click to select multiple options'},'ListString',Contrasts);

if ismember(1,Files_Acquired)
    vent2mrd(Subj_ID);
end
if ismember(2,Files_Acquired)
    diff2mrd(Subj_ID);
end
if ismember(3,Files_Acquired)
    gx2mrd(Subj_ID);
end