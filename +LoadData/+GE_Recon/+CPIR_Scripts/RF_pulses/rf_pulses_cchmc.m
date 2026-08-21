%% Hard pulses
% hard pulse for large flip angle/nutation
[rf(1).wf, rf(1).cvs] = design_hardpulse4fidall(195,true); %1ms
[rf(2).wf, rf(2).cvs] = design_hardpulse4fidall(176,true); %0.9ms
[rf(3).wf, rf(3).cvs] = design_hardpulse4fidall(156,true); %0.8ms
[rf(4).wf, rf(4).cvs] = design_hardpulse4fidall(137,true); %0.7ms
[rf(5).wf, rf(5).cvs] = design_hardpulse4fidall(116,true); %0.6ms