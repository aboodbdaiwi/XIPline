classdef input_params 
%% Input parameters for MATLAB VDP Analysis Script
% The goal of this function is to define the input settings for the Main
% VDP analysis. 
%
% Instructions:
% input_params = input_params((medfilter,...
%                 RFcorrection,...
%                 savedata,...
%                 calculateSNR,...
%                 original,...
%                 N4,...
%                 incomplete,...
%                 complete,...
%                 hyper)
%
% ...where all input variables are a numeric.
%
% Usage:
% example = Functions.input_params(1,0,1,1,0,60,15,200)
%
% Author: Joseph Plummer
% Date: 04/29/2021
	properties
        % These are the values that will be outputted by the classdef.
        
        % Settings:
        Median_Filter 
        RF_correction
        savedata
        calculate_SNR 
        
        % Analysis type settings:
        N4_bias_analysis 
        
        % Thresholds:
        Incomplete_threshold 
        Complete_threshold 
        Hyper_threshold
	end
	
	methods 
        % Default input values:
%         medfilter = 1
%         RFcorrection = 0
%         savedata = 1
%         calculateSNR = 1
%         N4 = 1
%         incomplete = 60
%         complete = 15
%         hyper = 200
        
		function input_params = input_params(medfilter,...
                RFcorrection,...
                savedata,...
                calculateSNR,...
                N4,...
                incomplete,...
                complete,...
                hyper)
            
            if medfilter == 1
                medfilter = "yes";
            elseif medfilter == 0
                medfilter = "no";
            else 
                error('Enter medfilter value of 1 (for yes) or 0 (for no)')
            end
            
            if RFcorrection == 1
                RFcorrection = "yes";
            elseif RFcorrection == 0
                RFcorrection = "no";
            else 
                error('Enter RF Correction value of 1 (for yes) or 0 (for no)')
            end
            
            if savedata == 1
                savedata = "yes";
            elseif savedata == 0
                savedata = "no";
            else 
                error('Enter savedata value of 1 (for yes) or 0 (for no)')
            end
            
            if calculateSNR == 1
                calculateSNR = "yes";
            elseif calculateSNR == 0
                calculateSNR = "no";
            else 
                error('Enter calculateSNR value of 1 (for yes) or 0 (for no)')
            end
            
            if N4 == 1
                N4 = "yes";
            elseif N4 == 0
                N4 = "no";
            else 
                error('Enter N4 value Enter 1 (for yes) or 0 (for no)')
            end
            
            if 0 <= incomplete <= 100
                incomplete = incomplete;
            else 
                error('Incomplete threshold must lie between 0 and 100 percent')
            end
            
            if 0 <= complete <= incomplete
                complete = complete;
            else 
                error('Complete threshold must lie between 0 and 100 percent')
            end
            
            if 1 <= hyper 
                hyper = hyper;
            else 
                error('Hyperventilated threshold must be greater than 1')
            end
            
            % Intialize the settings for each instance:
            input_params.Median_Filter = medfilter
            input_params.RF_correction = RFcorrection
            input_params.savedata = savedata
            input_params.calculate_SNR = calculateSNR
            input_params.N4_bias_analysis = N4
            input_params.Incomplete_threshold = incomplete
            input_params.Complete_threshold = complete
            input_params.Hyper_threshold = hyper

		end
	
	end
end