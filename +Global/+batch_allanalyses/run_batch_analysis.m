function run_batch_analysis(inputs, analysisTypes, previewFolder, analyst)
%RUN_BATCH_ANALYSIS Run selected batch-analysis pipelines.

    if isempty(inputs)
        error('No input datasets were provided for batch analysis.');
    end

    if ~istable(inputs)
        error('Batch-analysis inputs must be a MATLAB table.');
    end

    if isempty(analysisTypes)
        error('No analysis types were selected.');
    end

    analysisTypes = string(analysisTypes);

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Starting Batch Analysis\n');
    fprintf('========================================\n');
    fprintf('Total datasets : %d\n', height(inputs));
    fprintf('Analysis types : %s\n', strjoin(analysisTypes, ', '));
    fprintf('Analyst        : %s\n', string(analyst));
    fprintf('Preview folder : %s\n', previewFolder);
    fprintf('========================================\n\n');

    % Analysis names in the input table
    rowAnalysis = lower(strtrim(string(inputs.Analysis)));

    % ================================================================
    % Ventilation
    % ================================================================
    if any(strcmpi(analysisTypes, "Ventilation"))

        ventInputs = inputs(rowAnalysis == "vent", :);

        if ~isempty(ventInputs)

            fprintf('\n----------------------------------------\n');
            fprintf('VENTILATION ANALYSIS\n');
            fprintf('Datasets: %d\n', height(ventInputs));
            fprintf('----------------------------------------\n\n');

            Global.batch_allanalyses.run_ventilation_batch( ...
                ventInputs, ...
                previewFolder, ...
                analyst);

        else
            fprintf('No vent datasets found in the filtered data.\n');
        end

    end

    % ================================================================
    % Diffusion
    % ================================================================
    if any(strcmpi(analysisTypes, "Diffusion"))

        diffInputs = inputs(rowAnalysis == "diff", :);

        if ~isempty(diffInputs)

            fprintf('\n----------------------------------------\n');
            fprintf('DIFFUSION ANALYSIS\n');
            fprintf('Datasets: %d\n', height(diffInputs));
            fprintf('----------------------------------------\n\n');

            Global.batch_allanalyses.run_diffusion_batch( ...
                diffInputs, ...
                previewFolder, ...
                analyst);

        else
            fprintf('No diff datasets found in the filtered data.\n');
        end

    end

    % ================================================================
    % Gas Exchange
    % ================================================================
    if any(strcmpi(analysisTypes, "GasExchange"))

        gxInputs = inputs(rowAnalysis == "gx", :);

        if ~isempty(gxInputs)

            fprintf('\n----------------------------------------\n');
            fprintf('GAS EXCHANGE ANALYSIS\n');
            fprintf('Datasets: %d\n', height(gxInputs));
            fprintf('----------------------------------------\n\n');

            Global.batch_allanalyses.run_gasexchange_batch( ...
                gxInputs, ...
                previewFolder, ...
                analyst);

        else
            fprintf('No gx datasets found in the filtered data.\n');
        end

    end

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Batch Analysis Complete\n');
    fprintf('========================================\n\n');

end