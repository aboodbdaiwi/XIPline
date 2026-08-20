% RH: new file — shared by N4_bias_correction.m's Z/X/Y passes to avoid
% tripling the Windows-exe vs Python-script branch. Runs one
% N4BiasFieldCorrection pass with the given -s/-c/-b/-t arguments.
function run_N4_pass(N4Path, imagePath, maskPath, weightPath, shrinkFactor, convergence, bspline, histogram, correctedPath, biasPath)
    if ispc
        pathReg = fullfile(N4Path, 'N4BiasFieldCorrection.exe');
        cmd = sprintf('"%s" -d 3 -i "%s" -s %s -x "%s" -w "%s" -c %s -b %s -t %s -o ["%s","%s"]', ...
            pathReg, imagePath, shrinkFactor, maskPath, weightPath, convergence, bspline, histogram, correctedPath, biasPath);
    else
        XIPlineRoot = fullfile(getenv('HOME'), 'XIPline');
        registration_path = fullfile(XIPlineRoot, 'Registration');
        pythonExe = Registration.ANTs.getAntsPythonExe(registration_path);
        n4Script = fullfile(N4Path, 'n4_bias_correction.py');
        cmd = sprintf('"%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"', ...
            pythonExe, n4Script, imagePath, maskPath, weightPath, shrinkFactor, ...
            convergence, bspline, histogram, correctedPath, biasPath);
    end
    system(cmd);
end
