function checkForUpdates(currentVersion)
    % Function to check for new releases on GitHub
    %
    % Input:
    %   currentVersion - (string) the current version of the app (e.g., 'V1.1.0')
    %   repoOwner      - (string) the GitHub username or organization name
    %   repoName       - (string) the name of the GitHub repository
    %
    % Example usage:
    %   checkForUpdates('V1.1.0', 'myusername', 'myapp')
    %currentVersion = 'V1.1.0';
    repoOwner = 'aboodbdaiwi';
    repoName = 'XIPline';
    % GitHub API URL for releases
    % https://github.com/aboodbdaiwi/XIPline/tree/V1.1.0
    apiUrl = sprintf('https://github.com/%s/%s/releases', repoOwner, repoName);

    try
        % Fetch release data from GitHub API
        options = weboptions('Timeout', 10);
        textData = webread(apiUrl, options);

        % Regular expression pattern to match the desired text and extract the version
        pattern = '/aboodbdaiwi/XIPline/releases/tag/(V[0-9]+(?:\.[0-9]+)*)';
        
        % Use regexp to extract all matches
        matches = regexp(textData, pattern, 'tokens');
        
        % Flatten the nested cell array into a single cell array
        versionTags = [matches{:}];

        % Extract the latest release version
        latestVersion = versionTags{end}; % e.g., 'V1.2.0'
        pause(2);
        % Compare versions
        if ~strcmp(currentVersion, latestVersion)
            % Notify the user about the new version
            message = sprintf('A new version (%s) of the app is available. You are using version %s.', ...
                latestVersion, currentVersion);
              titleBarCaption = 'Updates';
              button = questdlg(message, titleBarCaption, 'Update Now', 'Later','yes');
              if strcmpi(button, 'Update Now')
                 web('https://github.com/aboodbdaiwi/XIPline', '-browser'); % Open the release page in the default browser
              else
                delete(message);
              end 
        end
    catch 
        % Handle errors (e.g., network issues or invalid repo)
        message = 'Could not check for updates. Please check your internet connection or try again later.';
          titleBarCaption = 'Updates';
          button = questdlg(message, titleBarCaption, 'Okay','yes');
          if strcmpi(button, 'Okay')
            delete(message);
          end 
    end
end
