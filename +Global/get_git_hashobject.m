function hash = get_git_hashobject
    git_folder = mfilename('fullpath');%get path to this function
    % idcs = strfind(git_folder,filesep);%determine location of file separators
    % git_folder = git_folder(1:idcs(end-2)-1);%remove file and move up two folders
    git_folder = getBaseXIPlinePath(git_folder);
    git_folder = [git_folder,filesep,'.git'];%add path to .git folder
    command = [ 'git --git-dir="', git_folder ,'" rev-parse HEAD '];%make git command
    [status,hash] = system(command);%run command
    if( status ~= 0 )
        error('Unable to get hash from file.');
    end
end

function basePath = getBaseXIPlinePath(fullPath)
    % Extracts base path "XIPline" folder from a full subpath
    parts = strsplit(fullPath, filesep);
    idx = find(strcmpi(parts, 'XIPline'), 1, 'last');
    if ~isempty(idx)
        basePath = fullfile(parts{1:idx});
    else
        error('XIPline folder not found in path.');
    end
end
