function hash = get_git_hashobject
    git_folder = mfilename('fullpath');%get path to this function
    idcs = strfind(git_folder,filesep);%determine location of file separators
    git_folder = git_folder(1:idcs(end-2)-1);%remove file and move up two folders
    git_folder = [git_folder,filesep,'.git'];%add path to .git folder
    command = [ 'git --git-dir="', git_folder ,'" rev-parse HEAD '];%make git command
    [status,hash] = system(command);%run command
    if( status ~= 0 )
        error('Unable to get hash from file.');
    end
end