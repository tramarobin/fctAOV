function [nFiles]=findcol(files,string)

allFiles=zeros(numel(files),1);

if ~isempty(string)
    for i=1:numel(files)
        allFiles(i)=contains(files(i),string);
    end
end

nFiles=find(allFiles);

end