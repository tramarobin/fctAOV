function [nFiles]=findcol(files,string)

allFiles=zeros(numel(files),1);

for i=1:numel(files)
    allFiles(i)=contains(files(i),string);
end


nFiles=find(allFiles);

end