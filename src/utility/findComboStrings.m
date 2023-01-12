function [nFiles]=findComboStrings(files,strings)

allFiles=zeros(numel(files),numel(strings));

for i=1:numel(files)
    for s=1:numel(strings)
        allFiles(i,s)=contains(files(i),strings{s});
    end
end


nFiles=find(all(allFiles')',1);

end