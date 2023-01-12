function getFunctionDirectory

cd=pwd;
idx=strfind(cd,'\');
if isempty(idx)
    idx=strfind(cd,'/');
end
addpath(genpath(cd(1:idx(end))));

end