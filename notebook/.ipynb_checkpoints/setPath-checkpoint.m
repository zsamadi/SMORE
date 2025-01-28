function setPath()

cd(fileparts(which(mfilename)));
cd('..\')
addpath(genpath(pwd))
