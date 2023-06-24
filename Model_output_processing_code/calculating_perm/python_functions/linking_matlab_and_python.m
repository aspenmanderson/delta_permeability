% Author: Aspen Anderson
% This code uses python functions to get the hydrualic conductivity and
% permeability of Delft3D delta models
% Last Modified: 4/05/2021

clc; close all; clear all;

% add function files
addpath("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\conductivity_paper\calculating_perm\python_functions"); %python functions

% define output destination 
ex_pth = ("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\Sed Modeling\matlab codes\hydrocon_outputs\fluvial_data\");

set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Times') 
set(0,'DefaultAxesFontWeight','Demi')
set(0,'DefaultAxesLineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load python functions (before loading data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cwd need to be the folder that the python functions are in
cd ("C:\Users\labuser10\Anaconda2\Scripts"); %python functions

% set python environment (https://blogs.mathworks.com/loren/2020/03/03/matlab-speaks-python/)
%pyversion('C:\Users\labuser10\Anaconda2\pythonw.exe')
pe = pyenv;
setenv('path',['C:\Users\labuser10\Anaconda2\Library\bin;', getenv('path')]);

if pe.Status == "NotLoaded"
    %[~,exepath] = system("where python");
    pypath = 'C:\Users\labuser10\Anaconda2\pythonw.exe'; %if where doesn't work, hard code path
        % found path using > import sys > sys.executable in Spyder
    pe = pyenv('Version',pypath);
end

% testing python
myPythonVersion = pe.Version;
py.print("Hello, Python!");

% testing my python script (in python functions folder
clear classes
mod_test = py.importlib.import_module('testing_python')
py.reload(mod_test)
out = py.testing_python.test(5);

% testing my python script 2
clear classes
mod = py.importlib.import_module('calc_hyd_params_for_matlab')
py.reload(mod)
out = py.calc_hyd_params_for_matlab.test(10)


% calculate permeability (doesn't work because the python to matlab data
% exchange puts the output in an odd format)
%perm_py = py.calc_hyd_params_for_matlab.calc_perm(x50,M,N)
%perm_py2 = double(py.array.array('d',py.numpy.nditer(perm_py)));
%perm = reshape(perm_py2,M,N);
%imagesc(perm)
