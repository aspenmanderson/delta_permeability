# Aspen Anderson
# Tests the functionality of calling python functions from matlab
# Last Modified: 4/5/2021

'''
To call from matlab:
    1) Change working directory to folder where function resides (addpath is not enough)
    2) Excute the following code:
    
% set python environment 
pe = pyenv;
if pe.Status == "NotLoaded"
    %[~,exepath] = system("where python");
    pypath = 'C:\Users\labuser10\Anaconda2\pythonw.exe'; %if where doesn't work, hard code path
        % found path using > import sys > sys.executable in Spyder
    pe = pyenv('Version',pypath);
end
 
clear mod_test
mod_test = py.importlib.import_module('testing_python')
py.reload(mod_test)
out = py.testing_python.test(5);  
'''

def test(input1):
    print("Test successful!") 
    
    output1 = input1 + input1
    
    return output1
    