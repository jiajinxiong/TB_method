fullPath = mfilename('fullpath');
[filePath, ~, ~] = fileparts(fullPath);
addpath(filePath);
fun1_path = 'slanCM';
fullPath1 = fullfile(filePath, fun1_path);
if exist(fullPath1,'dir') == 7
    addpath(fullPath1);
end