function [x,y,t] = importfile(filename, dataLines,variableNames,variableTypes)
%IMPORTFILE 从文本文件中导入数据
%  UNTITLED = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以表形式返回数据。
%
%  UNTITLED = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  Untitled = importfile("/Volumes/KODAK/data/tianchi/VIS/hy_round1_train_20200102/0.csv", [2, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2021-02-25 11:16:07 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [2, Inf];
end
if (nargin<3)
    variableNames = ["ID", "x", "y", "VarName4", "VarName5", "time", "type"];
end
if (nargin<4)
    variableTypes = ["double", "double", "double", "double", "double", "string", "categorical"];
end
xcol = find(variableNames=="x");
ycol = find(variableNames=="y");

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 7);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = variableNames;
opts.VariableTypes = variableTypes;

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, "time", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["time", "type"], "EmptyFieldRule", "auto");

% 导入数据
r = readtable(filename, opts);
x = table2array(r(:,xcol));
y = table2array(r(:,ycol));
t = [1:1:size(x)]';


end