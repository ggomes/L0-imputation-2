clear
close all

root = fileparts(fileparts(mfilename('fullpath')));

addpath(fullfile(root,'src','def'));

se = ScenarioPtr;
se.load(fullfile(root,'cfg','210W_v8.xml'));

fwy = Freeway(se);
