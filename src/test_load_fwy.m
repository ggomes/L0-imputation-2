clear
close all

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','def'));

% test load scenario
se = ScenarioPtr;
se.load(fullfile(root,'cfg','210W_v8.xml'));
fwy = Freeway(se);

% test get link id
L = fwy.get_link_for_id(196);

% test get sensor with vds
S = fwy.get_sensor_with_vds(717624);

% test get sensor with link id
S = fwy.get_sensor_with_linkid(196);




disp('done')