clear
close all

% root = fileparts(fileparts(mfilename('fullpath')));
% addpath(fullfile(root,'src','def'));
% 
% % test load scenario
% se = ScenarioPtr;
% se.load(fullfile(root,'cfg','210W_v8.xml'));
% 
% fwy = Freeway(se);
% 
% % test get link id
% L = fwy.get_link_for_id(196);
% 
% % test get sensor with vds
% S = fwy.get_sensor_with_vds(717624);
% 
% % test get sensor with link id
% S = fwy.get_sensor_with_linkid(196);
%
% % test load loop health for days
% fwy.load_loop_health_for_days(fullfile(root,'generated','data'),[735632 234532 735634 735637])
% 
% % test load data
% fwy.load_5min_data_for_days(fullfile(root,'generated','data'),[735632 234532 735634 735637])
% 
% % test has good sensor on day
% fwy.seg(1).has_good_sensor_on_day('ml',735632,90)

load aaa

% test
x=fwy.get_sensor_health_for_day('ml',735632,90);

% test build cell array
cells = build_cell_array(fwy,735632,20);

y=cells.get_num_seg_per_cell()

disp('done')