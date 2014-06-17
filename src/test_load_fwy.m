clear
close all

% root = fileparts(fileparts(mfilename('fullpath')));
% addpath(fullfile(root,'src','def'));
% 
% day = 735633;
% goodthresh = 0.7;
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
% fwy.load_loop_health_for_days(fullfile(root,'generated','data'),[735632 735633 234532 735634 735637])
% 
% % test load data
% fwy.load_5min_data_for_days(fullfile(root,'generated','data'),[735632 735633 234532 735634 735637])
% 
% % test has good sensor on day
% fwy.seg(1).has_good_sensor_on_day('ml',day,goodthresh)

load aaa

fwy.disable_sensor_vds([772918 717682 717673 761374 717657 717649 717642 717624]);

% test
x=fwy.get_sensor_health_for_day('ml',day,goodthresh);

x = fwy.get_all_sensor_health_for_day(day,goodthresh);

% test build cell array
cells = build_cell_array(fwy,day,goodthresh);

y=cells.get_num_seg_per_cell();

disp('done')