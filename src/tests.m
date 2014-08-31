clear
close all

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','def'));

day = 735633;
goodthresh = 0.7;
MajorOR = [31 96];
MajorFR = [25 92];
        
% test load scenario
se = ScenarioPtr;
se.load(fullfile(root,'cfg','210W_v8.xml'));

fwy = Freeway(se);

% test set major
fwy.set_major(MajorOR,'or');
fwy.set_major(MajorFR,'fr');
or_links = [fwy.seg.or_links];
find([or_links.is_major])

fr_links = [fwy.seg.fr_links];
find([fr_links.is_major])

% test get link id
L = fwy.get_link_for_id(196);

% test get sensor with vds
S = fwy.get_sensor_with_vds(717624);

% test get sensor with link id
S = fwy.get_sensor_with_linkid(196);

% test load loop health for days
fwy.load_loop_health_for_days(fullfile(root,'generated','data'),[735632 735633 234532 735634 735637])

% test load data
fwy.load_5min_data_for_days(fullfile(root,'generated','data'),[735632 735633 234532 735634 735637])

% test has good sensor on day
fwy.seg(1).has_good_sensor_on_day('ml',day,goodthresh)

fwy.disable_sensor_vds([772918 717682 717673 761374 717657 717649 717642 717624]);

% test
x=fwy.get_sensor_health_for_day('ml',day,goodthresh);

x = fwy.get_all_sensor_health_for_day(day,goodthresh);

time = (3600:30:7200);
data_type = 'flw';
% x=fwy.seg(1).ml_link.sensor.get_5min_data(day,time,data_type);

% test build cell array
cell_array = build_cell_array(fwy,day,goodthresh,'use_first_segment');

save(fullfile(root,'analysis','fwy'),'cell_array','fwy')

% test has_or_dwnstr(obj)
fwy.seg(1).has_or_dwnstr

% cells.get_num_seg_per_cell()

% load(fullfile(root,'generated','imputation','AfterFilledCellData_test1_03-Feb-2014'));

load('C:\Users\gomes\code\L0\L0-imputation-2\analysis\AfterFilledCellData_test1_03-Feb-2014')

%.ml_link.sensor.get_5min_data(day,time,'flw')

x=[cell_array.cells(1).segments(1:3).ml_link]
y=[x.sensor]

make_imputation_input(cell_array)

fwy=imputation_core_algorithm(cell_array);
disp('done')