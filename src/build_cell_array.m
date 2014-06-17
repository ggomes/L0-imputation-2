function [ cell_array ] = build_cell_array(fwy,day,threshold)

 cell_array=CellArray;
for i=1:length(fwy.seg)
    has_good_ml_sensor = fwy.seg(i).has_good_sensor_on_day('ml',day,threshold);
    make_new_cell = i==1 || has_good_ml_sensor;
    if(make_new_cell)
        cell_array.add_cell(Cell);
    end
    cell_array.add_segment_to_last_cell(fwy.seg(i));
end