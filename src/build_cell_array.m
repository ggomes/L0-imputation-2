function [ cell_array ] = build_cell_array(fwy,day,threshold,param_rule)

cell_array=CellArray;
for i=1:length(fwy.seg)
    has_good_ml_sensor = fwy.seg(i).has_good_sensor_on_day('ml',day,threshold);
    make_new_cell = i==1 || has_good_ml_sensor;
    if(make_new_cell)
        cell_array.add_cell(Cell);
    end
    cell_array.add_segment_to_last_cell(fwy.seg(i));
end

% prescribe an fd and lanes to each cell
for i=1:length(cell_array.cells)
    cell = cell_array.cells(i);
    switch(param_rule)
        case 'use_first_segment'
            ml_fd = cell.segments(1).ml_link.fd;        % reference
            ml_lanes = cell.segments(1).ml_link.lanes;  % copy
        case 'use_last_segment'
            error('not implemented')
        otherwise
            error('unknown parameter assignment rule')
    end
    set(cell,'ml_fd',ml_fd);
    set(cell,'ml_lanes',ml_lanes);
end