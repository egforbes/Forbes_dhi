% MDS_read_multisig.m
% Aug. 5, 2016
% Author:  Michael Ross
% Description:  Reads MDS data into MATLAB employing python.  Reads in data
% from multiple signals, finds a common timebase between all the signals,
% and resamples the signals onto that common timebase.  Using python
% enables calls from the MDS data server without usual memory issues
% encountered when calling for data directly in MATLAB.

function [timebase,signals] = MDS_read_multisig(shotnum,node_name_list)

% extract the data from the MDSplus tree:
for i = 1:size(node_name_list,1)
    node_name_spaces = isspace(node_name_list(i,:));
    node_name_tmp = node_name_list(i,not(node_name_spaces));
    eval(sprintf(['[time{i},',node_name_tmp,'_raw] = MDS_MATLAB_python(shotnum,''',node_name_tmp,''');']));
end

% compute a common timebase:
for i = 1:length(time)
    
    time_tmp = time{i};
    min_time(i) = min(time_tmp);
    max_time(i) = max(time_tmp);
    del_time(i) = time_tmp(3)-time_tmp(2);
    
end
timebase = [max(max(min_time)):max(max(del_time)):min(min(max_time))]';

% resample the data onto the common timebase:
for i = 1:size(node_name_list,1)
    node_name_spaces = isspace(node_name_list(i,:));
    node_name_tmp = node_name_list(i,not(node_name_spaces));
    eval(sprintf(['ts_',node_name_tmp,' = timeseries(',node_name_tmp,'_raw,time{i});']));
    eval(sprintf(['res_',node_name_tmp,' = resample(ts_',node_name_tmp,',timebase);']));
    eval(sprintf([node_name_tmp,' = get(res_',node_name_tmp,',''Data'');']));
    signals(i,:) = eval(node_name_tmp);
end
