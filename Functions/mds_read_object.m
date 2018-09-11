% mds_read_object.m

% Description:  This code extracts data from MDSplus data trees into
% MATLAB.  It interpolates to put all extracted signals on a common time
% base.  This common time base has the total extent of the shortest signal
% and takes the slowest sampling rate.  

% Inputs:
% shotnum:  A scalar value of the desired shot number.

% node_string:  An array of strings of the desired signal tag names.  Each row
% contains a tag name, and all rows must be of the same length (use extra
% spaces to fill in at the end of shorter tag names.  For instance:
% node_string = [
%     '\m_0_p0          ';
%     '\m_1_p0_norm     '
%     ];

% tree_string:  String identifying the tree.  Either 'zaphd' or 'zapmain'.

function [timebase,data] = mds_read_object(shotnum,node_string,tree_string)

import MDSplus.*

% mdsconnect('landau.hit');
% mdsconnect('zappa.zap');

node_name = cellstr([node_string]);   

for j = 1:length(node_name)
    for i = 1:length(shotnum)
        
        % open the experiment tree and extract the node name
        tree = Tree(tree_string,shotnum);
        n = tree.getNode(node_name(j));

        % extract each node's data and timebase
%         signal{j,i} = NATIVEvalue(n.getFloatArray());
        signal{j,i} = NATIVEvalue(n.getData().data());
        time{j,i} = NATIVEvalue(n.getDimensionAt(0).data());
        
        tree.close();

        % identify the min, max and time step for the timebase of each node
        time_tmp = time{j,i};
        min_time(j,i) = min(time_tmp);
        max_time(j,i) = max(time_tmp);
        del_time(j,i) = time_tmp(3)-time_tmp(2);
     
    end
end

% compute the common timebase
timebase = [max(max(min_time)):max(max(del_time)):min(min(max_time))];

m = length(node_name);
n = length(shotnum);
count = 1;

% resample each signal on the common timebase
for j = 1:length(node_name)
    for i = 1:length(shotnum)
        
        ts_node{j,i} = timeseries(signal{j,i},time{j,i});
        res_node{j,i} = resample(ts_node{j,i},timebase);
        data{j,i} = get(res_node{j,i},'Data');  %resampled data ready for plotting (j,i) = (node, shot)

%         node_final{j,i} = get(res_node{j,i},'Data');  %resampled data ready for plotting (j,i) = (node, shot)
%         data{j,i} = reshape(node_final{j,i},1,size(node_final{j,i},3));
        
        shot_number=num2str(shotnum(i));

    end
end

% mdsclose;
% mdsdisconnect;
