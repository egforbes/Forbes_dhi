function [time,signal] = mds_read_efficient(shotnum,node_string,tree_string)

import MDSplus.*;

% mdsconnect('landau.hit');
% mdsconnect('zappa.zap');

node_name = cellstr([node_string]);   
        
% open the experiment tree and extract the node name
tree = Tree(tree_string,shotnum);
n = tree.getNode(node_name);

% extract each node's data and timebase
%         signal{j,i} = NATIVEvalue(n.getFloatArray());
signal = NATIVEvalue(n.getData().data());
time = NATIVEvalue(n.getDimensionAt(0).data());

% identify the min, max and time step for the timebase of each node
time_tmp = time;
min_time = min(time_tmp);
max_time = max(time_tmp);
del_time = time_tmp(3)-time_tmp(2);

tree.close();


     
