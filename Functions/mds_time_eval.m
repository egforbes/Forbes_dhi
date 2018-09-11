% This function extracts the value of the signal identified in node_string
% at a time instant given in time_in.

function [time_out,value] = mds_time_eval(shotnum,time_in,node_string,tree_string)

% constants for debugging:
% shotnum = 160203041;
% time_in = 50e-6;
% node_string = '\m_0_p15';
% tree_string = 'zaphd';

node_name = cellstr([node_string]);

import MDSplus.*
tree = Tree(tree_string,shotnum);
n = tree.getNode(node_name);

signal = NATIVEvalue(n.getData().data());
time = NATIVEvalue(n.getDimensionAt(0).data());

tree.close();

ind_time = find((time_in - time)<0,1,'first');
value = signal(ind_time);
time_out = time(ind_time);