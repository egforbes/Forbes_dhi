% This function extracts the value of a signal identified in node_string.
% The identified signal cannot be a time base.  This function is used
% primarily to extract trigger times programmed into the DG11s and L4222s.

function [value] = mds_eval(shotnum,node_string,tree_string);

% constants for debugging:
% shotnum = 160203041;
% node_string = '\ZAPHD::TOP.DIGITIZERSHD.L4222:ZAP_4222_B:CHANNEL_1';
% tree_string = 'zaphd';

import MDSplus.*
tree = Tree(tree_string,shotnum);
n = tree.getNode(node_string);
value = NATIVEvalue(n.getData().data());

tree.close();
