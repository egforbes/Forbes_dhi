% For a range of shotnumbers, this code identifies the baseline and plasma
% shots.  It then, pairs each plasma shot with the nearest prior baseline
% and outputs the matrix "shots" for use in subsequent analysis codes.

function [shots] = shot_base_pair(shotnum,varargin)

% inspect the YAG laser trigger output
node_string = '\ZAPHD::TOP.DIGITIZERSHD.L4222:ZAP_4222_B:CHANNEL_1';
tree_string = 'zaphd';
node_name = node_string;
shotnum_def = [];
shotnum_base = [];

    % identify the deformed and baseline shots based on laser trigger time
    for i = 1:length(shotnum) 
        import MDSplus.*
        tree = Tree(tree_string,shotnum(i));
        n = tree.getNode(node_name);
        signal = NATIVEvalue(n.getData().data());
        tree.close();

        if signal > 0
            shotnum_def = [shotnum_def shotnum(i)];
        else
            shotnum_base = [shotnum_base shotnum(i)];
        end
        
    end
        
    % using varargin to accept input identifying if the user desires
    % deformed-baseline or baseline-baseline pairs
    switch(nargin),
        case 0,
            whichPair = 1;
        case 1,
            whichPair = 1;
        case 2,
            whichPair = varargin{1};
    end
    
    if whichPair == 1
        % pair the nearest deformed-baseline hologram pairs
        i = 1;
        shots = zeros(length(shotnum_def),2);
        for shotnum_def_it = shotnum_def
            ind_base_tmp = find((shotnum_def_it-shotnum_base)>0,1,'last');
            shots(i,1) = shotnum_def_it;
            shots(i,2) = shotnum_base(ind_base_tmp);
            i = i+1;
        end
        
    elseif whichPair == 0
        % pair the nearest baseline-baseline hologram pairs
        i = 1;
        shots = zeros(length(shotnum_base),2);
        for shotnum_base_it = shotnum_base
            ind_base_tmp = find((shotnum_base_it-shotnum_base)>0,1,'last');
            if isempty(ind_base_tmp)
                shots = zeros(length(shotnum_base)-1,2);
            else
                if isempty(shots)
                    shots = zeros(length(shotnum_base),2);
                else
                end
            shots(i,1) = shotnum_base_it;
            shots(i,2) = shotnum_base(ind_base_tmp);
            i = i+1;
            end
        end
        
    end