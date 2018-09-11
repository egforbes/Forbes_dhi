

function [B_prof,T_prof,v_drift,den_char,B_char,T_char,a_char] = ...
    DHI_profiles(rad,den,I,varargin)

[den_max,~] = max(den);
ind_tmp = length(den);

edge_frac = .5;

if isempty(varargin)
    % use the same integration bounds for both integrals:
    ind_int_d = [1:1:find(den<den_max*edge_frac,1,'first')];
%         ind_int_d = [1:1:length(rad)];
    ind_int = find(den<den_max*edge_frac,1,'first');
elseif nargin == 4
    % use the same integration bounds for both integrals, but allow for
    % specification of the upper integration boundary:
    integration_bounds_d = varargin{1};
    ind_int_d = find(rad<integration_bounds_d);
    ind_int = ind_int_d(end);
%     ind_int = find(den<den_max/2,1,'first');
elseif nargin == 5
    integration_bounds_d = varargin{1};
    ind_int_d = find(rad<integration_bounds_d);
    integration_bounds = varargin{2};
    ind_int = max(find(rad<integration_bounds));
end

% ind_tmp = length(den);
% ind_tmp = find(den<den_max/2,1,'first'); % integrate out to half density max
if isempty(ind_tmp) || ind_tmp == 1 || isempty(ind_int_d) || length(ind_int_d) <10% if the density doesn't drop to half the max value
    den_char = NaN;
    B_char = NaN;
    T_char = NaN;
    a_char = NaN;
    T_prof = NaN*ones(size(rad));
    B_prof = NaN*ones(size(rad));
    v_drift = NaN;
    
else 
    
    q = 1.6e-19;
    k_boltz =1.38e-23;
    mu_0 = (4*pi)*1e-7;
    
    % integrating the density:
%     try
        
if size(den(ind_int_d),1)~=size(rad(ind_int_d),1)
    rad = rad';
%     n_int = cumtrapz(rad(ind_int_d),den(ind_int_d)'.*rad(ind_int_d));
%     n_tot = n_int(end);
else
%     n_int = cumtrapz(rad(ind_int_d),den(ind_int_d).*rad(ind_int_d));
%     n_tot = n_int(end);
end

    n_int = cumtrapz(rad(ind_int_d),den(ind_int_d).*rad(ind_int_d));
    n_tot = n_int(end);
        
        % finding the drift velocity:
        v_drift = -I./(q*2*pi*n_tot);
        
        ind = find(rad<=rad(ind_int));
        %
        %     % integrating the density:
        %     n_int = cumtrapz(rad(ind),den(ind).*rad(ind));
        %     n_tot = n_int(end);
        %
        %     % finding the drift velocity:
        %     v_drift = -I./(q*2*pi*n_tot);
        
        % magnetic field profile:
        j_z = -v_drift.*q.*den(ind);
        B_prof_tmp = mu_0*cumtrapz(rad(ind),j_z.*rad(ind))./rad(ind);
        B_prof_tmp(1) = 0;
        [B_max, B_max_ind] = max(B_prof_tmp);
        B_decay = rad(B_max_ind)*B_max./rad;
        B_decay(1:B_max_ind-1) = 0;
        B_prof = zeros(size(rad));
        % B_prof = B_prof_tmp;
        B_prof(1:length(B_prof_tmp)) = B_prof_tmp;
        if length(B_prof_tmp) ~= length(B_prof)
            B_prof(length(B_prof_tmp)+1:end) = B_decay(length(B_prof_tmp)+1:end);
        end
        % B_prof(B_max_ind+1:end) = B_decay(B_max_ind+1:end);

        
        % temperature profile:
        % ind = find(rad<rad(B_max_ind)); % integrate out to max B
%         T = (1/11600)*cumtrapz(rad(ind),den(ind).*B_prof(ind)).*v_drift*q./(2*den(ind)*k_boltz);
        T = (1/11600)*cumtrapz(flipud(rad(ind)),flipud(den(ind)).*flipud(B_prof(ind))).*v_drift*q./(2*flipud(den(ind))*k_boltz);
        T_prof = zeros(length(rad),1);
        T_prof(1:length(T)) = flipud(T);
        T_char = max(T);
        
        % another means of doing the integration?:
%         n_edge_inx = max(ind);
%         rad_diff = (rad(2:n_edge_inx)-rad(1:n_edge_inx-1));
%         
%         T_prof = (1/11600)*(q*v_drift./(2*flipud(den(2:n_edge_inx))*k_boltz))...
%             .*cumsum(flipud(den(2:n_edge_inx)).*flipud(B_prof(2:n_edge_inx)).*(rad(2:n_edge_inx)-rad(1:n_edge_inx-1)));
        
%         T_prof = zeros(length(rad),1);
%         T_prof(1:length(T)) = (T-T(end));
        den_char = den_max;
        B_char = B_prof(ind_int);
%         T_char = max(T_prof);
        a_char = rad(ind_int);
        
%     catch
%         den_char = NaN;
%         B_char = NaN;
%         T_char = NaN;
%         a_char = NaN;
%         T_prof = NaN*ones(size(rad));
%         B_prof = NaN*ones(size(rad));
%         v_drift = NaN;
%     end
end






% q = 1.6e-19;
% k_boltz =1.38e-23;
% mu_0 = (4*pi)*1e-7;
% 
% % integrating the density:
% n_int = cumtrapz(rad,den.*rad);
% n_tot = n_int(end);
% 
% % finding the drift velocity:
% v_drift = -I./(q*2*pi*n_tot);
% 
% % magnetic field profile:
% j_z = -v_drift.*q.*den;
% B_prof = mu_0*cumtrapz(rad,j_z.*rad)./rad;
% B_prof(1) = 0;
% [B_max, B_max_ind] = max(B_prof);
% 
% % temperature profile:
% [den_max,~] = max(den);
% ind_tmp = find(den<den_max/2,1,'first'); % integrate out to half density max
% if isempty(ind_tmp) || ind_tmp == 1 % if the density doesn't drop to half the max value
%     den_char = NaN;
%     B_char = NaN;
%     T_char = NaN;
%     a_char = NaN;
%     T_prof = NaN;
%     B_prof = NaN;
%     
% else % if the density does drop to half the max value:
%     ind = find(rad<=rad(ind_tmp));
%     % ind = find(rad<rad(B_max_ind)); % integrate out to max B
%     T = (1/11600)*cumtrapz(rad(ind),den(ind).*B_prof(ind)).*v_drift*q./(2*den(ind)*k_boltz);
%     T_test = (1/11600)*cumtrapz(rad(ind),rad(ind).*den(ind).*B_prof(ind)).*v_drift*q./(2*den(ind)*k_boltz);
%     
%     T_prof = zeros(size(rad,2),1);
%     T_prof(1:length(T)) = (T-T(end));
%     den_char = den_max;
%     B_char = B_prof(ind_tmp);
%     T_char = max(T_prof);
%     a_char = rad(ind_tmp);
% end