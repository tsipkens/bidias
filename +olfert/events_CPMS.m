function [value,isterminal,direction] = events_couAPM(z,r_ode,v_flow,ro,ri,rc,vc,vi,B,m,q,V)
% Locate the time when radial distance passes through
%the inner or outer radius of the cylinder
rc = (ro+ri)/2;
value = abs(r_ode-rc)-(ro-ri)/2;     % Detect when it hits wall 
                                     %when value =0;
isterminal = 1;     % Stop the integration
direction = 0;      %don't care which direction the particle is going  