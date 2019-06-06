
function G = G_fun(min_fun,rL,rs,r1,r2,condit)

if nargin<6
    condit = ones(size(rs));
end

ns = length(rs); % number of equilirbrium radii to consider (one per particle mass considered)
nL = length(rL); % number of final radial positions to consider (usually only one)
offset = 1e-9; % offset of rs in which rL is considered equal to rs

G = zeros(nL,ns); % pre-allocate output variable
for jj=1:nL
    for ii=1:ns
        if rL(jj)<r1
        % particles cannot exist here, return r1
            G(jj,ii) = r1;
            
        elseif rL(jj)>r2
        % particles cannot exist here, return r2
            G(jj,ii) = r2;
            
        elseif and(rL(jj)<(rs(ii)+5*offset),rL(jj)>(rs(ii)-5*offset))
        % if rL = rs, then r0 = rs
            G(jj,ii) = rL(jj);
            
        elseif rL(jj)>rs(ii)
        % if rL > rs, then r0 > rL
        
            if condit(ii)
                if sign(min_fun(rL(jj),rL(jj),ii))==sign(min_fun(rL(jj),r2,ii))
                % if sign does not change in interval [rL,r2], return r2
                    G(jj,ii) = r2;
                    
                else
                % find zero in interval [rL,r2]
                    G(jj,ii) = fzero(@(r0) min_fun(rL(jj),r0,ii), [rL(jj),r2]);
                end
                
            else
                if sign(min_fun(rL(jj),rL(jj),ii))==sign(min_fun(rL(jj),max(r1,rs(ii)+offset),ii))
                % if sign does not change in interval [max(r1,rs),rL], return max(r1,rs)
                    G(jj,ii) = max(r1,rs(ii)+offset);
                    
                else
                % find zero in interval [max(r1,rs),rL]
                    G(jj,ii) = fzero(@(r0) min_fun(rL(jj),r0,ii), [max(r1,rs(ii)+offset),rL(jj)]);
                    
                end
            end
            
        elseif rL(jj)<rs(ii) % should be equivalent to else
        % if rL < rs, then r0 < rL
        
            if condit(ii)
                if sign(min_fun(rL(jj),rL(jj),ii))==sign(min_fun(rL(jj),r1,ii))
                % if sign does not change in interval [r1,rL], return r1
                    G(jj,ii) = r1;
                    
                else
                % find zero in interval [r1,rL]
                    G(jj,ii) = fzero(@(r0) min_fun(rL(jj),r0,ii), [r1,rL(jj)]);
                    
                end
                
            else
                if sign(min_fun(rL(jj),rL(jj),ii))==sign(min_fun(rL(jj),min(rs(ii)-offset,r2),ii))
                % if sign does not change in interval [rL,min(r2,rs)], return min(r2,rs)
                    G(jj,ii) = min(rs(ii)-offset,r2);
                    
                else
                % find zero in interval [rL,min(r2,rs)]
                    G(jj,ii) = fzero(@(r0) min_fun(rL(jj),r0,ii), [rL(jj),min(rs(ii)-offset,r2)]);
                    
                end
            end
            
        end
    end
end

end

