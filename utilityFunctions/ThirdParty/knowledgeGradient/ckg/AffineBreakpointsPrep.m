% Prepares vectors for passing to AffineEmaxBreakpoints, changing their
% order and removing elements with duplicate slope.
function [a,b] = AffineBreakpointsPrep(a,b)
    % Make sure a and b are column vectors.
    rows = size(a); if (rows == 1), a=a'; end
    rows = size(b); if (rows == 1), b=b'; end

    % 11/29/2008 PF: Experimental preprocessing step, which I hope will remove
    % a large number of the entries.
    [b1, i1] = min(b); % [a1,b1] is best at z=-infinity
    [a2, i2] = max(a); % [a2,b2] is best at z=0
    [b3, i3] = max(b); % [a3,b3] is best at z=+infinity
    a1 = a(i1);
    b2 = b(i2);
    a3 = a(i3);
    cleft = (a - a1)./(b1 - b); % intersection with leftmost line. 
    cright = (a - a3)./(b3 - b); % intersection with rightmost line.
    c2left = (a2 - a1)./(b1 - b2); % intersection with leftmost line. 
    c2right = (a2 - a3)./(b3 - b2); % intersection with rightmost line.
    keep = find(b==b1 | b==b3 | cleft <= c2left | cright >= c2right);
    %disp(sprintf('Preprocessing cut %d of %d entries', length(a)-length(keep), length(a)));
    a = a(keep);
    b = b(keep);
    clear keep cleft cright
   


        
    % Form a matrix for which ba(x,1) is the slope b(x) and ba(x,2) is the
    % y-intercept a(x).  Sort this matrix in ascending order of slope, 
    % breaking ties in slope with the y-intercept.  
    ba = [b, a];
    ba = sortrows(ba,[1,2]);
    a = ba(:,2);
    b = ba(:,1);
    
    % Then, from each pair of indices with the b component equal, remove
    % the one with smaller a component.  This code works because the sort
    % above enforced the condition: if b(i) == b(i+1), then a(i) <= a(i+1).
    keep = [find(diff(b)); length(b)];
    % This previous line is equivalent to:
    % keep = [];
    % for i=[1:length(b)-1]
    %    if b(i)~=b(i+1)
    %        keep = [keep, i];
    %    end
    %end 
    %keep = [keep, length(b)];  % We always keep the last one.
    
    % Note that the elements of keep are in ascending order.
    % This makes it so that b(keep) is still sorted in ascending order.
    a = a(keep);
    b = b(keep);
end
