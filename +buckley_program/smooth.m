function [smoothAnswer, SIGMA] = smooth(distribution, iter, response, Gtot1, Gtot2, dp, mp, SIGMA, error)
% inputs to this function are the distribution itself (matrix), a maximum
% number of smoothing steps (iter), the measured concentration matrix
% (resposne), the Kernel functions for both charge states, which include
% the DMA and APM transfer functions (Gtot1 and Gtot2, cells containing a
% matrix for each measurement channel combination), arrays for mobility
% diameter and mass description (dp and mp), chi-square values (SIGMA), and
% the error matrix (error)

[ll, mm, ~, ~, ~, ~, smoothFract] = Old.InversionParameters;
[dpStarLen, mpStarLen] = size(response);
smoothCount = 0;
y_1 = zeros(dpStarLen, mpStarLen);
y_2 = zeros(dpStarLen, mpStarLen);

SIGMA_new = SIGMA(end);
LIMIT = 1;
total = smoothFract*8 + 0.5;

SIGMA_0 = SIGMA(end);

if SIGMA_new < LIMIT
    while smoothCount < iter && SIGMA_new < LIMIT
        
        distribution(2:ll-1,2:mm-1) = 0.5*distribution(2:ll-1,2:mm-1) + smoothFract*distribution(3:ll,3:mm) + smoothFract*distribution(1:ll-2,1:mm-2) + ...
            smoothFract*distribution(3:ll,1:mm-2) + smoothFract*distribution(1:ll-2,3:mm) + smoothFract*distribution(1:ll-2,2:mm-1) + smoothFract*distribution(3:ll,2:mm-1)+...
            smoothFract*distribution(2:ll-1,1:mm-2) + smoothFract*distribution(2:ll-1,3:mm);
        distribution(2:ll-1,2:mm-1) = distribution(2:ll-1,2:mm-1)/total;
       
        distribution(1,:) = 3*distribution(1,:)/4 + distribution(2,:)/4;
        distribution(:,1) = 3*distribution(:,1)/4 + distribution(:,2)/4;
        distribution(end,:) = 0.25*(distribution(dpStarLen-1,:) + 3*distribution(dpStarLen,:)/4);
        distribution(:,end) = 0.25*(distribution(:,mpStarLen-1) + 3*distribution(:,mpStarLen)/4);
        
        smoothCount = smoothCount + 1;
        
        if 1 - SIGMA_new < 0.10 || rem(smoothCount,2) == 0 

            for i = 1:dpStarLen
                for j = 1:mpStarLen

                        G1 = Gtot1{i,j}.*distribution;
                        G2 = Gtot2{i,j}.*distribution;
                        y_1(i,j) = trapz(mp,trapz(dp,G1,1)); % totalInt
                        y_2(i,j) = trapz(mp,trapz(dp,G2,1));

                end
            end
            y = y_1 + y_2;

            SIGMA_new = 0;
            numErrorCalcs = 0;

            for i  = 1:dpStarLen
                for j = 1:mpStarLen
                    if response(i,j) > 0.01
                        SIGMA_new = SIGMA_new + ((response(i,j) - y(i,j))/(error(i,j)))^2;
                        numErrorCalcs = numErrorCalcs + 1;
                    end
                end
            end
            SIGMA_new = SIGMA_new/numErrorCalcs;
        end
        
        SIGMA = [SIGMA,SIGMA_new];

    end
else
    while smoothCount < 2 || SIGMA_new == 2*SIGMA_0

        distribution(2:ll-1,2:mm-1) = 0.5*distribution(2:ll-1,2:mm-1) + smoothFract*distribution(3:ll,3:mm) + smoothFract*distribution(1:ll-2,1:mm-2) + ...
            smoothFract*distribution(3:ll,1:mm-2) + smoothFract*distribution(1:ll-2,3:mm) + smoothFract*distribution(1:ll-2,2:mm-1) + smoothFract*distribution(3:ll,2:mm-1)+...
            smoothFract*distribution(2:ll-1,1:mm-2) + smoothFract*distribution(2:ll-1,3:mm);
        distribution(2:ll-1,2:mm-1) = distribution(2:ll-1,2:mm-1)/total;
            
        distribution(1,:) = 3*distribution(1,:)/4 + distribution(2,:)/4;
        distribution(:,1) = 3*distribution(:,1)/4 + distribution(:,2)/4;
        distribution(end,:) = 0.25*(distribution(dpStarLen-1,:) + 3*distribution(dpStarLen,:)/4);
        distribution(:,end) = 0.25*(distribution(:,mpStarLen-1) + 3*distribution(:,mpStarLen)/4);
        smoothCount = smoothCount + 1;
        
            %calc expected measurement
            for i = 1:dpStarLen
                for j = 1:mpStarLen

                        G1 = Gtot1{i,j}.*distribution;
                        G2 = Gtot2{i,j}.*distribution;
                        y_1(i,j) = trapz(mp,trapz(dp,G1,1)); % totalInt
                        y_2(i,j) = trapz(mp,trapz(dp,G2,1));

                end
            end
            y = y_1 + y_2;
        
            %calc chi-squared (SIGMA)
            SIGMA_new = 0;
            numErrorCalcs = 0;
            for i  = 1:dpStarLen
                for j = 1:mpStarLen
                    if response(i,j) > 0.01
                        SIGMA_new = SIGMA_new + ((response(i,j) - y(i,j))/(error(i,j)))^2;
                        numErrorCalcs = numErrorCalcs + 1;
                    end
                end
            end
            SIGMA_new = SIGMA_new/numErrorCalcs;
            SIGMA = [SIGMA,SIGMA_new];

    end
end

smoothAnswer = distribution;