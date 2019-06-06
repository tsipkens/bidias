% function [tferDT, roughness, SIGMA] = TwomeyCorrect(tferDT,Zp,DMAKernel, DTData, dn_dzp)
function [distribution, R_dp, R_mp, SIGMA, y] = Twomey(distribution, Gtot1, Gtot2, dp, mp, response, LIMIT, error, SIGMA)

%% Setup Parameters

[dpStarLen, mpStarLen] = size(response);
[ll, mm, factorGen, factorInit, intermediateLim, maxTwomeySteps, ~] = Old.InversionParameters;

y = zeros(dpStarLen, mpStarLen);
y_1 = zeros(dpStarLen, mpStarLen);
y_2 = zeros(dpStarLen, mpStarLen);


%% Algorithm

factor = factorInit;

for k = 2:intermediateLim
  

    for i = 1:dpStarLen
        for j = 1:mpStarLen

            G1 = Gtot1{i,j}.*distribution;
            G2 = Gtot2{i,j}.*distribution;
            y_1(i,j) = trapz(mp,trapz(dp,G1,1)); % totalInt
            y_2(i,j) = trapz(mp,trapz(dp,G2,1));

        end
    end
    y = y_1 + y_2;

    
    % Calc weight factor
    X = response./y; 
    
    distribution0 = distribution;

    %% NEW Transfer Function Calc
    count = 0;

    for i = 1:dpStarLen % Twomey loop
        for j = 1:mpStarLen

            if abs(y(i,j) - response(i,j)) > error(i,j)   % this is used so as not to operate on measured values of 0
                                                                                      % this value can be changed to be larger or smaller in # cm^-3
                                                                                      % to ignore low counts
                
                GtotBoth = Gtot1{i,j} + Gtot2{i,j};

                isRelevant = GtotBoth > 1E-8;
                GtotBoth = GtotBoth.*isRelevant;
                distribution = distribution.*(1 + (X(i,j)-1).*GtotBoth/factor);
                count = count + 1;
                
            end
        end
    end
    
    % SIGMA Calculation
    SIGMA_new = 0;

    numErrorCalcs = 0;
    for i  = 1:dpStarLen
        for j = 1:mpStarLen
            if response(i,j) > 0.01
                SIGMA_new = SIGMA_new + ((response(i,j) - y(i,j))/error(i,j))^2;
                numErrorCalcs = numErrorCalcs + 1;
            end
        end
    end
    
    SIGMA_new = SIGMA_new/numErrorCalcs;

    SIGMA = [SIGMA,SIGMA_new];

    if SIGMA_new < LIMIT
        display(k)
        display(SIGMA_new)
        distribution = distribution0;
        break
    end
end


%% Enter into WHILE Loop allowing for exit if certain error conditions are met

factor = factorGen;

factorOrig = factorGen;

SIGMA0 = SIGMA(k)*1.5;

if k == intermediateLim

    while SIGMA(end) > LIMIT && maxTwomeySteps > k
        if SIGMA(end) - SIGMA(end-1) > 0
            factor = factor*2;
        else
            factor = factorOrig;
        end

        for i = 1:dpStarLen
            for j = 1:mpStarLen
                
                G1 = Gtot1{i,j}.*distribution;
                G2 = Gtot2{i,j}.*distribution;
                y_1(i,j) = trapz(mp,trapz(dp,G1,1)); % totalInt
                y_2(i,j) = trapz(mp,trapz(dp,G2,1));

            end
        end
        y = y_1 + y_2;
        
        distribution0 = distribution;

        % Calc weight factor
        X = response./y;  
        count = 0;

        for i = 1:dpStarLen
            for j = 1:mpStarLen

                if abs(y(i,j) - response(i,j)) > error(i,j)

                    Gtot1_ij = Gtot1{i,j};
                    Gtot2_ij = Gtot2{i,j};

                    isRelevant = Gtot1_ij + Gtot2_ij > 1E-8; 
                    Gtot1_ij = Gtot1_ij.*isRelevant;
                    Gtot2_ij = Gtot2_ij.*isRelevant;

                    distribution = distribution.*(1 + (X(i,j)-1).*(Gtot1_ij + Gtot2_ij)/factor);
                    count = count + 1;
                end
            end
        end

        k = k + 1;

        SIGMA_new = 0;

        numErrorCalcs = 0;
        for i  = 1:dpStarLen
            for j = 1:mpStarLen
                if response(i,j) > 0.01
                    SIGMA_new = SIGMA_new + ((response(i,j) - y(i,j))/error(i,j))^2;
                    numErrorCalcs = numErrorCalcs + 1;
                end
            end
        end

        SIGMA_new = SIGMA_new/numErrorCalcs;

        SIGMA = [SIGMA,SIGMA_new];

    end
end

if k > intermediateLim
    display(k)
    display(SIGMA_new)
    distribution = distribution0;
end
    
%% Roughness Calculation

R_dp_test = sum(sum(abs(distribution(2:ll-1,3:mm) + distribution(2:ll-1,1:mm-2) - 2*distribution(2:ll-1,2:mm-1)),2)/(mm-2))/(ll-2);

R_mp_test = sum(sum(abs(distribution(3:ll,2:mm-1) + distribution(1:ll-2,2:mm-1) - 2*distribution(2:ll-1,2:mm-1)),1)/(ll-2))/(mm-2);


R_dp = 0;
R_mp = 0;

for ii = 2:ll-1
    R_dp_inner = 0;
    for jj = 2:mm-1
        R_dp_inner = R_dp_inner + abs(distribution(ii,jj+1) + distribution(ii,jj-1) - 2*distribution(ii,jj));
    end
    R_dp_inner = R_dp_inner/(mm-2);
    R_dp = R_dp + R_dp_inner;
end
R_dp = R_dp/(ll-2);

for jj = 2:mm-1
    R_mp_inner = 0;
    for ii = 2:ll-1
        R_mp_inner = R_mp_inner + abs(distribution(ii+1,jj) + distribution(ii-1,jj) - 2*distribution(ii,jj));
    end
    R_mp_inner = R_mp_inner/(ll-2);
    R_mp = R_mp + R_mp_inner;
end

R_mp = R_mp/(mm-2);