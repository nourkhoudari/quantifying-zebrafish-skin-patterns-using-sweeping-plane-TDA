%==========================================================================
% Striped and Spotted Binary Image Generation Function
% asociated with the manuscript titled:
%   "Quantifying Topological Features and Irregularities in Zebrafish 
%   Patterns Using the Sweeping-plane Filtration"
%
%   Version 1.0
%   (C) 2025/09/30
%    Nour Khoudari, John T. Nardini, Alexandria Volkening
%
%   A function that loads zebrafish skin pattern data collected from 
%   simulations of an agent-based model [A. Volkening and B. Sandstede, 
%   Iridophores as a source of robustness in zebrafish stripes and 
%   variability in Danio patterns, Nat. Commun., 9 (2018)],
%   to generate binary images based on a user defined pixel size.
%
%   Inputs: 
%   stepX: pixel size used to generate binary images from simulation data
%   dataset: specifies the data set to be loaded: wild-type or pfeffer
%
%   We use two data set of 1000 simulation each, the user can choose between
%   dataset = 1 for wild-type zebrafish skin patterns
%   dataset = 2 for pfeffer mutant zebrafish skin patterns
%   Both data set can be found on Figshare: M. McGuirl, A. Volkening, 
%   and B. Sandstede, Zebrafish simulation data, 2020.
%   https://figshare.com/projects/Zebrafish_simulation_data/72689
%   
%   Outputs:
%   Binary images generated from simulation data
%==========================================================================

function [] = generate_binary_image(stepX,dataset)

%Initialize processing parameters
numSim = 1000;                     % number of random simulations
numDays = 45;                      % includes day 0, which has index 1 (so goes to time = numDays - 1)
focusDays = 45;                    % days (corresponds to days post fertilization or dpf = 20 + focusDays) to analyze; note that 1 is the initial condition
numFocusDays = length(focusDays);  % number of days to consider
stepY = stepX;                     % dimensions of pixel in \mum

%Initialize matrices used in processing
numMelMaster = zeros(numSim,numFocusDays);


% Loop over simulations
for simI = 1:numSim
    
    % Load data
    if dataset == 1
        load(['./WT_rawImageData/Out_WT_default_' num2str(simI) '.mat'], 'cellsM', 'cellsXc', 'cellsXsn', 'cellsIl', 'cellsId', 'numMel',...
        'numIrid', 'numIril', 'numXanc', 'numXansn', 'boundaryX', 'boundaryY')
    elseif dataset == 2
        load(['./pfef_rawImageData/Out_pfef_default_' num2str(simI) '.mat'], 'cellsM', 'cellsXc', 'cellsXsn', 'cellsIl', 'cellsId', 'numMel',...
        'numIrid', 'numIril', 'numXanc', 'numXansn', 'boundaryX', 'boundaryY')
    end

    % Extracting final size of domain
    finalBoundaryY = boundaryY(numDays);    % final domain height (end of simulation)
    finalBoundaryX = boundaryX(numDays);    % final domain length (end of simulation)
    
    % Plot one example pattern (black melanophores only)
    if simI == 1
        figure
        plot(cellsM(1:numMel(45),1,45), cellsM(1:numMel(45),2,45), 'k*')
        axis equal
    end

    % Loop through focus days
    for dayStep = 1:numFocusDays
        
        dayI = focusDays(dayStep);
        
        % Define periodic functions that will be used to calculate cell-cell distances for this domain size
        myfunction = myboundary(dayI,boundaryX,boundaryY,'perX');
        
        %Collect cells in the current time step 
        tempXc = cellsXc(1:numXanc(dayI),:,dayI);
        tempXsn = cellsXsn(1:numXansn(dayI),:,dayI);
        tempId = cellsId(1:numIrid(dayI),:,dayI);
        tempIl = cellsIl(1:numIril(dayI),:,dayI);
        if dataset == 1
            tempM = cellsM(1:numMel(dayI),:,dayI);
        elseif dataset == 2
            tempM = tempIl;
        end
        
        % Shift cells up so that the pattern always appears in the middle of the final domain
        tempM(:,2) = tempM(:,2) + (finalBoundaryY - boundaryY(dayI))/2;
        tempXc(:,2) = tempXc(:,2) + (finalBoundaryY - boundaryY(dayI))/2;
        tempXsn(:,2) = tempXsn(:,2) + (finalBoundaryY - boundaryY(dayI))/2;
        tempId(:,2) = tempId(:,2) + (finalBoundaryY - boundaryY(dayI))/2;
        tempIl(:,2) = tempIl(:,2) + (finalBoundaryY - boundaryY(dayI))/2;

        % Shift cells to the middle of the domain
        tempM(:,1) = tempM(:,1) + (finalBoundaryX - boundaryX(dayI))/2;
        tempXc(:,1) = tempXc(:,1) + (finalBoundaryX - boundaryX(dayI))/2;
        tempXsn(:,1) = tempXsn(:,1)+ (finalBoundaryX - boundaryX(dayI))/2;
        tempId(:,1) = tempId(:,1) + (finalBoundaryX - boundaryX(dayI))/2;
        tempIl(:,1) = tempIl(:,1)+ (finalBoundaryX - boundaryX(dayI))/2;

        % Redefine boundaryY(dayI) to be finalBoundaryY
        boundaryY(dayI) = finalBoundaryY;
        

        % Discretize space and count the number of cells per grid square
        xgrid = 0:stepX:finalBoundaryX;%0:stepX:boundaryX(dayI);%
        ygrid = 0:stepY:finalBoundaryY;%0:stepY:boundaryY(dayI);%
        
        for xi = 2:length(xgrid)
            for yi = 2:length(ygrid)
                meanMM(xi-1, yi-1) = sum(logical((tempM(:,1) >= xgrid(xi-1)).*(tempM(:,1) < xgrid(xi)).*(tempM(:,2) >= ygrid(yi-1)).*(tempM(:,2) < ygrid(yi))));
            end
        end
        
        % Set max to 1 and transe and fill in holes
        meanMM = min(meanMM,1)';
        meanMMlarge = [meanMM, meanMM, meanMM];
        meanMMlargefilled = imfill(meanMMlarge);
        holes = meanMMlargefilled & ~meanMMlarge;
        bigholes = bwareaopen(holes,4);
        smallholes = holes & ~bigholes;
        meanMMlarge = meanMMlarge | smallholes;
        meanMM = meanMMlarge(:,(1+end/3):2*end/3);
       
        % Some image processing to show results and data storage
        sizes(dayStep,:) = [dayStep, size(meanMM)];
        
        if dataset == 1
            imwrite(meanMM, ['./WT_binnedImages/centered_wt_day' num2str(dayI) '_sim' num2str(simI) '_step' num2str(stepX) '.png'], 'PNG')
        elseif dataset == 2
            imwrite(meanMM, ['./WT_binnedImages_spotted/pfef_day' num2str(dayI) '_sim' num2str(simI) '_step' num2str(stepX) '.png'], 'PNG')
        end

        % Clear temporary variables to start again
        clear meanMM;
        
    end     % end of loop over time
    
end         % end of loop over simulations


end         % end of function


%==========================================================================
% Function myboundary defines the distance function used (depends on the
% boundary conditions). We use perX - periodic boundary conditions in the
% x-direction only
function[myfunction] = myboundary(time,boundaryX,boundaryY,option)

if strcmp(option, 'perX') == 1
    myfunction = @(x,y)(sqrt((x(:,2)-y(:,2)).^2 + min((x(:,1)-y(:,1)).^2, (boundaryX(time)-abs(x(:,1)-y(:,1))).^2)));
else
    myfunction = 'euclidean';
end

end
%==========================================================================

 

