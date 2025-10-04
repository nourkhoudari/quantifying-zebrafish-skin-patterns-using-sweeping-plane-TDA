%==========================================================================
% Striped and Spotted Binary Image Classification and Quantification Code
% asociated with the manuscript titled:
%   "Quantifying Topological Features and Irregularities in Zebrafish 
%   Patterns Using the Sweeping-plane Filtration"
%
%   Version 1.0
%   (C) 2025/09/30
%    Nour Khoudari, John T. Nardini, Alexandria Volkening
%
%   Inputs: 
%   persistent homology summaries for all binary images in a data set in
%   all sweeping directions: top(TB), bottom(BT), left(LR), right(RL)
%   for a choosen pixel size (stepX) used to generate those binary images
%
%   We use two data sets of 1000 image each and the user can choose between
%   dataset = 1 for wild-type zebrafish skin patterns
%   dataset = 2 for pfeffer mutant zebrafish skin patterns
%   
%
%   Outputs:
%   excel sheets summarizing all the classifications and quatifications
%   of all images for the loaded data set
%==========================================================================

close all
clear

% choose dataset=1 for striped patterns or dataset=2 for spotted patterns 
dataset = 1;

if dataset == 1
    filename = './WT_binnedImages/centered_wt_day45_sim_plane';
else
    filename = './pfef_binnedImages_spotted/pfef_day45_sim_plane';
end

stepX = 80; % choice of voxel width [\mum] used to generate the binary images

% Define variables to store measurements
sim = cell(1, 1000); % image or simulation number in the data set
class = cell(1, 1000); % classification: 0 for unbroken stripes, 1 for broken stripe(s), 2 for broken interstripe(s), 3 for breaks in both, 4 for spots
stripes = cell(1, 1000); % number of stripes
loops = cell(1, 1000); % number of loops or holes
max_stripes_width = cell(1, 1000); % maximum stripes widths
min_stripes_width = cell(1, 1000); % minimum stripes widths
max_interstripes_width = cell(1, 1000); % maximum interstripes widths
min_interstripes_width = cell(1, 1000); % minimum interstripes widths
bbridge = cell(1, 1000); % number of breaks in stripe(s)
wbridge = cell(1, 1000); % number of breaks in interstripe(s)
bbridge_width = cell(1, 1000); % widths of breaks in stripe(s)
wbridge_width = cell(1, 1000); % widths of breaks in interstripes
bbridge_center = cell(1, 1000); % position of center of breaks in stripe(s) in the x-direction
wbridge_center = cell(1, 1000); % position of center of breaks in interstripe(s) in the x-direction
bbridge_zone = cell(1, 1000); % position of breaks in stripe(s) in the y-direction, i.e. in which stripe 
wbridge_zone = cell(1, 1000); % position of breaks in interstripe(s) in the y-direction, i.e. in which interstripe
spots = cell(1, 1000); % number of spots
max_spot_width = cell(1, 1000); % vertical width of spots
max_spot_length = cell(1, 1000); % horizontal width of spots
rng(1);


for j = 1:1000 %loop over all binarized images 

    for i = stepX  %loop over replicas of voxel widths (if any) 
                       
        % load all stored persistent homology summaries for each image in
        % all sweeping directions: top(TB), bottom(BT), left(LR), right(RL)
        intervals_b_0 = load([filename,'_bottom_dim0_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_b_0), births_b_0 = intervals_b_0(:, 1); else, births_b_0 = [];end
        if ~isempty(intervals_b_0), deaths_b_0 = intervals_b_0(:, 2); else, deaths_b_0 = [];end        
        intervals_b_1 = load([filename,'_bottom_dim1_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_b_1), births_b_1 = intervals_b_1(:, 1); else, births_b_1 = [];end
        if ~isempty(intervals_b_1), deaths_b_1 = intervals_b_1(:, 2); else, deaths_b_1 = [];end            
        intervals_t_0 = load([filename,'_top_dim0_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_t_0), births_t_0 = intervals_t_0(:, 1); else, births_t_0 = [];end
        if ~isempty(intervals_t_0), deaths_t_0 = intervals_t_0(:, 2); else, deaths_t_0 = [];end            
        intervals_t_1 = load([filename,'_top_dim1_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_t_1), births_t_1 = intervals_t_1(:, 1); else, births_t_1 = [];end
        if ~isempty(intervals_t_1), deaths_t_1 = intervals_t_1(:, 2); else, deaths_t_1 = [];end            
        intervals_r_0 = load([filename,'_right_dim0_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_r_0), births_r_0 = intervals_r_0(:, 1); else, births_r_0 = [];end
        if ~isempty(intervals_r_0), deaths_r_0 = intervals_r_0(:, 2); else, deaths_r_0 = [];end    
        intervals_r_1 = load([filename,'_right_dim1_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_r_1), births_r_1 = intervals_r_1(:, 1); else, births_r_1 = [];end
        if ~isempty(intervals_r_1), deaths_r_1 = intervals_r_1(:, 2); else, deaths_r_1 = [];end             
        intervals_l_0 = load([filename,'_left_dim0_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_l_0), births_l_0 = intervals_l_0(:, 1); else, births_l_0 = [];end
        if ~isempty(intervals_l_0), deaths_l_0 = intervals_l_0(:, 2); else, deaths_l_0 = [];end    
        intervals_l_1 = load([filename,'_left_dim1_',num2str(j),'_step',num2str(i),'.txt']);
        if ~isempty(intervals_l_1), births_l_1 = intervals_l_1(:, 1); else, births_l_1 = [];end
        betti_t_0 = load([filename,'_top_Betti_',num2str(j),'_step',num2str(i),'_b0']);
        betti_r_0 = load([filename,'_right_Betti_',num2str(j),'_step',num2str(i),'_b0']);
                
        % Define some conditions used in classification
        any_dim1 = any(isinf(deaths_t_1) | isinf(deaths_b_1)); % at least one loop exists in TB or BT barcodes
        all_persistent_zeroborn_dim0_l = all(isinf(deaths_l_0)) & all(births_l_0==0); 
        all_persistent_zeroborn_dim0_r =all(isinf(deaths_r_0)) & all(births_r_0==0);
        any_persistent_nonzeroborn_dim0_l = any(isinf(deaths_l_0) & (births_l_0>0));
        any_persistent_nonzeroborn_dim0_r = any(isinf(deaths_r_0) & (births_r_0>0));
        any_nonpersistent_zeroborn_dim0_l = any(~isinf(deaths_l_0) & (births_l_0==0));
        any_nonpersistent_zeroborn_dim0_r = any(~isinf(deaths_r_0) & (births_r_0==0));
        any_nonpersistent_nonzeroborn_dim0_l = any(~isinf(deaths_l_0) & (births_l_0>0));
        any_nonpersistent_nonzeroborn_dim0_r = any(~isinf(deaths_r_0) & (births_r_0>0));

        sim{j} = j;
        
        if any_dim1 % striped pattern 
            % classify striped pattern as: unbroken, stripe(s) broken,interstripe(s) broken, or both broken
            if ((all_persistent_zeroborn_dim0_l & all_persistent_zeroborn_dim0_r)...
                | ((length(deaths_t_0(isinf(deaths_t_0)))==length(deaths_t_1(isinf(deaths_t_1)))) & all(any(deaths_t_0(~isinf(deaths_t_0))<= births_t_1(isinf(deaths_t_1))',2 )) & (length(deaths_b_0(isinf(deaths_b_0)))==length(deaths_b_1(isinf(deaths_b_1)))) & all(any(deaths_b_0(~isinf(deaths_b_0))<= births_b_1(isinf(deaths_b_1))',2)) ) )...
                & ~((any_persistent_nonzeroborn_dim0_l | any_persistent_nonzeroborn_dim0_r)   &   (~any_nonpersistent_zeroborn_dim0_l & ~any_nonpersistent_zeroborn_dim0_r))...
                & ~((any_nonpersistent_zeroborn_dim0_l | any_nonpersistent_zeroborn_dim0_r)  &  (~any_persistent_nonzeroborn_dim0_l & ~any(isinf(deaths_r_0) & (births_r_0>0)) & ~any(~isinf(deaths_l_0) & (births_l_0>0)) & ~any_persistent_nonzeroborn_dim0_r))
                % striped pattern with no breaks
                class{j} = 0;
                stripes{j} = max(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0==0))));
                loops{j} = max(length(deaths_t_1(isinf(deaths_t_1))),length(deaths_b_1(isinf(deaths_b_1))));
                bbridge{j} = 0;
                wbridge{j} = 0;
                % compute stripe widths
                max_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_0(isinf(deaths_t_0)),'ascend')+sort(births_b_0(isinf(deaths_b_0)),'descend')));
                min_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_1(isinf(deaths_t_1)),'ascend')+sort(births_b_1(isinf(deaths_b_1)),'descend')));
                % compute interstripe widths
                max_birth_t_0 = sort(births_t_0(isinf(deaths_t_0)),'ascend');
                max_birth_b_0 = (length(betti_t_0(:, 1))-1) - sort(births_b_0(isinf(deaths_b_0)),'descend');
                max_birth_t_1 = sort(births_t_1(isinf(deaths_t_1)),'ascend');
                max_birth_b_1 = (length(betti_t_0(:, 1))-1) - sort(births_b_1(isinf(deaths_b_1)),'descend');
                
                min_interstripes_width{j} = stepX*(max_birth_t_0(2:end)-max_birth_b_0(1:end-1));                
                max_interstripes_width{j} = stepX*(max_birth_t_1(2:end)-max_birth_b_1(1:end-1));
                clear max_birth_b_0; clear max_birth_t_0;clear max_birth_b_1; clear max_birth_t_1; % clear temporary variables
          
            elseif (any_persistent_nonzeroborn_dim0_l | any_persistent_nonzeroborn_dim0_r)...
                    &  (~any_nonpersistent_zeroborn_dim0_l & ~any_nonpersistent_zeroborn_dim0_r)
                % striped pattern with breaks in stripe(s)
                class{j} = 1;
                stripes{j} = max(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0==0))));
                loops{j} = max(length(deaths_t_1(isinf(deaths_t_1))),length(deaths_b_1(isinf(deaths_b_1))));
                bbridge{j} = max(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0>0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0>0))));
                wbridge{j} = 0;
   
                births_l_0_temp = births_l_0(isinf(deaths_l_0)&(births_l_0>0));
                births_r_0_temp = births_r_0(isinf(deaths_r_0)&(births_r_0>0));

                % find the widths of stripe breaks that are on the edges
                if length(births_l_0(isinf(deaths_l_0)&(births_l_0>0)))~=length(births_r_0(isinf(deaths_r_0)&(births_r_0>0)))
    
                    [m,ind] = max([length(births_l_0(isinf(deaths_l_0)&(births_l_0>0))),length(births_r_0(isinf(deaths_r_0)&(births_r_0>0)))]);
                    [n,~] = min([length(births_l_0(isinf(deaths_l_0)&(births_l_0>0))),length(births_r_0(isinf(deaths_r_0)&(births_r_0>0)))]);

                    for k=1:m-n
                        if ind==1
                            bbridge_width{j} = [bbridge_width{j}; -stepX*min(births_l_0_temp)];
                            bbridge_center{j} = [bbridge_center{j}; abs(bbridge_width{j}(end)/2)];
                            births_l_0_temp(births_l_0_temp == min(births_l_0_temp)) = [];
                            
                        elseif ind==2
                            bbridge_width{j} = [bbridge_width{j}; -stepX*min(births_r_0_temp)];
                            bbridge_center{j} = [bbridge_center{j}; stepX*(length(betti_r_0(:,1))-1) - abs(bbridge_width{j}(end)/2)];
                            births_r_0_temp(births_r_0_temp == min(births_r_0_temp)) = [];
                        end
                    end

                end
                
                % find the widths of all other stripe breaks (the ones that are not on the edges)
                if ~isempty(births_r_0_temp)
                    bbridge_width_new = stepX*((length(betti_r_0(:,1))-1) - (sort(births_l_0_temp,'ascend')+sort(births_r_0_temp,'descend')));
                    bbridge_width{j} = [bbridge_width{j}; bbridge_width_new]; 
                    bbridge_center{j} = [bbridge_center{j}; stepX*sort(births_l_0_temp,'ascend')+ abs(bbridge_width_new/2)];
                end      
                
                % stray blue cells in interstripes present and counted as a stripe break, amend break count, and classification (if necessary)
                if any(bbridge_width{j}>0) 
                    bbridge{j} = bbridge{j} - length(bbridge_width{j}(bbridge_width{j}>0));
                    bbridge_width{j} = bbridge_width{j}(bbridge_width{j}<0);
                    bbridge_center{j} = bbridge_center{j}(bbridge_width{j}<0);
                    if bbridge{j}==0 % striped pattern with no breaks, amend classification
                        class{j} = 0;           
                    end
                end
                
                bbridge_width{j} = abs(bbridge_width{j});
                
                clear births_l_0_temp births_r_0_temp m n ind bbridge_width_new % clear temporary variables

                % stripe break width is too big (more than 40% of image length), adjust classification
                if any(bbridge_width{j}>0.4*stepX*(length(betti_t_0(:, 1))-1)) 
               
                    if ~any_nonpersistent_nonzeroborn_dim0_l & ~any_nonpersistent_nonzeroborn_dim0_r % there is a stray blue impurity at the edge of an interstripe counted as a stripe break
                        bbridge{j} = bbridge{j} - 1;
                        if ~isempty(bbridge_width{j}(bbridge_width{j} <= 0.4*stepX*(length(betti_t_0(:, 1))-1)))
                            bbridge_width{j} = bbridge_width{j}(bbridge_width{j} <= 0.4*stepX*(length(betti_t_0(:, 1))-1));
                            bbridge_center{j} = bbridge_center{j}(bbridge_center{j} <= 0.4*stepX*(length(betti_t_0(:, 1))-1));
                        else
                            bbridge_width{j} = [];
                            bbridge_center{j} = [];
                        end
                        if bbridge{j}==0 % striped pattern with no breaks, amend classification and find measurements of stripe and interstripe widths
                            class{j} = 0;
                            max_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_0(isinf(deaths_t_0)),'ascend')+sort(births_b_0(isinf(deaths_b_0)),'descend')));
                            min_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_1(isinf(deaths_t_1)),'ascend')+sort(births_b_1(isinf(deaths_b_1)),'descend')));
                            
                            max_birth_t_0 = sort(births_t_0(isinf(deaths_t_0)),'ascend');
                            max_birth_b_0 = (length(betti_t_0(:, 1))-1) - sort(births_b_0(isinf(deaths_b_0)),'descend');
                            max_birth_t_1 = sort(births_t_1(isinf(deaths_t_1)),'ascend');
                            max_birth_b_1 = (length(betti_t_0(:, 1))-1) - sort(births_b_1(isinf(deaths_b_1)),'descend');
                            
                            min_interstripes_width{j} = stepX*(max_birth_t_0(2:end)-max_birth_b_0(1:end-1));                
                            max_interstripes_width{j} = stepX*(max_birth_t_1(2:end)-max_birth_b_1(1:end-1));
                            clear max_birth_b_0; clear max_birth_t_0;clear max_birth_b_1; clear max_birth_t_1; % clear temporary variables

                        end
                        stripes{j} = stripes{j} - 1;
                    else % there is an uncaptured stripe and interstripe breaks next to each other
                        bbridge{j} = bbridge{j} + 1;
                        wbridge{j} = wbridge{j} + 1;
                        class{j} = 3;
                    end

                else % no change in clasification or stripe break count, find the break zones

                % use kmeans to sort features in TB and BT filtration (based on birth time) into (number of stripes) groups    
                births_t = sort([births_t_0;births_t_1]);
                births_b = sort([births_b_0;births_b_1]);
                [idx_tt, centroids_t] = kmeans(births_t,stripes{j}); 
                [~, sortorder_t] = sort(centroids_t, 'ascend'); idx_t = 0*idx_tt;
                for k = 1 : length(idx_tt), idx_t(k) = find(sortorder_t == idx_tt(k));end
                [idx_bb, centroids_b] = kmeans(births_b,stripes{j}); 
                [~, sortorder_b] = sort(centroids_b, 'ascend'); idx_b = 0*idx_bb;
                for k = 1 : length(idx_bb), idx_b(k) = find(sortorder_b == idx_bb(k));end
                
                % find the group that is missing a persistent feature and match with the stripe number
                for k=1:stripes{j}
                    intersect_t(k) = length(intersect(births_t(idx_t == k),births_t_1));
                    intersect_b(k) = length(intersect(births_b(idx_b == stripes{j}-k+1),births_b_1));
                    if (intersect_t(k)==0 | intersect_b(k)==0) 
                        bbridge_zone{j} = [bbridge_zone{j},k]; 
                    end
                end          
                
                clear intersect_b intersect_t births_t births_b idx_tt centroids_t sortorder_t idx_t idx_bb centroids_b sortorder_b idx_b
                
                end        
    
            elseif (any_nonpersistent_zeroborn_dim0_l | any_nonpersistent_zeroborn_dim0_r)  &  (~any_persistent_nonzeroborn_dim0_l & ~any_persistent_nonzeroborn_dim0_r & ~any_nonpersistent_nonzeroborn_dim0_l & ~any_nonpersistent_nonzeroborn_dim0_r)
                % striped pattern with breaks in interstripe(s)
                class{j} = 2;
                loops{j} = max(length(deaths_t_1(isinf(deaths_t_1))),length(deaths_b_1(isinf(deaths_b_1))));
                wbridge{j} = max(length(deaths_l_0(~isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(~isinf(deaths_r_0) & (births_r_0==0))));
                stripes{j} = wbridge{j} + min(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0==0))));
                bbridge{j} = 0;
                
                if length(births_l_0(~isinf(deaths_l_0)&(births_l_0==0)))==length(births_r_0(~isinf(deaths_r_0)&(births_r_0==0)))
                    wbridge_width{j} = stepX*((length(betti_r_0(:, 1))-1) - (sort(deaths_r_0(~isinf(deaths_r_0)&(births_r_0==0)),'ascend')+sort(deaths_l_0(~isinf(deaths_l_0)&(births_l_0==0)),'descend')));
                    wbridge_center{j} = stepX*sort(deaths_l_0(~isinf(deaths_l_0)&(births_l_0==0)),'ascend')+wbridge_width{j}/2;
                
                else % stray gold cell on edge of an unbroken pattern or a pattern with breaks in interstripes

                    if length(births_l_0(~isinf(deaths_l_0)&(births_l_0==0)))+length(births_r_0(~isinf(deaths_r_0)&(births_r_0==0)))==1 % there is a stray gold cell on the edge of a stripe in an unbroken pattern, reclassify
                        class{j} = 0;
                        max_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_0(isinf(deaths_t_0)),'ascend')+sort(births_b_0(isinf(deaths_b_0)),'descend')));
                        min_stripes_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_1(isinf(deaths_t_1)),'ascend')+sort(births_b_1(isinf(deaths_b_1)),'descend')));
                        
                        max_birth_t_0 = sort(births_t_0(isinf(deaths_t_0)),'ascend');
                        max_birth_b_0 = (length(betti_t_0(:, 1))-1) - sort(births_b_0(isinf(deaths_b_0)),'descend');
                        max_birth_t_1 = sort(births_t_1(isinf(deaths_t_1)),'ascend');
                        max_birth_b_1 = (length(betti_t_0(:, 1))-1) - sort(births_b_1(isinf(deaths_b_1)),'descend');
                        
                        min_interstripes_width{j} = stepX*(max_birth_t_0(2:end)-max_birth_b_0(1:end-1));                
                        max_interstripes_width{j} = stepX*(max_birth_t_1(2:end)-max_birth_b_1(1:end-1));
                        clear max_birth_b_0; clear max_birth_t_0;clear max_birth_b_1; clear max_birth_t_1; % clear temporary variables

                    else % there is a stray gold cell on the edge of a stripe and break in interstripe(s)
                        deaths_l_0_temp = deaths_l_0(~isinf(deaths_l_0)&(births_l_0==0));
                        deaths_r_0_temp = deaths_r_0(~isinf(deaths_r_0)&(births_r_0==0));
                        % remove all the connected components that correspond to all the stray edge gold pixel
                        [m,ind] = max([length(deaths_l_0_temp),length(deaths_r_0_temp)]);
                        [n,~] = min([length(deaths_l_0_temp),length(deaths_r_0_temp)]);
                        for k=1:m-n
                            if ind == 1
                               deaths_l_0_temp(deaths_l_0_temp == min(deaths_l_0_temp)) = [];
                            else
                               deaths_r_0_temp(deaths_r_0_temp == min(deaths_r_0_temp)) = [];
                            end
                        end
        
                        wbridge_width{j} = stepX*((length(betti_r_0(:, 1))-1) - (sort(deaths_r_0_temp,'ascend')+sort(deaths_l_0_temp,'descend')));
                        wbridge_center{j} = stepX*sort(deaths_l_0_temp,'ascend')+wbridge_width{j}/2;

                    end
                end

                % use kmeans to sort features in TB and BT filtration (based on birth time) into (number of stripes) groups    
                births_t = sort([births_t_0;births_t_1]);
                births_b = sort([births_b_0;births_b_1]);
                [idx_tt, centroids_t] = kmeans(births_t,stripes{j}); 
                [~, sortorder_t] = sort(centroids_t, 'ascend'); idx_t = 0*idx_tt;
                for k = 1 : length(idx_tt), idx_t(k) = find(sortorder_t == idx_tt(k));end
                [idx_bb, centroids_b] = kmeans(births_b,stripes{j}); 
                [~, sortorder_b] = sort(centroids_b, 'ascend'); idx_b = 0*idx_bb;
                for k = 1 : length(idx_bb), idx_b(k) = find(sortorder_b == idx_bb(k));end
                
                % find the group that is missing a persistent feature and match with the interstripe number
                for k=1:stripes{j}
                    intersect_t(k) = length(intersect(births_t(idx_t == k),births_t_0(isinf(deaths_t_0))));
                    intersect_b(k) = length(intersect(births_b(idx_b == stripes{j}-k+1),births_b_0(isinf(deaths_b_0))));
                    if k>1 & (intersect_t(k)==0 | intersect_b(k-1)==0) 
                        wbridge_zone{j} = [wbridge_zone{j},k-1];
                    end
                end
        
                clear intersect_b intersect_t births_t births_b idx_tt centroids_t sortorder_t idx_t idx_bb centroids_b sortorder_b idx_b
               
            else % striped pattern with breaks in both stripe(s) and interstripe(s)
                
                class{j} = 3;
                loops{j} = max(length(deaths_t_1(isinf(deaths_t_1))),length(deaths_b_1(isinf(deaths_b_1))));
                wbridge{j} = max(length(deaths_l_0(~isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(~isinf(deaths_r_0) & (births_r_0==0))));
                bbridge{j} = max(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0>0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0>0)))) + max(length(deaths_l_0(~isinf(deaths_l_0) & (births_l_0>0))),length(deaths_r_0(~isinf(deaths_r_0) & (births_r_0>0))));
                stripes{j} = wbridge{j} + min(length(deaths_l_0(isinf(deaths_l_0) & (births_l_0==0))),length(deaths_r_0(isinf(deaths_r_0) & (births_r_0==0))));
                
                if (length(births_l_0(~isinf(deaths_l_0)&(births_l_0==0)))==length(births_r_0(~isinf(deaths_r_0)&(births_r_0==0)))) & (length(births_l_0(isinf(deaths_l_0)&(births_l_0>0)))==length(births_r_0(isinf(deaths_r_0)&(births_r_0>0))))
                    wbridge_width{j} = stepX*((length(betti_r_0(:, 1))-1) - (sort(deaths_r_0(~isinf(deaths_r_0)&(births_r_0==0)),'ascend')+sort(deaths_l_0(~isinf(deaths_l_0)&(births_l_0==0)),'descend')));
                    wbridge_center{j} = stepX*sort(deaths_l_0(~isinf(deaths_l_0)&(births_l_0==0)),'ascend')+wbridge_width{j}/2;
                    bbridge_width{j} = abs(stepX*((length(betti_r_0(:,1))-1) - (sort(births_l_0(isinf(deaths_l_0)&(births_l_0>0)),'ascend')+sort(births_r_0(isinf(deaths_r_0)&(births_r_0>0)),'descend'))));
                    bbridge_center{j} = stepX*sort(births_l_0(isinf(deaths_l_0)&(births_l_0>0)),'ascend')+ bbridge_width{j}/2;
                end

                % use kmeans to sort features in TB and BT filtration (based on birth time) into (number of stripes) groups    
                births_t = sort([births_t_0;births_t_1]);
                births_b = sort([births_b_0;births_b_1]);
                [idx_tt, centroids_t] = kmeans(births_t,stripes{j}); 
                [~, sortorder_t] = sort(centroids_t, 'ascend'); idx_t = 0*idx_tt;
                for k = 1 : length(idx_tt), idx_t(k) = find(sortorder_t == idx_tt(k));end
                [idx_bb, centroids_b] = kmeans(births_b,stripes{j}); 
                [~, sortorder_b] = sort(centroids_b, 'ascend'); idx_b = 0*idx_bb;
                for k = 1 : length(idx_bb), idx_b(k) = find(sortorder_b == idx_bb(k));end
                
                % find the group that is missing a persistent feature and match with the stripe or interstripe number
                for k=1:stripes{j}
                    intersect_t_loops(k) = length(intersect(births_t(idx_t == k),births_t_1));
                    intersect_b_loops(k) = length(intersect(births_b(idx_b == stripes{j}-k+1),births_b_1));
                    intersect_t_connectedcomponent(k) = length(intersect(births_t(idx_t == k),births_t_0(isinf(deaths_t_0))));
                    intersect_b_connectedcomponent(k) = length(intersect(births_b(idx_b == stripes{j}-k+1),births_b_0(isinf(deaths_b_0))));
                    if (intersect_t_loops(k)==0 | intersect_b_loops(k)==0) 
                        bbridge_zone{j} = [bbridge_zone{j},k]; 
                    end  
                    if k>1 & ((intersect_t_connectedcomponent(k)==0 | intersect_b_connectedcomponent(k-1)==0)) 
                        wbridge_zone{j} = [wbridge_zone{j},k-1];
                    end
                end          
                
                clear intersect_b_loops intersect_t_connectedcomponent intersect_b_connectedcomponent intersect_t_loops births_t births_b idx_tt centroids_t sortorder_t idx_t idx_bb centroids_b sortorder_b idx_b
                


            end


        else % spotted pattern
            class{j} = 4;
            spots{j} = max(length(deaths_t_0(isinf(deaths_t_0))),length(deaths_b_0(isinf(deaths_b_0))));        
            max_spot_width{j} = stepX*((length(betti_t_0(:, 1))-1) - (sort(births_t_0(isinf(deaths_t_0)),'ascend')+sort(births_b_0(isinf(deaths_b_0)),'descend')));     
            max_spot_length{j} = stepX*((length(betti_r_0(:, 1))-1) - (sort(births_l_0(isinf(deaths_l_0)),'ascend')+sort(births_r_0(isinf(deaths_r_0)),'descend')));
  
        end

    end

end

% Save classifications and quatifications in excel sheet
if dataset == 1

    A = [sim(:),class(:),stripes(:),max_stripes_width(:),min_stripes_width(:), max_interstripes_width(:), min_interstripes_width(:), loops(:),bbridge(:),bbridge_zone(:),bbridge_center(:),bbridge_width(:),wbridge(:),wbridge_zone(:),wbridge_center(:),wbridge_width(:)];
    T = cell2table(A,'VariableNames',{'Simulation','Classification','#Stripes','max_stripes_width','min_stripes_width','max_interstripes_width','min_interstripes_width','#Loops','#black bridges', 'black_bridge_zone','black_bridge_center','black_bridge_width','#white bridges','white_bridge_zone','white_bridge_center','white_bridge_width'});
    writetable(T,'output_classification_images_stripes.xlsx')

else

    A = [sim(:),spots(:),max_spot_width(:),max_spot_length(:)];
    T = cell2table(A,'VariableNames',{'Simulation','#Spots','max_spot_width','max_spot_length'});
    writetable(T,'classification_images_spotted.xlsx')

end
