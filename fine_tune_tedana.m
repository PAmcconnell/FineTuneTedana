% CJL; cjl2007@med.cornell.edu
% 4.26.2019

%   Fine-tuning designed for manual rejection and/or 
%   acceptance of ICA components after an initial run of Tedana. 
%   note: steps 1-2 consider temporal characteristics of the component
%   time-series , steps 3-4 consider spatial information. 
%   
%   Note: the last part of this script creates a movie summarizing the 
%   overall denoising performance and each individual 
%   component classification.
%
%   Disclaimer: just a first attempt; some steps may be prove to be
%   more useful than others. 

%% variables to define

rp = []; % this is a path to a head position traces [assumes they are in the format that is produced by FSL's MCFLIRT; nTR x dimension (x,y,z,yaw,pitch,roll)]
tdir = []; % origin tedana directory; 
mdir = []; % mask directory; assumes some tissues masks with specific naming conventions are contained within (look down at lines 132-139); see prep_fine_tune.sh for example of how these masks can be made 

% define densities
densities = 90:1:99; % This is the range of thresholds that we have used in the past to match FC maps to this set of functional brain network templates. Seems to work well here too. 

% define the artifact template names; used for notes
artifact_names = {'veins','brain_edge','csf'};

% load network templates file;
t = niftiRead([tdir '/fine_tune/network_templates.nii.gz']); % see prep_fine_tune.sh to see how this file is created; other templates could be used

% network labels; order corresponding to 4th dim. in variable "t"
networks = {'DMN','V2','FP','V1','DAN','Pmot','VAN','SAL',...
    'CON','SMh','SMfa','AUD','aMTL','pMTL','PMN','CAN','SMfo'};

% load ICA components classifications; pre fine-tuning 
ds = dataset('File',[tdir '/comp_table_ica.txt']); 

% preallocate some new columns;
ds.networks = cell(size(ds,1),1); % network matches
ds.man_acc = zeros(size(ds,1),1); % manually accepted 
ds.man_rej = zeros(size(ds,1),1); % manually rejected 
ds.notes = cell(size(ds,1),1); % fine tuning notes; rationale for any manual rejection/acceptance will be listed here 

% extract original classifications;
orig_classifications = ds.classification; % preserve original classifications; 

%% Step 1) evaluate relationships between component time-course & head motion traces

rp = load(rp); % load head motion rps.; 

% load ICA component time-courses
c_ts = load([tdir '/meica_mix.1D']);

% normalize
% head motion
for i = 1:size(rp,2)
    rp(:,i) = normalize01(detrend(rp(:,i),'linear','constant'));
end

% normalize
% component
% time-courses
for i = 1:size(c_ts,2)
    c_ts(:,i) = normalize01(detrend(c_ts(:,i),'linear','constant'));
end

% note: component time-courses and head position traces are normalized 
% (range 0-1) and linearly detrended. This makes the traces easier to plot 
% and visualize later on. 

% calc. correlations
for i = 1:size(c_ts,2)
    for ii = 1:size(rp,2)
        [mot_r(i,ii),mot_p(i,ii)] = corr(c_ts(:,i),rp(:,ii)); % mot_r = Pearson correlations; mot_p = p-values
    end
end

% correct for multiple comparisons. 
% note: this step does not influence which components are rejected.
% This is because we are identifying the max correlations, and not every sig. correlation. 
% We set non-sig.   correlations to zero here because later on, we will overlay head position 
% traces over component time-courses, but only if they were sig.
mot_r(mot_p>(0.01 / size(c_ts,2))) = 0;
mot_r = abs(mot_r); % absolute correlation

% find max
% abs corrs.
for i = 1:6
    mot_comp(i) = find(mot_r(:,i)==max(mot_r(:,i)));
end

% find "motion" components
mot_comp = unique(nonzeros(mot_comp)); 

% if motion
% comps. were found... 
if ~isempty(mot_comp)
    
    % sweep motion components
    for i = 1:length(mot_comp)
        
        % mark component as (manually) rejected;
        if ~strcmp(ds.classification{mot_comp(i)},'rejected')
            ds.classification{mot_comp(i)} = 'rejected'; % change classification 
            ds.man_rej(mot_comp(i)) = 1; % mark as manually rejected 
        end
        
        % cite motion criteria;
        ds.notes{mot_comp(i)} = 'mot_corr'; % log reason for manual rejection 
        
    end
    
end

%% Step 2) evaluate component frequency content 

% sweep ica components
for i = 1:size(c_ts,2)
    
    % calculate power spectra
    [a,f] = calc_ps(c_ts(:,i),TR); % a=amplitude; f=frequency

    % if high freq. > low freq.
    if max(a(f>.1)) > max(a(f<.1)) % probably can find a more elegant way of doing this; for example, using Matlab's findpeaks.m function; this seems to work well in my data 
        
        % mark component as (manually) rejected;
        if strcmp(ds.classification{i},'accepted')
            ds.classification{i} = 'rejected';
            ds.man_rej(i) = 1;
        end
        
        % note reason for manual exclusion
        ds.notes{i}= ['high_freq;'...
            ds.notes{i}];
        
    end
    
end

%% Step 3) Evaluate component overlap with artifact templates 

% unzip; if needed
if exist([tdir '/betas_OC.nii.gz'],'file')
    system(['gunzip ' tdir '/betas_OC.nii.gz']); % the niftiRead function complains sometimes about loading larger files that are gzipped; 
end

% load component betas
c = niftiRead([tdir '/betas_OC.nii']); 

% load tissue masks 
wm = niftiRead([mdir '/white_ero_0.nii.gz']); % white matter
csf = niftiRead([mdir '/ventricles_ero_0.nii.gz']); % ventricle pulsation; also seems to catch subependymal/transmedullary veins.
cr = niftiRead([mdir '/grey_ero_0.nii.gz']); % cortical ribbon
gm = niftiRead([mdir '/grey_subcort+cerebellum_ero_0.nii.gz']); % grey matter + subcortical structures + cerebellum

% load artifact templates 
veins = niftiRead([mdir '/veins.nii.gz']); % veins
brain_edge = niftiRead([mdir '/brain_edge.nii.gz']); % brain edge (proxy for head motion); 

% adjust artifact templates; 
idx_a = find(gm.data==1 | veins.data==1);
idx_b = find(gm.data==1 | brain_edge.data==1);
brain_edge.data(idx_a)=0; % remove grey matter & veins voxels 
veins.data(idx_b) = 0; % remove grey matter & brain edge voxels 
clear idx_a idx_b ; % clear some variables 

% preallocate artifact similarities (a_s) matrix 
a_s = zeros(size(c.data,4),2,length(densities));  

% sweep artifacts
for i = 1:size(c.data,4)
    
    % sweep densities
    for d = 1:length(densities)
        a = c.data(:,:,:,i); % tedana ica component "i"
        gm_j = jaccard(a(:)>=prctile(nonzeros(a(:)),densities(d)),gm.data(:)==1); % grey matter control; the idea is that we want to reject components that are more strongly represented in either veins, edge, or csf more than gm
        a_s(i,1,d) = jaccard(a(:)>=prctile(nonzeros(a(:)),densities(d)),veins.data(:)==1) - gm_j; % veins 
        a_s(i,2,d) = jaccard(a(:)>=prctile(nonzeros(a(:)),densities(d)),brain_edge.data(:)==1) - gm_j; % brain edge
        a_s(i,3,d) = jaccard(a(:)>=prctile(nonzeros(a(:)),densities(d)),csf.data(:)==1) - gm_j; % csf 
    end
    
end

% average densities
a_s = mean(a_s,3);

% sweep artifact templates
for i = 1:length(artifact_names)
    
    % artifact i components
    idx = find(a_s(:,i)>0);
    
    % check
    % for matches
    if ~isempty(idx)
        
        % sweep through indices
        for ii = 1:length(idx)
            
            % mark component as (manually) rejected; if needed
            if strcmp(ds.classification{idx(ii)},'accepted')
                ds.classification{idx(ii)} = 'rejected';
                ds.man_rej(idx(ii)) = 1;
            end
            
            % note reason for manual exclusion
            ds.notes{idx(ii)} = [artifact_names{i} ';'...
                ds.notes{idx(ii)}];
            
        end
        
    end
    
end

a_s(a_s < 0) = 0; % set "good" (gm > artifact) components to zero; 
% similar to when we set non-sig. head motion correlations to zero, this is
% done to emphasize rejected components later in the visualization steps. 

% tissue overlap - Not used currently for screening out bad/good components, 
% but useful for visualization and QA. 

% preallocate tissue overlap
t_p = zeros(size(c.data,4),3,...
    length(densities));

% sweep ica components
for i = 1:size(c.data,4)
    
    % sweep densities
    for d = 1:length(densities)
        comp_i = c.data(:,:,:,i); % ica component "i"
        idx = find(csf.data(:)==1 | wm.data(:)==1 | gm.data(:)==1); % indices of csf, wm, and gm voxels
        t_p(i,1,d) = mean(ismember(find(comp_i(idx)>=prctile(nonzeros(comp_i(idx)),densities(d))),find(csf.data(idx)==1))); % csf
        t_p(i,2,d) = mean(ismember(find(comp_i(idx)>=prctile(nonzeros(comp_i(idx)),densities(d))),find(wm.data(idx)==1))); % wm
        t_p(i,3,d) = mean(ismember(find(comp_i(idx)>=prctile(nonzeros(comp_i(idx)),densities(d))),find(gm.data(idx)==1))); % grey matter
    end
    
end

% average densities
t_p = mean(t_p,3) * 100; % convert to 0-100 range

%% step 4) Evaluate component overlap with brain network templates 
% calculate the overlap between components mapped thresholded using 
% the specified percentile and template networks. Components that are
% succesfully matched to a brain network are manually accepted; regardless
% of any manual rejection made to this point. Idea is to avoid rejecting
% components containing signal of interest. 

% preallocate matrix of network similarities
n_s = zeros(size(c.data,4),size(t.data,4),...
    length(densities));

% sweep ica components
for i = 1:size(c.data,4)
    
    % sweep through templates
    for ii = 1:size(t.data,4)
        
        % sweep through densities
        for d = 1:length(densities)
            a = c.data(:,:,:,i); % tedana ica component "i"
            b = t.data(:,:,:,ii); % network template "ii"
            n_s(i,ii,d) = jaccard(a(cr.data==1)>=prctile(nonzeros(a(cr.data==1)),densities(d)),...
                b(cr.data==1)>=prctile(nonzeros(b(cr.data==1)),densities(d))); % calc. network similarity
        end
        
    end
    
end

n_s = mean(n_s,3); % average densities
n_s(n_s<0.1) = 0; % "good" matches only

% sweep components
for i = 1:size(ds,1)
    
    % if match exists
    if sum(n_s(i,:))~=0
        
        % if rejected previously & rationale & kappa > rho;
        if ~strcmp(ds.classification(i),'accepted') && ds.kappa(i) > ds.rho(i)
            
            % mark component as (manually) accepted
            ds.classification{i} = 'accepted';
            ds.man_acc(i) = 1;
            ds.man_rej(i) = 0;
            
        end
        
        % index of sig. matches
        idx = find(n_s(i,:)~=0);
        
        % sort idx by match quality
        [~,b] = sort(n_s(i,idx),'Descend');
        idx = idx(b); % apply reordering
        
        nets = []; % preallocate
        
        % sweep network matches
        for ii = 1:length(idx)
            if ii == 1
                nets = [nets networks{idx(ii)}];
            else
                nets = [nets ';' networks{idx(ii)}];
            end
        end
        
        % log network identities
        ds.networks{i} = nets;
        
    end
    
end

% write out updated component classifications;
export(ds,'file',[tdir '/comp_table_ica.txt']);

% at this point "fine-tuning" is complete and tedana can be re-run using the
% updated comp_table_ica.txt file; the steps below are purely for
% visualization of component classifications; good for quality control. 

%% movie creation

% read some anatomicals for subplot 1
T1 = niftiRead([tdir '/fine_tune/T1_func.nii.gz']); % anatomical in functional space
T1_brain = niftiRead([tdir '/fine_tune/T1_bet_func.nii.gz']); % brain extracted anatomical in functional space
c_s = niftiRead([tdir '/fine_tune/betas_OC_s2.55.nii.gz']); % tedana components; smoothed for visualization ; this needs to be done prior; can overlay non-smoothed data, but this looks better, I think 

% preallocate 
top = size(T1_brain.data,3);
b = 0;

% find top of brain
while b < mean(nonzeros(T1_brain.data(:)))
    b = max(T1_brain.data(:,:,top))';
    top = top - 1;
end

% preallocate 
bottom = 1;
b = 0;

% find bottom of brain
while b < mean(nonzeros(T1_brain.data(:)))
    b = max(T1_brain.data(:,:,bottom))';
    bottom = bottom + 1;
end

% define 
slices = round(linspace(bottom*2,top,16)); % 16 slices seems to fit this figure size best; can automate in future (i.e., take into consideration figure size and relative dims. of anatomical image)

% preallocate 
lh_edge = 1;
b = 0;

% find bottom of brain
while b < mean(nonzeros(T1.data(:)))
    b = max(squeeze(T1.data(lh_edge,:,:)));
    lh_edge = lh_edge + 1;
end

% preallocate 
rh_edge = size(T1.data,1);
b = 0;

% find bottom of brain
while b < mean(nonzeros(T1.data(:)))
    b = max(squeeze(T1.data(rh_edge,:,:)));
    rh_edge = rh_edge - 1;
end

% add some padding
lh_edge = lh_edge - 2; % left hemisphere "edge"
rh_edge = rh_edge + 2; % right hemisphere "edge" 

close all; % close all windows; just in case. 

% sweep through comps.
for i = 1:size(c_s.data,4)
    
    anat = []; 
    func = [];
    
    % sweep slices
    for ii = 1:length(slices)
        anat = [anat T1.data(lh_edge:rh_edge,:,slices(ii))'];
        func = [func c_s.data(lh_edge:rh_edge,:,slices(ii),i)'];
    end

    func(func==0 | anat==0) = nan; % facilitates overlaying func. on anat.
    func = normalize01(func); % set range to 0-1
    func(func < prctile(nonzeros(func(:)),95)) = nan; % 95 %tile represents the "middle" density;  
    imoverlay(flip(anat,1),flip(func,1),[],[],'hot',.9); % overlay thresholded comp. on anatomical
    saveas(gcf,[tdir '/fine_tune/imgs/comp' num2str(i-1) '_axial.jpg']); % save image 
    close; % close image 

    H = figure; % prellocate parent figure
    set(H,'position',[1 596 1680 359]);
    [~] = evalc('hold');
    set(H,'Color','w');
    
    % sub plot
    h = subplot(2,6,1:6);
    img = imread([tdir '/fine_tune/imgs/comp' num2str(i-1) '_axial.jpg']); imagesc(img); % display image 
    title(['ICA Component ' num2str(i-1) ' (' num2str(round(ds.normalizedVarianceExplained(i),2)) '% of total variance explained); Kappa: ' num2str(round(ds.kappa(i),2)) ' & Rho: ' num2str(round(ds.rho(i),2))]);
    set(gca,'FontName','Arial','FontSize',12);
    axis 'off';
    
    ax = h; % minimize white space
    ax.Position = [0.0067 0.5338 0.9851 0.3842]; 
    set(gca,'Yticklabel',[],'Xticklabel',[],'TickLength',[0 0]);
    
    h = subplot(2,6,7); % second sub-plot; network templates
    bar(n_s(i,:),'Facecolor',[.5 .5 .5],'Edgecolor',[0 0 0]);
    [~] = evalc('hold');
    box 'off';
    
    % set figure properties; set x axis limits to reasonable range
    set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],'Xticklabel',...
        networks,'XTick',1:17,'XTickLabelRotation',90,'TickLength',[0 0]);
    xlim([0 18]);
    
    % format y axis and title
    ylabel('Overlap (Jaccard)');
    title('Network Templates');
    
    ax = h; % min. white space
    ax.Position = [0.0350 0.1100 0.2079 0.3412];
    
    % if no match;
    if sum(n_s(i,:))==0
        ylim([0 1]); % change ylim; 
    end
    
    h = subplot(2,6,8); % third sub-plot; tissue overlap
    bar(t_p(i,:),'Facecolor',[.5 .5 .5],'Edgecolor',[0 0 0]);
    [~] = evalc('hold');
    box 'off';
    
    % set format figure;
    set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],...
        'Xtick',1:size(t_p,2),'Xticklabel',{'CSF','WM','GM'},...
        'XTickLabelRotation',45,'TickLength',[0 0]);
    title('Tissue Templates');
    
    % assign axis labels
    ylabel('% in Tissue Template');
    xlim([0 size(t_p,2)+1]);
    
    ax = h; % min. white space
    ax.Position = [0.28 0.1100 0.0928 0.3376];
    
    h = subplot(2,6,9); % fourth sub-plot; artifact templates
    barh(a_s(i,:),'Facecolor',[.5 .5 .5],'Edgecolor',[0 0 0]);
    [~] = evalc('hold');
    box 'off';
    
    % set format figure;
    set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],'Ytick',1:size(a_s,2),...
        'Yticklabel',{'Veins','Edge','CSF'},'TickLength',[0 0]);
    
    % assign axis/title labels 
    xlabel('Overlap (Jaccard) > GM');
    title('Artifact Templates');
    
    % if no match;
    if sum(a_s(i,:)>0)==0
        xlim([0 .1]); % change ylim
    end
    
    ax = h; % min. white space
    ax.Position = [0.4041 0.1100 0.0928 0.3376];
    
    h = subplot(2,6,10); % fifth sub-plot; power spectrum
    [a,f] = calc_ps(c_ts(:,i),TR);
    plot(f,a,'LineWidth',.6,'Color',[.5 .5 .5]);
    [~] = evalc('hold');
    box 'off';
    
    % format the axis;
    set(gca,'FontName','Arial','FontSize',12,...
        'TickLength',[0 0],'TickLength',[0 0]);
    
    % assign axis/title labels 
    xlabel('Frequency [Hz]');
    ylabel('Amplitude');
    title('Power Spectrum');
    
    ax = h; % min. white space
    ax.Position = [0.54 0.1100 0.1123 0.3376];
    
    % sixth sub-plot;
    h = subplot(2,6,11); % head motion correlations
    [~] = evalc('hold');
    
    % find sig. correlations
    idx = find(mot_r(i,:)>0);
    
    % sweep position traces
    for r = 1:length(idx)
        ta = nan(1,6); % temp array 
        ta(idx(r)) = mot_r(i,idx(r));
        bar(ta,'Facecolor',[.5 .5 .5],'Edgecolor',[0 0 0]);
    end

    % remove axes; set y & x limits
    box 'off'; ylim([0 1]); xlim([0 7]);
    
    ax = h; % min. white space
    ax.Position = [0.6932 0.1100 0.1028 0.3429];
    
    % set format figure;
    set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],...
        'Xtick',1:6,'Xticklabel',{'X','Y','Z','y','p','r'},...
        'XTickLabelRotation',45,'TickLength',[0 0]);
    
    % set axis and title labels
    ylabel('Absolute Correlation');
    xlabel('Direction');
    title('Head Motion');
    
    % final subplot;
    h = subplot(2,6,12); % component time-courses & head motion
    
    % plot component time-course;
    if strcmp(ds.classification{i},'accepted')
        if ds.man_acc(i)==1
            plot(moving(c_ts(:,i),3),'LineWidth',.6,'Color',[0 .75 0]); % bright green for manually accepted
        else
            plot(moving(c_ts(:,i),3),'LineWidth',.6,'Color',[0 .5 0]); % green for accepted
        end
    end
    
    % plot component time-course;
    if strcmp(ds.classification{i},'rejected')
        if ds.man_rej(i)==1
            plot(moving(c_ts(:,i),3),'LineWidth',.6,'Color',[.75 0 0]); % bright red for manually rejected
        else
            plot(moving(c_ts(:,i),3),'LineWidth',.6,'Color',[.5 0 0]); % red for rejected
        end
    end
    
    % plot component time-course;
    if strcmp(ds.classification{i},'ignored')
        plot(moving(c_ts(:,i),3),'LineWidth',.6,'Color',[.5 .5 .5]); % gray for ignored
    end
    
    % hold current plot;
    [~] = evalc('hold'); 
    
    % set format figure;
    set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0])
    h.YAxis.Visible='off'; % turn off Y axis
    xlim([0 size(c_ts,1)]);
    xlabel('TRs');
    ylim([0 3]);
    box 'off';
    
    % find sig. correlations
    idx = find(mot_r(i,:)>0);
    
    % sweep position traces
    for r = 1:length(idx)
        l = plot(moving(rp(:,idx(r)),3)+(r/6)+.5,'LineWidth',.6,'Color',[0.5 0.5 0.5]);
        l.Color(4) = mot_r(i,idx(r)); % note: line transparency will scale with correlation strength; stronger correlations are darker / more visible
    end

    ax = h; % min. white space
    ax.Position = [0.8027 0.1100 0.175 0.3429];
    
    % log frame for movie
    m(i+1) = getframe(H);
    close all;
    
end

% now create a summary figure that serves as a "cover sheet" for 
% all the individual frames ; summarizing overall denoising 

% preallocate a
% classifications index
c_idx = ones(size(ds,1),1)*3; % add one to account for "cover sheet"

% sweep components
for i = 1:size(ds,1)
    if strcmp(ds.classification{i},'accepted')
        c_idx(i) = 1;
    end
    if strcmp(ds.classification{i},'rejected')
        c_idx(i) = 2;
    end
end

H = figure; % prellocate denoising summary plot; serves as "cover sheet"
set(H,'position',[1 596 1680 359]);
set(H,'Color','w');
[~] = evalc('hold');

% extract normalized variance explained 
var_x = ds.normalizedVarianceExplained;

% sub plot
h = subplot(4,6,[1:2 7:8 13:14 19:20]);  % rho vs. kappa plot
[~] = evalc('hold'); % apply hold 

% index of rejected components 
idx = find(strcmp(ds.classification,'rejected'));
[~,b] = sort(var_x(idx),'Descend'); % sort index by variance explained; which avoids having larger components covering up smaller ones in plot 
idx = idx(b);

% sweep 
% rejected
% components 
for i = 1:length(idx)
    
    % check if component 
    % was manually rejected
    if ds.man_rej(idx(i))==1
        s = scatter(ds.kappa(idx(i)),ds.rho(idx(i)),100,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.75 0 0]);
    else
        s = scatter(ds.kappa(idx(i)),ds.rho(idx(i)),100,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.5 0 0]);
    end
    
    % size scales w/
    % normalized variance explained
    s.SizeData = var_x(idx(i))*100;  
    
end

% index of accepted components 
idx = find(strcmp(ds.classification,'accepted'));
[~,b] = sort(var_x(idx),'Descend'); % sort index by variance explained; avoid large components covering up smaller ones in plot 
idx = idx(b);

% sweep 
% accepted 
% components 
for i = 1:length(idx)
    
    % check if component 
    % was manually accepted 
    if ds.man_acc(idx(i))==1
        s = scatter(ds.kappa(idx(i)),ds.rho(idx(i)),100,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .75 0]);
    else
        s = scatter(ds.kappa(idx(i)),ds.rho(idx(i)),100,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .5 0]);
    end
    
    % size scales w/
    % normalized variance explained
    s.SizeData = var_x(idx(i))*100;  
 
end

% index of ignored components 
idx = find(strcmp(ds.classification,'ignored'));
[~,b] = sort(var_x(idx),'Descend'); % sort index by variance explained; avoid large components covering up smaller ones in plot 
idx = idx(b);

% sweep 
% ignored 
% components 
for i = 1:length(idx)
    
    % plot ignored component 
    s = scatter(ds.kappa(idx(i)),ds.rho(idx(i)),100,'MarkerEdgeColor',...
        [0 0 0],'MarkerFaceColor',[.5 .5 .5]);
  
    % size scales w/
    % normalized variance explained
    s.SizeData = var_x(idx(i))*100;  
 
end

% format the figure and axes
set(gca,'FontName','Arial','FontSize',12,...
    'TickLength',[0 0],'TickLength',[0 0]);

% assign 
% axis labels 
xlabel('Kappa (T2*)');
ylabel('Rho (S0)');

ax = h; % minimize white space
ax.Position = [0.0351 0.15 0.3617 0.7366];

% secondsubplot
h = subplot(4,6,[10:11 16:17 22:23]); % grayplots
[c_idx_sort,b] = sort(c_idx); % sort components 
imagesc(c_ts(:,b)'); hold; % plot ICA component time-courses 
colormap(gray); % set colormap 

% format the figure and axes
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],...
    'TickLength',[0 0],'Yticklabel','');

% add colorbar on y-axis
plot(ones(length(find(c_idx_sort==1))+1,1),[find(c_idx_sort==1)-1 ; max(find(c_idx_sort==1))+1],'Color',[0 .5 0],'LineWidth',5);
plot(ones(length(find(c_idx_sort==2))+1,1),[find(c_idx_sort==2)-1 ; max(find(c_idx_sort==2))+1],'Color',[.5 0 0],'LineWidth',5);
plot(ones(length(find(c_idx_sort==3))+1,1),[find(c_idx_sort==3)-1 ; max(find(c_idx_sort==3))+1],'Color',[.5 .5 .5],'LineWidth',5);

ax = h; % min. white space
ax.Position = [0.42 0.15 0.36 0.566];
axis off;

% plot motion
% traces on top 
h = subplot(4,6,[4:5]);

% define some steps 
% for position traces
steps = .4:.10:.9; % selected manually when I first wrote this code; may not generalize well to other data
hold;

% sweep through 
% position estimates 
for i = 1:size(rp,2)
    rp(:,i) = normalize01(rp(:,i));
    l = plot((rp(:,i))+steps(i),'Color',[.5 .5 .5]);
    l.Color(4)=.3;
end

fd = rp; % preallocate FD
fd(1:2,:) = 0; % by convention

% calc. backward
% differences (2 TRs)
for i = 1:size(rp,2)
    for ii = 3:size(fd,1)
        fd(ii,i) = abs(rp(ii,i)-rp(ii-2,i));
    end
end

% convert rotation columns into angular displacement
fd_ang = fd(:,[1:3]); fd_ang = fd_ang / 360; 
fd_ang = fd_ang * 100 * pi;

fd(:,1:3) = []; % delete rotation columns,
fd = [fd_ang fd]; % add back in as angular displacement
fd = sum(fd,2); % sum

% plot frame-
% wise displacement 
plot(fd,'r','LineWidth',1);

% set x-axis
xlim([0 size(rp,1)]);

ax = h; % min. white space
ax.Position = [0.42 0.7473 0.36 0.1711];
axis off;

%%

% final subplot
h = subplot(1,6,6); % relationships with head motion

% sweep 
% directions 
for i = 1:6
    acc{i} = nonzeros(mot_r(c_idx==1,i)');
    rej{i} = nonzeros(mot_r(c_idx==2,i)');
end

try % plot relationships with head motion (rejected)
    plotSpread(rej,'distributionColors',[.5 0 0]);
catch
end

try % plot relationships with head motion (accepted)
    plotSpread(acc,'distributionColors',[0 .5 0]); 
catch
end

% set format figure;
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0],...
    'Xtick',1:size(mot_r,2),'Xticklabel',{'X','Y','Z','y','p','r',}...
    ,'TickLength',[0 0]); xlim([0 size(mot_r,2)+1]); 

% assign axis labels 
ylabel('Absolute Correlation with Head Motion');
xlabel('Direction');

ax = h; % min. white space
ax.Position = [0.83 0.15 0.145 0.7366];

% log frame for movie
m(1) = getframe(H);
close all;

% create quality control vid.; summary + all components
v = VideoWriter([tdir '/fine_tune/all_comps.avi']);
v.FrameRate = 1;
v.Quality=100; % max quality
open(v);
writeVideo(v,m);
close(v);

m(1) = []; % delete "cover sheet"

tm = m; % preallocate temporary movie 
tm(c_idx~=1) = [];

% create quality control vid.; accepted components
v = VideoWriter([tdir '/fine_tune/acc_comps.avi']);
v.FrameRate = .5;
v.Quality = 100; % max quality
open(v);
writeVideo(v,tm);
close(v);

% check if comps.
% were manaually acc.
if sum(ds.man_acc)>0
    
    tm = m; % preallocate
    tm(ds.man_acc~=1) = []; % manually accepted components only
    
    % create quality control vid.; manually accepted components
    v = VideoWriter([tdir '/fine_tune/man_acc_comps.avi']);
    v.FrameRate = 1;
    v.Quality=100; % max quality
    open(v);
    writeVideo(v,tm);
    close(v);
    
end

tm = m; % preallocate
tm(c_idx~=2) = []; % rejected components only

% create quality control vid.; rejected components
v = VideoWriter([tdir '/fine_tune/rej_comps.avi']);
v.FrameRate = .5;
v.Quality = 100; % max quality
open(v);
writeVideo(v,tm);
close(v);

% check if comps.
% were manaually acc.
if sum(ds.man_rej)>0
    
    tm = m; % preallocate
    tm(ds.man_rej~=1) = []; % manually rejected components only
    
    % create quality control vid.; manually rejected components
    v = VideoWriter([tdir '/fine_tune/man_rej_comps.avi']);
    v.Quality = 100; % max quality
    v.FrameRate = 1;
    open(v);
    writeVideo(v,tm);
    close(v);
    
end