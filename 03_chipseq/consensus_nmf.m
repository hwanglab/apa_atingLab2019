clc; close all; clear

%--------------------------------------------------------------------------
% Consensus clustering (or clustering ensemble) 
%
% Input: mm_Xfea (input data: samples X features)
% Ouput: mm_ClusterMat (final clustering reults: samples X K =[2,3,...]), 
%        where each column inclules clustering results of samples 
%        with the fixed number of clusters (i.e., K is set to 2)
%
% The concept of consensus clustering is to aggregate clustering results
% from multiple runs of simple clustering methods. In this script, we use
% a variant of nonnegative matrix factorization (NMF) as a single clustering 
% method. Our approach to generate multiple clustering results with a single
% clustering method (i.e., the NMF) can be understood as doing bootstrapping: 
% for each K (number of clusters) we run the NMF with the fixed latent 
% dimensionality (set to K) on randomly selected 80% samples from the input 
% data (without replacement) and repeat this experiment 100 times.
% To aggregate clustering results, we follow the approach proposed in        
% S. Monti, et al., "Consensus Clustering: A Resampling-Based Method for 
% Class Discovery and Visualization of Gene Expression Microarray Data",
% Machine Learning 2003. (For each K) We construct a consensus matrix 
% (samples X samples) from clustering results, where each element in the 
% matrix represents the probability of the corresponding two samples being   
% in the same cluster. (For each K) We then apply hierarchical clustering on
% the consensus matrix in order to obtain final cluster labels.
% 
%                                                               Mar/15/2019
%                                   programmed by Sunho Park (parks@ccf.org)
%--------------------------------------------------------------------------

%- load the data from a formatted text file
fid = fopen('./fastq/cluster_analy/chipseq_xing_ublocks_ublocks_filt10a.tsv', 'r');

tline = fgets(fid);
m_cheaders_en = strsplit(tline,'\t');
m_nLen = length(m_cheaders_en);

m_strpattern = ['%s', repmat('%f',[1,m_nLen-1]), '%[^\n\r]'];
m_mData = textscan(fid,m_strpattern,'Delimiter','\t');

fclose(fid);  

mc_Feat_Names = strtrim(m_cheaders_en(2:9));
mc_Sample_Names = m_mData{1,1};

mm_Xfea = NaN*ones(length(mc_Sample_Names), length(mc_Feat_Names));

for mn_i = 1:(length(mc_Feat_Names))
    mm_Xfea(:,mn_i) = m_mData{1,1+mn_i};
end



%- set the path for the NMF implementation
% please read the following paper for the detailed information 
% about the NMF implementation used in this scirpt 
%
% For using this software, please cite:
%     Jingu Kim and Haesun Park, "Toward Faster Nonnegative Matrix Factorization: 
%     A New Algorithm and Comparisons, ICDM 2008, 353-362, 2008.
% also visit the atuhors's website: https://www.cc.gatech.edu/~hpark/nmfbpas.html

addpath('../resource/libs/nmf_bpas');

m_nMaxIters = 500;
m_rResRate = 0.8;

m_nSamples = length(mc_Sample_Names);

m_strMethod = 'NMF';
disp(['Sub Clutering: ', m_strMethod]);

m_str_savefile = [pwd(), '/fastq/cluster_analy/figures/chipseq_xing_ublocks_ublocks_filt10a_', m_strMethod, '_'];

%- input 
m_strNMFtrans = '';
 
if (sum(sum(mm_Xfea<0)) ==  0)
    X = mm_Xfea;
else
    % To handle negative values in the input data matrix: we divide the matrix 
    % into positive and negative parts and take absolute values from the negative part.    
    X = [max(mm_Xfea,0) max(-mm_Xfea,0)];
    
    m_strNMFtrans = [m_strNMFtrans, ' Duplication'];
end

disp(['Tranformation: ', m_strNMFtrans]);

m_nSubSamples = round(m_nSamples*m_rResRate);

%- K: the number of clusters [2,3,...]
m_vKs = 2:size(X,2);

m_cResults.featur_names = mc_Feat_Names;
m_cResults.discription = './fastq/cluster_analy/chipseq_xing_ublocks_ublocks_filt10a.tsv';
m_cResults.Ks = m_vKs;

%- A(K): the area under the CDF (please see Eq. (6) in the cited paper,
%  "Consensus Clustering:~", Machine Learning 2003
m_vAk = zeros(length(m_vKs), 1);

%- A matrix which contains all the clustering results (each column == each K)  
mm_ClusterMat = zeros(m_nSamples, length(m_vKs));

for m_nselK = 1:length(m_vKs)
    K = m_vKs(m_nselK);
    
    fprintf('K = %d \n', K);
    
    %- run the NMF 500 times
    m_mGroup = inf*ones(m_nSamples, m_nMaxIters);
    for m_niter = 1:m_nMaxIters
        %- resampling: randomly select 80% samples without replacement
        m_vRandIDX = randperm(m_nSamples);
        m_vDataIDX = m_vRandIDX(1:m_nSubSamples);

        m_mX = X(m_vDataIDX,:);
        [n, m1] = size(m_mX);
        
        %------------------------------------------------------------------
        %- NMF
        % random initialization for W and H (X \approx WH)
        W_init = rand(n, K);
        H_init = rand(K, m1);
        
        m_cResults.W_init{m_niter,m_nselK} = W_init;
        m_cResults.H_init{m_niter,m_nselK} = H_init;
        m_cResults.DataID{m_niter,m_nselK} = m_vRandIDX;
        
        % call the NMF function
        [W, H] = nmf(m_mX, K, 'W_INIT',W_init, 'H_INIT',H_init, ...
            'type','plain', 'tol',1e-4, 'NNLS_SOLVER','bp');
        
        m_cResults.W{m_niter,m_nselK} = W;
        m_cResults.H{m_niter,m_nselK} = H;
        
        %- clustering results from the NMF
        [val, m_vidx] = max(W, [], 2);
        %------------------------------------------------------------------
        
        m_mGroup(m_vDataIDX,m_niter) = m_vidx;
    end
    
    %- Construct the consensus matrix  m_mSim (samples X samples)
    m_mSim = zeros(m_nSamples, m_nSamples);
    for m_ni = 1:m_nSamples
        m_vCurIDs = m_mGroup(m_ni,:);
        
        m_mChk = (m_mGroup)./repmat(m_vCurIDs,[m_nSamples,1]);
        m_mChk = m_mChk == 1.0;
        
        m_mInc = (m_mGroup<Inf) & repmat(m_vCurIDs<Inf,[m_nSamples,1]);
        
        m_mSim(m_ni,:) = sum(m_mChk,2)./sum(m_mInc,2);
    end
    
    % diagonal elements should be 1  
    m_mSim(1:m_nSamples+1:m_nSamples^2) = 1;
        
    m_cResults.Group{m_nselK,1} = m_mGroup;
    m_cResults.SimMat{m_nselK,1} = m_mSim;
        
    %- final clustering results 
    mv_pdist = squareform(1-m_mSim);
    Z = linkage(mv_pdist, 'average');
    m_vgroupIDX = cluster(Z, 'maxclust', K); % cluster assignment
    
    mm_ClusterMat(:,m_nselK) = m_vgroupIDX;
    
    %- Draw a Heatmap 
    [m_vdummy, m_vReorderIDX] = sort(m_vgroupIDX);
        
    fig_open = figure('visible','off');
    colormap('jet');   % set colormap
    h = imagesc((m_mSim(m_vReorderIDX, m_vReorderIDX)));        % draw image and scale colormap to values range
    colorbar; 
    
    m_strname = [m_str_savefile, 'Cluster_K_', num2str(K), '.jpg'];
    saveas(fig_open, m_strname)
    close all hidden
    
    %- Calculate A(k)  
    m_vupperidx = find(~tril(ones(size(m_mSim))));
    m_vUpperVals = sort(m_mSim(m_vupperidx), 'ascend');
    
    m_vuniquevals = unique(m_vUpperVals);
    
    N = hist(m_vUpperVals, m_vuniquevals);
    m_vdf = N'/length(m_vUpperVals);
    m_vCDF = cumsum(m_vdf);
    
    m_vdiff = m_vuniquevals(2:end) - m_vuniquevals(1:end-1);

    m_vAk(m_nselK) = sum(m_vdiff.*m_vCDF(2:end));
end
 
save([m_str_savefile, 'Results'], 'm_cResults', '-V7.3');

%- Delta (K): change between A(K) and A(K-1)
m_vDeltaK_num2 = zeros(length(m_vKs),1);

m_vAhatk = zeros(length(m_vKs),1);
for m_ni = 1:length(m_vKs)
    m_vAhatk(m_ni) = max(m_vAk(1:m_ni));
end

m_vDeltaK_num2(1) = m_vAhatk(1);
m_vDeltaK_num2(2:end-1) = ( m_vAhatk(3:end)-m_vAhatk(2:end-1) )./m_vAhatk(3:end);
m_vDeltaK_num2(end) = NaN;

% Display some statistics 
disp('>>>>> K ----- Ak ----- Delta')
disp([m_vKs', m_vAk, m_vDeltaK_num2])

fig_open = figure('visible','off');

plot(m_vKs(1:end-1), m_vDeltaK_num2(1:end-1), '-bs',...
    'LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',10)
title('Relative Change in Area (\Delta) under the CDF plot');
ylabel('\Delta(K)');
xlabel('Number of subtypes(K)');

m_strname = [m_str_savefile, '_Delta.jpg'];
saveas(fig_open, m_strname)
close(fig_open)

%- save the clustering results into a text file
fid = fopen(['./fastq/cluster_analy/figures/chipseq_xing_ublocks_ublocks_filt10a_', m_strMethod, '_groupInfo.txt'], 'w');

m_strpattern = repmat('%s\t',[1,1+size(mm_ClusterMat,2)]);
m_strpattern(end-1:end) = '\n';

mc_header = cellfun(@(x) sprintf('(K=%d)', x), num2cell(m_vKs),'UniformOutput',false);
fprintf(fid, m_strpattern, 'ID', mc_header{:});

m_strpattern = ['%s\t', repmat('%d\t',[1,size(mm_ClusterMat,2)])];
m_strpattern(end-1:end) = '\n';

for m_ni = 1:length(mc_Sample_Names)
        fprintf(fid, m_strpattern, mc_Sample_Names{m_ni}, mm_ClusterMat(m_ni,:));
end

fclose(fid);
