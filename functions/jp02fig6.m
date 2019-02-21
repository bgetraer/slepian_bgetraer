function jp02fig6(originalimage, wname, level, ptl)
%JP02FIG6 Creates Figure 6 from "REGIONAL FORCING OF GREENLAND ICE LOSS 
%   2002-2017" Spring 2018 Junior Paper, Princeton Department of Geosciences
%
% SEE ALSO:
%   WAVEFUN, WAVEDEC2, WTHCOEF2, WAVEREC2
%
% last modified by: bgetraer@princeton.edu, 7/25/2018


% throw errors and warnings as necessary
if length(wname)>3
    error('Maximum of 3 wavelets for testing!')
end
if length(wname)<3
    warning('Expected 3 wavelets for testing, but can run with fewer.')
end


figure(1)
% Plot the wavelet shapes
for i = 1:length(wname)
    [~,psi,xval] = wavefun(wname{i});
    
    subplot(3,3,i*3-2)
    plot(xval,psi,'linewidth',2)
    title(sprintf('%s wavelet',wname{i}))
    axis tight
end

% Plot the different effects of data compression
subplot(3,3,[2,3,5,6,8,9])
cla
hold on
% plotting color scheme
c = {'r','b','k'};

% Goal of 90% invariance (our somewhat arbitrarily 
    %   chosen baseline for "good"):
pthresh = 0.9;

% initialize variables for outer loop
hb = cell(1,length(wname));
h = cell(1,length(wname));
qc = zeros(1,length(wname));
np = zeros(1,length(wname));
Q = cell(1,length(wname));
txth = cell(1,length(wname));
txtb = cell(1,length(wname));


% Loop over each wavelet
for j = 1:length(wname)
    % perform the wavelet deconstruction
    [wdiff,sdiff]=wavedec2(originalimage,level,wname{j});
    % absolute values of wavelet coefficients
    abwdiff = abs(wdiff);
    N = 1:level;
    
    % initialize variables for inner loop
    T = zeros(1,length(ptl)); 
    NC = cell(1,length(ptl));
    wD = cell(1,length(wname));
    er = zeros(1,length(ptl));
    r2 = zeros(1,length(ptl));
    b = zeros(1,length(ptl));
    
    % Evaluate image reconstruction after throwing away a certain
    %   percentile of the wavelet coefficients
    for i = 1:length(ptl)
        % Percentile threshold level 
        T(i) = prctile(abwdiff,ptl(i));
        % Hard thresholding of wavelet coefficients at the threshold level
        NC{i} = wthcoef2('t',wdiff,sdiff,N,repmat(T(i),size(N)),'h');
        % Reconstruction of the image after thresholding
        wD{i} = waverec2(NC{i},sdiff,wname{j});
        % mean squared error between the original and the thresholded
        %   reconstructions
        er(i) = immse(originalimage,wD{i});
        % pixel-to-pixel invariance between the original and the 
        %   thresholded reconstructions
        r2(i) = iminvar(wD{i}, originalimage);
        % total energy loss (bias) between the original and the 
        %   thresholded reconstructions
        b(i) =  imbias(wD{i}, originalimage);
    end
    
    % plotting routines
    hb{j} = plot(ptl,b,':','color',c{j},'linewidth',2);
    h{j} = plot(ptl,r2,'color',c{j},'linewidth',2);
    
    % where does the curve reach 90% invariance (our somewhat arbitrarily 
    %   chosen baseline for "good"):
    [~, qc(j)] = min(abs(r2-pthresh));  
    QC = wthcoef2('t',wdiff,sdiff,N,repmat(T(qc(j)),size(N)),'h');
    np(j) = sum(QC~=0);
    Q{j} = waverec2(QC,sdiff,wname{j});
end
% figure text
for j = 1:length(wname)
    txth{j}=strcat('\begin{tabular}{lr} \textbf{',wname{j},'} inv. &',...
        num2str(np(j)), sprintf(' coefficients at %0.2f',pthresh), '\end{tabular}');
    txtb{j}=strcat('\begin{tabular}{l} \textbf{',wname{j},'} bias \end{tabular}');
end
% axis plotting and legend
axis([96 100 0 1])
plot(xlim,[pthresh pthresh],'--','color',0.6*[1 1 1],'linewidth',1)
ylabel('image invariance or bias')
xlabel('percentile threshold')
lgd = legend([h{1} h{2} h{3} ,...
    hb{1} hb{2} hb{3}],[txth txtb]);
set(lgd,'interpreter','latex')
grid on

end

