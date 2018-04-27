addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

order = 10;
load(strcat('ptsGL',num2str(order)))
load(strcat('im_tools',num2str(order)))
load(strcat('im_endsSH',num2str(order)))
%%
level = 8;
% choose orthoganal wavelet for enery preservation
% we choose one that is blocky - we don't want to overfit our test case.
% choose fejer-korovkin for simplicity of representation (in minimum number
% of wavelets needed to cover the image, orthogonality,
% and improvement over haar in invariance.
wname = {'haar','fk4','bior3.7'};
level = 10;
ptl = linspace(96,100,50);
Ddiff =  D(:,:,end)-D(:,:,1);


figure(1)
for i = 1:length(wname)
    [phi,psi,xval] = wavefun(wname{i});
    
    subplot(3,3,i*3-2)
    plot(xval,psi)
    title(sprintf('%s wavelet',wname{i}))
    axis tight
end
%
pthresh = 0.95;

subplot(3,3,[2,3,5,6,8,9])
cla
hold on

c = {'r','b','k'};
for j = 1:length(wname)
    [wdiff,sdiff]=wavedec2(Ddiff,level,wname{j});
    abwdiff = abs(wdiff);
    N = 1:level;
        clear T er r2 wD b NC b
    for i = 1:length(ptl)
        T(i) = prctile(abwdiff,ptl(i));
        NC{i} = wthcoef2('t',wdiff,sdiff,N,repmat(T(i),size(N)),'h');
        wD{i} = waverec2(NC{i},sdiff,wname{j});
        er(i) = immse(Ddiff,wD{i});
        r2(i) = 1-var(Ddiff(:)-wD{i}(:))/var(Ddiff(:));
        b(i) =  abs((sum(Ddiff(:))-sum(wD{i}(:)))/sum(Ddiff(:)));
    end
    hb{j} = plot(ptl,b,':','color',c{j},'linewidth',2);
    h{j} = plot(ptl,r2,'color',c{j},'linewidth',2);
    [~, qc(j)] = min(abs(r2-pthresh));  % where does the curve reach 90% invariance
    QC = wthcoef2('t',wdiff,sdiff,N,repmat(T(qc(j)),size(N)),'h');
    np(j) = sum(QC~=0);
    Q{j} = waverec2(QC,sdiff,wname{j});
end
%
for j = 1:length(wname)
    txth{j}=strcat('\begin{tabular}{lr} \textbf{',wname{j},'} inv. &',...
        num2str(np(j)), sprintf(' coefficients at %0.2f',pthresh), '\end{tabular}');
    txtb{j}=strcat('\begin{tabular}{l} \textbf{',wname{j},'} bias \end{tabular}');
end
axis([96 100 ylim])
plot(xlim,[pthresh pthresh],'--','color',0.6*[1 1 1],'linewidth',1)
ylabel('image invariance or bias')
xlabel('percentile threshold')
lgd = legend([h{1} h{2} h{3} ,...
    hb{1} hb{2} hb{3}],[txth txtb]);
set(lgd,'interpreter','latex')
grid on