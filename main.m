close all
clear all

ubcdat = readtable([pwd filesep 'ubc_raw_expression.csv']);
load('umap_mock.mat')
%%
mcount = table2array(ubcdat(:,2:end));
mnames = ubcdat.Properties.VariableNames(2:end);
[n, k] = size(mcount);
mlogcount = log(1+mcount);
mncount = mcount./repmat(sum(mcount,2),[1 k])*mean(sum(mcount,2));
mlogncount = log(1+mncount);
mnncount = mncount./repmat(mean(mncount,1),[n 1]);
mlognncount = log(mnncount);
mc = mean(mlogcount);
[sm sidxm] = sort(mc,'descend');
msmall = mlogcount(:,sidxm(1:5000));
[mM mN] = nnmf(msmall,3);

[smM, psuid] = sort(mM(:,1));
[cc, pvals] = corr(smM,mncount(psuid,:),'Type','Pearson');
cc(isnan(cc)) = 0;

topn = 400;
ssidxm = sidxm(1:topn);
[scc, scidx] = sort(cc(sidxm(1:topn)),'descend');
dispcnt = mncount(flipud(psuid),ssidxm(scidx));
dispncnt = dispcnt./repmat(nanmean(dispcnt),[numel(psuid) 1]);
%% correlation of 1d sorting with expression
figure('Position',[1458 442 300 500]); 
subplot(8,1,1:2)
scatter(Y(:,1),Y(:,2),5,mM(:,1),'filled','MarkerFaceAlpha',.5)
axis equal off; 
colormap('fire')
subplot(8,1,3:8)
imagesc(log(medfilt2(dispncnt,[10 2]))')
% top 1 Pearson R over all genes
sigcbar = ones(size(scc))*n; sigcbar(abs(scc)<.2537) = nan;
hold on; plot(sigcbar,1:topn,'g-')
caxis([-.5 .5]); colormap('fire')
set(gca,'YAxisLocation','right','TickDir','out','YDir','reverse')
lidx = find(contains(mnames(ssidxm(scidx)),{'Grm','Dgk','Dagl','Plcb','Gria','Grin','Grid','Kcnj','Trpc'}));
yticks(lidx);
yticklabels(mnames(ssidxm(scidx(lidx))));
xlabel('Factor sorted cells')
xticks([1 numel(psuid)])
box off
colorbar('XTick', [-.5 0 .5])

%% plot the umap gradient maps for all genes of interest
mnfcount = [mncount(psuid,:) mM(psuid,:) [1:numel(psuid)]'];
mfnames = [mnames {'NNMF1','NNMF2','NNMF3','PsuID'}];
%cplots = {'NNMF1','NNMF2','NNMF3','PsuID','Gria2','Grin1','Grid2','Grm2','Kcnj6','Grm1','Gnaq','Gna11','Plcb4','Dgkb','Dgkg','Trpc1','Plpp4','Plppr2','Plppr4','Daglb','Dagla'};
cplots = ['NNMF1','NNMF2','NNMF3','PsuID', sort(mnames(ssidxm(scidx(find(contains(mnames(ssidxm(scidx)),{'Grm','Dgk','Dagl','Plcb','Gria','Grin','Grid','Kcnj','Trpc'}))))))];
%cplots = mfnames(find(contains(mfnames,'Prkc')));
fig = upanel(mnfcount,Y(psuid,:),mfnames,cplots,[1:numel(psuid)]','tsne');

%% for PC UBC paper
cplots = ['NNMF1','NNMF2','NNMF3','PsuID',mfnames(find(contains(mfnames,{'Glr'})))];
%cplots = mfnames(find(contains(mfnames,'Prkc')));
fig = upanel(mnfcount,Y(psuid,:),mfnames,cplots,[1:numel(psuid)]','tsne');