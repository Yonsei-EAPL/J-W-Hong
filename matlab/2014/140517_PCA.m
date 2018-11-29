data = load('sediments.txt');

for i = 1:10
    sample(i,:) = ['sample',sprintf('%02.0f',i)];
end
clear i

minerals = ['amp';'pyr';'pla';'ksp';'qtz';'cla';'flu';'sph';'gal'];

corrmatrix = corrcoef(data);
corrmatrix = flipud(corrmatrix);

imagesc(corrmatrix), colormap(hot)
title('Correlation Matrix')
axis square, colorbar, hold

set(gca,'XTickLabel',minerals,'YTickLabel',flipud(minerals))

[pcs,newdata,variances,t2] = princomp(data);
[pcs,newdata] = princomp(data);
pcs(:,1:5)