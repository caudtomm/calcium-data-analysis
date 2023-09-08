data = experiment.series{1}.data;

A = traceFormat(data.tracesdn,data.L);
A = transpose(traceFormat(permute(A,[1 3 2])));

Aztime = zscore(A,[],1);
Azcells = zscore(A,[],2);

%%
[W,H] = nnmf(A,100);
[Wztime,Hztime] = nnmf(Aztime,100);
[Wzcells,Hzcells] = nnmf(Azcells,100);

figure;
subplot(231); imagesc(W); title('raw 100 params'); subplot(234); imagesc(H);
subplot(232); imagesc(Wztime);title('z-time 100 params'); subplot(235); imagesc(Hztime);
subplot(233); imagesc(Wzcells);title('z-cells 100 params'); subplot(236); imagesc(Hzcells);

%%
[W20,H20] = nnmf(A,20);
[Wzcells20,Hzcells20] = nnmf(Azcells,20);
[Wztime20,Hztime20] = nnmf(Aztime,20);

figure;
subplot(231); imagesc(W20); title('raw 20 params'); subplot(234); imagesc(H20);
subplot(232); imagesc(Wztime20);title('z-time 20 params'); subplot(235); imagesc(Hztime20);
subplot(233); imagesc(Wzcells20);title('z-cells 20 params'); subplot(236); imagesc(Hzcells20);

%%

contributions = Wztime;

fs = data.meta.framerate;
interval = transpose(floor(data.stim_on_sec*fs) : 1+floor(data.stim_off_sec*fs));
distTraceData(data,[],[],[],'cosine',data.common_units,interval');
