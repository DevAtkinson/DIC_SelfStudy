function testGPU()
V=rand(50);
X=1:1:50;
Y=1:1:50;
xq = linspace(-3, 3, 100);
yq=linspace(-3,3,100);
tic
Vq = interp2(X, Y, V, xq, yq);
toc
tic
Vg=gpuArray(V);
Xg=gpuArray(X);
Yg=gpuArray(Y);
xqg=gpuArray(xq);
yqg=gpuArray(yq);
toc
tic
Vq = interp2(X, Y, V, xq, yq);
toc

end