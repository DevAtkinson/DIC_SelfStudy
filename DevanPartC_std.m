function DevanPartC()
	close all
	clc
	% f='C:\Users\devan\Downloads\varsity\DIC selfstudy\assessmentImages';
	% f='I:\PC install\DIC selfstudy\assessmentImages';
	f='D:\Work\DIC selfstudy\assessmentImages';
	im='calimg';
	for i=1:9
		ims{i}=sprintf('0%d',i);
	end
	calib=DevanCalibration('folder',f,'ImageName',im,'images',ims)
	subsize=41;
	subpos=[425,230];


	% for i=1:max(size(ims))
	% 	images=strcat(strcat(im,num2str(ims{i})),'.tif');
	% 	ImageFileName{i}=fullfile(f,images);
	% 	I{i} = imread(ImageFileName{i});
	% end

	% F=undistort(im2double(I{1}),calib);
	% figure()
	% imagesc(I{1})
	% figure
	% imagesc(F)

	imageFileName2{1}=fullfile(f,'img00.tif');	%original image
	imageFileName2{2}=fullfile(f,'img01.tif');	%camera distortion image
	imageFileName2{3}=fullfile(f,'img02.tif');	%deformed image
	imageFileName2{4}=fullfile(f,'img03.tif');	%deformed and distorted image
	for i=1:4
		II{i}=imread(imageFileName2{i});
	end
	FF=undistort(im2double(II{4}),calib); %for correlation



	F_in=im2double(II{1});
	% G_in=imnoise(FF,'gaussian',0);
	G_in=FF;
	F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);


	figure(1)
	subplot(2,2,1);
	imagesc(F_in);
	hold on;
	plot(subpos(1),subpos(2),'rx')
	subplot(2,2,2);
	imagesc(F);
	hold on;
	[xf,yf]=ginput(1);
	plot(xf,yf,'rx');
	subplot(2,2,3);
	imagesc(G_in);
	key=0;
	while(key==0)
		[xg1,yg1]=ginput(1);
		subplot(2,2,4);
		imagesc(G_in(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		key = waitforbuttonpress;
	end
	subplot(2,2,4);
	imagesc(G_in(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
	[xg,yg]=ginput(1);
	xguess=(xf+subpos(1))-(xg+floor(xg1)-floor(subsize/2));
	yguess=(yf+subpos(2))-(yg+floor(yg1)-floor(subsize/2));
	% guess=[xguess;0;0;yguess;0;0]
	tic
	guess=[(xguess);0;0;(yguess);0;0];
	parfor i=1:1000
		[P{i},Corr(i)]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess);
		Corr2(i)=(1-Corr(i)/2);
	end
	result.P=P;
	result.Corr=Corr;
	result.Corr2=Corr2;
	save('std_results_partc.mat','result');
	Final_corrected_image=undeform(G_in,P,subsize,subpos);

	% meshcompare(F_in,Final_corrected_image);

	figure()
	imagesc(F_in)
	figure()
	imagesc(Final_corrected_image)










	% figure()
	% imagesc(II{1})
	% figure
	% imagesc(FF)
	% figure
	% imagesc(II{2})

	% meshcompare(im2double(II{1}),FF)
	% meshcompare(II{1},II{2})

	% [r,c]=size(FF);
	% for i=1:r
	% 	for j=1:c
	% 		thing(i,j)=(im2double(II{1}(i,j))-FF(i,j))/im2double(II{1}(i,j));
	% 	end
	% end
	% figure
	% surf(thing)
	% figure
	% surf(II{1})

end

function I = undistort(Idistorted, params)
	% https://stackoverflow.com/questions/12117825/how-can-i-undistort-an-image-in-matlab-using-the-known-camera-parameters
	fx = params(3);
	fy = params(4);
	cx = params(1);
	cy = params(2);
	fs = params(5);
	k1 = params(end-1);
	k2 = params(end);
	% k3 = params(1);
	% p1 = params(1);
	% p2 = params(1);

	K = [fx fs cx; 0 fy cy; 0 0 1];

	I = zeros(size(Idistorted));
	[i j] = find(~isnan(I));

	% Xp = the xyz vals of points on the z plane
	Xp = inv(K)*[j i ones(length(i),1)]';

	% Now we calculate how those points distort i.e forward map them through the distortion
	r2 = Xp(1,:).^2+Xp(2,:).^2;
	x = Xp(1,:);
	y = Xp(2,:);

	x = x.*(1-k1*r2 - k2*r2.^2);
	y = y.*(1-k1*r2 - k2*r2.^2);

	% u and v are now the distorted cooridnates
	u = reshape(fx*x + cx,size(I));
	v = reshape(fy*y + cy,size(I));

	% Now we perform a backward mapping in order to undistort the warped image coordinates
	I = interp2(Idistorted, u, v);
end

function send_this_out=undeform(Ideformed,P,subsize,subpos)
	G=Ideformed;
	[r,c]=size(G);
	X=1:1:c;										% x values over subset
	Y=1:1:r;										% y values over subset
	x0=subsize/2+subpos(1);													% x value for subset centre
	y0=subsize/2+subpos(2);													% y value for subset centre

	[Xmesh,Ymesh]=meshgrid(1:1:c,1:1:r);

	
	for l=1:r
		for j=1:c
			dx=X(j)-x0;
			dy=Y(l)-y0;
			xp(l,j)=x0+dx*(1+P(2))+P(3)*dy+P(1);
			yp(l,j)=y0+dy*(1+P(6))+P(5)*dx+P(4);
			
		end
	end
	send_this_out=interp2(Xmesh,Ymesh,G,xp,yp,'linear');

end

function out=Warp(p,d)
	a=[(1+p(2)), p(3), p(1);
		p(5), (p(6)+1), p(4);
		0 0 1];
	b=[(1+d(2)), d(3), d(1);
		d(5), (d(6)+1), d(4);
		0 0 1];
	out=a*inv(b);
end

function CorrelationStation(F,G)
	[r,c]=size(F);
	for i=1:r
		for j=1:c
			temp(i,j)=abs(F(i,j)-G(i,j));
		end
	end
	figure
	bar3(temp)
	xlabel('x')
	ylabel('y')
end
