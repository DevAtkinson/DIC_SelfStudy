function DevanPartC()
	close all
	clc
	% find folder containing the images
	current_folder=pwd;
	if ispc
		image_folder=strcat(current_folder,'\Images');
	elseif ismac||isunix
		image_folder=strcat(current_folder,'/Images');
	end

	im='calimg';
	for i=1:9 
		if i<10
			ims{i}=sprintf('0%d',i);
		else
			ims{i}=sprintf('%d',i);
		end
	end
	% perform calibration
	calib=DevanCalibration('folder',image_folder,'ImageName',im,'images',ims)
	% define subset size
	subsize=41;
	% define subset position
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

	% import images
	imageFileName2{1}=fullfile(image_folder,'img00.tif');	%original image
	imageFileName2{2}=fullfile(image_folder,'img01.tif');	%camera distortion image
	imageFileName2{3}=fullfile(image_folder,'img02.tif');	%deformed image
	imageFileName2{4}=fullfile(image_folder,'img03.tif');	%deformed and distorted image
	for i=1:4
		II{i}=imread(imageFileName2{i});
	end
	% undistort image 'img03.tif'
	FF=undistort(im2double(II{4}),calib); %for correlation
	F_in=im2double(II{1});
	G_in=FF;
	%define reference subset
	F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);

	% this code allows the user to select a seed point in the reference and
	% investigated image that serve as an initial guess for the correlation function
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
	fprintf('\nSelect the point in the investigated image where the subset will be located\n');
	fprintf('This function will loop until a keyboard key is pressed at which point you can select the seed point in the subset\n')
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
	tic
	guess=[(xguess);0;0;(yguess);0;0]
	% perform correlation
	[P,Corr]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess)

	% undeform the investigated image
	Final_corrected_image=undeform(G_in,P,subsize,subpos);

	% meshcompare(F_in,Final_corrected_image);

	figure()
	imagesc(F_in)
	title('Reference image')
	figure()
	imagesc(Final_corrected_image)
	title('Investigated image undistorted and undeformed (in that order)')
end

function I = undistort(Idistorted, params)
	% function to undistort an image using the calibration parameters. This was sourced from the following website
	% https://stackoverflow.com/questions/12117825/how-can-i-undistort-an-image-in-matlab-using-the-known-camera-parameters
	fx = params(3);
	fy = params(4);
	cx = params(1);
	cy = params(2);
	fs = params(5);
	k1 = params(end-1);
	k2 = params(end);

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
	% function to undeform the investigated image
	G=Ideformed;
	[r,c]=size(G);
	X=1:1:c;							% x values over subset
	Y=1:1:r;							% y values over subset
	x0=subsize/2+subpos(1);				% x value for subset centre
	y0=subsize/2+subpos(2);				% y value for subset centre

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