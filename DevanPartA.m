function DevanPartA()
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

	for i=1:max(size(ims))
		images=strcat(strcat(im,num2str(ims{i})),'.tif');
		ImageFileName{i}=fullfile(f,images);
		I{i} = imread(ImageFileName{i});
	end

	% F=undistort(im2double(I{1}),calib);
	% figure()
	% imagesc(I{1})
	% figure
	% imagesc(F)

	imageFileName2{1}=fullfile(f,'img00.tif');
	imageFileName2{2}=fullfile(f,'img01.tif');
	for i=1:2
		II{i}=imread(imageFileName2{i});
	end
	FF=undistort(im2double(II{2}),calib);

	figure()
	imagesc(II{1})
	figure
	imagesc(FF)
	figure
	imagesc(II{2})

	meshcompare(im2double(II{1}),FF)
	meshcompare(II{1},II{2})

	[r,c]=size(FF);
	for i=1:r
		for j=1:c
			thing(i,j)=(im2double(II{1}(i,j))-FF(i,j))/im2double(II{1}(i,j));
		end
	end
	figure
	surf(thing)
	figure
	surf(II{1})
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