function DevanPartA()
	close all
	clc
	% find folder containing the images
	current_folder=pwd;
	if ispc
		image_folder=strcat(current_folder,'\Images');
	elseif ismac||isunix
		image_folder=strcat(current_folder,'/Images');
	end
	% define the names of the images
	im='calimg';
	for i=1:9
		if i<10
			ims{i}=sprintf('0%d',i);
		else
			ims{i}=sprintf('%d',i);
		end
	end
	% call the calibration function to perform the calibration
	calib=DevanCalibration('folder',image_folder,'ImageName',im,'images',ims)

	% Import the images of Einstein (distorted and undistorted)
	imageFileName2{1}=fullfile(image_folder,'img00.tif');
	imageFileName2{2}=fullfile(image_folder,'img01.tif');
	for i=1:2
		II{i}=imread(imageFileName2{i});
	end
	% undistort the distorted image of Einstein
	FF=undistort(im2double(II{2}),calib);

	% display images
	figure()
	imagesc(II{1})
	title('Reference image of Einstein')
	figure
	imagesc(FF)
	title('Undistorted image of Einstein')
	figure
	imagesc(II{2})
	title('Original distorted image of Einstein')

	% meshcompare the original images and the reference image and the undistorted image
	% meshcompare(II{1},II{2})
	% meshcompare(im2double(II{1}),FF)
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