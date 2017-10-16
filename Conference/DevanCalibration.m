function actualout=MatlabCalibration(varargin)
	% format long
	%handle the input variables
	for i=1:nargin/2
		switch varargin{i*2-1}
		case 'folder'
			folder=varargin{i*2};
		case 'ImageName'
			ImageName=varargin{i*2};
		case 'images'
			imageNumbers=varargin{i*2};
		end
	end

	%create path to specific images
	for i=1:max(size(imageNumbers))
		images=strcat(strcat(ImageName,num2str(imageNumbers{i})),'.tif');
		ImageFileName{i}=fullfile(folder,images);
		I{i} = imread(ImageFileName{i});
	end
	%detect the location of checker pattern intersections
	[imagePoints, boardSize] = detectCheckerboardPoints(ImageFileName);

	% assign world coordinates of the calibration targets
	worldPoints=worldPointsvalues;
	% determine the homography of each calibration image
	Homo=Homographies(worldPoints,imagePoints);
	% calculate the vij vectors of each image
	V=calculateV(Homo);
	% calculate the B vector
	B=calculateB(V);
	% determine the instrinsic parameters
	A=intrinsicParam(B);
	% determine the extrinsic parameters
	[rot,trans]=extrinsicParam(Homo,A);

	% convert the 9 element rotation matrix to 3 rotation angles about each axis
	% this is done to reduce the amount of variables that need to be optimised by
	% the iterative method
	[rr,cc,kk]=size(rot);
	for i=1:kk
		Wtemp=R2W(rot(:,:,i));
		W(i,1)=Wtemp(1);
		W(i,2)=Wtemp(2);
		W(i,3)=Wtemp(3);
	end

	% assign necessary variables
	[r,c,k]=size(imagePoints);
	cx=A(1,3);
	cy=A(2,3);
	fx=A(1,1);
	fy=A(2,2);
	fs=A(1,2);
	tz=trans(:,3);
	tx=trans(:,1);
	ty=trans(:,2);
	R11=rot(1,1,:);
	R12=rot(1,2,:);
	R13=rot(1,3,:);
	R21=rot(2,1,:);
	R22=rot(2,2,:);
	R23=rot(2,3,:);
	R31=rot(3,1,:);
	R32=rot(3,2,:);
	R33=rot(3,3,:);
	k1=0;
	k2=0;

	% convert form world coordinates to coordinates on the image plane (plane points)
	PlanePoints=Step1(tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints);
	% determine the coordinates on the sensor plane for the plane points determined above
	sensorPoints=affine_transformation(cx,cy,fx,fy,fs,PlanePoints);
	% estimate the radial distortion values
	[k1,k2]=est_radial_distortion(k1,k2,cx,cy,imagePoints,sensorPoints,PlanePoints)

	% assign necessary variables for non-linear least squares optimisation
	[r,c,k]=size(imagePoints);
	lsqin2(1)=A(1,3);
	lsqin2(2)=A(2,3);
	lsqin2(3)=A(1,1);
	lsqin2(4)=A(2,2);
	lsqin2(5)=A(1,2);
	lsqin2(6:5+k)=trans(:,3);
	lsqin2(6+k:5+2*k)=trans(:,1);
	lsqin2(6+2*k:5+3*k)=trans(:,2);
	lsqin2(6+3*k:5+4*k)=W(:,1);
	lsqin2(6+4*k:5+5*k)=W(:,2);
	lsqin2(6+5*k:5+6*k)=W(:,3);
	lsqin2(6+6*k)=k1;
	lsqin2(7+6*k)=k2;

	% define a matrix which controls which values are optimised in the non-linear least squares optimisation
	% Ones indicate that the variable is optimised and zero means it is not
	% (here the extrinsic parameters are fist optimised before the intrinsic parameters are optimised)
	improve=[0 0 0 0 0 1 1 1 1 1 1 0 0;
			1 1 1 1 1 1 1 1 1 1 1 1 1];
	[improve_r,improve_c]=size(improve);

	% perform non-linear least squares optimisation
	lsqin=lsqin2;
	options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','FunctionTolerance',1e-40,'StepTolerance',1e-40,'MaxFunctionEvaluations',1e6,'MaxIterations',100,'UseParallel',false);
	tic
	for i=1:improve_r
		disp(improve(i,:))
		out=lsqnonlin(@(lsqin) lsqcontrol(improve(i,:),lsqin,worldPoints,imagePoints,lsqin2),lsqin2,[],[],options);
		lsqin=out;
		lsqin2=out;
	end
	toc
	indexess=[1 2 3 4 5 (6+6*k) (7+6*k)];
	actualout=out(indexess);
end

function out=Step1(tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints)
	% function to determine the image plane coordinates of the world coordinates (using rotation matrix)
	[r,c,k]=size(imagePoints);
	for j=1:k
		T=[R11(j) R12(j) R13(j) tx(j);R21(j) R22(j) R23(j) ty(j);R31(j) R32(j) R33(j) tz(j)];
		for i=1:r
			Z=[worldPoints(i,1);worldPoints(i,2);0;1]; 
			answer=T*Z; %#2.21

			%project to normalised ideal image plane
			out(i,1,j)=answer(1)/answer(3);
			out(i,2,j)=answer(2)/answer(3);
		end
	end
end

function out=Step1W(tz,tx,ty,W1,W2,W3,worldPoints,imagePoints)
	% function to determine the image plane coordinates of the world coordinates (using rotation angles)
	[r,c,k]=size(imagePoints);
	for j=1:k
		R=W2R([W1(j); W2(j); W3(j)]);
		T=[R(1,1) R(1,2) R(1,3) tx(j);R(2,1) R(2,2) R(2,3) ty(j);R(3,1) R(3,2) R(3,3) tz(j)];
		for i=1:r
			Z=[worldPoints(i,1);worldPoints(i,2);0;1];

			answer=T*Z; %#2.21

			%project to normalised ideal image plane
			out(i,1,j)=answer(1)/answer(3);
			out(i,2,j)=answer(2)/answer(3);
		end
	end
end

function out=radial_distortion(points,k1,k2)
	% function to apply radial distortion to points
	[r,c,k]=size(points);
	for j=1:k
		for i=1:r
			r1=points(i,1,j)^2+points(i,2,j)^2; % #3.47
			Dapply=1+k1*r1+k2*r1^2;	% #3.47
			out(i,1,j)=points(i,1,j)*Dapply; % #3.46
			out(i,2,j)=points(i,2,j)*Dapply; % #3.46
		end
	end
end

function out=affine_transformation(cx,cy,fx,fy,fs,points)
	% function to convert coordinates on the image plane to coordinates on the sensor plane
	A=[fx, fs, cx; 0 fy, cy];
	[r,c,k]=size(points);
	for j=1:k
		for i =1:r
			X=[points(i,1,j); points(i,2,j); 1];
			temp=A*X; % #2.25
			out(i,1,j)=temp(1);
			out(i,2,j)=temp(2);
		end
	end
end

function [k1,k2]=est_radial_distortion(k1,k2,cx,cy,imagePoints,modelPoints,posPoints)
	% function to estimate radial distortion
	[rr,c,k]=size(imagePoints);
	count=1;
	for j=1:k
		for i=1:rr
			r=(posPoints(i,1,j))^2+(posPoints(i,2,j))^2; % #3.47
			D2(2*count-1,1)=imagePoints(i,1,j)-modelPoints(i,1,j); % #3.48
			D2(2*count,1)=imagePoints(i,2,j)-modelPoints(i,2,j); % #3.48
			D(2*count-1,1)=(imagePoints(i,1,j)-cx)*r; % creating matrix of equation #3.55
			D(2*count-1,2)=(imagePoints(i,1,j)-cx)*r^2; % creating matrix of equation #3.55
			D(2*count,1)=(imagePoints(i,2,j)-cy)*r; % creating matrix of equation #3.55
			D(2*count,2)=(imagePoints(i,2,j)-cy)*r^2; % creating matrix of equation #3.55
			count=count+1;
		end
	end
	ktemp=D\D2; % solving equation #3.55
	k1=ktemp(1);
	k2=ktemp(2);
end


function errorval=lsqcontrol(vec,other1,worldPoints,imagePoints,other2)
	% function on which non-linear least squares optimisation is applied
	[r,c,k]=size(imagePoints);
	% This for loop containing if statements controls which parameters are optimised and which
	% aren't according to the vector passed to the function
	for i=1:max(size(vec))
		if vec(i)==1
			if i==1
				cx=other1(1);
			elseif i==2
				cy=other1(2);
			elseif i==3
				fx=other1(3);
			elseif i==4
				fy=other1(4);
			elseif i==5
				fs=other1(5);
			elseif i==6
				tz=other1(6:5+k);
			elseif i==7
				tx=other1(6+k:5+2*k);
			elseif i==8
				ty=other1(6+2*k:5+3*k);
			elseif i==9
				W1=other1(6+3*k:5+4*k);
			elseif i==10
				W2=other1(6+4*k:5+5*k);
			elseif i==11
				W3=other1(6+5*k:5+6*k);
			elseif i==12
				k1=other1(6+6*k);
			elseif i==13
				k2=other1(7+6*k);
			end
		else
			if i==1
				cx=other2(1);
			elseif i==2
				cy=other2(2);
			elseif i==3
				fx=other2(3);
			elseif i==4
				fy=other2(4);
			elseif i==5
				fs=other2(5);
			elseif i==6
				tz=other2(6:5+k);
			elseif i==7
				tx=other2(6+k:5+2*k);
			elseif i==8
				ty=other2(6+2*k:5+3*k);
			elseif i==9
				W1=other2(6+3*k:5+4*k);
			elseif i==10
				W2=other2(6+4*k:5+5*k);
			elseif i==11
				W3=other2(6+5*k:5+6*k);
			elseif i==12
				k1=other2(6+6*k);
			elseif i==13
				k2=other2(7+6*k);
			end
		end
	end

	% convert from world coordinate points to image plane coordinates
	PlanePoints=Step1W(tz,tx,ty,W1,W2,W3,worldPoints,imagePoints);
	% determine the distorted image plane coordinates using the distortion parameters
	DistortedPoints=radial_distortion(PlanePoints,k1,k2);
	% determine the sensor plane coordinates from the distorted image plane coordinates
	sensorPoints=affine_transformation(cx,cy,fx,fy,fs,DistortedPoints);

	% calculate a vector that represents the error between the calculated and actual sensor coordinate points
	% this vector is used by the optimiser to improve on the initial guesses
	count=1;
	for i=1:k
		for j=1:r
			for tipples=1:c
				errorval(count)=abs(imagePoints(j,tipples,i)-sensorPoints(j,tipples,i));
				count=count+1;
			end
		end
	end
end

function out=R2W(R)
	% function to convert form the rotation matrix to rotation angles about each axis
	theta=acos((trace(R)-1)/2);
	out=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end

function out=W2R(W)
	% function to convert from the rotation angles about each axis to the rotation matrix
	theta=sqrt(W(1)^2 + W(2)^2 + W(3)^2);
	omega=[0 -W(3) W(2); W(3) 0 -W(1); -W(2) W(1) 0];
	out=eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);
end

function [Homo]=Homographies(worldPoints,imagePoints)
	% function to determine the homographies of all the calibration images
	img_count=size(imagePoints,3);
	Homo=zeros(3,3,img_count);
	for i=1:img_count
		Homo(:,:,i)=Homography(worldPoints,imagePoints(:,:,i));
	end
end

function [T]=Homography(worldPoints,imagePoints)
	% function to determine the homography of a calibration image

	% normalise points for numerical stability
	[worldPoints,NormMat1]=NormalisePoints(worldPoints);
	[imagePoints,NormMat2]=NormalisePoints(imagePoints);

	% create necessary variables
	M=size(imagePoints,1);
	x=imagePoints(:,1);
	y=imagePoints(:,2);
	vector1=ones(M,1);
	vector0=zeros(M,1);
	u=worldPoints(:,1);
	v=worldPoints(:,2);
	% create the matrix of equation #3.12
	X=[x y vector1 vector0 vector0 vector0 -u.*x -u.*y;
		vector0 vector0 vector0 x y vector1 -v.*x -v.*y];
	% create the RHS vector of equation #3.12
	U=[u;v];
	% check for rank deficiency and if it is not rank deficient determine the homography matrix (H')
	if rank(X)>=8
		T_temp=X\U; % solve #3.12
	else
		fprintf('\nWent wrong at Homography');
	end
	% apply constraint on homography matrix to avoid trivial solution
	T_temp(9)=1;
	Tinv=reshape(T_temp,3,3);
	% undo normalization
	Tinv=NormMat2\(Tinv*NormMat1); % #3.15
	T=inv(Tinv);
	% continue to apply constant on the homography matrix
	T=T./T(3,3);
	T=T';
	T=T./T(3,3);
end

function [pointNorm,NormMatInv]=NormalisePoints(points)
	% function to normalise the coordinate points prior to determining the homographies
	N=size(points,1);
	% determine the centroid of the points
	centre=mean(points,1);
	% shift centroid of points (translation step)
	pointNorm = bsxfun(@minus,points,centre);
	% determine the sum of the squared differences of the points from the origin
	SumOfSquareDistances= sum( hypot(pointNorm(:,1),pointNorm(:,2)).^2 );

	% determine the scale factor to be applied to the points to bring their average distance 
	% from the origin is equal to the square root of 2
	if SumOfSquareDistances > 0
	    scaleFactor = sqrt(2*N) / sqrt(SumOfSquareDistances);
	else
	    % If all input control points are at the same location, the denominator
	    % of the scale factor goes to 0. Don't rescale in this case.
	    if isa(points,'single')
	        scaleFactor = single(1);
	    else
	        scaleFactor = 1.0;
	    end
	end
	% apply scale factor
	pointNorm = pointNorm .* scaleFactor;
	% determine the normalisation matrix (#3.13 & #3.14)
	NormMatInv= [...
    1/scaleFactor,     0,            0;...
    0,            1/scaleFactor,     0;...
    centre(1),      centre(2),      1];
end

function [V]=calculateV(H)
	% function to calculate V in equation #3.35
	img_count=size(H,3);
	V = zeros(2 * img_count, 6);
	for i=1:img_count
		Homo=H(:,:,i)';
		V(i*2-1,:) = getv(Homo, 1, 2);
    	V(i*2, :) = getv(Homo, 1, 1) - getv(Homo, 2, 2);
	end
end

function [v]=getv(H,i,j)
	% #3.33
	v=[H(i,1)*H(j,1), H(i,1)*H(j,2)+H(i,2)*H(j,1), H(i,2)*H(j,2),...
    H(i,3)*H(j,1)+H(i,1)*H(j,3), H(i,3)*H(j,2)+H(i,2)*H(j,3), H(i,3)*H(j,3)];
end

function [B]=calculateB(V)
	% function to calculate B by solving #3.35
	[~, ~, U] = svd(V);
	b = U(:, end);
	B = [b(1), b(2), b(4); b(2), b(3), b(5); b(4), b(5), b(6)];
end

function A=intrinsicParam(B)
	% function to determine the intrinsic parameters 
	cy = (B(1,2)*B(1,3) - B(1,1)*B(2,3)) / (B(1,1)*B(2,2)-B(1,2)^2); % #3.36
	lambda = B(3,3) - (B(1,3)^2 + cy * (B(1,2)*B(1,3) - B(1,1)*B(2,3))) / B(1,1); % #3.37
	fx = sqrt(lambda / B(1,1)); % #3.38
	fy = sqrt(lambda * B(1,1) / (B(1,1) * B(2,2) - B(1,2)^2)); % #3.39
	skew = -B(1,2) * fx^2 * fy / lambda; % #3.40
	cx = skew * cy / fx - B(1,3) * fx^2 / lambda; % #3.41
	A = [fx, skew, cx; ...
     0, fy, cy; ...
     0, 0, 1];
end

function [rot,trans]=extrinsicParam(Homo,A)
	% function to determine the extrinsic parameters
	img_count=size(Homo, 3);
	rot=zeros(3,3,img_count);
	trans=zeros(3,img_count); 
	Ainv = inv(A);
	for i = 1:img_count;
	    H = Homo(:, :, i);
	    h1 = H(:, 1);
	    h2 = H(:, 2);
	    h3 = H(:, 3);
	    lambda = 1 / norm(Ainv * h1); 
	    
	    % 3D rotation matrix
	    r1 = lambda * Ainv * h1;  % #3.42
	    r2 = lambda * Ainv * h2;  % #3.43
	    r3 = cross(r1, r2); % #3.44
	    R = [r1,r2,r3];

	    rot(:,:,i)=R;

	    t = lambda * Ainv * h3;  % #3.45
	    trans(:, i) = t;
	end
	trans = trans';
end

function out=worldPointsvalues
	% calibration target coordinates in world coordinate system
	out=[110 0;
	100 0;
	90 0;
	80 0;
	70 0;
	60 0;
	50 0;
	40 0;
	30 0;
	20 0;
	10 0;
	0 0;
	110 10;
	100 10;
	90 10;
	80 10;
	70 10;
	60 10;
	50 10;
	40 10;
	30 10;
	20 10;
	10 10;
	0 10;
	110 20;
	100 20;
	90 20;
	80 20;
	70 20;
	60 20;
	50 20;
	40 20;
	30 20;
	20 20;
	10 20;
	0 20;
	110 30;
	100 30;
	90 30;
	80 30;
	70 30;
	60 30;
	50 30;
	40 30;
	30 30;
	20 30;
	10 30;
	0 30;
	110 40;
	100 40;
	90 40;
	80 40;
	70 40;
	60 40;
	50 40;
	40 40;
	30 40;
	20 40;
	10 40;
	0 40;
	110 50;
	100 50;
	90 50;
	80 50;
	70 50;
	60 50;
	50 50;
	40 50;
	30 50;
	20 50;
	10 50;
	0 50;
	110 60;
	100 60;
	90 60;
	80 60;
	70 60;
	60 60;
	50 60;
	40 60;
	30 60;
	20 60;
	10 60;
	0 60;
	110 70;
	100 70;
	90 70;
	80 70;
	70 70;
	60 70;
	50 70;
	40 70;
	30 70;
	20 70;
	10 70;
	0 70;
	110 80;
	100 80;
	90 80;
	80 80;
	70 80;
	60 80;
	50 80;
	40 80;
	30 80;
	20 80;
	10 80;
	0 80;
	110 90;
	100 90;
	90 90;
	80 90;
	70 90;
	60 90;
	50 90;
	40 90;
	30 90;
	20 90;
	10 90;
	0 90;
	110 100;
	100 100;
	90 100;
	80 100;
	70 100;
	60 100;
	50 100;
	40 100;
	30 100;
	20 100;
	10 100;
	0 100;
	110 110;
	100 110;
	90 110;
	80 110;
	70 110;
	60 110;
	50 110;
	40 110;
	30 110;
	20 110;
	10 110;
	0 110;
	110 120;
	100 120;
	90 120;
	80 120;
	70 120;
	60 120;
	50 120;
	40 120;
	30 120;
	20 120;
	10 120;
	0 120];
end