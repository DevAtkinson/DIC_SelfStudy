function actualout=MatlabCalibration(varargin)
	format long
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

	% display calibration images with calibration targets marked
	% for i=1:max(size(imageNumbers))
	% 	figure
	% 	% Display detected points
	% 	J = insertText(I{i}, imagePoints(:,:,i), 1:size(imagePoints(:,:,i), 1));
	% 	J = insertMarker(J, imagePoints(:,:,i), 'o', 'Color', 'red', 'Size', 5);
	% 	imshow(J);
	% 	title(sprintf('Detected a %d x %d Checkerboard', boardSize));
	% end

	% assign world coordinates of the calibration targets
	worldPoints=worldPointsvalues;
	% params = estimateCameraParameters(imagePoints, worldPoints);
	% params=estimateCameraParametersCustom(imagePoints, worldPoints);
	Homo=Homographies(worldPoints,imagePoints);
	V=calculateV(Homo);
	B=calculateB(V);
	A=intrinsicParam(B);
	[rot,trans]=extrinsicParam(Homo,A);
	[v0,u0,lambda,a,b,g]=intrinsicParam2(B);
	% params.IntrinsicMatrix
	% params


	[rr,cc,kk]=size(rot);
	for i=1:kk
		Wtemp=R2W(rot(:,:,i));
		W(i,1)=Wtemp(1);
		W(i,2)=Wtemp(2);
		W(i,3)=Wtemp(3);
	end


	

	[r,c,k]=size(imagePoints);
	lsqin(1)=A(1,3);
	lsqin(2)=A(2,3);
	lsqin(3)=A(1,1);
	lsqin(4)=A(2,2);
	lsqin(5)=A(1,2); %fs
	lsqin(6:5+k)=trans(:,3);
	lsqin(6+k:5+2*k)=trans(:,1);
	lsqin(6+2*k:5+3*k)=trans(:,2);
	lsqin(6+3*k:5+4*k)=rot(1,1,:);
	lsqin(6+4*k:5+5*k)=rot(1,2,:);
	lsqin(6+5*k:5+6*k)=rot(1,3,:);
	lsqin(6+6*k:5+7*k)=rot(2,1,:);
	lsqin(6+7*k:5+8*k)=rot(2,2,:);
	lsqin(6+8*k:5+9*k)=rot(2,3,:);
	lsqin(6+9*k:5+10*k)=rot(3,1,:);
	lsqin(6+10*k:5+11*k)=rot(3,2,:);
	lsqin(6+11*k:5+12*k)=rot(3,3,:);
	lsqin(6+12*k)=0;%0.01;
	lsqin(7+12*k)=0;%0.000000000000839;
	lsqin(8+12*k)=0;%0.01;
	lsqin(9+12*k)=0;%0.01;


	cx=lsqin(1);
	cy=lsqin(2);
	fx=lsqin(3);
	fy=lsqin(4);
	fs=lsqin(5);
	tz=lsqin(6:5+k);
	tx=lsqin(6+k:5+2*k);
	ty=lsqin(6+2*k:5+3*k);
	R11=lsqin(6+3*k:5+4*k);
	R12=lsqin(6+4*k:5+5*k);
	R13=lsqin(6+5*k:5+6*k);
	R21=lsqin(6+6*k:5+7*k);
	R22=lsqin(6+7*k:5+8*k);
	R23=lsqin(6+8*k:5+9*k);
	R31=lsqin(6+9*k:5+10*k);
	R32=lsqin(6+10*k:5+11*k);
	R33=lsqin(6+11*k:5+12*k);
	k1=lsqin(6+12*k);
	k2=lsqin(7+12*k);
	k3=lsqin(8+12*k);
	k4=lsqin(9+12*k);



	PlanePoints=Step1(tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints);
	sensorPoints=affine_transformation(cx,cy,fx,fy,fs,PlanePoints);
	[k1,k2]=est_radial_distortion(k1,k2,cx,cy,imagePoints,sensorPoints,PlanePoints)


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

	improve=[0 0 0 0 0 1 1 1 1 1 1 0 0;
			1 1 1 1 1 1 1 1 1 1 1 1 1];
	[improve_r,improve_c]=size(improve);

	clear lsqin
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

	% checkthis=fminsearch(@particleswarnfun,lsqin,opts)
	% lsqin(6+12*k:9+12*k)=out(6+12*k:9+12*k);




	% out=lsqnonlin(@(lsqin) lsqfun2(lsqin,worldPoints,imagePoints),lsqin,[],[],options)

	% out=lsqnonlin(@(cx,cy,fx,fy,fs,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33) lsqfun(cx,cy,fx,fy,fs,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints),[cx,cy,fx,fy,fs,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33])
% 	% m_undist=Distorted2StraightLines(imagePoints)

% 	for i=1:max(size(imageNumbers))
% 		figure
% 		% Display detected points
% 		J = insertText(I{i}, m_undist(:,:,i), 1);
% 		J = insertMarker(J, m_undist(:,:,i), 'o', 'Color', 'red', 'Size', 5);
% 		imshow(J);
% 		title(sprintf('Detected a %d x %d Checkerboard', boardSize));
% 	end

% 	[cameraParams, imagesUsed, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints)
% 	CamParam.cameraParams=cameraParams;
% 	CamParam.imagesUsed=imagesUsed;
% 	CamParam.estimationErrors=estimationErrors;
% 	save('CamParam','CamParam')
end

function out=Step1(tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints)
	[r,c,k]=size(imagePoints);
	for j=1:k
		T=[R11(j) R12(j) R13(j) tx(j);R21(j) R22(j) R23(j) ty(j);R31(j) R32(j) R33(j) tz(j)];
		for i=1:r
			Z=[worldPoints(i,1);worldPoints(i,2);0;1];
			% size(T)
			% size(Z)
			answer=T*Z;
			% out(i,:)=answer';
			%project to normalised ideal image plane
			out(i,1,j)=answer(1)/answer(3);
			out(i,2,j)=answer(2)/answer(3);
		end
	end
end

function out=Step1W(tz,tx,ty,W1,W2,W3,worldPoints,imagePoints)
	[r,c,k]=size(imagePoints);
	for j=1:k
		R=W2R([W1(j); W2(j); W3(j)]);
		T=[R(1,1) R(1,2) R(1,3) tx(j);R(2,1) R(2,2) R(2,3) ty(j);R(3,1) R(3,2) R(3,3) tz(j)];
		for i=1:r
			Z=[worldPoints(i,1);worldPoints(i,2);0;1];
			% size(T)
			% size(Z)
			answer=T*Z;
			% out(i,:)=answer';
			%project to normalised ideal image plane
			out(i,1,j)=answer(1)/answer(3);
			out(i,2,j)=answer(2)/answer(3);
		end
	end
end

function out=radial_distortion(points,k1,k2)
	[r,c,k]=size(points);
	for j=1:k
		for i=1:r
			r1=points(i,1,j)^2+points(i,2,j)^2;
			Dapply=1+k1*r1+k2*r1^2;
			out(i,1,j)=points(i,1,j)*Dapply;
			out(i,2,j)=points(i,2,j)*Dapply;
		end
	end
end

function out=affine_transformation(cx,cy,fx,fy,fs,points)
	A=[fx, fs, cx; 0 fy, cy];
	[r,c,k]=size(points);
	for j=1:k
		for i =1:r
			X=[points(i,1,j); points(i,2,j); 1];
			temp=A*X;
			out(i,1,j)=temp(1);
			out(i,2,j)=temp(2);
		end
	end
end

function [k1,k2]=est_radial_distortion(k1,k2,cx,cy,imagePoints,modelPoints,posPoints)
	[rr,c,k]=size(imagePoints);
	count=1;
	for j=1:k
		for i=1:rr
			r=(posPoints(i,1,j))^2+(posPoints(i,2,j))^2;
			D2(2*count-1,1)=imagePoints(i,1,j)-modelPoints(i,1,j);
			D2(2*count,1)=imagePoints(i,2,j)-modelPoints(i,2,j);
			D(2*count-1,1)=(imagePoints(i,1,j)-cx)*r;
			D(2*count-1,2)=(imagePoints(i,1,j)-cx)*r^2;
			D(2*count,1)=(imagePoints(i,2,j)-cy)*r;
			D(2*count,2)=(imagePoints(i,2,j)-cy)*r^2;
			count=count+1;
		end
	end
	ktemp=D\D2;
	k1=ktemp(1);
	k2=ktemp(2);
end


function errorval=lsqcontrol(vec,other1,worldPoints,imagePoints,other2)
	[r,c,k]=size(imagePoints);
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

	PlanePoints=Step1W(tz,tx,ty,W1,W2,W3,worldPoints,imagePoints);
	% size(PlanePoints)
	DistortedPoints=radial_distortion(PlanePoints,k1,k2);
	% size(DistortedPoints)
	sensorPoints=affine_transformation(cx,cy,fx,fy,fs,DistortedPoints);
	% size(sensorPoints)
	% size(imagePoints)

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
	theta=acos((trace(R)-1)/2);
	out=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end

function out=W2R(W)
	theta=sqrt(W(1)^2 + W(2)^2 + W(3)^2);
	omega=[0 -W(3) W(2); W(3) 0 -W(1); -W(2) W(1) 0];
	out=eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);
end

function [Homo]=Homographies(worldPoints,imagePoints)
	img_count=size(imagePoints,3);
	Homo=zeros(3,3,img_count);
	for i=1:img_count
		Homo(:,:,i)=Homography(worldPoints,imagePoints(:,:,i));
	end
end

function [T]=Homography(worldPoints,imagePoints)
	%normalise points for numerical stability
	[worldPoints,NormMat1]=NormalisePoints(worldPoints);
	[imagePoints,NormMat2]=NormalisePoints(imagePoints);

	M=size(imagePoints,1);
	x=imagePoints(:,1);
	y=imagePoints(:,2);
	vector1=ones(M,1);
	vector0=zeros(M,1);
	u=worldPoints(:,1);
	v=worldPoints(:,2);
	X=[x y vector1 vector0 vector0 vector0 -u.*x -u.*y;
		vector0 vector0 vector0 x y vector1 -v.*x -v.*y];
	U=[u;v];
	if rank(X)>=8
		T_temp=X\U;
	else
		fprintf('\nWent wrong at Homography');
	end
	T_temp(9)=1;
	Tinv=reshape(T_temp,3,3);
	Tinv=NormMat2\(Tinv*NormMat1); %undo normalization
	T=inv(Tinv);
	T=T./T(3,3);

	T=T';
	T=T./T(3,3);

	% tform = projective2d(T);
end

function [pointNorm,NormMatInv]=NormalisePoints(points)
	N=size(points,1);
	centre=mean(points,1);
	%shift centroid of points to 
	% pointNorm(:,1)=points(:,1)-centre(1);
	% pointNorm(:,2)=points(:,2)-centre(2);
	pointNorm = bsxfun(@minus,points,centre);
	SumOfSquareDistances= sum( hypot(pointNorm(:,1),pointNorm(:,2)).^2 );

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

	pointNorm = pointNorm .* scaleFactor;
	NormMatInv= [...
    1/scaleFactor,     0,            0;...
    0,            1/scaleFactor,     0;...
    centre(1),      centre(2),      1];
end

function [V]=calculateV(H)
	img_count=size(H,3);
	V = zeros(2 * img_count, 6);
	for i=1:img_count
		Homo=H(:,:,i)';
		V(i*2-1,:) = getv(Homo, 1, 2);
    	V(i*2, :) = getv(Homo, 1, 1) - getv(Homo, 2, 2);
	end
end

function [v]=getv(H,i,j)
	v=[H(i,1)*H(j,1), H(i,1)*H(j,2)+H(i,2)*H(j,1), H(i,2)*H(j,2),...
    H(i,3)*H(j,1)+H(i,1)*H(j,3), H(i,3)*H(j,2)+H(i,2)*H(j,3), H(i,3)*H(j,3)];
end

function [B]=calculateB(V)
	[~, ~, U] = svd(V);
	b = U(:, end);
	B = [b(1), b(2), b(4); b(2), b(3), b(5); b(4), b(5), b(6)];
end

function A=intrinsicParam(B)
	cy = (B(1,2)*B(1,3) - B(1,1)*B(2,3)) / (B(1,1)*B(2,2)-B(1,2)^2);
	lambda = B(3,3) - (B(1,3)^2 + cy * (B(1,2)*B(1,3) - B(1,1)*B(2,3))) / B(1,1);
	fx = sqrt(lambda / B(1,1));
	fy = sqrt(lambda * B(1,1) / (B(1,1) * B(2,2) - B(1,2)^2));
	skew = -B(1,2) * fx^2 * fy / lambda;
	cx = skew * cy / fx - B(1,3) * fx^2 / lambda;
	A = [fx, skew, cx; ...
     0, fy, cy; ...
     0, 0, 1];
end

function [rot,trans]=extrinsicParam(Homo,A)
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
	    r1 = lambda * Ainv * h1; 
	    r2 = lambda * Ainv * h2; 
	    r3 = cross(r1, r2);
	    R = [r1,r2,r3];
	    % rot(:,:,i)=R';
	    rot(:,:,i)=R;
	    
	    % rotationVectors(:, i) = vision.internal.calibration.rodriguesMatrixToVector(R);
	    % rot(:,i)= vision.internal.calibration.rodriguesMatrixToVector(R);
	    % translation vector
	    t = lambda * Ainv * h3;  
	    trans(:, i) = t;
	end

	% rot = rot';
	trans = trans';
end

function [v0,u0,lambda,a,b,g]=intrinsicParam2(B)
	v0=(B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
	lambda=B(3,3)-(B(1,3)^2 + v0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
	a=sqrt(lambda/B(1,1));
	b=sqrt(lambda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
	g=-B(1,2)*a^2*b/lambda;
	u0=g*v0/b-B(1,3)*a^2/lambda;
end

function varargout=makereal(varargin)
	for i=1:nargin
		if ~isreal(varargin{i})
			varargout{i}=0;
		else 
			varargout{i}=varargin{i};
		end
	end
end

function out=worldPointsvalues
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