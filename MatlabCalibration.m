function varargout=MatlabCalibration(varargin)
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

	% for i=1:max(size(imageNumbers))
	% 	figure
	% 	% Display detected points
	% 	J = insertText(I{i}, imagePoints(:,:,i), 1:size(imagePoints(:,:,i), 1));
	% 	J = insertMarker(J, imagePoints(:,:,i), 'o', 'Color', 'red', 'Size', 5);
	% 	imshow(J);
	% 	title(sprintf('Detected a %d x %d Checkerboard', boardSize));
	% end

	% [r,c,k]=size(imagePoints);
	% worldPoints=zeros(r,c)
	% vals=0:10:120;
	% for i=1:r
	% 	if r/12<1
	% 		adx=adx+10
	% 	end
	% 	worldPoints(j,i)=
		worldPoints=worldPointsvalues;
		% params = estimateCameraParameters(imagePoints, worldPoints);
		params=estimateCameraParametersCustom(imagePoints, worldPoints);
		Homo=Homographies(worldPoints,imagePoints);
		V=calculateV(Homo);
		B=calculateB(V);
		A=intrinsicParam(B);
		[rot,trans]=extrinsicParam(Homo,A);
		[v0,u0,lambda,a,b,g]=intrinsicParam2(B);
		params.IntrinsicMatrix
		% params



	% 	m_undist=Distorted2StraightLines(imagePoints);
	% 	for i=1:max(size(imageNumbers))
	% 		[cx(i),cy(i),fx(i),fy(i),tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i)]=world2camera(m_undist(:,:,i),worldPoints)
	% 	end
	% 	% for i=1:max(size(imageNumbers))
	% 	% 	[cx(i),cy(i),fx(i),fy(i),tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i)]=world2camera(imagePoints(:,:,i),worldPoints)
	% 	% end
	% 	fx=mean(nonzeros(fx));
	% 	fy=mean(nonzeros(fy));
	% 	cx=mean(nonzeros(cx));
	% 	cy=mean(nonzeros(cy));
	% 	fs=0.0001;

	% [r,c,k]=size(imagePoints);
	% 	lsqin(1)=cx;
	% 	lsqin(2)=cy;
	% 	lsqin(3)=fx;
	% 	lsqin(4)=fy;
	% 	lsqin(5)=0; %fs
	% 	lsqin(6:5+k)=tz;
	% 	lsqin(6+k:5+2*k)=tx;
	% 	lsqin(6+2*k:5+3*k)=ty;
	% 	lsqin(6+3*k:5+4*k)=R11;
	% 	lsqin(6+4*k:5+5*k)=R12;
	% 	lsqin(6+5*k:5+6*k)=R13;
	% 	lsqin(6+6*k:5+7*k)=R21;
	% 	lsqin(6+7*k:5+8*k)=R22;
	% 	lsqin(6+8*k:5+9*k)=R23;
	% 	lsqin(6+9*k:5+10*k)=R31;
	% 	lsqin(6+10*k:5+11*k)=R32;
	% 	lsqin(6+11*k:5+12*k)=R33;
	% 	lsqin(6+12*k)=0.00001;
	% 	lsqin(7+12*k)=0;
	% 	lsqin(8+12*k)=0;
	% 	lsqin(9+12*k)=0;
	% 	lsqin(10+12*k)=0;

	% 	% lsqin(1)=cx;
	% 	% lsqin(2)=cy;
	% 	% lsqin(3)=fx;
	% 	% lsqin(4)=fy;
	% 	% lsqin(5)=0.001;
	% 	% lsqin(6:5+k)=zeros(1,k);
	% 	% lsqin(6+k:5+2*k)=zeros(1,k);
	% 	% lsqin(6+2*k:5+3*k)=zeros(1,k);
	% 	% lsqin(6+3*k:5+4*k)=zeros(1,k);
	% 	% lsqin(6+4*k:5+5*k)=zeros(1,k);
	% 	% lsqin(6+5*k:5+6*k)=zeros(1,k);
	% 	% lsqin(6+6*k:5+7*k)=zeros(1,k);
	% 	% lsqin(6+7*k:5+8*k)=zeros(1,k);
	% 	% lsqin(6+8*k:5+9*k)=zeros(1,k);
	% 	% lsqin(6+9*k:5+10*k)=zeros(1,k);
	% 	% lsqin(6+10*k:5+11*k)=zeros(1,k);
	% 	% lsqin(6+11*k:5+12*k)=zeros(1,k);
	% 	% lsqin(6+12*k)=0;
	% 	% lsqin(7+12*k)=0;
	% 	% lsqin(8+12*k)=0;
	% 	% lsqin(9+12*k)=0;
	% 	% lsqin(10+12*k)=0;

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
		% lsqin(10+12*k)=0.0001;


		


		function errorval=particleswarnfun(other1)
			[r,c,k]=size(imagePoints);
			for i=1:21
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
					R11=other1(6+3*k:5+4*k);
				elseif i==10
					R12=other1(6+4*k:5+5*k);
				elseif i==11
					R13=other1(6+5*k:5+6*k);
				elseif i==12
					R21=other1(6+6*k:5+7*k);
				elseif i==13
					R22=other1(6+7*k:5+8*k);
				elseif i==14
					R23=other1(6+8*k:5+9*k);
				elseif i==15
					R31=other1(6+9*k:5+10*k);
				elseif i==16
					R32=other1(6+10*k:5+11*k);
				elseif i==17
					R33=other1(6+11*k:5+12*k);
				elseif i==18
					k1=other1(1);
				elseif i==19
					k2=other1(2);
				elseif i==20
					k3=other1(3);
				elseif i==21
					k4=other1(4);
				end
			end

			for i=1:k
				out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
			end

			for i=1:k
				m(:,:,i)=out{i};
			end

			m=Distortion2(m,cx,cy,k1,k2,k3,k4);

			for i=1:k
				out{i}=m(:,:,i);
			end

			% errorval=0;
			% for i=1:k
			% 	for j=1:r
			% 		errorval=errorval+(imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2;
			% 	end
			% end
			errorval=0;
			for i=1:k
				for j=1:r
					for tipples=1:c
						errorval=errorval+(abs(imagePoints(j,tipples,i)-out{i}(j,tipples)))^2;
					end
				end
			end
		end

		function errorval=lsqagain(other)
			[r,c,k]=size(imagePoints);
			cx=other(1);
			cy=other(2);
			fx=other(3);
			fy=other(4);
			fs=other(5);
			tz=other(6:5+k);
			tx=other(6+k:5+2*k);
			ty=other(6+2*k:5+3*k);
			R11=other(6+3*k:5+4*k);
			R12=other(6+4*k:5+5*k);
			R13=other(6+5*k:5+6*k);
			R21=other(6+6*k:5+7*k);
			R22=other(6+7*k:5+8*k);
			R23=other(6+8*k:5+9*k);
			R31=other(6+9*k:5+10*k);
			R32=other(6+10*k:5+11*k);
			R33=other(6+11*k:5+12*k);
			k1=other(6+12*k);
			k2=other(7+12*k);
			k3=other(8+12*k);
			k4=other(9+12*k);


			for i=1:k
				out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
			end

			for i=1:k
				m(:,:,i)=out{i};
			end

			m=Distortion2(m,cx,cy,k1,k2,k3,k4);

			for i=1:k
				out{i}=m(:,:,i);
			end

			errorval=0;
			for i=1:k
				for j=1:r
					errorval=errorval+(imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2;
				end
			end

			% count=1;
			% for i=1:k
			% 	for j=1:r
			% 		for tipples=1:c
			% 			errorval(count,1)=abs(imagePoints(j,tipples,i)-out{i}(j,tipples));
			% 			count=count+1;
			% 		end
			% 	end
			% end

		end

		current_path=pwd;
		addpath(strcat(current_path,'\levenberg-Marquardt toolbox'));
		addpath(strcat(current_path,'\Jacobian'));
		% [x,S,cnt]=LMFsolve(@lsqagain,lsqin)
		% opts3.TolX=1e-40;
		% opts3.RelTolX=1e-25;
		% opts3.Display='iter';
		checkthis=LevenbergMarquardt(@lsqagain,lsqin)


		% opts2=psooptimset('MaxIterations',300,'UseParallel',true);

		% current_path=pwd;
		% addpath(strcat(current_path,'\PSO'));
		% checkthis=pso(@particleswarnfun,4)
		% lsqin(end-3:end)=checkthis;

		% opts=optimoptions('particleswarm','Display','iter','UseParallel',true)
		% checkthis=particleswarm(@particleswarnfun,(9+12*k))
		% lsqin=checkthis;
		opts=optimset('Display','iter');
		% checkthis=fminsearch(@particleswarnfun,lsqin,opts)


		improve=[	%0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1;
					% 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0;
					% 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
					% 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1;
					% % 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
					% % 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0;
					% 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
					% % 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
					% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1;
					1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
		[improve_r,improve_c]=size(improve);
		lsqin2=lsqin;
		options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','FunctionTolerance',1e-40,'StepTolerance',1e-999,'MaxFunctionEvaluations',1e6,'MaxIterations',100,'UseParallel',true);
		tic
		for i=1:improve_r
			disp(improve(i,:))
			out=lsqnonlin(@(lsqin) lsqcontrol(improve(i,:),lsqin,worldPoints,imagePoints,lsqin2),lsqin,[],[],options)
			lsqin=out;
			lsqin2=out;
		end
		toc

		% checkthis=fminsearch(@particleswarnfun,lsqin,opts)
		lsqin(6+12*k:9+12*k)=out(6+12*k:9+12*k);




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

function [cx,cy,fx,fy,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33]=world2camera(imagePoints,worldPoints)
	[r,c,k]=size(imagePoints);
	for i=1:r
		M(i,:)=[1 0 worldPoints(i,1) worldPoints(i,2) 0 0 -imagePoints(i,1)*worldPoints(i,1) -imagePoints(i,1)*worldPoints(i,2)];
		M(i+1,:)=[0 1 0 0 worldPoints(i,1) worldPoints(i,2) -imagePoints(i,2)*worldPoints(i,1) -imagePoints(i,2)*worldPoints(i,2)];
		b(i,1)=imagePoints(i,1);
		b(i+1,1)=imagePoints(i,2);
	end
	% n=pinv(M)*b
	n=M\b

	R=[-n(3)*n(4);(n(4)^2-n(3)^2)];
	L=[(n(3)*n(8)+n(4)*n(7)), n(5)*n(6), (n(5)*n(8)+n(6)*n(7)), n(7)*n(8);
	2*(n(3)*n(7)-n(4)*n(8)), (n(5)^2-n(6)^2),2*(n(5)*n(7)-n(6)*n(8)), (n(7)^2-n(8)^2)];
	epsilon=L\R;
	cx=epsilon(1);
	cy=epsilon(2);
	fx=sqrt(epsilon(4)-epsilon(3)^2/epsilon(2)-epsilon(1)^2);
	fy=sqrt((epsilon(4)-epsilon(3)^2/epsilon(2)-epsilon(1)^2)/epsilon(2));

	tz=sqrt(2*(((n(3)-n(7)*cx)^2+(n(4)-n(8)*cx)^2)/(fx^2)+((n(5)-n(7)*cy)^2+(n(6)-n(8)*cy)^2)/(fy^2)+n(7)^2+n(8)^2)^-1);
	tx=(n(1)-cx)*tz/fx;
	ty=(n(2)-cy)*tz/fy;

	R11=tz*(n(3)-n(7)*cx)/fx;
	R12=tz*(n(4)-n(8)*cx)/fx;
	R21=tz*(n(5)-n(7)*cy)/fy;
	R22=tz*(n(6)-n(8)*cy)/fy;
	R31=tz*n(7);
	R32=tz*n(8);
	R33=-R11*R22+R21*R12;
	R23=R11*R32-R31*R12;
	R13=-R21*R32+R31*R22;
	[cx,cy,fx,fy,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33]=makereal(cx,cy,fx,fy,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33);
end

function out=world2sensor(cx,cy,fx,fy,fs,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints)
	% convert from the world coordinate system to the sensor plane coordinate system (M->m)
	K=[fx,fs,cx,0;0,fy,cy,0;0,0,1,0];
	T=[R11 R12 R13 tx;R21 R22 R23 ty;R31 R32 R33 tz; 0 0 0 1];
	for i=1:max(size(worldPoints))
		Z=[worldPoints(i,1);worldPoints(i,2);0;1];
		answer=K*T*Z;
		out(i,:)=answer';
	end
end

function out=Distortion(m,cx,cy,k1,k2,k3,k4,k5)
	%m->m~
	[r,c,k]=size(m);
	for i=1:k
		for j=1:r
			row=sqrt((m(j,1,i)-cx)^2+(m(j,2,i)-cy)^2);
			ruv=[(m(j,1,i)-cx);(m(j,2,i)-cy)];	%radial unit vector
			ruv=norm(ruv);
			tuv=null(ruv(:).'); % (.') - transpose even for complex numbers; otherwise would give the complex conjugate
			tuv=norm(tuv);
			%radial
			Drad=k1*row^3*ruv+k2*row^5*ruv+k3*row^7*ruv;
			%de-centering
			theta=atan2(m(j,1,i),m(j,2,i));
			Dspher=k4*row^2*(3*sin(theta-k5)*ruv+cos(theta-k5)*tuv);

			

			answer=m(j,:,i)'+Dspher+Drad;
			out(j,1,i)=answer(1);
			out(j,2,i)=answer(2);
		end
	end

end

function out=Distortion2(m,cx,cy,k1,k2,k3,k4)
	%m->m~
	[r,c,k]=size(m);
	for i=1:k
		for j=1:r
			row=sqrt((m(j,1,i)-cx)^2+(m(j,2,i)-cy)^2);
			ruv=[(m(j,1,i)-cx);(m(j,2,i)-cy)];	%radial unit vector
			ruv=norm(ruv);
			tuv=null(ruv(:).'); % (.') - transpose even for complex numbers; otherwise would give the complex conjugate
			tuv=norm(tuv);
			%radial
			Drad=k1*row^3*ruv+k2*row^5*ruv;
			%de-centering
			theta=atan2(m(j,1,i),m(j,2,i));
			Dspher=k3*row^2*(3*sin(theta-k4)*ruv+cos(theta-k4)*tuv);

			answer=m(j,:,i)'+Dspher+Drad;
			out(j,1,i)=answer(1);
			out(j,2,i)=answer(2);
		end
	end

end

function errorval=lsqfun(var,worldPoints,imagePoints)
	[r,c,k]=size(imagePoints);
	cx=var(1);
	cy=var(2);
	fx=var(3);
	fy=var(4);
	fs=var(5);
	tz=var(6:5+k);
	tx=var(6+k:5+2*k);
	ty=var(6+2*k:5+3*k);
	R11=var(6+3*k:5+4*k);
	R12=var(6+4*k:5+5*k);
	R13=var(6+5*k:5+6*k);
	R21=var(6+6*k:5+7*k);
	R22=var(6+7*k:5+8*k);
	R23=var(6+8*k:5+9*k);
	R31=var(6+9*k:5+10*k);
	R32=var(6+10*k:5+11*k);
	R33=var(6+11*k:5+12*k);
	k1=var(6+12*k);
	k2=var(7+12*k);
	k3=var(8+12*k);
	k4=var(9+12*k);
	k5=var(10+12*k);


	for i=1:k
		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
	end

	for i=1:k
		m(:,:,i)=out{i};
	end

	m=Distortion(m,cx,cy,k1,k2,k3,k4,k5);

	for i=1:k
		out{i}=m(:,:,i);
	end

	errorval=0;
	for i=1:k
		for j=1:r
			errorval=errorval+norm((imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2);
		end
	end

end

function errorval=lsqfun2(var,worldPoints,imagePoints)
	[r,c,k]=size(imagePoints);
	cx=var(1);
	cy=var(2);
	fx=var(3);
	fy=var(4);
	fs=var(5);
	tz=var(6:5+k);
	tx=var(6+k:5+2*k);
	ty=var(6+2*k:5+3*k);
	R11=var(6+3*k:5+4*k);
	R12=var(6+4*k:5+5*k);
	R13=var(6+5*k:5+6*k);
	R21=var(6+6*k:5+7*k);
	R22=var(6+7*k:5+8*k);
	R23=var(6+8*k:5+9*k);
	R31=var(6+9*k:5+10*k);
	R32=var(6+10*k:5+11*k);
	R33=var(6+11*k:5+12*k);
	k1=var(6+12*k);
	k2=var(7+12*k);
	k3=var(8+12*k);
	k4=var(9+12*k);


	for i=1:k
		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
	end

	for i=1:k
		m(:,:,i)=out{i};
	end

	m=Distortion2(m,cx,cy,k1,k2,k3,k4);

	for i=1:k
		out{i}=m(:,:,i);
	end

	errorval=0;
	for i=1:k
		for j=1:r
			errorval=errorval+norm((imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2);
		end
	end

end

function errorval=lsqfundist(var,worldPoints,imagePoints,other)
	[r,c,k]=size(imagePoints);
	cx=other(1);
	cy=other(2);
	fx=other(3);
	fy=other(4);
	fs=other(5);
	tz=other(6:5+k);
	tx=other(6+k:5+2*k);
	ty=other(6+2*k:5+3*k);
	R11=other(6+3*k:5+4*k);
	R12=other(6+4*k:5+5*k);
	R13=other(6+5*k:5+6*k);
	R21=other(6+6*k:5+7*k);
	R22=other(6+7*k:5+8*k);
	R23=other(6+8*k:5+9*k);
	R31=other(6+9*k:5+10*k);
	R32=other(6+10*k:5+11*k);
	R33=other(6+11*k:5+12*k);
	k1=var(6+12*k);
	k2=var(7+12*k);
	k3=var(8+12*k);
	k4=var(9+12*k);


	for i=1:k %parfor
		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
	end

	for i=1:k
		m(:,:,i)=out{i};
	end

	m=Distortion2(m,cx,cy,k1,k2,k3,k4);

	for i=1:k
		out{i}=m(:,:,i);
	end

	% errorval=0;
	% for i=1:k
	% 	for j=1:r
	% 		errorval=errorval+(imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2;
	% 	end
	% end

	count=1;
	for i=1:k
		for j=1:r
			for tipples=1:c
				errorval(count)=abs(imagePoints(j,tipples,i)-out{i}(j,tipples));
				count=count+1;
			end
		end
	end

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
				R11=other1(6+3*k:5+4*k);
			elseif i==10
				R12=other1(6+4*k:5+5*k);
			elseif i==11
				R13=other1(6+5*k:5+6*k);
			elseif i==12
				R21=other1(6+6*k:5+7*k);
			elseif i==13
				R22=other1(6+7*k:5+8*k);
			elseif i==14
				R23=other1(6+8*k:5+9*k);
			elseif i==15
				R31=other1(6+9*k:5+10*k);
			elseif i==16
				R32=other1(6+10*k:5+11*k);
			elseif i==17
				R33=other1(6+11*k:5+12*k);
			elseif i==18
				k1=other1(6+12*k);
			elseif i==19
				k2=other1(7+12*k);
			elseif i==20
				k3=other1(8+12*k);
			elseif i==21
				k4=other1(9+12*k);
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
				R11=other2(6+3*k:5+4*k);
			elseif i==10
				R12=other2(6+4*k:5+5*k);
			elseif i==11
				R13=other2(6+5*k:5+6*k);
			elseif i==12
				R21=other2(6+6*k:5+7*k);
			elseif i==13
				R22=other2(6+7*k:5+8*k);
			elseif i==14
				R23=other2(6+8*k:5+9*k);
			elseif i==15
				R31=other2(6+9*k:5+10*k);
			elseif i==16
				R32=other2(6+10*k:5+11*k);
			elseif i==17
				R33=other1(6+11*k:5+12*k);
			elseif i==18
				k1=other2(6+12*k);
			elseif i==19
				k2=other2(7+12*k);
			elseif i==20
				k3=other2(8+12*k);
			elseif i==21
				k4=other2(9+12*k);
			end
		end
	end

	parfor i=1:k
		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
	end

	for i=1:k
		m(:,:,i)=out{i};
	end

	m=Distortion2(m,cx,cy,k1,k2,k3,k4);

	for i=1:k
		out{i}=m(:,:,i);
	end

	% errorval=0;
	% for i=1:k
	% 	for j=1:r
	% 		errorval=errorval+(imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2;
	% 	end
	% end
	count=1;
	for i=1:k
		for j=1:r
			for tipples=1:c
				errorval(count)=abs(imagePoints(j,tipples,i)-out{i}(j,tipples));
				count=count+1;
			end
		end
	end


end

% function errorval=particleswarnfun(vec,other1,worldPoints,imagePoints,other2)
% 	[r,c,k]=size(imagePoints);
% 	for i=1:max(size(vec))
% 		if vec(i)==1
% 			if i==1
% 				cx=other1(1);
% 			elseif i==2
% 				cy=other1(2);
% 			elseif i==3
% 				fx=other1(3);
% 			elseif i==4
% 				fy=other1(4);
% 			elseif i==5
% 				fs=other1(5);
% 			elseif i==6
% 				tz=other1(6:5+k);
% 			elseif i==7
% 				tx=other1(6+k:5+2*k);
% 			elseif i==8
% 				ty=other1(6+2*k:5+3*k);
% 			elseif i==9
% 				R11=other1(6+3*k:5+4*k);
% 			elseif i==10
% 				R12=other1(6+4*k:5+5*k);
% 			elseif i==11
% 				R13=other1(6+5*k:5+6*k);
% 			elseif i==12
% 				R21=other1(6+6*k:5+7*k);
% 			elseif i==13
% 				R22=other1(6+7*k:5+8*k);
% 			elseif i==14
% 				R23=other1(6+8*k:5+9*k);
% 			elseif i==15
% 				R31=other1(6+9*k:5+10*k);
% 			elseif i==16
% 				R32=other1(6+10*k:5+11*k);
% 			elseif i==17
% 				R33=other1(6+11*k:5+12*k);
% 			elseif i==18
% 				k1=other1(6+12*k);
% 			elseif i==19
% 				k2=other1(7+12*k);
% 			elseif i==20
% 				k3=other1(8+12*k);
% 			elseif i==21
% 				k4=other1(9+12*k);
% 			end
% 		else
% 			if i==1
% 				cx=other2(1);
% 			elseif i==2
% 				cy=other2(2);
% 			elseif i==3
% 				fx=other2(3);
% 			elseif i==4
% 				fy=other2(4);
% 			elseif i==5
% 				fs=other2(5);
% 			elseif i==6
% 				tz=other2(6:5+k);
% 			elseif i==7
% 				tx=other2(6+k:5+2*k);
% 			elseif i==8
% 				ty=other2(6+2*k:5+3*k);
% 			elseif i==9
% 				R11=other2(6+3*k:5+4*k);
% 			elseif i==10
% 				R12=other2(6+4*k:5+5*k);
% 			elseif i==11
% 				R13=other2(6+5*k:5+6*k);
% 			elseif i==12
% 				R21=other2(6+6*k:5+7*k);
% 			elseif i==13
% 				R22=other2(6+7*k:5+8*k);
% 			elseif i==14
% 				R23=other2(6+8*k:5+9*k);
% 			elseif i==15
% 				R31=other2(6+9*k:5+10*k);
% 			elseif i==16
% 				R32=other2(6+10*k:5+11*k);
% 			elseif i==17
% 				R33=other1(6+11*k:5+12*k);
% 			elseif i==18
% 				k1=other2(6+12*k);
% 			elseif i==19
% 				k2=other2(7+12*k);
% 			elseif i==20
% 				k3=other2(8+12*k);
% 			elseif i==21
% 				k4=other2(9+12*k);
% 			end
% 		end
% 	end

% 	parfor i=1:k
% 		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
% 	end

% 	for i=1:k
% 		m(:,:,i)=out{i};
% 	end

% 	m=Distortion2(m,cx,cy,k1,k2,k3,k4);

% 	for i=1:k
% 		out{i}=m(:,:,i);
% 	end

% 	% errorval=0;
% 	% for i=1:k
% 	% 	for j=1:r
% 	% 		errorval=errorval+(imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2;
% 	% 	end
% 	% end
% 	count=1;
% 	for i=1:k
% 		for j=1:r
% 			for tipples=1:c
% 				errorval(count)=abs(imagePoints(j,tipples,i)-out{i}(j,tipples));
% 				count=count+1;
% 			end
% 		end
% 	end


% end

function m_undist=Distorted2StraightLines(imagePoints)
	warning off
	[r,c,k]=size(imagePoints);
	maximumx=max(max(imagePoints(:,1,:)));
	minimumx=min(min(imagePoints(:,1,:)));
	maximumy=max(max(imagePoints(:,2,:)));
	minimumy=min(min(imagePoints(:,2,:)));
	for i=1:k
		clear X
		clear Y
		count=1;
		for j=1:13
			X(:,1)=imagePoints(count:count+11,1,i);
			Y=imagePoints(count:count+11,2,i);
			X(:,2)=ones(12,1);
			A(j,:,i)=X\Y;
			count=count+12;
		end
		clear X
		clear Y
		coord=[1 13 25 37 49 61 73 85 97 109 121 133 145];
		for j=1:12
			X(:,1)=imagePoints(coord,1,i);
			Y=imagePoints(coord,2,i);
			X(:,2)=ones(13,1);
			B(j,:,i)=X\Y;
			coord=coord+1;
		end
		count=1;
		for j=1:13
			curves2(1,:)=[(minimumx-100), (maximumx+100)];
			curves2(2,:)=[(A(j,1,i)*(minimumx-100)+A(j,2,i)), (A(j,1,i)*(maximumx+100)+A(j,2,i))];
			for l=1:12
				curves1(1,:)=[(minimumx-100), (maximumx+100)];
				curves1(2,:)=[(B(l,1,i)*(minimumx-100)+B(l,2,i)), (B(l,1,i)*(maximumx+100)+B(l,2,i))];
				out=InterX(curves1,curves2);
				m_undist(count,1,i)=out(1,1);
				m_undist(count,2,i)=out(2,1);
				count=count+1;
			end
		end
	end
	warning on
end

function [Homo]=Homographies(worldPoints,imagePoints)
	img_count=size(imagePoints,3);
	Homo=zeros(3,3,img_count);
	for i=1:img_count
		Homo(:,:,i)=Homography(worldPoints,imagePoints(:,:,i));
	end
end

function [T]=Homography(worldPoints,imagePoints)
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
	Tinv=NormMat2\(Tinv*NormMat1);
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

% function errorval=lsqfun(cx,cy,fx,fy,fs,tz,tx,ty,R11,R12,R13,R21,R22,R23,R31,R32,R33,worldPoints,imagePoints)
% 	[r,c,k]=size(imagePoints)
% 	cx=var(1);
% 	cy=var(2);
% 	fx=var(3);
% 	fy=var(4);
% 	fs=var(5);
% 	tz=var(6:5+k);
% 	tx=var(6+k:5+2*k);
% 	ty=var(6+2*k:5+3*k);
% 	R11=var(6+3*k:5+4*k);
% 	R12=var(6+4*k:5+5*k);
% 	R13=var(6+5*k:5+6*k);
% 	R21=var(6+6*k:5+7*k);
% 	R22=var(6+7*k:5+8*k);
% 	R23=var(6+8*k:5+9*k);
% 	R31=var(6+9*k:5+10*k);
% 	R32=var(6+10*k:5+11*k);
% 	R33=var(6+11*k:5+12*k);


% 	for i=1:k
% 		out{i}=world2sensor(cx,cy,fx,fy,fs,tz(i),tx(i),ty(i),R11(i),R12(i),R13(i),R21(i),R22(i),R23(i),R31(i),R32(i),R33(i),worldPoints);
% 	end

% 	errorval=0;
% 	for i=1:k
% 		for j=1:r
% 			errorval=errorval+norm((imagePoints(j,1,i)-out{i}(j,1))^2+(imagePoints(j,2,i)-out{i}(j,2))^2);
% 		end
% 	end

% end

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