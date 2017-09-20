function DevanDICtracking(varargin)
	%handle the input variables
	for i=1:nargin/2
		switch varargin{i*2-1}
		case 'folder'
			folder=varargin{i*2};
		case 'ImageName'
			ImageName=varargin{i*2};
		case 'images'
			imageNumbers=varargin{i*2};
		case 'subset size'
			subsize=varargin{i*2};
		case 'subset position'
			subpos=varargin{i*2};
		end
	end
	%create path to specific images
	for i=1:max(size(imageNumbers))
		images=strcat(strcat(ImageName,num2str(imageNumbers{i})),'.tif');
		ImageFileName{i}=fullfile(folder,images);
		I{i} = imread(ImageFileName{i});
	end

	[r,c]=size(I{1});
	F=im2double(I{1}(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1));	%assign subset from reference image
	% figure
	% imagesc(I{1}(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1))
	% [x,y]=ginput(1)
	% figure
	% imagesc(I{2}(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1))
	% [x,y]=ginput(1)
	
	F_interp=im2double(I{1}(subpos(2)-1:subpos(2)+subsize,subpos(1)-1:subpos(1)+subsize));	%assign subset from reference image for interpolation
	
	X=subpos(1):subpos(1)+subsize-1;										% x values over subset
	Y=subpos(2):subpos(2)+subsize-1;										% y values over subset
	x0=subsize/2+subpos(1);													% x value for subset centre
	y0=subsize/2+subpos(2);													% y value for subset centre
	[r_int,c_int]=size(F_interp);

	% for i=2:r_int-1
	% 	for j=2:c_int-1
	% 		Interp_coef{i-1,j-1}=bicubicInterp(F_interp,i,j);
	% 	end
	% end

	[r,c]=size(F);
	sumF=0;
	sumG2=0;
	Fmean=mean(mean(F));
	H=0;
	dF_temp=0;
	[gradx,grady]=gradient(F_interp);
	% gradx=conv2(F_interp,[1,0,-1;0,0,0;1,0,-1]/4,'same');
	% grady=conv2(F_interp,[-1,0,-1;0,0,0;1,0,1]/4,'same');

	%sobel edge operator
	% gradx=conv2(F_interp,[1,0,-1;-2,0,2;1,0,-1]/4,'same');
	% grady=conv2(F_interp,[1,2,1;0,0,0;-1,-2,-1]/4,'same');

	%convolution from paper
	% for i=1:r_int
	% 	gradx(i,:)=conv(F_interp(i,:),[1/12, -8/12,0,8/12,-1/12],'same');
	% end
	% for i=1:c_int
	% 	grady(:,i)=-conv(F_interp(:,i)',[1/12,-8/12,0,8/12,-1/12],'same')';
	% end

	for i=1:r
		for j=1:c
			xp=X(j)-x0;
			yp=Y(i)-y0;
			% dfdx=Interp_coef{i,j}(2,1)+2*Interp_coef{i,j}(3,1)*xp+3*Interp_coef{i,j}(4,1)*xp^2+Interp_coef{i,j}(2,2)*yp+Interp_coef{i,j}(2,3)*yp^2+...
			% Interp_coef{i,j}(2,4)*yp^3+2*Interp_coef{i,j}(3,2)*xp*yp+2*Interp_coef{i,j}(3,3)*xp*yp^2+2*Interp_coef{i,j}(3,4)*xp*yp^3+3*Interp_coef{i,j}(4,2)*xp^2*yp+...
			% 3*Interp_coef{i,j}(4,3)*xp^2*yp^2+3*Interp_coef{i,j}(4,4)*xp^2*yp^3;

			% dfdy=Interp_coef{i,j}(1,2)+2*Interp_coef{i,j}(1,3)*yp+3*Interp_coef{i,j}(1,4)*yp^2+Interp_coef{i,j}(2,2)*xp+2*Interp_coef{i,j}(2,3)*xp*yp+...
			% 3*Interp_coef{i,j}(2,4)*xp*yp^2+Interp_coef{i,j}(3,2)*xp^2+2*Interp_coef{i,j}(3,3)*xp^2*yp+3*Interp_coef{i,j}(3,4)*xp^2*yp^2+Interp_coef{i,j}(4,2)*xp^3+...
			% 2*Interp_coef{i,j}(4,3)*xp^3*yp+3*Interp_coef{i,j}(4,4)*xp^3+yp^2;

			
			dfdx=gradx(i+1,j+1);
			dfdy=grady(i+1,j+1);

			dWdp=[1 xp yp 0 0 0; 0 0 0 1 xp yp];

			dfdw{i,j}=[dfdx, dfdy]*dWdp;
			H=H+dfdw{i,j}'*dfdw{i,j};
			% dF_temp=dF_temp+F(i,j)-Fmean;
			dF_temp=dF_temp+(F(i,j)-Fmean)^2;
		end
	end
	% dF=sqrt(dF_temp^2);
	dF=sqrt(dF_temp);
	% for i=1:r
	% 	for j=1:c
	% 		xp=X(i)-x0;
	% 		yp=Y(j)-y0;
	% 		sumdF=sumdF+Interp_coef{i,j}(2,1)+2*Interp_coef{i,j}(3,1)*xp+3*Interp_coef{i,j}(4,1)*xp^2+Interp_coef{i,j}(2,2)*yp+Interp_coef{i,j}(2,3)*yp^2+...
	% 		Interp_coef{i,j}(2,4)*yp^3+2*Interp_coef{i,j}(3,2)*xp*yp+2*Interp_coef{i,j}(3,3)*xp*yp^2+2*Interp_coef{i,j}(3,4)*xp*yp^3+3*Interp_coef{i,j}(4,2)*xp^2*yp+...
	% 		3*Interp_coef{i,j}(4,3)*xp^2*yp^2+3*Interp_coef{i,j}(4,4)*xp^2*yp^3;

	% 		sumG2=sumG2+(G)
	% 	end
	% end


	for i=2:max(size(imageNumbers))
		% G=I{i}(subpos(1):subpos(1)+subsize-1,subpos(2):subpos(2)+subsize-1);	%assign subset from deformed image
		% G=im2double(I{i});	%assign subset from deformed image
		G=im2double(I{1}(1:end,6:end));
		% G_interp=I{i}(subpos(1)-1:subpos(1)+subsize,subpos(2)-1:subpos(2)+subsize);	%assign subset from reference image for interpolation
		% [r_int,c_int]=size(G_interp);
		P(:,1)=[-4.5;0;0;0;0;0];

		% for i=2:r_int-1
		% 	for j=2:c_int-1
		% 		Interp_coef{i-1,j-1}=bicubicInterp(G_interp,i,j);
		% 	end
		% end
		[r_g,c_g]=size(G);
		[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
		[r,c]=size(F);
		for converge=1:100
			tic
			for l=1:r
				parfor j=1:c
					dx=X(j)-x0;
					dy=Y(l)-y0;
					xp=x0+dx*(1+P(2,converge))+P(3,converge)*dy+P(1,converge);
					yp=y0+dy*(1+P(6,converge))+P(5,converge)*dx+P(4,converge);
					G_defromed(l,j)=interp2(Xmesh,Ymesh,G,xp,yp,'cubic');
				end
			end
			G_def_mean=mean(mean(G_defromed));
			dG_temp=0;
			for l=1:r
				for j=1:c
					dG_temp=dG_temp+(G_defromed(l,j)-G_def_mean)^2;
					% dG_temp=dG_temp+G_defromed(l,j)-G_def_mean;
				end
			end
			dG=sqrt(dG_temp^2);
			% dG=sqrt(dG_temp^2);
			q=0;
			for l=1:r
				for j=1:c
					q=q+dfdw{l,j}'*(F(l,j)-Fmean-dF/dG*(G_defromed(l,j)-G_def_mean));
				end
			end
			deltaP=-H\q
			W=Warp(P(:,converge),deltaP)
			P(:,converge+1)=[W(1,3);(W(1,1)-1);W(1,2);W(2,3);W(2,1);(W(2,2)-1)]
			toc
		end
	end
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

% function deltaP(F,G,Interp)
% 	[r,c]=size(F);
% 	for i=1:r
% 		for j=1:c
% 			sumdF=sumdF+Interp{i,j}(2,1)+2*Interp{i,j}(3,1)
% 		end
% 	end
% end

function Out=bicubicInterp(I,n,m)
	%get bicubic interpolation coefficients
	A=[1 0 0 0;
	0 0 1 0;
	-3 3 -2 -1;
	2 -2 1 1];
	C=[1 0 -3 2;
	0 0 3 -2;
	0 1 -2 1;
	0 0 -1 1];
	% [Ix,Iy]=imgradient(I);
	% [Ixy,Iyy]=imgradient(Iy);
	[Ix,Iy]=gradient(I);
	[Ixx,Ixy]=gradient(Ix);
	% B=[double(I(n,m)) double(I(n,m+1)) Iy(n,m) Iy(n,m+1);
	% double(I(n+1,m)) double(I(n+1,m+1)) Iy(n+1,m) Iy(n+1,m+1);
	% Ix(n,m) Ix(n,m+1) Ixy(n,m) Ixy(n,m+1);
	% Ix(n+1,m) Ix(n+1,m+1) Ixy(n+1,m) Ixy(n+1,m+1)];
	% B=[I(n,m) I(n,m+1) Iy(n,m) Iy(n,m+1);
	% I(n+1,m) I(n+1,m+1) Iy(n+1,m) Iy(n+1,m+1);
	% Ix(n,m) Ix(n,m+1) Ixy(n,m) Ixy(n,m+1);
	% Ix(n+1,m) Ix(n+1,m+1) Ixy(n+1,m) Ixy(n+1,m+1)];
	B=[I(n,m) I(n+1,m) Iy(n,m) Iy(n+1,m);
	I(n,m+1) I(n+1,m+1) Iy(n,m+1) Iy(n+1,m+1);
	Ix(n,m) Ix(n+1,m) Ixy(n,m) Ixy(n+1,m);
	Ix(n,m+1) Ix(n+1,m+1) Ixy(n,m+1) Ixy(n+1,m+1)];
	Out=A*B*C;
end

% do we treat fbar as a constant in the correlation function
% convolution for the gradient
% if use bicubic interpolation coefficients to find gradient what is x' and y'

%resources
%Advancement of Optical Methods in Experimental Mechanics, Volume 3
%http://www.ncorr.com/index.php/dic-algorithms#3_5

% http://geoserver.ing.puc.cl/info/docencia/ice1603/DIC/Correlation_Tracking_Guide_2010.htm  - DIC tracking program for matlab