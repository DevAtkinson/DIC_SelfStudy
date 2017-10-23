function DIC
	close all
	clc
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	%define subset size
	subsize=71;
	B=[(1), 0, P1;
		0, (1), P2;
		0 0 1]
	X=[dx;dy;1];


	% B=[(1+P2), P3, P1;
	% 	P5, (P6+1), P4;
	% 	0 0 1]
	% X=[dx;dy;1];
	% B=[(1+P2), P3,P4,0, P1;
	% P6, (P7+1),0,P8, P5;
	% 0 0 1 0 0;
	% 0 0 0 1 0;
	% 0 0 0 0 1]

	% B=[(1), 0,P2,0, P1;
	% 0, (1),0,0, P3;
	% 0 0 1 0 0;
	% 0 0 0 1 0;
	% 0 0 0 0 1]
	% X=[dx;dy;sin(abs(dx)/(subsize/2));sin(abs(dy)/(subsize/2));1];
	symbolic_warp(B,X)

	% find folder containing the images
	current_folder=pwd;
	if ispc
		image_folder=strcat(current_folder,'\Images');
	elseif ismac||isunix
		image_folder=strcat(current_folder,'/Images');
	end

	% define the names of the reference and investigated images
	im='img';
	ims{1}=sprintf('0%d',0);
	ims{2}=sprintf('0%d',2);
	
	%define subset position
	subposin=[1907,1517];
	selecting=0;

	%import images
	% for i=1:max(size(ims))
	% 	images=strcat(strcat(im,num2str(ims{i})),'.tif');
	% 	ImageFileName{i}=fullfile(image_folder,images);
	% 	I{i} = imread(ImageFileName{i});
	% end
	addpath(strcat(current_folder,'\readimxstuff'));
	% addpath('\readimx') 
	I{1}=readimx('G:\Work\Old data\20160527_Sergey_PMMA_Test1_CalGlibble_75mmLens\Test1\B00039.im7');
	I{2}=readimx('G:\Work\Old data\20160527_Sergey_PMMA_Test1_CalGlibble_75mmLens\Test1\B00060.im7');
	% imagesc(im2double(images.Frames{1,1}.Components{1,1}.Planes{1,1}))
	% colormap(gray)
	% imagesc(im2double(images2.Frames{1,1}.Components{1,1}.Planes{1,1}))
	% colormap(gray)
	% save('imm.mat','images')
	% disp(images.Frames)
	% images.Attributes

	% convert images to double
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
	
	% P_def=[-2    -6   5     -4]
	P_def=[-2    9   -4]
	G_deformed=warp_image(F_in,P_def,subposin,subsize);
	figure(1000)
	imagesc(G_deformed)
	colormap(gray)

	% G_in=im2double(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1});
	G_in=G_deformed;
	% save('G_in_final.mat','G_in');
	clear G_in;
	load('G_in_final.mat')
	% F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);

	% this code allows the user to select a seed point in the reference and
	% investigated image that serve as an initial guess for the correlation function
	if selecting==1
		figure(1)
		subplot(2,2,1);
		imagesc(F_in);
		colormap(gray)
		hold on;
		plot(subposin(1),subposin(2),'cx')
		[subx,suby]=ginput(1)
		
		subpos(1)=floor(subx);
		subpos(2)=floor(suby);
		plot(subpos(1),subpos(2),'rx')
		subplot(2,2,2);
		F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);
		imagesc(F);
		colormap(gray)
		hold on;
		[xf,yf]=ginput(1);
		plot(xf,yf,'rx');
		subplot(2,2,3);
		imagesc(G_in);
		colormap(gray)
		key=0;
		fprintf('\nSelect the point in the investigated image where the subset will be located\n');
		fprintf('This function will loop until a keyboard key is pressed at which point you can select the seed point in the subset\n')
		while(key==0)
			[xg1,yg1]=ginput(1);
			subplot(2,2,4);
			imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
			colormap(gray)
			key = waitforbuttonpress;
		end
		subplot(2,2,4);
		imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		colormap(gray)
		[xg,yg]=ginput(1);
		% determine the guess for the displacements
		xguess=(xf+subpos(1))-(xg+floor(xg1)-floor(subsize/2));
		yguess=(yf+subpos(2))-(yg+floor(yg1)-floor(subsize/2));

		% guess=[xguess,0,0,yguess,0,0]
		% guess=[xguess,0,0,0,yguess,0,0,0]
		guess=[xguess,0,yguess]
	else
		subpos=[1911,1511]
		F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);
		% guess=[0.5412,0,0,0.5111,0,0]
		guess=[-0.5474,0.7269]
		% guess=[0.4557,-0.0309]
		% guess=[-5.05,0,0,1,0,0]
	end
	% call the correlation algorithm
	% syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	
	% B=[(1+P2), P3, P1;
	% 	P5, (P6+1), P4;
	% 	0 0 1]
	% X=[dx;dy;1];

	% B=[(1+P2), P3,P4,0, P1;
	% P6, (P7+1),0,P8, P5;
	% 0 0 1 0 0;
	% 0 0 0 1 0;
	% 0 0 0 0 1]
	% X=[dx;dy;sin(dx);sin(dy);1];
	% symbolic_warp(B,X)
	% [P,Corr]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess)

	[P,Corr,G_store]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess)
	ZNCC=1-Corr/2
	save('G_final_zero.mat','G_store');
	% save('F_in_final.mat','F');



	% correct the investigated image to be identical to the reference image 
	% Final_corrected_image=undeform(G_in,P,subsize,subpos);
	% % plot images
	% figure()
	% imagesc(F_in)
	% title('Reference image')
	% figure
	% imagesc(Final_corrected_image)
	% title('investigated image')
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

function G_deformed=warp_image(G,P,subpos,subsize)
	[r_g,c_g]=size(G);
	XX=1:c_g;										% x values over subset
	YY=1:r_g;										% y values over subset
	X=repmat(XX,r_g,1);
	Y=repmat(YY',1,c_g);
	x0=subsize/2+subpos(1);													% x value for subset centre
	y0=subsize/2+subpos(2);
	[r_g,c_g]=size(G);
	[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	dx=X-x0;
	dy=Y-y0;
	% size(X)
	% size(dx)
	dx=reshape(dx,[r_g*c_g,1]);
	dy=reshape(dy,[r_g*c_g,1]);
	[temp]=WarpFuncFinal(dx,dy,P);
	% size(temp)
	xpp=x0+temp(:,1);
	ypp=y0+temp(:,2);
	xp=reshape(xpp,[r_g,c_g]);
	yp=reshape(ypp,[r_g,c_g]);

	% interpolate the investigated subset to obtain the subset used for comparison purposes
	G_deformed=interp2(Xmesh,Ymesh,G,xp,yp,'cubic');
end