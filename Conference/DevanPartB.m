function DevanPartB
	close all
	clc
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
	%define subset size
	subsize=41;
	%define subset position
	subpos=[425,230];

	%import images
	for i=1:max(size(ims))
		images=strcat(strcat(im,num2str(ims{i})),'.tif');
		ImageFileName{i}=fullfile(image_folder,images);
		I{i} = imread(ImageFileName{i});
	end

	% convert images to double
	F_in=im2double(I{1});
	G_in=im2double(I{2});
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
		imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		key = waitforbuttonpress;
	end
	subplot(2,2,4);
	imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
	[xg,yg]=ginput(1);
	% determine the guess for the displacements
	xguess=(xf+subpos(1))-(xg+floor(xg1)-floor(subsize/2));
	yguess=(yf+subpos(2))-(yg+floor(yg1)-floor(subsize/2));

	% guess=[xguess;0;0;yguess;0;0]
	guess=[xguess;yguess]
	% call the correlation algorithm
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8
	B=[(1+P2), P3,P7,0, P1;
	P5, (P6+1),0,P8, P4]
	X=[dx;dy;sin(dx);sin(dy);1];
	[a,b,c]=symbolic_warp(B,X)
	[P,Corr]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess)

	% [P,Corr]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess,'warp',a,'warp2',b,'warp3',c)
	ZNCC=1-Corr/2




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