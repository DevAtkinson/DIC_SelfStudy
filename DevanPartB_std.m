function DevanPartB_std
	close all
	clc
	% f='C:\Users\devan\Downloads\varsity\DIC selfstudy\assessmentImages';
	f='D:\Work\DIC selfstudy\assessmentImages';
	im='img';
	ims{1}=sprintf('0%d',0);
	ims{2}=sprintf('0%d',2);
	subsize=41;
	subpos=[425,230];

	for i=1:max(size(ims))
		images=strcat(strcat(im,num2str(ims{i})),'.tif');
		ImageFileName{i}=fullfile(f,images);
		I{i} = imread(ImageFileName{i});
	end

	F_in=im2double(I{1});
	G_in=im2double(I{2});
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
		imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		key = waitforbuttonpress;
	end
	subplot(2,2,4);
	imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
	[xg,yg]=ginput(1);
	xguess=(xf+subpos(1))-(xg+floor(xg1)-floor(subsize/2));
	yguess=(yf+subpos(2))-(yg+floor(yg1)-floor(subsize/2));

	guess=[xguess;0;0;yguess;0;0]
	parfor i=1:1000
		[P{i},Corr{i}]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess);
		Corr2{i}=(1-Corr{i}/2);
	end
	result.P=P;
	result.Corr=Corr;
	result.Corr2=Corr2;
	save('std_results_partb','result');
	Final_corrected_image=undeform(G_in,P,subsize,subpos);
	figure()
	imagesc(F_in)
	figure
	imagesc(Final_corrected_image)
	% DevanDICtracking('folder',f,'ImageName',im,'images',ims,'subset size',subsize,'subset position',subpos)

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