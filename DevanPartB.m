function DevanPartB
	close all
	clc
	% f='C:\Users\devan\Downloads\varsity\DIC selfstudy\assessmentImages';
	f='D:\Work\DIC selfstudy\assessmentImages';
	im='img';
	ims{1}=sprintf('0%d',0);
	ims{2}=sprintf('0%d',2);
	subsize=71;
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
	[xf,yf]=ginput(1)
	plot(xf,yf,'rx');
	subplot(2,2,3);
	imagesc(G_in);
	key=0;
	while(key==0)
		[xg1,yg1]=ginput(1)
		subplot(2,2,4);
		imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		key = waitforbuttonpress;
	end
	subplot(2,2,4);
	imagesc(I{2}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
	[xg,yg]=ginput(1)
	xguess=(xf+subpos(1))-(xg+floor(xg1)-floor(subsize/2));
	yguess=(yf+subpos(2))-(yg+floor(yg1)-floor(subsize/2));
	% guess=[xguess;0;0;yguess;0;0]
	tic
	count=1;
	for i=-50:10:50
		for j=-50:10:50
			fprintf('i:%d j:%d\n',i,j);
			guess=[(xguess+i);0;0;(yguess+j);0;0]
			[P,Cond,figs]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess);
			result{count}.P=P;
			result{count}.guess=guess;
			result{count}.Cond=Cond;
			result{count}.figs=figs;
			result{count}.original_fig=F;
			count=count+1;
			save('stability_results2.mat','result');
		end
	end
	% save('stability_results.mat','results');
	% total_time=toc

	% DevanDICtracking('folder',f,'ImageName',im,'images',ims,'subset size',subsize,'subset position',subpos)

end