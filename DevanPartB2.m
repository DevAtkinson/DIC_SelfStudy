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
	multiplier=5;
	n=21;
	adapter=ceil(n/2);
	for i=1:n
		parfor j=1:n
			fprintf('i:%d j:%d\n',i,j);
			guess{i,j}=[(xguess+(i-adapter)*multiplier);0;0;(yguess+(j-adapter)*multiplier);0;0];
			[P{i,j},Cond{i,j},figs{i,j}]=DevanDICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos,'guess',guess{i,j});
			% result{i,j}.P=P;
			% result{i,j}.guess=guess;
			% result{i,j}.Cond=Cond;
			% result{i,j}.figs=figs;
			% result{i,j}.original_fig=F;
			% count=count+1;
			
		end
	end
	result=getback(P,Cond,figs,guess,F,n);
	tic
	save('stability_results_subset_71_2.mat','result','-v7.3');
	toc
	% save('stability_results.mat','results');
	% total_time=toc

	% DevanDICtracking('folder',f,'ImageName',im,'images',ims,'subset size',subsize,'subset position',subpos)

end

function result=getback(P,Cond,figs,guess,F,n)
	for i=1:n
		for j=1:n
			result{i,j}.P=P{i,j};
			result{i,j}.Cond=Cond{i,j};
			result{i,j}.figs=figs{i,j};
			result{i,j}.guess=guess{i,j};
			result{i,j}.original_fig=F;
		end
	end
end