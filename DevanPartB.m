function DevanPartB
	close all
	clc
	% f='C:\Users\devan\Downloads\varsity\DIC selfstudy\assessmentImages';
	f='D:\Work\DIC selfstudy\assessmentImages';
	im='img';
	ims{1}=sprintf('0%d',0);
	ims{2}=sprintf('0%d',2);
	subsize=75;
	subpos=[425,230];
	DevanDICtracking('folder',f,'ImageName',im,'images',ims,'subset size',subsize,'subset position',subpos)
end