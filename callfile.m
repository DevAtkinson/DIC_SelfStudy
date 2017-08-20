function callfile()
	close all
	clc
	f='C:\Users\devan\Downloads\varsity\DIC selfstudy\assessmentImages';
	im='calimg';
	for i=1:9
		ims{i}=sprintf('0%d',i);
	end
	DevanCalibration('folder',f,'ImageName',im,'images',ims)


end