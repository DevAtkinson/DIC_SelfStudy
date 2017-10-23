function fulldic
	close all
	clc
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	inc=5;
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	% [FileName,PathName] = uigetfile('*.im7','Select the images','MultiSelect','on');
	% image_count=max(size(FileName));



	% for i=1:2
	% 	% image_folder=strcat(PathName,FileName(i))
	% 	image_folder = fullfile( PathName , FileName{i} );
	% 	I{i}=readimx(image_folder);
	% end

	% image_folder = fullfile( PathName , FileName{1} );
	% I{1}=readimx(image_folder);
	% image_folder = fullfile( PathName , FileName{1+inc} );
	% I{2}=readimx(image_folder);

	%define subset size
	subsize=61;
	stepsize=29;
	% B=[(1), 0, P1;
	% 	0, (1), P2;
	% 	0 0 1]
	% X=[dx;dy;1];

	B=[(1+P2), P3, P1;
		P5, (P6+1), P4;
		0 0 1]
	X=[dx;dy;1];
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



	% symbolic_warp(B,X)



	% clear P;
	%define subset position
	% subposin=[1907,1517];
	selecting=0;

	
	% % addpath('\readimx') 
	% I{1}=readimx('G:\Work\Old data\20160527_Sergey_PMMA_Test1_CalGlibble_75mmLens\Test1\B00039.im7');
	% I{2}=readimx('G:\Work\Old data\20160527_Sergey_PMMA_Test1_CalGlibble_75mmLens\Test1\B00060.im7');
	% imagesc(im2double(images.Frames{1,1}.Components{1,1}.Planes{1,1}))
	% colormap(gray)
	% imagesc(im2double(images2.Frames{1,1}.Components{1,1}.Planes{1,1}))
	% colormap(gray)
	% save('imm.mat','images')
	% disp(images.Frames)
	% images.Attributes

	% convert images to double

	% F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});


	
	% imagesc(F_in)
	% polygon=imrect();
	% [subpos,horizon_count,vert_cont,xtick,ytick]=Mask2Subsets(polygon,subsize,stepsize);



	% mask=createMask(polygon);
	% figure
	% imagesc(mask)
	% BB=bwboundaries(mask)
	% figure
	% imagesc(F_in)
	% hold on
	% visboundaries(BB)
	% save('BB.mat','BB')
	[r_F,c_F]=size(F_in);
	
	% P_def=[-2    -6   5     -4]
	% P_def=[-2    9   -4]
	% G_deformed=warp_image(F_in,P_def,subposin,subsize);
	% figure(1000)
	% imagesc(G_deformed)
	% colormap(gray)

	% G_in=im2double(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1});


	% G_in=im2double(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1});


	% save('G_in_final.mat','G_in');

	% F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);

	% this code allows the user to select a seed point in the reference and
	% investigated image that serve as an initial guess for the correlation function
	if (selecting==1)
		figure(1)
		subplot(2,2,1);
		imagesc(F_in);
		colormap(gray)
		hold on;
		for i=1:max(size(xtick))
			plot([xtick(i) xtick(i)],[ytick(1) ytick(end)],'c')
		end
		for i=1:max(size(ytick))
			plot([xtick(1) xtick(end)],[ytick(i) ytick(i)],'c')
		end
		key=0;
		while(key==0)
			[subx,suby]=ginput(1);
			out_temp=whichSubpos(subpos,stepsize,subx,suby);
			subx=subpos{out_temp(1),out_temp(2)}(1)
			suby=subpos{out_temp(1),out_temp(2)}(2)
			subplot(2,2,2);
			imagesc(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)));
			colormap(gray)
			key = waitforbuttonpress;
		end
		subplot(2,2,1);
		plot(subx,suby,'rx')
		subplot(2,2,2);
		imagesc(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)));
		colormap(gray)
		hold on
		[xf,yf]=ginput(1)
		plot(xf,yf,'rx');
		subplot(2,2,3);
		imagesc(G_in);
		colormap(gray)
		hold on
		key=0;
		while(key==0)
			[xg1,yg1]=ginput(1)
			subplot(2,2,4);
			imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
			colormap(gray)
			key = waitforbuttonpress;
		end
		subplot(2,2,3);
		plot(xg1,yg1,'rx')
		subplot(2,2,4);
		imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		colormap(gray)
		[xg,yg]=ginput(1)
		hold on
		plot(xg,yg,'cx')
		xguess=((xf+subx)-(xg+floor(xg1)))
		yguess=((yf+suby)-(yg+floor(yg1)))
		figure
		imagesc(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)))
		figure
		imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby+yguess)-floor(subsize/2):floor(suby+yguess)+floor(subsize/2),floor(subx+xguess)-floor(subsize/2):floor(subx+xguess)+floor(subsize/2)))

		figure
		meshcompare(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)),I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby+yguess)-floor(subsize/2):floor(suby+yguess)+floor(subsize/2),floor(subx+xguess)-floor(subsize/2):floor(subx+xguess)+floor(subsize/2)))
		guess=[xguess,0,0,yguess,0,0]
	elseif (selecting==0)
		% subpos=[1911,1511]
		% F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);
		% guess=[0.5412,0,0,0.5111,0,0]
		guess=[0.0949,0,0,0,0,0]
		subx=188;
		suby=654;
		xf=47.8149;
		yf=29.9574

		% guess=[0.4557,-0.0309]
		% guess=[-5.05,0,0,1,0,0]
	end
	% out=whichSubpos(subpos,stepsize,subx+xf,suby+yf)
	% subpos{out(1),out(2)}
	% [PP,CCorr]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{out(1),out(2)},'guess',guess)
	% % ZNCC=1-Corr/2
	% result{out(1),out(2)}.P=PP;
	% result{out(1),out(2)}.Corr=CCorr;
	% [process,dismax,dissing]=correlation_order(subpos,out);
	% dissing
	% for i=1:dismax
	% 	for j=1:dissing(i)
	% 		P_in(j,:)=result{process{i}(j,3),process{i}(j,4)}.P;
	% 		subposindex(j,1)=subpos{process{i}(j,1),process{i}(j,2)}(1);
	% 		subposindex(j,2)=subpos{process{i}(j,1),process{i}(j,2)}(2);
	% 		endcondition=dissing(i);
	% 	end
	% 	parfor j=1:endcondition
	% 		% [P{process{i}(j,1),process{i}(j,2)},Corr{process{i}(j,1),process{i}(j,2)}]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process{i}(j,1),process{i}(j,2)},'guess',P{process{i}(j,3),process{i}(j,4)});
	% 	% [PP{j},CCorr{j}]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',[subposindex(j,1),subposindex(j,2)],'guess',P_in(j,:));
	% 	[PP{j},CCorr{j}]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',[6,6],'guess',guess);

	% 	end
	% 	for j=1:dissing(i)
	% 		result{process{i}(j,1),process{i}(j,2)}.P=PP{j};
	% 		result{process{i}(j,1),process{i}(j,2)}.Corr=CCor{j};
	% 	end
	% end
	% save('P_save.mat','P');
	% save('Corr_save.mat','Corr');





	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	image_folder = fullfile( PathName , FileName{1+inc} );
	I{2}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
	subpos=Proc.subpos;
	out=Proc.starting_subset;
	% Proc.starting_pos=[(subx+xf),(suby+yf)];
	% Proc.polygon=polygon;
	% Proc.stepsize=stepsize;
	% Proc.subsize=subsize;
	% Proc.Warp=B;
	% Proc.WarpVec=X;
	FileName=Proc.FileName;
	PathName=Proc.PathName;
	inc=Proc.inc;
	current_image=Proc.correlated_to;
	image_count=max(size(FileName));

	out=whichSubpos(subpos,stepsize,subx+xf,suby+yf)
	subpos{out(1),out(2)}
	Proc.subpos=subpos;
	Proc.starting_subset=out;
	Proc.starting_pos=[(subx+xf),(suby+yf)];
	Proc.polygon=polygon;
	Proc.stepsize=stepsize;
	Proc.subsize=subsize;
	Proc.Warp=B;
	Proc.WarpVec=X;
	Proc.FileName=FileName;
	Proc.PathName=PathName;
	Proc.inc=inc;
	[process,dismax,dissing,D,Index]=correlation_order(subpos,out);
	Proc.process=process;
	Proc.dismax=dismax;
	proc.dissing=dissing;
	Proc.Index=Index;
	counter_image=1;
	for k=(current_image+inc):inc:image_count
		fprintf('image %d',k);
		
		if (k==(1+inc))
			G_in=im2double(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1});
			[PP,CCorr]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{out(1),out(2)},'guess',guess)
		else
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
			[PP,CCorr]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{out(1),out(2)},'guess',Proc.im{k-1-inc}.D(Index(out(1),out(2)),6:11))
		end
		% ZNCC=1-Corr/2
		% P{out(1),out(2)}=PP;
		% Corr{out(1),out(2)}=CCorr;
		D(1,6:11)=PP;
		D(1,12)=CCorr;
		D(1,1)=1;
		D(1,2)=out(1);
		D(1,3)=out(2);
		D(1,4:5)=[0 0];


		tic
		for i=1:dismax
			% for j=1:dissing(i)
			% 	% process{i}(j,3)
			% 	% process{i}(j,4)
			% 	% [rr,cc]=size(P{process{i}(j,3),process{i}(j,4)})
				% P_in(j,:)=P{process{i}(j,3),process{i}(j,4)};
				% process{i}(j,3)
				% process{i}(j,4)
				% P_in(j,:)=D(Index(process{i}(j,3),process{i}(j,4)),6:11);
			% 	subposindex(j,1)=subpos{process{i}(j,1),process{i}(j,2)}(1);
			% 	subposindex(j,2)=subpos{process{i}(j,1),process{i}(j,2)}(2);
				
			% end
			% if i~=1
			% 	prevous_count=sum(dissing(1:i-1));
			% else
			% 	prevous_count=0;
			% end
			endcondition=dissing(i);
			parfor j=1:endcondition
			[PP(j,:),CCorr(j)]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process{i}(j,1),process{i}(j,2)},'guess',D(Index(process{i}(j,3),process{i}(j,4)),6:11));
				% [P{process{i}(j,1),process{i}(j,2)},Corr{process{i}(j,1),process{i}(j,2)}]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process{i}(j,1),process{i}(j,2)},'guess',P{process{i}(j,3),process{i}(j,4)});
				% [PP(j,:),CCorr(j)]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',[subposindex(j,1),subposindex(j,2)],'guess',P_in(j,:));
			% [PP{j},CCorr{j}]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',[6,6],'guess',guess);
			end
			for j=1:dissing(i)
				% P{process{i}(j,1),process{i}(j,2)}=PP(j,:);
				% Corr{process{i}(j,1),process{i}(j,2)}=CCorr(j);
				D(Index(process{i}(j,1),process{i}(j,2)),6:11)=PP(j,:);
				D(Index(process{i}(j,1),process{i}(j,2)),12)=CCorr(j);
			end
		end
		Proc.im{k-1}.D=D;
		Proc.correlated_to=k;
		save('Richard_CTL_05_2.mat','Proc');
		% result{k-1}.P=P;
		% result{k-1}.Corr=Corr;
		
		% save('Corr_save.mat','Corr');
		
	end
	toc
	% save('Result2.mat','result');

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

function [subpos,horizontal,vertical,xtick,ytick]=Mask2Subsets(rect,subsize,stepsize)
	coords=getPosition(rect)
	xvals=[coords(1),(coords(1)+coords(3))];
	yvals=[coords(2),(coords(2)+coords(4))];
	xmin=ceil(min(xvals));
	ymin=ceil(min(yvals));
	xmax=floor(max(xvals));
	ymax=floor(max(yvals));
	wid=xmax-xmin;
	height=ymax-ymin;
	vertical=floor((height-subsize)/stepsize);
	horizontal=floor((wid-subsize)/stepsize);
	for i=1:vertical
		for j=1:horizontal
			subpos{i,j}(1)=xmin+(j-1)*stepsize;
			subpos{i,j}(2)=ymin+(i-1)*stepsize;

		end
	end
	for i=1:vertical+1
		ytick(i)=ymin+(i-1)*stepsize;
	end
	for i=1:horizontal+1
		xtick(i)=xmin+(i-1)*stepsize;
	end

end

function out=whichSubpos(subpos,stepsize,x,y)
	[r,c]=size(subpos);
	for i=1:r
		for j=1:c
			if (x>=subpos{i,j}(1))&(x<subpos{i,j}(1)+stepsize)&(y>=subpos{i,j}(2))&(y<subpos{i,j}(2)+stepsize)
				out=[i,j];
			end
		end
	end
end

function [processing,dismax,dissing,processout,Mat]=correlation_order(subpos,selected)
	[r,c]=size(subpos);
	% count=zeros(max([r,c]));
	count=1;
	for i=1:r
		for j=1:c
			dis=max([abs(selected(2)-i), abs(selected(1)-j)]);
			% process{dis}(count(dis),1:2)=[i,j,0,0];
			process1(count,:)=[dis,i,j,0,0];
			% count(dis)=count(dis)+1;
			A(i,j)=dis;
			count=count+1;
		end
	end
	count=1;
	for i=1:r
		for j=1:c
			dis=max([abs(selected(1)-i), abs(selected(2)-j)]);
			if (i~=1)&(A(i-1,j)<A(i,j))
				% fprintf('%f %f',dis,count(dis))
				process(count,:)=[dis,i,j,(i-1),j];
				% process{dis}(count(dis))=[i,j,(i-1),j];
			elseif (j~=1)&(A(i,j-1)<A(i,j))
				% fprintf('%f %f',dis,count(dis))
				process(count,:)=[dis,i,j,i,(j-1)];
				% process{dis}(count(dis))=[i,j,i,(j-1)];
			elseif (i~=r)&(A(i+1,j)<A(i,j))
				% fprintf('%f %f',dis,count(dis))
				process(count,:)=[dis,i,j,(i+1),j];
				% process{dis}(count(dis))=[i,j,(i+1),j];
			elseif (j~=c)&(A(i,j+1)<A(i,j))
				% fprintf('%f %f',dis,count(dis))
				process(count,:)=[dis,i,j,i,(j+1)];
				% process{dis}(count(dis))=[i,j,i,(j+1)];
			elseif (j~=c)&(i~=r)&(A(i+1,j+1)<A(i,j))
				process(count,:)=[dis,i,j,(i+1),(j+1)];
			elseif (j~=c)&(i~=1)&(A(i-1,j+1)<A(i,j))
				process(count,:)=[dis,i,j,(i-1),(j+1)];
			elseif (j~=1)&(i~=r)&(A(i+1,j-1)<A(i,j))
				process(count,:)=[dis,i,j,(i+1),(j-1)];
			elseif (j~=1)&(i~=1)&(A(i-1,j-1)<A(i,j))
				process(count,:)=[dis,i,j,(i-1),(j-1)];
			end
			count=count+1;
		end
	end
	processout1=sortrows(process)

	
	dismax=max(processout1(:,1))
	for j=1:dismax
		count=1;
		for i=1:max(size(processout1))
			if processout1(i,1)==j
				processing{j}(count,:)=processout1(i,2:5);
				count=count+1;
			end
		end
		dissing(j)=count-1;
	end
	processout(1,:)=[processout1(1,1), processout1(1,2), processout1(1,3), processout1(1,4), processout1(1,5), 0, 0, 0, 0, 0, 0, 0];
	Mat(selected(1),selected(2))=1;
	for i=2:max(size(processout1))
		Mat(processout1(i,2),processout1(i,3))=i;
		processout(i,:)=[processout1(i,1), processout1(i,2), processout1(i,3), processout1(i,4), processout1(i,5), 0, 0, 0, 0, 0, 0, 0];
	end
end