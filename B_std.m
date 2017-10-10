load('std_results_partb');
for i=1:1000
	P1(i)=result.P{i}(1);
	P2(i)=result.P{i}(2);
	P3(i)=result.P{i}(3);
	P4(i)=result.P{i}(4);
	P5(i)=result.P{i}(5);
	P6(i)=result.P{i}(6);
	% C(i)=result.Corr{i}(1);
	% C2(i)=result.Corr2{i}(1);
end

P=[std(P1);std(P2);std(P3);std(P4);std(P5);std(P6)]
% Corr=std(C);
% Corr2=std(C2);