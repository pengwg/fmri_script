function  [Number_of_Frame_To_Remove]= check_motionDisplacement_NHM_v1(path_to_func,threshold,session,subjects)
%This function add a column of framwise displacement (using bramila
%function (ref % Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also
% Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048) and determine a
% threshold in which the motion was too high , in that case the frame will
% be labeld as 1, otherwise 0

%input are: path of the subjects , name of the output text file and
%treshold


v = importdata([path_to_func ls([path_to_func 'rp_c*'])]);
cfg.motionparam = [path_to_func ls([path_to_func 'rp_c*'])];
cfg.prepro_suite =  'spm' ;
cfg.radius = 50 ; %default value

[fwd,rms]=bramila_framewiseDisplacement(cfg)
Frame2remove = zeros(length(fwd),1);
Frame2remove(find(fwd > threshold)) = 1;

Number_of_Frame_To_Remove = length(find(fwd > threshold));
cd ('C:\Users\natha\Downloads\FUS_transferToNIDA(1)\FUS_transferToNIDA\')
a = figure(1)
plot(v(:,1))
hold on
plot(v(:,2))
hold on
plot(v(:,3))
ylim([-2,2])
saveas(a,['Figure_' char(subjects) '_' char(session) '.png']);
close(a)