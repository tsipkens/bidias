
for kk=1:length(F.n_mat)
    
    kernel.plot_FD_solutions;
    fr(kk) = getframe(gcf);
    drawnow;
    
end

writerObj = VideoWriter('n(m)_APM.avi');
writerObj.FrameRate = 20; % set the seconds per image

% open the video writer
open(writerObj);

% write the frames to the video
writeVideo(writerObj, fr);

% close the writer object
close(writerObj);




