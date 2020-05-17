clear all
close all
clc

N = 100;
writer = VideoWriter('simulation.avi');
writer.FrameRate = 10;

i = 0;
fname = 'frame0';

scatter([], [])
set(gca, 'nextplot', 'replacechildren')
axis tight manual 

open(writer);

while isfile(fname)
	fdata = readmatrix(fname);

	scatter(fdata(:,1), fdata(:,2));

	frame = getframe(gcf);
	writeVideo(writer, frame);	
	
	i = i+1;
	fname = 'frame' + string(i);

	if mod(i, 100) == 0
		display('Frame: ' + string(i))
	end
end

close(writer);
