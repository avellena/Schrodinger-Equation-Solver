   tmax = 0.01;
   level = 8;
   lambda = 0.02;
   idtype = 1;
   idpar = [0.5 0.2 0.05 0.05 0 8];
   vtype = 1;
   vpar = [0.3 0.6 0.3 0.6 -exp(10)];
   
   [x, y, t, ~, ~, ~, psimod, v] = ...
sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
   nt = length(t);
   prob = psimod;
   % Initialize and open the video object ...
   v = VideoWriter('well', 'MPEG-4');
   open(v)

   % For each time step ...
   for it = 1 : nt
      figure(1);
      clf
      contourf(x,y,squeeze(prob(it, :, :)), 'edgecolor', 'none');
      colormap bone;
      drawnow
      title(sprintf('Time step %d of %d: t = %.3g', it, nt, t(it)));
         % Record the frame ...
      writeVideo(v, getframe(gcf));
   end
   % Close the video object ...
   close(v);

