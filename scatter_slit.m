function scatter_slit
   
   tmax = 0.005;
   level = 8;
   lambda = 0.005;
   idtype = 1;
   idpar = [0.5 0.2 0.08 0.01 0 10];
   vtype = 2;
   vpar = [0.42 0.46 0.54 0.58 exp(11)];
   % if for single slit ...
   % vpar = [0.48 0.5 0.5 0.52 exp(11)];
   
   [x, y, t, ~, ~, ~, psimod, vol] = ...
sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
   nt = length(t);
   prob = psimod;
   % Initialize and open the video object ...
   v = VideoWriter('slit', 'MPEG-4');
   open(v)

   % For each time step ...
   for it = 1 : nt
      z = squeeze(prob(it, :, :));
      figure(1);
      clf
      contourf(x,y,z, 'edgecolor', 'none');
      pcolor(x,y,z)
      shading interp
      colormap bone
      drawnow
      writeVideo(v, getframe(gcf));
   end
   % Close the video object ...
   close(v);
end
