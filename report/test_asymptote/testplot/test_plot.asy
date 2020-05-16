settings.outformat="pdf";
settings.prc=true;
settings.render=32;

import three;
import fontsize;
import graph;

// font settings and pen settings
TeXHead2.size=new real(pen p=currentpen) {return 0.2mm;};
defaultpen(fontsize(4pt));
defaultpen(0.5);

pen dashed = linetype(new real[] {2,2});
pen thin = black + 0.3pt;
pen redpen  = dashed + 0.3pt + red;
pen bluepen = blue + 0.3pt;

// size settings
real azm_r = 3.5;
real zen = pi/5;
real azm = pi/3;

// point
triple core_point = (azm_r*sin(zen)*cos(azm), azm_r*sin(zen)*sin(azm), azm_r*cos(zen));
triple z_proj = (0., 0., azm_r*cos(zen));
triple xy_proj = (azm_r*sin(zen)*cos(azm), azm_r*sin(zen)*sin(azm), 0.);

// output settings
size(4.00cm, 4.00cm);
currentprojection = orthographic(6,4,4);

// draw azimuthal plane
path3 azm_circle = O -- arc(c=O, v1=azm_r*Z, v2=cos(azm)*X+sin(azm)*Y) -- cycle;
draw(azm_circle, p=thin);

// draw zenith line
path3 zen_line = O -- core_point;
draw(zen_line, p=thin);

// draw the red box
path3 azi_plane = plane(O=O, z_proj, xy_proj);
draw(azi_plane, p=redpen);
// draw(azi_plane, surfacepen=red);

// draw axis system
draw(O--4X, arrow=Arrow3(TeXHead2), L=Label("$x$", position=EndPoint, align=4N));
draw(O--4Y, arrow=Arrow3(TeXHead2), L=Label("$y$", position=EndPoint));
draw(O--4Z, arrow=Arrow3(TeXHead2), L=Label("$z$", position=EndPoint));

// shipout(scale(4.0) * currentpicture.fit());
