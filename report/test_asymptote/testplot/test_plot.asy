settings.outformat="html";
settings.prc=false;
settings.render=32;

import three;
import fontsize;
import graph3;

// font settings and pen settings
usepackage("fourier");
defaultpen(font("T1","fut\textfamilyextension","m","n"));

defaultpen(fontsize(4pt));
defaultpen(0.5);

pen dashed = linetype(new real[] {5,3});
pen thin = black + 0.3pt;
pen redpen  = dashed + 0.3pt + red;
pen bluepen = blue + 0.3pt + fontsize(3pt);
pen blackpen = black + 0.3pt + fontsize(3pt);
material redsurface = material(red+opacity(0.3));

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
path3 azm_circle = arc(c=O, v1=azm_r*Z, v2=cos(azm)*X+sin(azm)*Y) -- xy_proj;
draw(azm_circle, p=thin);

// draw zenith line
path3 zen_line = O -- core_point;
draw(zen_line, p=thin);

// draw the red box
path3 azi_plane = plane(O=O, z_proj, xy_proj);
draw(azi_plane, p=redpen);
draw(surface(azi_plane), surfacepen=redsurface, light=nolight);

// draw unit vector at core vector
TeXHead2.size=new real(pen p=currentpen) {return 0.1mm;};

triple psi = unit((sin(zen)*cos(azm), sin(zen)*sin(azm), cos(zen)));
triple phi = unit((-sin(azm), cos(azm), 0.));
triple theta = unit((cos(zen)*cos(azm), cos(zen)*sin(azm), -sin(zen)));

draw(shift(core_point) * scale3(1.0) * (O -- psi), p=bluepen, 
L=Label("$\hat{\psi}_s$", position=EndPoint), 
arrow=Arrow3(TeXHead2(normal=phi), emissive(blue))
);

draw(shift(core_point) * scale3(0.4) * (O -- phi), p=bluepen, 
L=Label("$\hat{\phi}_s = \hat{h}_s$", position=EndPoint), 
arrow=Arrow3(TeXHead2(normal=psi), emissive(blue))
);

draw(shift(core_point) * scale3(0.4) * (O -- theta), p=bluepen, 
L=Label("$\hat{\theta}_s = \hat{v}_s$", position=EndPoint), 
arrow=Arrow3(TeXHead2(normal=psi), emissive(blue))
);

draw(shift(xy_proj) * scale3(0.4) * (O -- phi), p=blackpen, 
L=Label("$\hat{\phi}_s = \hat{h}_s$", position=EndPoint), 
arrow=Arrow3(TeXHead2(normal=Z), emissive(black))
);

// draw the angles
real azm_angle_r = 0.8;
real zen_angle_r = 1.7;

path3 azm_ang = arc(c=O, v1=azm_angle_r*X, v2=xy_proj);
path3 zen_ang = arc(c=O, v1=zen_angle_r*Z, v2=core_point);

draw(zen_ang, p=blackpen,
L=Label("$\theta_s$", position=MidPoint),
arrow=Arrow3(TeXHead2(normal=phi), emissive(black))
);

draw(azm_ang, p=blackpen,
L=Label("$\phi_s$", position=MidPoint),
arrow=Arrow3(TeXHead2(normal=Z), emissive(black))
);


// draw axis system
TeXHead2.size=new real(pen p=currentpen) {return 0.2mm;};
draw(O--4X, arrow=Arrow3(TeXHead2), L=Label("$x$", position=EndPoint, align=4N));
draw(O--4Y, arrow=Arrow3(TeXHead2), L=Label("$y$", position=EndPoint));
draw(O--4Z, arrow=Arrow3(TeXHead2), L=Label("$z$", position=EndPoint));

// shipout(scale(4.0) * currentpicture.fit());
