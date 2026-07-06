kill
filename='3typ.out'
mol = G09_structure(filename);
G09_draw_molecule(mol, 'ShowLabels', false);
%print -dpdf 3typ.pdf
sp0  = G09_spectra(filename, 'FWHM', 10, 'plot', true);
%print -dpdf 3typ_spectra.pdf
nm  = G09_nmodes(filename);
G09_charges(filename);
%print -dpdf 3typ_charge.pdf
G09_draw_mode(mol, nm, 44,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf 3tpy_44.pdf
G09_draw_mode(mol, nm, 45,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf 3tpy_45.pdf
G09_draw_mode(mol, nm, 46,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf 3tpy_46.pdf
G09_draw_mode(mol, nm, 58,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf 3tpy_58.pdf
G09_draw_mode(mol, nm, 59,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf 3tpy_58.pdf

%[m,i]=max(sp.Raman);
%G09_draw_mode(mol, nm, 79,'scale',1,'flipsign', true);
save 3typ mol sp0


kill
filename='cd-3typ.out'
mol = G09_structure(filename);
G09_draw_molecule(mol, 'ShowLabels', false);
%print -dpdf cd-3typ.pdf
spcd  = G09_spectra(filename, 'FWHM', 10, 'plot', true);
%print -dpdf cd-3typ_spectra.pdf
nm  = G09_nmodes(filename);
G09_charges(filename);
%print -dpdf cd-3typ_charge.pdf
G09_draw_mode(mol, nm, 47,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf cd-3tpy_47.pdf
G09_draw_mode(mol, nm, 48,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf cd-3tpy_48.pdf
G09_draw_mode(mol, nm, 49,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf cd-3tpy_49.pdf
G09_draw_mode(mol, nm, 61,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf cd-3tpy_61.pdf
G09_draw_mode(mol, nm, 62,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf cd-3tpy_62.pdf

% [m,i]=max(sp.Raman);
% G09_draw_mode(mol, nm, 79,'scale',1,'flipsign', true);
save cd-3typ mol spcd
kill
filename='zn-3typ.out'
mol = G09_structure(filename);
G09_draw_molecule(mol, 'ShowLabels', false);
%print -dpdf zn-3typ.pdf
spzn  = G09_spectra(filename, 'FWHM', 10, 'plot', true);
%print -dpdf zn-3typ_spectra.pdf
nm  = G09_nmodes(filename);
G09_charges(filename);
save zn-3typ mol spzn
print -dpdf zn-3typ_charge.pdf
% G09_draw_mode(mol, nm, 43);
% [m,i]=max(sp.Raman);
G09_draw_mode(mol, nm, 47,'scale',1,'flipsign', true);
print -dpdf zn-3tpy_47.pdf
G09_draw_mode(mol, nm, 48,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf zn-3tpy_48.pdf
G09_draw_mode(mol, nm, 49,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf zn-3tpy_49.pdf
G09_draw_mode(mol, nm, 61,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf zn-3tpy_61.pdf
G09_draw_mode(mol, nm, 62,'flipsign',true);
set(gca,'cameraposition',[-17.4431 -54.1421 49.0335]);
print -dpdf zn-3tpy_62.pdf

