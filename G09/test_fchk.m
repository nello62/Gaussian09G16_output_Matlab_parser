kill
% ── Carica tutto dal .fchk ────────────────────────────────────────────────
data = G09_fchk_read('3typ.fchk');
mol=G09_structure('3typ.out');
nm=G09_nmodes('3typ.out');
% ── Molecola ──────────────────────────────────────────────────────────────
G09_draw_molecule(data.mol)

% ── Cariche Mulliken con visualizzazione ──────────────────────────────────
ch = G09_charges_fchk(data.mol, data.ch);

% ── Con opzioni ───────────────────────────────────────────────────────────
% ch = G09_charges_fchk(data.mol, data.ch, 'threshold', 0.05);
% ch = G09_charges_fchk(data.mol, data.ch, 'plot', false);
% ch = G09_charges_fchk(data.mol, data.ch, 'ColorScale', 'none');

% ── Modi normali ──────────────────────────────────────────────────────────
G09_draw_mode(data.mol, data.nm, 85)
G09_draw_mode(mol, nm, 85) % modo Raman più intenso