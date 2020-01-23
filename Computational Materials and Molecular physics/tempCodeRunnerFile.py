h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.7]])
h2.calc = NWChem(xc='PBE')
opt = BFGS(h2)
opt.run(fmax=0.02)