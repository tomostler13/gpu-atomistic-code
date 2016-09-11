To calculate M(T) the simulation flags for MvT should be set. The essential
ones are the starting, ending and change in temperatures. The only files
needed are:
-config file
-unit cell file
-exchange file

The damping in the unit cell file should be set to 1.0 (critical damping) so
that the magnetisation relaxes to it's equilibrium quickly.

For the three simple ferromagnetic structures (sc, bcc and fcc) the Curie
temperature should follow the following relations.

sc - 1.44*Jij = kB*Tc
bcc - 2.05*Jij = kB*Tc
fcc - 3.18*Jij = kB*Tc
