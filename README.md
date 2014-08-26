Welcome!

This is GIZMO (beta version: likely to be Google-style and stay in beta for quite some time).

The simulation code GIZMO is a flexible, multi-method mhd+gravity code. The code includes a hydro solver using the Lagrangian meshless finite-mass method, a meshless finite volume method, the modern "pressure-sph" method, and "traditional" SPH. Self-gravity is solved fast, with a BH-Tree (optionally a hybrid PM-Tree for periodic boundaries), with adaptive gravitational softenings. Hydrodynamics and gravity are both optional. 

The code is descended from P-SPH/P-GADGET, itself descended from GADGET-3 (so a huge debt owes to the work of Volker Springel), and many of the naming conventions and routines remain (for the sake of compatibility with the large library of GADGET work and analysis software). You should see the source code for appropriate attribution of the code elements. Currently available modules include things like: hydrodynamics, MHD, cosmological integrations, galaxy/star/black hole formation with feedback from stars and black holes (both explicit, detailed models and sub-grid models), self-interacting dark matter, adaptive gravitational softening lengths for all particle types, conduction, sub-grid turbulent diffusion, the ability to insert arbitrary external gravitational fields, integration in non-standard cosmologies, sink particles, "dust fluids" (particulate-gas interactions), cosmic rays (in progress, partially implemented), nuclear+degenerate equations of state (in progress, partially implemented), and radiation hydrodynamics (in progress, partially implemented). Most of these are not in the public release of the code, but in the private (code development) branch; see the users guide for details.

No, the code title is not an acronym, I just liked it. It refers both to the code's multi-purpose applications and to its historical relation ship to GADGET.

The bitbucket site is where I will post code updates, news, and documentation, so check back regularly. The main reference for the numerical methods, setting up the code, code policies, branching etc, is the user's guide, available through download on the bitbucket site or at my website: 

http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html

Read it!

The code is written in standard ANSI C, and should run on all parallel platforms that support MPI. The portability of the code has been confirmed on a large number of systems -- if it can run GADGET (and it almost certainly can), it can run GIZMO.

