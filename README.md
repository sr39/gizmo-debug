Welcome!

This is GIZMO: a flexible, multi-method mhd+gravity+radiation code. The code includes a hydro solver using Lagrangian mesh-free finite-volume Godunov methods (or SPH), and self-gravity solved fast with hybrid PM-Tree methods and fully-adaptive gravitational softenings.

Other physics modules include: magnetic fields (ideal and non-ideal), cosmological integrations, anisotropic conduction and viscosity, sub-grid turbulent diffusion, the ability to insert arbitrary external gravitational fields, radiative cooling (optically thick and thin), integration in non-standard cosmologies, sink particles, "dust-gas mixtures" (with drag, Lorentz forces, back-reaction, collisions), cosmic rays (with advection, diffusion, streaming, heating/cooling, and injection by SNe), nuclear+degenerate equations of state, radiation-hydrodynamics (with several different solvers available), galaxy/star/black hole formation and feedback, self-interacting dark matter, on-the-fly halo finding, and more. Many of these are not in the public code, but in the private (code development) branch; see the User Guide for details.

The code is descended from P-SPH/P-GADGET, itself descended from GADGET-3 (so a huge debt owes to the work of Volker Springel), and many of the GADGET conventions remain (for compatibility with the large library of GADGET analysis software). See the source code for appropriate attribution of the code elements. 

The code is written in standard ANSI C, and should run on all parallel platforms that support MPI. The portability of the code has been confirmed on a large number of systems -- if it can run GADGET (and it almost certainly can), it can run GIZMO.

No, the code title is not an acronym, I just liked it. It refers both to the code's multi-purpose applications and to its historical relationship to GADGET.

The BitBucket site is where I will post code updates, news, and documentation, so check back regularly. If you have code issues, feature requests, bugs, or just questions, please post to the GIZMO Google Group: https://groups.google.com/d/forum/gizmo-code

Basic Rules: 

1. The main reference for the numerical methods, setting up the code, use policies, citations, etc, is the User Guide, available as the file "gizmo_documentation.html" in the "scripts" folder as part of this repository, or via download on the bitbucket site or at my website: 

http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html

Read it! The code is extensively documented. Important things you must know before running are there. Most of the questions people ask are already answered there.

2. Access to the development code does not imply permission to use any modules identified as proprietary either in the User Guide or Template\_Config.sh file or elsewhere. The private (development) version of the code can only be used or distributed with explicit permission from the code authors. Please note that many of the non-public modules are proprietary and developed by active students/postdocs for their ongoing research - it is not acceptable to use or share these routines without first obtaining the explicit permission of both the lead code author and the author(s) of the relevant routines. If in doubt, ask. Anyone who violates these terms will have their access and permissions immediately revoked.

The public version of the code is free software, distributed under the GNU General Public License (http://www.gnu.org/copyleft/gpl.html). This implies that you may freely distribute and copy the (public version) software. You may also modify it as you wish, and distribute these modified versions as long as you indicate prominently any changes you made in the original code, and as long as you leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. Note that the authors retain their copyright on the code. The public code is available at: http://www.tapir.caltech.edu/~phopkins/public/gizmo_public.tgz

If you use any version of the code, please reference the code paper at: http://arxiv.org/abs/1409.7395 (Hopkins 2015); you should also reference Volker Springel's GADGET paper (Springel, 2005, MNRAS, 364, 1105) for the domain decomposition and N-body algorithms. Appropriate citations for specific modules are described in the Users Guide.
