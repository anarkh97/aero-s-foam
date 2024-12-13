AERO-S is a general-purpose finite element structural analyzer developed by the Farhat Research Group (www.stanford.edu/group/frg).
It features capabilities for coupled multi-physics thermoelastic, aeroelastic and aerothermoelastic simulations,
and parallel implementation using domain-decomposition methods. The source code is published online and placed in the "public domain".
Stanford University makes no representations or warranties, expressed or implied, nor assumes any liability for the use of this software.

# Updates incorporated in this repository

## Non-linear crushable foam material model

This fork introduces a crushable foam non-linear material model, which has been validated against `*MAT_154` for uniaxial compression 
loads. The model requires the computation of velocity gradients to determine the deformation tensor. As a result, the following components 
have been updated: `StrainEvaluator`, `GaussIntgElem`, and `MDDynam`. Additionally, the source files for the foam model, `CrushableFoam.h` 
and `CrushableFoam.C`, have been added to the repository.

To use the crushable foam model, the following inputs must be specified in the Aero-S input file:
+-------------+---------------+--------+---+-------+------------+------------+----------+---------+--------------+-----------+-----+------------------+
| MATERIAL_ID | CrushableFoam | $\rho$ | E | $\nu$ | $\sigma_P$ | $\alpha_2$ | $\gamma$ | $\beta$ | $\epsilon_D$ | [$\alpha$ | Tol | $\epsilon_{cr}$] |
+-------------+---------------+--------+---+-------+------------+------------+----------+---------+--------------+-----------+-----+------------------+
|             |               |        |   |       |            |            |          |         |              |           |     |                  |
+-------------+---------------+--------+---+-------+------------+------------+----------+---------+--------------+-----------+-----+------------------+

## Other implementations

This fork also includes the following updates to the base `Aero-S` project:

    1. Fix for the mass computation bug in `FellipaShell`.
    2. Added plastic energy dissipation calculation for `BelytschkoTsayShell`.
    3. Introduced a user-defined parameter for `BlastLoading`, which explicitly defines the metric system for ConWep.
    4. Enabled output of plastic energy dissipation based on a user-defined attribute group from the input file.