# two-body-orbit-propagator
A MATLAB tool that computes classical orbital elements and propagates position and velocity vectors between ECI and perifocal frames using two-body Keplerian dynamics.

This MATLAB script calculates and propagates an object's orbital state based on two-body Keplerian motion. Given position and velocity vectors in the Earth-Centered Inertial (ECI) frame, it computes classical orbital elements, transforms states into the perifocal frame, and propagates the orbit to a future time.

## Features
- Computes:
  - Classical Orbital Elements (COEs)
  - Direction Cosine Matrix (DCM) from ECI to Perifocal
  - True and Eccentric Anomalies at initial and final time
  - Time of Periapsis Passage
- Converts position and velocity vectors between:
  - ECI frame
  - Perifocal frame
- Propagates orbital state vectors from `t1` to `t2`
- All outputs are labeled and include physical units

## Inputs
- Initial position and velocity vectors in ECI coordinates
- Initial and future time (`t1` and `t2`)
- Absolute and relative tolerances for iterative procedures

## Outputs
- Classical orbital elements (a, e, i, Ω, ω, ν)
- DCM converting ECI → Perifocal at `t1`
- Position and velocity vectors in:
  - ECI at `t1` and `t2`
  - Perifocal frame at `t1` and `t2`
- True and Eccentric Anomalies at `t1` and `t2`
- Time of Periapsis Passage (Tₚ)

## File Structure
- `main_script.m` — Main execution file
- `functions/` — Folder containing helper functions (if applicable)

## Related Topics
Orbital Mechanics • Keplerian Motion • ECI Frames • DCMs • Aerospace Engineering
