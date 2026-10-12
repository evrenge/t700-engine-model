"""T700 fuel control system constants.

Two groups:

* Table C.1 [TM-100991 pdf pp.83-84], 57 rows over two pages, in the report's own
  alphabetical order so this module checks line by line against the printed table.
* Sixteen constants printed only inside block diagrams, in no table. A transcription that
  reads only the tables loses them silently.

Read from page images; the OCR layer gives `CORR 2010` for 20.0 and `T10 0.050` for 0.060.

Table C.1 labels `CORR`, `CR`, `YHILIM` and `YLOLIM` "nondimensional". They are not: `CR`
is ft*lbf and the two limits ft*lbf*sec (open question #11). Use the printed numbers.

Appendix C has no equations, so every constant cites a figure; the traced signal paths
are in `docs/notes/inventory-appendix-c.md`.
"""

from __future__ import annotations

from typing import Final

# =========================================================================== Table C.1
# First printed page [pdf p.83, printed 69] -- 28 rows, AWFP through T17.

AWFP: Final = 0.05909
"""Deceleration schedule slope, lbm*in^2/(lbf*hr*%). [Table C.1]  Fig. C22 gain on PCNGHL."""

B6: Final = 1.0
"""Nonlinear NP-loop-gain enable, nondimensional. [Table C.1]  Fig. C3."""

BWFP: Final = -3.927
"""Deceleration schedule intercept, lbm*in^2/(lbf*hr). NEGATIVE. [Table C.1]  Fig. C22."""

CB: Final = 0.75
"""Load-share speed trim deadband, % NP. [Table C.1]  Fig. C2."""

CE: Final = 13.21
"""Load-share speed trim gain, % NP. [Table C.1]  Fig. C2."""

CH: Final = 0.05
"""NG sensor hysteresis width, % NG. [Table C.1]  One of the four memory elements."""

CLMV: Final = 0.03
"""Metering valve lag time constant, sec. [Table C.1]  Fig. C19."""

CLLDS: Final = 0.2
"""Load demand spindle lag time constant, sec. [Table C.1]  Fig. C12."""

CNTL: Final = 0.025
"""Load demand rate-limit time constant, sec. [Table C.1]  Fig. C12."""

CORR: Final = 20.0
"""Threshold of the B4 relay after the engine-torque integrator. [Table C.1]  Fig. C3.

Labelled "nondimensional" but gates an integral of (TRQL - CR) dt; see the module docstring."""

CR: Final = 20.0
"""Engine torque-level threshold for the nonlinear NP loop gain circuit, ft*lbf.
[Table C.1; the Appendix C nomenclature, pdf p.77, names it]  Fig. C3.

About 79.6 hp at NP_des, a sub-idle "engine is loaded" threshold. Engine shaft torque, not
the rotor-hub torque of Table B.1 past the 81:1 gearbox (open question #11)."""

CT2: Final = 0.088
"""ECU lag time constant, sec. [Table C.1]"""

CT7: Final = 1.0
"""Load-share lag time constant, sec. [Table C.1]  Fig. C2, as 1.1/(CT7*s + 1)."""

CT9: Final = 0.010
"""ZK3 path lag time constant, sec. [Table C.1]"""

CT12: Final = 0.088
"""ECU lag time constant, sec. [Table C.1]  Numerically equal to CT2, printed separately."""

CT13: Final = 1.0
"""ECU lag time constant, sec. [Table C.1]"""

CT14: Final = 0.10
"""ECU lag time constant, sec. [Table C.1]"""

CT16: Final = 1.0
"""ECU lag time constant, sec. [Table C.1]"""

CTPL: Final = 0.010
"""Power turbine speed sensor lag time constant, sec. [Table C.1]  Fig. C1."""

CTPS3: Final = 0.010
"""Compressor discharge static pressure sensor lag time constant, sec. [Table C.1]  Fig. C11."""

DBIAS: Final = 3.25
"""Droop bias, % NP. [Table C.1]"""

KNDRP: Final = 0.25
"""NG droop line slope, lbm*in^2/(lbf*hr*%). [Table C.1]

Used in the main fuel demand [Fig. C9], the load demand schedule [Fig. C17] and the idle
schedule [Fig. C20] (open question #10)."""

NGREF: Final = 101.0
"""NG droop line reference speed, % NG. [Table C.1]  Fig. C9."""

PS3HYS: Final = 0.375
"""Compressor discharge pressure sensor hysteresis, lbf/in^2. [Table C.1]  A memory element."""

T8: Final = 0.40
"""ECU lag time constant, sec. [Table C.1]"""

T10: Final = 0.060
"""ECU lag time constant, sec. [Table C.1]  The OCR layer reads `0.050`; the image says 0.060."""

T11: Final = 0.87
"""ECU lead-lag time constant, sec. [Table C.1]"""

T17: Final = 0.010
"""ECU lag time constant, sec. [Table C.1]"""

# Second printed page [pdf p.84, printed 70] -- 29 rows, T45COR through ZLOLIM.

T45COR: Final = 11.0
"""T4.5 correction offset, deg R. [Table C.1]  Fig. C6."""

T45REF: Final = 2004.0
"""T4.5 limiting reference temperature, deg R. [Table C.1]  Fig. C1.

Competes with the speed error through the maximum error selector."""

TL1: Final = 0.010
"""Torque sensor first lag time constant, sec. [Table C.1]  Fig. C2."""

TL2: Final = 0.010
"""Torque sensor second lag time constant, sec. [Table C.1]  Fig. C2.

TORQ45 -> 1/(TL1*s+1) -> 1/(TL2*s+1) -> TRQL has unity DC gain: no normalisation on that
path (open question #11)."""

TLGE: Final = 0.077
"""ECU lead-lag time constant, sec. [Table C.1]"""

TMDB: Final = 2.0
"""Torque motor deadband, ma. [Table C.1]  Fig. C10."""

TMGN: Final = 0.0159
"""Torque motor gain, in/(ma*sec). [Table C.1]  Fig. C10."""

TMLG: Final = 84.0
"""Torque motor linkage gain, lbm*in^2/(lbf*hr*in). [Table C.1]  Fig. C10."""

TMLVG: Final = 55.0
"""Torque motor voltage gain, volts/in. [Table C.1]  Fig. C10."""

WFMAX: Final = 785.0
"""Maximum metered fuel flow, lbm/hr. [Table C.1]  Fig. C19.

The Fig. 9 open-loop step [pdf p.45] runs 400 -> 775 lbm/hr, just under this ceiling."""

WFMIN: Final = 65.0
"""Minimum metered fuel flow, lbm/hr. [Table C.1]  Fig. C19."""

WFPDCH: Final = 2.10
"""Deceleration schedule upper limit, lbm*in^2/(lbf*hr). [Table C.1]  Fig. C22."""

WFPDCL: Final = 1.45
"""Deceleration schedule lower limit, lbm*in^2/(lbf*hr). [Table C.1]  Fig. C22."""

XHILIM: Final = 0.0787
"""Torque motor integrator rate limit, upper, in/sec. [Table C.1]  Fig. C10."""

XKINTG: Final = 0.18
"""ECU P+I integral gain, nondimensional. [Table C.1]  Fig. C7."""

XKPROP: Final = 0.20
"""ECU P+I proportional gain, nondimensional. [Table C.1]  Fig. C7."""

XLDHYS: Final = 2.5
"""Load demand spindle hysteresis, deg. [Table C.1]  A memory element; see open question #12
on whether the lookup abscissa is XLDSA or XLDSH."""

XLOLIM: Final = -0.0427
"""Torque motor integrator rate limit, lower, in/sec. NEGATIVE. [Table C.1]  Fig. C10."""

YHILIM: Final = 40.0
"""Engine torque integrator upper limit. [Table C.1]  Fig. C3.

Labelled "nondimensional"; it limits an integral of torque, so ft*lbf*sec."""

YLOLIM: Final = 0.0
"""Engine torque integrator lower limit. [Table C.1]  Fig. C3.  Units as YHILIM."""

ZHILIM: Final = 3.5
"""ECU P+I integrator upper limit, nondimensional. [Table C.1]  Figs. C7, C8."""

ZK1: Final = 1.7
"""NP loop additional proportional gain used during high power operation, nondimensional.
[Table C.1; named in the Appendix C nomenclature, pdf p.81]  Gated by B4; see CR."""

ZK3: Final = 0.045
"""ECU gain, nondimensional. [Table C.1]"""

ZK5: Final = 0.40
"""ECU gain, nondimensional. [Table C.1]"""

ZK7: Final = 0.30
"""ECU gain, nondimensional. [Table C.1]"""

ZK8: Final = 0.231
"""ECU gain, nondimensional. [Table C.1]"""

ZK9: Final = 0.625
"""ECU gain, nondimensional. [Table C.1]"""

ZK10: Final = 0.375
"""ECU gain, nondimensional. [Table C.1]"""

ZLOLIM: Final = -1.0
"""ECU P+I integrator lower limit, nondimensional. [Table C.1]  Figs. C7, C8.

Switched: -1.0 when engine 2 torque is below 180 ft*lbf, -0.3 above [Fig. C8]. With the
one-engine switch open TRQL(engine 2) = 0, so -1.0 is the value in force."""

# ================================================== constants printed only on the figures
# Sixteen, in no table (57 + 16 = 73). Nineteen numbers appear only on block diagrams, but
# three repeat Table C.1 values; `tests/test_control_constants.py` reconciles the count.

LOAD_SHARE_LAG_GAIN: Final = 1.1
"""Numerator gain of the load-share lag 1.1/(CT7*s + 1). [Fig. C2, pdf p.85]"""

NP_LOOP_RELAY_LOW: Final = -1.0
"""Lower input threshold of the nonlinear NP-loop-gain relay, % NP. [Fig. C3, pdf p.86]"""

NP_LOOP_RELAY_HIGH: Final = 4.0
"""Upper input threshold of the nonlinear NP-loop-gain relay, % NP. [Fig. C3, pdf p.86]"""

RELAY_OUTPUT: Final = 1.0
"""Output level of the NP-loop relay and of the CORR relay, nondimensional. [Fig. C3, pdf p.86]"""

NP_LOOP_RELAY_GAIN: Final = 1000.0
"""Gain between the NP-loop relay and the B6 gain. [Fig. C3, pdf p.86]"""

TM_INPUT_BIAS: Final = -0.44
"""Bias summed with SPDG at the torque motor input, volts. [Fig. C10, pdf p.89]"""

TM_LEAD: Final = 0.04
"""Torque motor compensation lead time constant, sec. [Fig. C10, pdf p.89]"""

TM_LAG: Final = 0.2
"""Torque motor compensation lag time constant, sec. [Fig. C10, pdf p.89]"""

TM_FORWARD_GAIN: Final = 564.0
"""Torque motor forward gain. [Fig. C10, pdf p.89]"""

TM_CURRENT_BIAS: Final = -31.0
"""Bias summed after the 564.0 forward gain, ma. [Fig. C10, pdf p.89]"""

COLLECTIVE_RIGGING_GAIN: Final = 0.914
"""UH-60A collective pitch to load demand spindle rigging gain, deg/%. [Fig. C13, pdf p.90]

The one place the UH-60A airframe reaches into the engine control: the pilot's input."""

COLLECTIVE_RIGGING_BIAS: Final = 5.34
"""UH-60A collective pitch to load demand spindle rigging bias, deg. [Fig. C13, pdf p.90]"""

LOAD_DEMAND_BIAS: Final = 3.26
"""Bias from which WFQPS3 is subtracted in the load-demand path. [Fig. C17, pdf p.91]"""

FUEL_TRANSPORT_DELAY: Final = 0.015
"""Fuel transport delay, sec, as exp(-0.015 s). [Fig. C19, pdf p.92]

The only pure transport delay in the model, and the fourth memory element."""

ENGINE2_TORQUE_THRESHOLD: Final = 180.0
"""Engine-2 torque threshold in the P+I integrator limit logic, ft*lbf. [Fig. C8, pdf p.88]

About 716 hp at NP_des, between the 80 kt trim (552.6 hp) and hover (911.1 hp)."""

ZLOLIM_HIGH_TORQUE: Final = -0.3
"""ECU P+I integrator lower limit above the 180 ft*lbf threshold. [Fig. C8, pdf p.88]"""

# Two figure constants are unit conversions and live in t700.units: the metering valve's
# 1/3600 lbm/hr -> lbm/sec [Fig. C19, pdf p.92] and the NG sensor's 100/NGDES [Fig. C15,
# pdf p.91], NGDES being Table A.1's NG_DES (open question #9).
