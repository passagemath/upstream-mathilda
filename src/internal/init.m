(* D, Dt and Derivative are now implemented in C (src/deriv.c) and
   registered during core_init(). The rule-based bootstrap that used
   to live in src/internal/deriv.m is no longer loaded. *)

(* The CRC integral table (CRCMathTablesIntegrals.m) is loaded lazily
   the first time Integrate's CRCTable stage is invoked — see
   src/integrate.c's try_crctable().  We do not Get it here so startup
   cost stays low for users who never call Integrate. *)

(* The Bessel J/Y/I/K symbolic rules (half-integer order -> elementary,
   negative-integer-order reflection) that used to live in bessel.m are in
   the C builtins (src/special_functions/bessel.c): as DownValues they were
   tried before the builtin on every numeric call. *)

(* HeavisideTheta / DiracDelta value rules on numeric arguments, used by
   DSolve's Green's-function / impulse solutions (see distributions.m). *)
LoadModule["distributions.m"];
