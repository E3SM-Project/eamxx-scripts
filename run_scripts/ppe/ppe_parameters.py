#!/usr/bin/env python3
"""
Registry of every parameter the PPE can vary. Shared by make_json_files_multires_append.py
(which samples from it) and create_ensembles.py (which validates JSON files against it).

To make a parameter available to a PPE, add it here. Which parameters an ensemble actually
varies is decided by the JSON file (its "parameters" list), not by this file.

Each entry: namelist key (namespace::name), sampling range, scale, kind, code default.
  scale   "linear" or "log" (log needs minimum > 0)
  kind    "real" or "int" (int: equal-probability levels minimum..maximum)
  default code default at the ZM checkout (d5090c84), for reference and metadata only
"""
from dataclasses import dataclass
from typing import Literal, Optional


@dataclass(frozen=True)
class Parameter:
    key: str
    minimum: float
    maximum: float
    scale: Literal["linear", "log"] = "linear"
    kind: Literal["real", "int"] = "real"
    default: Optional[float] = None
    legacy_env: Optional[str] = None   # env var used by the pre-registry run scripts

    @property
    def scheme(self) -> str:
        return self.key.split("::", 1)[0]

    @property
    def name(self) -> str:
        return self.key.split("::", 1)[1]


PARAMETERS = [
    # SHOC
    Parameter("shoc::thl2tune",        0.1,    10.0,  "linear", default=1.0,   legacy_env="SHOC_THL2TUNE"),
    Parameter("shoc::qw2tune",         0.1,    10.0,  "linear", default=1.0,   legacy_env="SHOC_QW2TUNE"),
    Parameter("shoc::length_fac",      0.1,    1.15,  "linear", default=0.5,   legacy_env="SHOC_LENGTH_FAC"),
    Parameter("shoc::c_diag_3rd_mom",  0.01,   10.0,  "log",    default=7.0,   legacy_env="SHOC_C_DIAG_3RD_MOM"),
    Parameter("shoc::coeff_kh",        0.01,   0.2,   "log",    default=0.1,   legacy_env="SHOC_COEFF_KH"),
    Parameter("shoc::coeff_km",        0.03,   0.2,   "log",    default=0.1,   legacy_env="SHOC_COEFF_KM"),
    Parameter("shoc::lambda_low",      0.0001, 0.1,   "linear", default=0.001, legacy_env="SHOC_LAMBDA_LOW"),
    Parameter("shoc::lambda_high",     0.0001, 0.1,   "linear", default=0.08,  legacy_env="SHOC_LAMBDA_HIGH"),
    # P3
    Parameter("p3::spa_ccn_to_nc_factor",                 100.0,  4000.0, "linear", default=2000.0,  legacy_env="P3_SPA_CCN_TO_NC_FACTOR"),
    Parameter("p3::cldliq_to_ice_collection_factor",      0.1,    1.0,    "linear", default=0.5,     legacy_env="P3_CLDLIQ_TO_ICE_COLLECTION_FACTOR"),
    Parameter("p3::rain_to_ice_collection_factor",        0.1,    1.0,    "linear", default=1.0,     legacy_env="P3_RAIN_TO_ICE_COLLECTION_FACTOR"),
    Parameter("p3::accretion_prefactor",                  0.01,   100.0,  "log",    default=67.0,    legacy_env="P3_ACCRETION_PREFACTOR"),
    Parameter("p3::deposition_nucleation_exponent",       0.2,    0.304,  "linear", default=0.304,   legacy_env="P3_DEPOSITION_NUCLEATION_EXPONENT"),
    Parameter("p3::max_total_ni",                         5.0e5,  1.8e6,  "log",    default=740.0e3, legacy_env="P3_MAX_TOTAL_NI"),
    Parameter("p3::ice_sedimentation_factor",             0.7,    2.0,    "linear", default=1.0,     legacy_env="P3_ICE_SEDIMENTATION_FACTOR"),
    Parameter("p3::rain_selfcollection_breakup_diameter", 0.0,    500e-6, "linear", default=0.00028, legacy_env="P3_RAIN_SELFCOLLECTION_BREAKUP_DIAMETER"),
    Parameter("p3::autoconversion_prefactor",             10.0,   20000.0,"log",    default=1350.0,  legacy_env="P3_AUTOCONVERSION_PREFACTOR"),
    Parameter("p3::autoconversion_qc_exponent",           2.0,    4.0,    "linear", default=2.47,    legacy_env="P3_AUTOCONVERSION_QC_EXPONENT"),
    Parameter("p3::autoconversion_radius",                25e-6,  50e-6,  "linear", default=25.0e-6, legacy_env="P3_AUTOCONVERSION_RADIUS"),
    # ZM
    Parameter("zm::tau",             1800.0, 7200.0, "linear", default=3600.0),
    Parameter("zm::alfa",            0.05,   0.30,   "linear", default=0.14),
    Parameter("zm::ke",              0.5e-6, 4.0e-6, "log",    default=2.5e-6),
    Parameter("zm::tpert_fac",       0.5,    3.0,    "linear", default=2.0),
    Parameter("zm::tiedke_add",      0.0,    3.0,    "linear", default=0.8),
    Parameter("zm::c0_lnd",          0.0005, 0.005,  "log",    default=0.002),
    Parameter("zm::c0_ocn",          0.0005, 0.005,  "log",    default=0.002),
    Parameter("zm::cape_threshold",  0.0,    500.0,  "linear", default=0.0),
    Parameter("zm::dcape_threshold", 0.0,    0.1,    "linear", default=0.0),
    Parameter("zm::mcsp_t_coeff",    0.05,   0.5,    "linear", default=0.3),
    Parameter("zm::mcsp_q_coeff",    0.0,    0.5,    "linear", default=0.0),
    # num_cin=0 makes CAPE zero everywhere (ZM off) and writes out of bounds in the CAPE code, so start at 1
    Parameter("zm::num_cin",         1,      10,     "linear", kind="int", default=1),
]

BY_KEY = {p.key: p for p in PARAMETERS}
SCHEMES = sorted({p.scheme for p in PARAMETERS})

# Column order of the original list-of-19 JSON files (no parameter names in the file)
LEGACY_KEYS = [p.key for p in PARAMETERS if p.legacy_env]


def select(schemes=None, exclude=()):
    """Registry parameters in registry order, limited to `schemes` and minus `exclude` keys."""
    exclude = set(exclude)
    unknown = exclude - set(BY_KEY)
    if unknown:
        raise ValueError(f"Unknown parameter(s) to exclude: {sorted(unknown)}")
    if schemes is not None:
        bad = set(schemes) - set(SCHEMES)
        if bad:
            raise ValueError(f"Unknown scheme(s) {sorted(bad)}; known: {SCHEMES}")
    return [p for p in PARAMETERS
            if (schemes is None or p.scheme in schemes) and p.key not in exclude]


def _check_registry():
    seen = set()
    for p in PARAMETERS:
        assert "::" in p.key, p.key
        assert p.key not in seen, f"duplicate {p.key}"
        seen.add(p.key)
        assert p.minimum < p.maximum, p.key
        assert p.scale in ("linear", "log"), p.key
        assert p.kind in ("real", "int"), p.key
        if p.scale == "log":
            assert p.minimum > 0, f"{p.key}: log scale needs minimum > 0"
        if p.kind == "int":
            assert p.scale == "linear" and float(p.minimum).is_integer() and float(p.maximum).is_integer(), p.key
    assert len(LEGACY_KEYS) == 19


_check_registry()


if __name__ == "__main__":
    # python ppe_parameters.py  -> print the registry
    for p in PARAMETERS:
        dflt = "" if p.default is None else f"default={p.default:g}"
        inside = "" if p.default is None or p.minimum <= p.default <= p.maximum else "  (default outside range)"
        print(f"{p.key:42s} {p.kind:4s} {p.scale:6s} {p.minimum:<10g} {p.maximum:<10g} {dflt}{inside}")
