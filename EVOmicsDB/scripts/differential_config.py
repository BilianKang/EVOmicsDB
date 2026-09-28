"""Dataset-specific differential-analysis configuration.

The configuration is generated from Supplementary Data 1, Table S1.  It is
kept separate from the public dataset model because these fields control the
server-side statistical workflow and are not user-editable plot parameters.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "2.0.0"
SUPPORTED_PROFILES = {
    "rnaseq_raw_counts",
    "rnaseq_fractional_counts",
    "transcriptomics_normalized",
    "microarray",
    "proteomics_label_free",
    "proteomics_tmt",
    "protein_spectral_counts",
    "metabolomics_abundance",
}
ALLOWED_INPUT_SCALES = {"auto", "count", "linear", "log2"}
ALLOWED_NORMALIZATION_TODO = {
    "auto", "none", "median_center", "quantile", "irs", "tmmwsp"
}
ALLOWED_TRANSFORM_TODO = {"auto", "none", "log2", "log2p1"}
ALLOWED_FEATURE_FILTER_MODES = {
    "auto", "none", "mean_detection", "group_valid", "filterByExpr"
}
ALLOWED_VISUAL_IMPUTATION = {"auto", "none", "row_median", "minprob", "knn", "zero"}


class DifferentialConfigError(RuntimeError):
    """Base exception for missing or invalid dataset analysis configuration."""


class DifferentialAnalysisNotSupported(DifferentialConfigError):
    """Raised before R is called when a dataset is not supported."""


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _config_path() -> Path:
    configured = os.getenv("DIFFERENTIAL_ANALYSIS_CONFIG", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    return _project_root() / "config" / "differential_analysis.json"


def _resolve_project_path(value: str | None) -> str | None:
    if value is None or not str(value).strip():
        return None
    path = Path(str(value).strip()).expanduser()
    if not path.is_absolute():
        path = _project_root() / path
    return str(path.resolve())


def _required_text(data: dict[str, Any], key: str) -> str:
    value = str(data.get(key, "")).strip()
    if not value:
        raise DifferentialConfigError(f"Differential configuration is missing '{key}'.")
    return value


def _optional_text(data: dict[str, Any], key: str) -> str | None:
    value = str(data.get(key, "")).strip()
    return value or None


def _text_with_default(data: dict[str, Any], key: str, default: str) -> str:
    """Return a trimmed text value, falling back to a schema default.

    Derivable Schema v2 fields (normalization_todo, transform_todo,
    feature_filter_mode, visual_imputation, ...) may be omitted from the
    registry; the R script resolves 'auto' deterministically per profile.
    """
    raw = data.get(key)
    value = "" if raw is None else str(raw).strip()
    return value or default


def _bool_value(data: dict[str, Any], key: str, default: bool | None = None) -> bool:
    value = data.get(key)
    if isinstance(value, bool):
        return value
    if (value is None or str(value).strip() == "") and default is not None:
        return default
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise DifferentialConfigError(f"Differential configuration '{key}' must be TRUE or FALSE.")


def _float_value(data: dict[str, Any], key: str, default: float) -> float:
    value = data.get(key)
    if value is None or str(value).strip() == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise DifferentialConfigError(
            f"Differential configuration '{key}' must be numeric."
        ) from exc


def _int_value(data: dict[str, Any], key: str, default: int) -> int:
    value = data.get(key)
    if value is None or str(value).strip() == "":
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise DifferentialConfigError(
            f"Differential configuration '{key}' must be an integer."
        ) from exc
    if parsed < 1:
        raise DifferentialConfigError(
            f"Differential configuration '{key}' must be a positive integer."
        )
    return parsed


@dataclass(frozen=True)
class DifferentialDatasetConfig:
    dataset_id: int
    supported: bool
    unsupported_reason: str
    analysis_profile: str
    input_scale: str
    normalization_applied: str
    normalization_todo: str
    transform_applied: str
    transform_todo: str
    missing_value_encoding: str
    zero_is_missing: bool
    feature_filter_mode: str
    min_valid_prop_per_group: float
    min_valid_n_per_group: int
    sample_max_missing_prop: float
    model_imputation: str
    visual_imputation: str
    feature_id_type: str
    output_feature_id_type: str
    organism: str
    organism_taxid: str
    group_source: str
    group_file: str | None
    sample_column: str | None
    group_column: str | None
    case_code: str
    control_code: str
    case_prefix: str
    control_prefix: str
    paired_design: bool
    pair_column: str | None
    batch_column: str | None
    plex_column: str | None
    probe_map_file: str | None
    probe_collapse_method: str | None
    quantification_method: str
    omics_layer: str
    molecule_type: str
    detection_method: str

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DifferentialDatasetConfig":
        try:
            dataset_id = int(data["dataset_id"])
        except (KeyError, TypeError, ValueError) as exc:
            raise DifferentialConfigError("Each configuration requires an integer dataset_id.") from exc

        supported = _bool_value(data, "supported")
        input_scale = _required_text(data, "input_scale")
        feature_id_type = _required_text(data, "feature_id_type")
        profile = _required_text(data, "analysis_profile")
        if supported and profile not in SUPPORTED_PROFILES:
            raise DifferentialConfigError(
                f"Dataset {dataset_id} is marked supported but has invalid analysis_profile={profile!r}."
            )

        # Fail invalid Table S1 values while loading the registry, rather than
        # allowing a malformed deployment row to reach R/argparse at request
        # time.  These are the same enumerations enforced by Schema v2 R.
        # Derivable fields may be omitted; they default to 'auto' and the R
        # script resolves them deterministically from the analysis profile.
        if supported:
            enumerations = {
                "input_scale": (ALLOWED_INPUT_SCALES, "input scale", None),
                "normalization_todo": (ALLOWED_NORMALIZATION_TODO, "normalization_todo", "auto"),
                "transform_todo": (ALLOWED_TRANSFORM_TODO, "transform_todo", "auto"),
                "feature_filter_mode": (ALLOWED_FEATURE_FILTER_MODES, "feature_filter_mode", "auto"),
                "visual_imputation": (ALLOWED_VISUAL_IMPUTATION, "visual_imputation", "auto"),
            }
            for key, (allowed, label, default) in enumerations.items():
                value = (
                    _text_with_default(data, key, default)
                    if default is not None
                    else _required_text(data, key)
                )
                if value not in allowed:
                    raise DifferentialConfigError(
                        f"Dataset {dataset_id} has invalid {label}={value!r}; "
                        f"allowed values are {sorted(allowed)}."
                    )
            if _text_with_default(data, "model_imputation", "none") != "none":
                raise DifferentialConfigError(
                    f"Dataset {dataset_id} requests model imputation, but Schema v2 production analysis forbids it."
                )

        group_source = _required_text(data, "group_source").lower()
        if group_source not in {"metadata", "prefix"}:
            raise DifferentialConfigError(
                f"Dataset {dataset_id} has invalid group_source={group_source!r}."
            )

        group_file = _resolve_project_path(_optional_text(data, "group_file"))
        sample_column = _optional_text(data, "sample_column")
        group_column = _optional_text(data, "group_column")
        if supported and group_source == "metadata":
            if not group_file or not sample_column or not group_column:
                raise DifferentialConfigError(
                    f"Dataset {dataset_id} uses metadata grouping but group_file/sample_column/group_column is incomplete."
                )

        paired_design = _bool_value(data, "paired_design", default=False)
        pair_column = _optional_text(data, "pair_column")
        if supported and paired_design and (group_source != "metadata" or not pair_column):
            raise DifferentialConfigError(
                f"Dataset {dataset_id} paired design requires metadata grouping and pair_column."
            )

        return cls(
            dataset_id=dataset_id,
            supported=supported,
            unsupported_reason=str(data.get("unsupported_reason", "")).strip(),
            analysis_profile=profile,
            input_scale=input_scale,
            normalization_applied=_text_with_default(data, "normalization_applied", "unknown"),
            normalization_todo=_text_with_default(data, "normalization_todo", "auto"),
            transform_applied=_text_with_default(
                data,
                "transform_applied",
                "log2" if input_scale == "log2" else "none",
            ),
            transform_todo=_text_with_default(data, "transform_todo", "auto"),
            missing_value_encoding=_text_with_default(data, "missing_value_encoding", "NA"),
            zero_is_missing=_bool_value(data, "zero_is_missing", default=False),
            feature_filter_mode=_text_with_default(data, "feature_filter_mode", "auto"),
            min_valid_prop_per_group=_float_value(
                data, "min_valid_prop_per_group", 0.5
            ),
            min_valid_n_per_group=_int_value(data, "min_valid_n_per_group", 3),
            sample_max_missing_prop=_float_value(
                data, "sample_max_missing_prop", 1.0
            ),
            model_imputation=_text_with_default(data, "model_imputation", "none"),
            visual_imputation=_text_with_default(data, "visual_imputation", "auto"),
            feature_id_type=feature_id_type,
            output_feature_id_type=_text_with_default(
                data, "output_feature_id_type", feature_id_type
            ),
            organism=_text_with_default(data, "organism", "Homo sapiens"),
            organism_taxid=_text_with_default(data, "organism_taxid", "9606"),
            group_source=group_source,
            group_file=group_file,
            sample_column=sample_column,
            group_column=group_column,
            case_code=_text_with_default(data, "case_code", "C"),
            control_code=_text_with_default(data, "control_code", "N"),
            case_prefix=_text_with_default(data, "case_prefix", "^C"),
            control_prefix=_text_with_default(data, "control_prefix", "^N"),
            paired_design=paired_design,
            pair_column=pair_column,
            batch_column=_optional_text(data, "batch_column"),
            plex_column=_optional_text(data, "plex_column"),
            probe_map_file=_resolve_project_path(_optional_text(data, "probe_map_file")),
            probe_collapse_method=_optional_text(data, "probe_collapse_method"),
            quantification_method=str(data.get("quantification_method", "")).strip(),
            omics_layer=str(data.get("omics_layer", "")).strip(),
            molecule_type=str(data.get("molecule_type", "")).strip(),
            detection_method=str(data.get("detection_method", "")).strip(),
        )

    def require_supported(self) -> None:
        if self.supported:
            return
        reason = self.unsupported_reason or "The dataset is not enabled for differential analysis."
        raise DifferentialAnalysisNotSupported(
            f"Dataset {self.dataset_id} does not support differential analysis: {reason}"
        )

    def validate_runtime_files(self, counts_file: str) -> None:
        counts_path = Path(counts_file).expanduser().resolve()
        if not counts_path.is_file():
            raise DifferentialConfigError(f"Expression matrix does not exist: {counts_path}")
        if self.group_source == "metadata":
            if not self.group_file or not Path(self.group_file).is_file():
                raise DifferentialConfigError(
                    f"Dataset {self.dataset_id} metadata file does not exist: {self.group_file}"
                )
        if self.probe_map_file and not Path(self.probe_map_file).is_file():
            raise DifferentialConfigError(
                f"Dataset {self.dataset_id} probe map does not exist: {self.probe_map_file}"
            )

    def r_arguments(self) -> dict[str, Any]:
        args: dict[str, Any] = {
            "--analysis_profile": self.analysis_profile,
            "--quantification_method": self.quantification_method,
            "--omics_layer": self.omics_layer,
            "--molecule_type": self.molecule_type,
            "--detection_method": self.detection_method,
            "--input_scale": self.input_scale,
            "--normalization_applied": self.normalization_applied,
            "--normalization_todo": self.normalization_todo,
            "--transform_applied": self.transform_applied,
            "--transform_todo": self.transform_todo,
            "--missing_value_encoding": self.missing_value_encoding,
            "--zero_is_missing": self.zero_is_missing,
            "--feature_filter_mode": self.feature_filter_mode,
            "--min_valid_prop_per_group": self.min_valid_prop_per_group,
            "--min_valid_n_per_group": self.min_valid_n_per_group,
            "--sample_max_missing_prop": self.sample_max_missing_prop,
            "--model_imputation": self.model_imputation,
            "--visual_imputation": self.visual_imputation,
            "--feature_id_type": self.feature_id_type,
            "--output_feature_id_type": self.output_feature_id_type,
            "--organism": self.organism,
            "--organism_taxid": self.organism_taxid,
            "--case_code": self.case_code,
            "--control_code": self.control_code,
            "--case_prefix": self.case_prefix,
            "--control_prefix": self.control_prefix,
            "--paired_design": self.paired_design,
        }
        if self.group_source == "metadata":
            args.update(
                {
                    "--group_file": self.group_file,
                    "--sample_column": self.sample_column,
                    "--group_column": self.group_column,
                }
            )
        if self.paired_design:
            args["--pair_column"] = self.pair_column
        if self.batch_column:
            args["--batch_column"] = self.batch_column
        if self.plex_column:
            args["--plex_column"] = self.plex_column
        if self.probe_map_file:
            args["--probe_map_file"] = self.probe_map_file
        if self.probe_collapse_method:
            args["--probe_collapse_method"] = self.probe_collapse_method
        return args

    def fingerprint(self, counts_file: str, script_path: str) -> str:
        script = Path(script_path)
        payload = {
            "config": self.__dict__,
            "counts": _file_signature(Path(counts_file)),
            "script": _file_signature(script),
            "rda_utils": _optional_file_signature(script.with_name("evomics_rda_utils.R")),
            "group_file": _optional_file_signature(Path(self.group_file)) if self.group_file else None,
            "probe_map_file": _optional_file_signature(Path(self.probe_map_file)) if self.probe_map_file else None,
        }
        serialized = json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
        return hashlib.sha256(serialized).hexdigest()[:16]


def _file_signature(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _optional_file_signature(path: Path) -> dict[str, Any] | None:
    return _file_signature(path) if path.expanduser().exists() else None


@lru_cache(maxsize=1)
def _load_registry() -> dict[int, DifferentialDatasetConfig]:
    path = _config_path()
    if not path.is_file():
        raise DifferentialConfigError(
            f"Differential configuration file was not found: {path}"
        )
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if str(payload.get("schema_version", "")) != SCHEMA_VERSION:
        raise DifferentialConfigError(
            f"Differential configuration schema must be {SCHEMA_VERSION}; "
            f"found {payload.get('schema_version')!r}."
        )
    rows = payload.get("datasets")
    if not isinstance(rows, list):
        raise DifferentialConfigError("Differential configuration 'datasets' must be a list.")
    registry: dict[int, DifferentialDatasetConfig] = {}
    for row in rows:
        if not isinstance(row, dict):
            raise DifferentialConfigError("Each differential configuration row must be an object.")
        config = DifferentialDatasetConfig.from_dict(row)
        if config.dataset_id in registry:
            raise DifferentialConfigError(
                f"Duplicate differential configuration for dataset {config.dataset_id}."
            )
        registry[config.dataset_id] = config
    return registry


def get_differential_config(
    dataset_id: int | str, *, require_supported: bool = True
) -> DifferentialDatasetConfig:
    try:
        key = int(dataset_id)
    except (TypeError, ValueError) as exc:
        raise DifferentialConfigError(
            f"Dataset ID must be an integer; received {dataset_id!r}."
        ) from exc
    config = _load_registry().get(key)
    if config is None:
        raise DifferentialConfigError(
            f"Dataset {key} has no schema {SCHEMA_VERSION} differential configuration."
        )
    if require_supported:
        config.require_supported()
    return config


def iter_differential_configs() -> tuple[DifferentialDatasetConfig, ...]:
    """Return the validated Schema v2 registry in dataset order."""

    return tuple(_load_registry()[key] for key in sorted(_load_registry()))


def clear_differential_config_cache() -> None:
    """Test/deployment helper for reloading a replaced configuration file."""

    _load_registry.cache_clear()
