from pathlib import Path
import sys


PIPELINE_ROOT = Path(workflow.basedir).resolve().parent
sys.path.insert(0, str(PIPELINE_ROOT / "src"))

from run_metadata import record_run_metadata_safely


def _is_placeholder(value):
    value = str(value or "")
    return (
        value == ""
        or value.startswith("/path/to")
        or value.startswith("/opt/cgi-tools")
    )


def _set_default_param(name, default):
    if _is_placeholder(config["params"].get(name)):
        config["params"][name] = default


config.setdefault("params", {})

if _is_placeholder(config["params"].get("src_dir")):
    config["params"]["src_dir"] = str(PIPELINE_ROOT) + "/"

_set_default_param(
    "adapter_ref",
    str(Path(config["params"]["src_dir"]) / "config" / "adapters" / "mgi_dnbseq_adapters.fa"),
)

for _name, _default in {
    "gatk_install": "gatk",
    "calc_frag_python": "python",
    "gcbias_python": "python",
    "general_python": "python",
    "rtg_install": "rtg",
    "tabix": "tabix",
    "bcftools": "bcftools",
    "star": "STAR",
    "hisat2": "hisat2",
    "minimap2": "minimap2",
    "megahit": "megahit",
    "bbduk": "bbduk.sh",
    "bgzip": "bgzip",
    "featurecounts": "featureCounts",
}.items():
    _set_default_param(_name, _default)


onstart:
    record_run_metadata_safely(PIPELINE_ROOT, config)
