from pathlib import Path
import sys


PIPELINE_ROOT = Path(workflow.basedir).resolve().parent


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
# adapter_ref_5p is a plain per-batch on/off switch (no default -> off).
# The 5' contaminant (adapter183) is confirmed only for the cLFR/stLFR mRNA
# batches investigated so far (bare/chr18/E250041951 -- see
# LFR_Pipeline/memory/project_lfr_chr18_troubleshoot.md). Opt in per batch
# by setting `adapter_ref_5p: true` in that run's own config.yaml -- do not
# flip it on for every sample/protocol this pipeline runs. The reference
# file path always defaults to the packaged copy below; override
# adapter_ref_5p_file only if some other batch needs a different sequence.
_set_default_param(
    "adapter_ref_5p_file",
    str(Path(config["params"]["src_dir"]) / "config" / "adapters" / "mgi_dnbseq_5p_contam.fa"),
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
