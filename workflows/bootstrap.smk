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

# adapter_ref_3p is a plain per-batch on/off switch (no default -> off) for the cLFR SE600 3' linker
# (GTATCTGAGTCC). Opt in per batch with `adapter_ref_3p: true` in that run's config.yaml;
# adapter_ref_3p_file defaults to the packaged copy below.
_set_default_param(
    "adapter_ref_3p_file",
    str(Path(config["params"]["src_dir"]) / "config" / "adapters" / "clfr_se600_3p_linker.fa"),
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


# ---- run identity lines in .snakemake/log/*.snakemake.log -------------------
# The Snakemake log never says which OS process wrote it. These one-line,
# key=value records name the master process (pid/host/user/cwd/cores), so a log
# reader -- a person, or the LLM triage in Agentic_Triage/lfr_watch.py -- can tell
# "this process finished" from "this process vanished". A log that has the start
# line but no end line means the master was killed (kill -9, OOM, server down).
import os
import socket
import getpass
from snakemake.logging import logger


def _run_identity():
    try:
        user = getpass.getuser()
    except Exception:
        user = os.environ.get("USER", "?")
    return "pid=%d ppid=%d host=%s user=%s cwd=%s cores=%s" % (
        os.getpid(), os.getppid(), socket.gethostname(), user, os.getcwd(),
        getattr(workflow, "cores", "?"))


onstart:
    logger.info("[LFR_RUN] start %s cmd=%s" % (_run_identity(), " ".join(sys.argv)))


onsuccess:
    logger.info("[LFR_RUN] end status=success pid=%d" % os.getpid())


onerror:
    logger.info("[LFR_RUN] end status=error pid=%d" % os.getpid())
