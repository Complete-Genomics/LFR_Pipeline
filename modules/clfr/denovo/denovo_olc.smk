import os

SEQUENCE_TYPE = config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()

# ---------------------------------------------------------------------------
# QC presets.
#
# The underlying QC has several thresholds, but asking an end user to pick
# them is asking the wrong question -- they cannot know what "min_local span
# ratio 0.2" means. What they CAN answer is the scientific one: do you want
# breadth (keep as many fragments as possible) or purity (only fragments you
# would deposit or build a tree from)? So a single `qc_preset` drives
# everything, and the numbers below are what each preset actually produced on
# the ZymoBIOMICS control (3000 UMI, identity vs the known references --
# denovo.md sec 31/32/33):
#
#   auto (DEFAULT) picks between these from the probe's measurement of this
#   sample's cross-barcode identity, so the user never has to know or declare
#   whether they handed the pipeline soil, a mock, or a gut sample. The choice
#   made is written to denovo/qc_decision.tsv.
#
#   Caveat worth knowing: because auto reads the data, two samples in one study
#   can end up under different QC. That is correct per sample but makes their
#   fragment counts less directly comparable, so for a batch meant to be
#   compared, pin the preset explicitly (e.g. qc_preset: high_div for a soil
#   series) rather than letting each sample choose.
#
#   preset      kept    mean id   >=97%   <90%   use when
#   none       100.0%    94.15    30.3%  14.3%   debugging / want raw output
#   high_div    n/a       n/a      n/a     n/a    high-diversity environmental
#                                                 samples (soil, sediment) --
#                                                 see note below
#   sensitive   91.6%    95.52    44.9%   7.1%   rare taxa matter, losing a
#                                                 fragment is worse than a
#                                                 few bad ones
#   balanced    84.2%    95.81    48.0%   5.3%   default; general community
#                                                 profiling
#   strict      23.0%    97.11    71.4%   3.0%   depositing sequences,
#                                                 phylogenetics, anything
#                                                 where a wrong sequence is
#                                                 expensive
#
# high_div exists because sample diversity changes the PRICE of QC, not
# whether it works. Controlled same-subset A/B, 2000 barcodes each
# (denovo.md sec 37):
#
#                              Zymo mock (8 spp)   soil
#   reads dropped by filter          13.6%          32.1%
#   junction_suspect                26.4 -> 15.7    33.5 -> 20.5
#   total assembled bases              -3.9%         -15.5%
#
# The chimera benefit is comparable; the cost is ~4x higher on the diverse
# sample, because far more of its reads look mutually dissimilar. So high_div
# still filters, but keeps post-QC lenient so breadth is not cut twice.
#
# The drop from balanced to strict is a real cliff, not a smooth dial -- the
# extra gate strict turns on is a different metric with much lower yield but
# much higher purity. There is no setting that gives both.
#
# Expert escape hatch: setting any of these keys explicitly in config
# overrides the preset for that key alone.
# ---------------------------------------------------------------------------
QC_PRESETS = {
    "none":      {"read_filter": False, "max_span_ratio": 0.0,
                  "min_local_span_ratio": 0.0},
    # high_div keeps the read filter ON. An earlier version turned it off on
    # the strength of a comparison that turned out to be uncontrolled (the two
    # arms had been assembled with different min_ctg_len); a proper same-subset
    # A/B showed the filter cuts soil's junction_suspect from 33.5% to 20.5%.
    # What is genuinely different for a diverse sample is the PRICE: -15.5% of
    # assembled bases versus -3.9% on the mock. So high_div pays that price but
    # keeps post-QC lenient to protect breadth (denovo.md sec 35/37).
    "high_div":  {"read_filter": True,  "max_span_ratio": 0.15,
                  "min_local_span_ratio": 0.0},
    "sensitive": {"read_filter": True,  "max_span_ratio": 0.15,
                  "min_local_span_ratio": 0.0},
    "balanced":  {"read_filter": True,  "max_span_ratio": 0.25,
                  "min_local_span_ratio": 0.0},
    "strict":    {"read_filter": True,  "max_span_ratio": 0.25,
                  "min_local_span_ratio": 0.20},
}

QC_PRESET_NAME = config['frag_de_novo'].get('qc_preset', 'auto')
# YAML 1.1 turns a bare `off`/`no` into the boolean False (and `on`/`yes` into
# True), so a user writing the obvious `qc_preset: off` never gets a string.
# Accept those spellings rather than failing on something that looks correct.
if QC_PRESET_NAME is False:
    QC_PRESET_NAME = "none"
elif QC_PRESET_NAME is True:
    raise ValueError(
        "frag_de_novo.qc_preset was parsed as the boolean True (YAML reads a "
        "bare on/yes that way); write one of {} instead".format(sorted(QC_PRESETS)))
if isinstance(QC_PRESET_NAME, str):
    QC_PRESET_NAME = {"off": "none", "raw": "none",
                       "soil": "high_div", "environmental": "high_div"}.get(
        QC_PRESET_NAME.strip().lower(), QC_PRESET_NAME.strip().lower())
if QC_PRESET_NAME not in QC_PRESETS and QC_PRESET_NAME != "auto":
    raise ValueError(
        "frag_de_novo.qc_preset={!r} is not a known preset; choose 'auto' or one of {}"
        .format(QC_PRESET_NAME, sorted(QC_PRESETS)))


# ---------------------------------------------------------------------------
# Candidate selection: which of a UMI's assembled candidates (k41_0, k41_1,
# ...) gets delivered. Default 'longest' is the historical, unchanged
# behaviour (always k41_0). Opt-in 'gated_switch' is denovo.md sec 109-120's
# rule-only switch -- NOT flipped as the new default, matching this file's
# existing precedent for the read-filter ML tie-break (sec above,
# frag_de_novo.read_filter_ml_model): validated on a one-time held-out eval
# but with a severe-loss caveat (sec 116/117/120) that needs to be understood
# before turning it on. See denovo_candidate_select.py's module docstring for
# the full rationale.
# ---------------------------------------------------------------------------
CANDIDATE_SELECT_MODE = config['frag_de_novo'].get('candidate_select', 'longest')
if CANDIDATE_SELECT_MODE not in ('longest', 'gated_switch'):
    raise ValueError(
        "frag_de_novo.candidate_select={!r} is not a known mode; choose "
        "'longest' (default) or 'gated_switch'".format(CANDIDATE_SELECT_MODE))
def resolve_qc(probe_path=None):
    """Return (preset_name, settings dict).

    `auto` (the default) takes the preset the probe recommends from this
    sample's own cross-barcode identity distribution, which is why the user
    does not have to know or declare their sample type. It has to be resolved
    here, at rule run time, rather than in `params`: Snakemake evaluates params
    while building the DAG, and the probe output does not exist yet at that
    point.

    Explicitly setting any underlying key in config still wins, so `auto` is a
    default rather than a lock-in.
    """
    name = QC_PRESET_NAME
    if name == "auto":
        chosen = "balanced"
        if probe_path and os.path.exists(probe_path):
            with open(probe_path) as fh:
                probe = dict(line.rstrip("\n").split("\t")[:2]
                              for line in fh if "\t" in line)
            chosen = probe.get("suggested_qc_preset", "balanced")
            if chosen not in QC_PRESETS:
                chosen = "balanced"
        name = chosen
    settings = dict(QC_PRESETS[name])
    for key in list(settings):
        if key in config['frag_de_novo']:
            settings[key] = config['frag_de_novo'][key]
    return name, settings


def qc_setting(key):
    """Preset value for `key`, unless config sets that key explicitly."""
    if key in config['frag_de_novo']:
        return config['frag_de_novo'][key]
    if QC_PRESET_NAME == "auto":
        return QC_PRESETS["balanced"][key]
    return QC_PRESETS[QC_PRESET_NAME][key]


# Preprocessing (barcode selection -> read filtering -> sgrep TSV reformat) lives
# in denovo_preprocess.smk, included unconditionally by workflows/clfr.smk so both
# the megahit (denovo_clfr.smk) and seedext/OLC (this file) assembler branches
# share it -- this file only covers denovo_seed_olc.py onward.

## Pre-assembly read-level contamination filter. Drops the minority
## contaminating reads inside a barcode while keeping the barcode itself.
## Measured on the ZymoBIOMICS control (denovo.md sec 32): mean identity
## 94.15 -> 95.04 at zero yield cost, and it roughly halves how many chimeras
## form at all (junction_suspect 26.5% -> 15.8%). Removing whole low-quality
## UMIs instead was tested and rejected -- at matched yield it was strictly
## worse than just filtering after assembly.
## Cheap probe on a small, artifact-free barcode subset: predicts what the
## full read filter would do before paying for it. Its drop-rate estimate is
## calibrated -- on the Zymo control the probe said 13.46% and the full run
## did 13.59% (denovo.md sec 33). Two uses:
##   - skip the full filter entirely when the library shows no real mixing
##   - surface an unusually high predicted drop rate for a human to look at
##     BEFORE an hour of filtering and 13h of assembly
## It reports the sample's cross-barcode identity distribution as a diagnostic
## but deliberately does NOT auto-tune the conflict threshold from it: that was
## implemented, measured against ground truth, and rejected (sec 33).


## --- Learned read-quality score (denovo.md sec 51/59/62-67): validated as a
## conflict-graph tie-break (mean identity +0.03..+0.07 vs the plain filter,
## significant across 4 independent samples + a 10k-UMI scale check).
## denovo_read_filter.py now has the --ml-model/--ml-features code path
## (previously only exercised as one-off analysis scripts, never as a
## reusable pipeline rule); readFilter_olc below wires it in as an opt-in
## config switch (frag_de_novo.read_filter_ml_model), NOT a new default --
## the reproduced model (subprojects/olc/mlpf/model_identity.lgb, since the
## original was lost with /tmp) has not yet had its own downstream assembly
## validation redone, so the default DAG is unchanged until that happens.
##
## sortFastq_olc keeps quality, which reformat_fasta2 (rule above) discards --
## the score model needs 5 quality-derived features reformat_fasta2's TSV
## cannot supply. samtools sort -t BX rather than sort(1): ~7x faster
## (denovo.md sec 65), and PE is out of scope for now (both rules are SE-only,
## matching denovo_read_features.py and --r2_format fastq).
rule sortFastq_olc:
    input:
        "data/split_read_2_trimmed.fastq.gz"
    output:
        "denovo/data_R2_sorted_qual.fastq.gz"
    benchmark:
        "Benchmarks/denovo.sortFastq_olc.txt"
    params:
        samtools = config['params'].get('samtools', 'samtools'),
        src_dir = config['params']['src_dir'],
        threads = config['frag_de_novo']['num_processes']
    shell:
        "bash {params.src_dir}/modules/clfr/denovo/denovo_sort_fastq.sh "
        "{input} {output} {params.threads} {params.samtools}"


## Per-read features for the score model (denovo_read_features.py). O(n)
## minimizer + per-barcode process pool: 75.5s -> 9.8s on 132k reads (sec 64).
rule extractReadFeatures_olc:
    input:
        "denovo/data_R2_sorted_qual.fastq.gz"
    output:
        "denovo/read_features.tsv"
    benchmark:
        "Benchmarks/denovo.extractReadFeatures_olc.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        num_processes = config['frag_de_novo']['num_processes']
    shell:
        "{params.python} {params.src_dir}/modules/clfr/denovo/denovo_read_features.py "
        "--fastq {input} --out {output} --num_processes {params.num_processes}"


rule readFilterProbe_olc:
    input:
        noisy_preprocess_reads
    output:
        "denovo/read_filter_probe.tsv"
    benchmark:
        "Benchmarks/denovo.readFilterProbe_olc.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        same_molecule_id = config['frag_de_novo'].get('read_filter_identity', 0.90),
        n_barcodes = config['frag_de_novo'].get('probe_barcodes', 600),
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_qc_probe.py",
                    "--r2 {input}",
                    "--same-molecule-id {params.same_molecule_id}",
                    "--n-barcodes {params.n_barcodes}",
                    "--out {output}"]
        shell(" ".join(command))


rule readFilter_olc:
    input:
        reads=noisy_preprocess_reads,
        probe="denovo/read_filter_probe.tsv",
        noisy_qc="denovo/noisy_preprocess_decision.tsv",
        features=("denovo/read_features.tsv"
                  if config['frag_de_novo'].get('read_filter_ml_model') else [])
    output:
        reads="denovo/data_R2_readfilt.tsv",
        dropped="denovo/data_R2_readfilt.dropped.tsv",
        report="denovo/read_filter_report.tsv",
        decision="denovo/qc_decision.tsv"
    benchmark:
        "Benchmarks/denovo.readFilter_olc.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        num_processes = config['frag_de_novo']['num_processes'],
        same_molecule_id = config['frag_de_novo'].get('read_filter_identity', 0.90),
        run_parallel = config['frag_de_novo'].get('run_parallel', False),
        enabled = qc_setting('read_filter'),
        ml_model = config['frag_de_novo'].get('read_filter_ml_model', '')
    run:
        preset_name, qc = resolve_qc(input.probe)
        enabled = qc["read_filter"]
        with open("denovo/qc_decision.tsv", "w") as fh:
            fh.write("setting\tvalue\n")
            fh.write("qc_preset_requested\t{}\n".format(QC_PRESET_NAME))
            fh.write("qc_preset_effective\t{}\n".format(preset_name))
            for k in sorted(qc):
                fh.write("{}\t{}\n".format(k, qc[k]))
        print("QC preset: requested={} effective={} -> {}".format(
            QC_PRESET_NAME, preset_name, qc))

        skip_by_probe = False
        if params.run_parallel and enabled:
            with open(input.probe) as fh:
                probe = dict(line.rstrip("\n").split("\t")[:2]
                              for line in fh if "\t" in line)
            skip_by_probe = probe.get("run_read_filter") == "False"
            if skip_by_probe:
                print("read filter skipped: probe drop rate {} below threshold"
                      .format(probe.get("probe_drop_rate")))

        if not params.run_parallel or not enabled or skip_by_probe:
            # pass through unchanged so downstream rules do not branch
            shell("cp {input.reads} {output.reads} && touch {output.dropped} {output.report}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_read_filter.py",
                    "--r2 {input.reads}",
                    "--same-molecule-id {params.same_molecule_id}",
                    "--num_processes {params.num_processes}",
                    "--out {output.reads}",
                    "--dropped-out {output.dropped}",
                    "--report {output.report}"]
        if params.ml_model:
            command += ["--ml-model {params.ml_model}",
                        "--ml-features denovo/read_features.tsv"]
        shell(" ".join(command))


rule run_denovoOLC_parallel:
    input:
        "denovo/data_R2_readfilt.tsv"
    output:
        "denovo/frag_denovo_done"
    benchmark:
        "Benchmarks/denovo.run_denovoOLC_parallel.txt"
    params:
        num_processes = config['frag_de_novo']['num_processes'],
        sequence_type = config['params']['sequence_type'],
        python = config['params']['general_python'],
        min_ctg_len = config['frag_de_novo']['min_ctg_len'],
        n_umi = config['frag_de_novo'].get('assembly_N_umi'),
        run_parallel = config['frag_de_novo'].get('run_parallel', False),
        max_contigs = config['frag_de_novo'].get('max_contigs', 8),
        src_dir = config['params']['src_dir']
    run:
        if not params.run_parallel:
            shell("mkdir -p denovo && touch {output}")
            return

        # chunk to run on multiple nodes, or a dummy number 300000000 to run on single node
        params.end_idx = config['frag_de_novo'].get('end_idx')
        params.start_idx = config['frag_de_novo'].get('start_idx', 0)

        # pure-Python OLC assembler, no megahit / no subprocess fork per UMI
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_seed_olc.py",
                    "--num_processes {params.num_processes} ",
                    "--sequence_type {params.sequence_type} ",
                    "--r2 denovo/data_R2_readfilt.tsv ",
                    "--n_line_chunk 2000000 ",
                    "--start_idx {params.start_idx} ",
                    "--min_ctg_len {params.min_ctg_len} ",
                    "--max_contigs {params.max_contigs} ",
                    "--nth_of_nodes 0"]
        if params.n_umi not in (None, "", "all"):
            command.append("--n {params.n_umi} ")

        if params.end_idx not in (None, "", "all"):
            command.append("--end_idx {params.end_idx} ")
        shell(" ".join(command))


# Candidate-level shadow scoring is an optional sidecar. It shares this run's
# OLC candidates and raw read pool but never affects delivery or QC; see
# denovo_shadow.smk.
include: config['params']['src_dir'] + "/modules/clfr/denovo/denovo_shadow.smk"


## denovo_seed_olc.py already writes each barcode's contigs longest-first
## (k41_0 == longest). By default this rule still just keeps k41_0
## unconditionally ("keep only the longest per UMI") -- unchanged from
## before candidate selection existed as a choice. Opt-in
## frag_de_novo.candidate_select: gated_switch switches to the first
## k41_rank-ordered candidate passing span_cov_ratio>=0.25 &&
## placed_reads>=2, falling back to k41_0 otherwise (denovo.md sec 109-112);
## read denovo_candidate_select.py's module docstring (severe-loss caveat,
## sec 116/117/120) before turning it on.
## The output filename stays denovo.longest.fasta for every downstream rule
## regardless of mode -- which candidate actually got delivered per barcode
## is always in denovo/candidate_select_report.tsv, never implied by the name.
rule filterOLC_longest:
    input:
        contigs_done="denovo/frag_denovo_done",
        candidate_qc="denovo/junction_qc_candidates.tsv",
        probe="denovo/read_filter_probe.tsv"
    output:
        fasta="denovo/denovo.longest.fasta",
        report="denovo/candidate_select_report.tsv",
        decision="denovo/candidate_select_decision.tsv"
    benchmark:
        "Benchmarks/denovo.filterOLC_longest.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        mode = CANDIDATE_SELECT_MODE,
        max_span_ratio = qc_setting('max_span_ratio'),
        min_placed_reads = config['frag_de_novo'].get('gated_switch_min_placed_reads', 2),
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output.fasta} {output.report} {output.decision}")
            return

        # Re-resolve here (rather than trusting the params default above)
        # so a gated_switch gate threshold under qc_preset=auto matches the
        # same probe-informed value junctionQCAllCandidates_olc used to
        # compute denovo/junction_qc_candidates.tsv -- the params default is
        # only a DAG-build-time placeholder (qc_preset=auto resolves to
        # 'balanced' there, since the probe output does not exist yet).
        _name, qc = resolve_qc(input.probe)
        params.max_span_ratio = qc["max_span_ratio"]

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_candidate_select.py",
                    "--contigs denovo/final_contigs_0.fa",
                    "--mode {params.mode}",
                    "--out-fasta {output.fasta}",
                    "--out-report {output.report}",
                    "--out-decision {output.decision}"]
        if params.mode == "gated_switch":
            command += ["--candidate-qc {input.candidate_qc}",
                        "--max-span-ratio {params.max_span_ratio}",
                        "--min-placed-reads {params.min_placed_reads}"]
        shell(" ".join(command))


rule plotOLC_frag_len_distribution:
    input:
        "denovo/denovo.longest.highconf.fasta"
    output:
        "denovo/frag_length_distribution.pdf"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py",
                    "--fasta denovo/denovo.longest.fasta",
                    "--outdir denovo/",
                    "--module metrics_olc"]
        shell(" ".join(command))

## Everything below is reference-free: production samples are environmental,
## with no curated reference database to check against (the ZymoBIOMICS mock
## community was only ever a control to validate these detectors -- denovo.md
## sec 28/29/31). Three independent signals, all computed from the run's own
## reads and contigs:
##   - libraryQC   (pre-assembly)  is the read pool mixed at all?
##   - junctionQC  (post-assembly) verified-spanning-depth collapse == chimera
##   - readbackQC  (post-assembly) is each contig covered/bracketed by its reads
##   - chimeraQC   (post-assembly) vsearch --uchime_denovo
## Note their measured power differs a lot (denovo.md sec 30/31): junctionQC
## carries the chimera signal (AUC 0.827), while readback/uchime_denovo were at
## chance for chimeras and are kept for the different failure mode they do
## catch (contigs with poor or missing read support).

## Pre-assembly gate: does this library have a read-pool mixing problem at all?
## Reference-free, and cheap enough to run before committing to a full assembly.
## Reports a library-level indicator only -- per-barcode mixing calls were tried
## and did not survive validation (denovo.md sec 30).
rule libraryQC_olc:
    input:
        noisy_preprocess_reads
    output:
        "denovo/library_qc.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        n_barcodes = config['frag_de_novo'].get('library_qc_barcodes', 3000),
        pairs = config['frag_de_novo'].get('library_qc_pairs', 6000),
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_library_qc.py",
                    "--r2 {input}",
                    "--n-barcodes {params.n_barcodes}",
                    "--pairs {params.pairs}",
                    "--out {output}"]
        shell(" ".join(command))


## Post-assembly chimera detection by verified-spanning-depth collapse.
## This is the one that actually discriminates: on the ZymoBIOMICS control it
## reaches AUC 0.827, where plain read-back QC and reference-free uchime_denovo
## were both at chance (denovo.md sec 30/31).
##
## r2 MUST be the pre-read_filter pool, not data_R2_readfilt.tsv: the 0.827
## AUC / max_span_ratio=0.25 default were both calibrated on the raw pool
## (denovo.md sec 91), and re-measuring on the filtered pool is not merely a
## different operating point -- it is a strictly worse one (AUC 0.827->0.771,
## clean-contig false-flag rate 14.2%->17.4%, kept yield 73.5%->71.1%, kept
## mean identity 95.35->95.22), because read_filter drops indel-rich reads
## unevenly and can hollow out the one true-molecule contig's already-thin
## local coverage into a spurious spanning-depth dip.
rule junctionQC_olc:
    input:
        contigs="denovo/denovo.longest.fasta",
        r2=noisy_preprocess_reads,
        probe="denovo/read_filter_probe.tsv"
    output:
        "denovo/junction_qc.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        max_span_ratio = qc_setting('max_span_ratio'),
        min_local_span_ratio = qc_setting('min_local_span_ratio'),
        num_processes = config['frag_de_novo']['num_processes'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        _name, qc = resolve_qc(input.probe)
        params.max_span_ratio = qc["max_span_ratio"]
        params.min_local_span_ratio = qc["min_local_span_ratio"]

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_junction_qc.py",
                    "--r2 {input.r2}",
                    "--contigs {input.contigs}",
                    "--max-span-ratio {params.max_span_ratio}",
                    "--min-local-span-ratio {params.min_local_span_ratio}",
                    "--num_processes {params.num_processes}",
                    "--out {output}"]
        shell(" ".join(command))


rule readbackQC_olc:
    input:
        contigs="denovo/denovo.longest.fasta",
        r2="denovo/data_R2_readfilt.tsv"
    output:
        "denovo/readback_qc.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        num_processes = config['frag_de_novo']['num_processes'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        # --n uncapped: no reference DB to subsample against in production,
        # check read-back support for every assembled UMI
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/test/readback_qc_single.py",
                    "--r2 {input.r2}",
                    "--contigs {input.contigs}",
                    "--n 100000000",
                    "--num_processes {params.num_processes}",
                    "--out {output}"]
        shell(" ".join(command))


rule chimeraQC_olc:
    input:
        "denovo/denovo.longest.fasta"
    output:
        uchimeout="denovo/uchime_denovo_report.tsv",
        chimeras="denovo/chimeras.fasta",
        nonchimeras="denovo/nonchimeras.fasta"
    params:
        vsearch = config['params'].get('vsearch', 'vsearch'),
        # Off by default: measured 0 chimeras on every set tried (0/3000 Zymo,
        # 0/20000 and 0/80000 real soil contigs) while costing an extrapolated
        # 7h at 1.5M and 25h at 3M contigs -- its runtime grows ~O(n^1.76), so
        # it is both the single most expensive QC step and the only one with no
        # measurable benefit (denovo.md sec 30/31/33). Set frag_de_novo.uchime:
        # True to re-enable.
        enabled = config['frag_de_novo'].get('uchime', False),
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel or not params.enabled:
            shell("touch {output.uchimeout} {output.chimeras} {output.nonchimeras}")
            return

        shell("""
            {params.vsearch} --uchime_denovo {input} \
                --uchimeout {output.uchimeout} \
                --chimeras {output.chimeras} \
                --nonchimeras {output.nonchimeras}
        """)


rule combineQC_olc:
    input:
        contigs="denovo/denovo.longest.fasta",
        readback="denovo/readback_qc.tsv",
        uchimeout="denovo/uchime_denovo_report.tsv",
        junction="denovo/junction_qc.tsv"
    output:
        report="denovo/qc_report.tsv",
        highconf="denovo/denovo.longest.highconf.fasta",
        flagged="denovo/denovo.longest.flagged.fasta"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output.report} {output.highconf} {output.flagged}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_qc_combine.py",
                    "--contigs {input.contigs}",
                    "--readback {input.readback}",
                    "--uchimeout {input.uchimeout}",
                    "--junction {input.junction}",
                    "--out-report {output.report}",
                    "--out-highconf-fasta {output.highconf}",
                    "--out-flagged-fasta {output.flagged}"]
        shell(" ".join(command))


# rule filterOLC_fasta_1k:
#     input:
#         inputfile="denovo/denovo.longest.fasta"
#     output:
#         outputfile="denovo/denovo.longest.1k.fasta"
#     params:
#         run_parallel = config['frag_de_novo'].get('run_parallel', False)
#     run:
#         if not params.run_parallel:
#             shell("touch {output.outputfile}")
#             return

#         from Bio import SeqIO
#         outfile = open(output.outputfile, 'w')

#         with open(input.inputfile, "r") as input_fasta:
#             for record in SeqIO.parse(input_fasta, "fasta"):
#                 if len(record.seq) >= 1000:
#                     SeqIO.write(record, outfile, "fasta")
#         outfile.close()

rule olc_done:
    input:
        "denovo/denovo.longest.fasta",
        "denovo/candidate_select_report.tsv",
        "denovo/frag_length_distribution.pdf",
        "denovo/qc_report.tsv",
        "denovo/denovo.longest.highconf.fasta",
        "denovo/library_qc.tsv",
        "denovo/data_R2_readfilt.dropped.tsv"
    output:
        touch("denovo/done.fq")
