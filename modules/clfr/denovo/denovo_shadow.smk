"""Optional, monitoring-only candidate chimera scoring sidecar.

This module is included by denovo_olc.smk so it scores the exact candidates
and pre-filter read pool from the production OLC run. Its outputs are never
inputs to delivery, QC, or olc_done; request denovo/shadow/umi_summary.tsv
explicitly after setting frag_de_novo.shadow_score_model.
"""

SHADOW_SCORE_MODEL = config['frag_de_novo'].get('shadow_score_model', '')


## Per-candidate QC is shared by rule-only gated_switch selection and shadow
## scoring. It is skipped by default because it repeats read placement for
## every candidate rather than only the delivered primary contig.
rule junctionQCAllCandidates_olc:
    input:
        contigs_done="denovo/frag_denovo_done",
        r2=noisy_preprocess_reads,
        probe="denovo/read_filter_probe.tsv"
    output:
        "denovo/junction_qc_candidates.tsv"
    benchmark:
        "Benchmarks/denovo.junctionQCAllCandidates_olc.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        max_span_ratio = qc_setting('max_span_ratio'),
        min_local_span_ratio = qc_setting('min_local_span_ratio'),
        num_processes = config['frag_de_novo']['num_processes'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False),
        needed = CANDIDATE_SELECT_MODE == 'gated_switch' or bool(SHADOW_SCORE_MODEL)
    run:
        if not params.run_parallel or not params.needed:
            shell("touch {output}")
            return

        _name, qc = resolve_qc(input.probe)
        params.max_span_ratio = qc["max_span_ratio"]
        params.min_local_span_ratio = qc["min_local_span_ratio"]

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_junction_qc.py",
                    "--r2 {input.r2}",
                    "--contigs denovo/final_contigs_0.fa",
                    "--all-candidates",
                    "--max-span-ratio {params.max_span_ratio}",
                    "--min-local-span-ratio {params.min_local_span_ratio}",
                    "--num-processes {params.num_processes}",
                    "--out {output}"]
        shell(" ".join(command))


## Scores all candidates for monitoring and compares the rule-selected output
## with the model's lowest-risk candidate. No production rule consumes either
## output.
rule shadowScoreChimera_olc:
    input:
        candidate_qc="denovo/junction_qc_candidates.tsv",
        selection="denovo/candidate_select_report.tsv"
    output:
        candidates="denovo/shadow/candidate_scores.tsv",
        summary="denovo/shadow/umi_summary.tsv"
    benchmark:
        "Benchmarks/denovo.shadowScoreChimera_olc.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        model = SHADOW_SCORE_MODEL,
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel or not params.model:
            shell("mkdir -p denovo/shadow && touch {output.candidates} {output.summary}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_shadow_score.py",
                    "--candidate-qc {input.candidate_qc}",
                    "--selection {input.selection}",
                    "--model {params.model}",
                    "--out-candidates {output.candidates}",
                    "--out-summary {output.summary}"]
        shell(" ".join(command))
