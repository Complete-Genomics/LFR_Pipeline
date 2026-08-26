'''
vc_polish: molecule-linkage SNP confidence re-scoring (cLFR_eval Claim A model).

CANARY / optional module -- runs downstream of variant_calling's SNP output and
produces an ADDITIONAL polished VCF + per-SNP decisions. It does not replace,
overwrite, or gate any existing variant_calling output; enabling it only adds
new targets under Make_Vcf/vc_polish/. Toggle with modules.vc_polish (default
False) -- requires modules.variant_calling to also be True.

Model training (01_make_candidates.py / 03_train_eval.py, GIAB-truth-labeled)
lives in cLFR_eval/src, not here -- this module only ships the trained model
(models/) plus the runtime inference path: VCF -> candidates -> features ->
re-score. See models/README.md for the shipped model's provenance and
validation numbers.
'''

rule vc_polish_candidates:
    input:
        vcf = "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz",
        tbi = "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz.tbi"
    output:
        "Make_Vcf/vc_polish/{id}_candidates.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir']
    benchmark:
        "Benchmarks/vc_polish.candidates.{id}.txt"
    shell:
        "{params.python} {params.src_dir}/modules/clfr/vc_polish/vcf_to_candidates.py "
        "--vcf {input.vcf} --out {output}"

rule vc_polish_features:
    input:
        bam = "keep/Align/{id}.sort.markdup.bam",
        bai = "keep/Align/{id}.sort.markdup.bam.bai",
        candidates = "Make_Vcf/vc_polish/{id}_candidates.tsv"
    output:
        "Make_Vcf/vc_polish/{id}_features.tsv"
    threads:
        config['threads'].get('vc_polish', 4)
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        ref = config['params']['ref_fa']
    benchmark:
        "Benchmarks/vc_polish.features.{id}.txt"
    shell:
        "{params.python} {params.src_dir}/modules/clfr/vc_polish/02_extract_features.py "
        "--bam {input.bam} --ref {params.ref} --candidates {input.candidates} "
        "--out {output} --threads {threads} "
        "--molecule-source readname_regex --readname-regex '#([ACGTN]+)'"

rule vc_polish_rescore:
    input:
        features = "Make_Vcf/vc_polish/{id}_features.tsv"
    output:
        vcf = "Make_Vcf/vc_polish/{id}.polished.pass.vcf",
        tsv = "Make_Vcf/vc_polish/{id}.polished.rescored.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        model = config['vc_polish'].get('model') or
                "{}/modules/clfr/vc_polish/models/claimA_snv_confidence.hg002_chr20holdout.2026-08-25.lgb.txt"
                .format(config['params']['src_dir']),
        threshold = config['vc_polish'].get('threshold', 0.5),
        editing_bed = config['vc_polish'].get('editing_bed')
    benchmark:
        "Benchmarks/vc_polish.rescore.{id}.txt"
    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/vc_polish/04_apply_rescore.py",
                   "--features {input.features}",
                   "--model {params.model}",
                   "--out-prefix Make_Vcf/vc_polish/{wildcards.id}.polished",
                   "--threshold {params.threshold}"]
        if params.editing_bed:
            command.append("--editing-bed {params.editing_bed}")
        shell(" ".join(command))
