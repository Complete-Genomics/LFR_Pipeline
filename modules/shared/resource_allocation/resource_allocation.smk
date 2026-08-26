"""Shadow-only resource recommendations from completed cLFR benchmark history."""

RESOURCE_ALLOCATION = config.get("resource_allocation", {})

if config["modules"].get("resource_allocation", False):
    history = RESOURCE_ALLOCATION.get("history", "")
    candidates = RESOURCE_ALLOCATION.get("candidates", "")
    if not history or not candidates:
        raise ValueError(
            "resource_allocation.history and resource_allocation.candidates are required "
            "when modules.resource_allocation is true"
        )

    rule resource_allocation_shadow:
        input:
            history = history,
            candidates = candidates
        output:
            RESOURCE_ALLOCATION.get("shadow_output", "resource_allocation/shadow.tsv")
        params:
            python = config["params"]["general_python"],
            src_dir = config["params"]["src_dir"],
            quantile = RESOURCE_ALLOCATION.get("quantile", 0.95),
            memory_margin = RESOURCE_ALLOCATION.get("memory_margin", 1.30),
            runtime_margin = RESOURCE_ALLOCATION.get("runtime_margin", 1.25),
            min_rule_samples = RESOURCE_ALLOCATION.get("min_rule_samples", 20),
            min_bucket_samples = RESOURCE_ALLOCATION.get("min_bucket_samples", 10)
        shell:
            "{params.python} {params.src_dir}/modules/shared/resource_allocation/resource_allocation.py "
            "shadow --history {input.history} --candidates {input.candidates} --output {output} "
            "--quantile {params.quantile} --memory-margin {params.memory_margin} "
            "--runtime-margin {params.runtime_margin} --min-rule-samples {params.min_rule_samples} "
            "--min-bucket-samples {params.min_bucket_samples}"
