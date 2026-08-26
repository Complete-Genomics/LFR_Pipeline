#!/usr/bin/env python3
"""Offline counterfactual learning-to-rank benchmark for read-filter A1.

The production A1 policy uses predicted single-read identity only when several
nodes tie for the highest conflict degree.  This benchmark replaces that proxy
target with the downstream decision target:

1. follow the index-order production policy to a tie state;
2. intervene once on every tied candidate in turn;
3. finish that graph with the fixed index-order policy;
4. reassemble the remaining reads and score the primary contig against a known
   control reference;
5. train a LambdaRank model to rank candidates within that decision state.

Training states and evaluation barcodes are disjoint.  A final held-out policy
run compares index-order, the current identity tie-break, and the learned
counterfactual ranker.  Nothing in this file is connected to the production
workflow or changes its default output.
"""
import argparse
import csv
import json
import math
import multiprocessing as mp
import os
import random
import subprocess
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np

import denovo_read_features as drf
import denovo_read_filter as read_filter
import denovo_seed_olc as olc


BASE_FEATURES = drf.FIELDS[2:]
CONTEXT_FEATURES = [
    "identity_prediction",
    "active_neighbor_identity_mean",
    "active_neighbor_identity_min",
    "active_neighbor_identity_max",
    "active_neighbor_degree_mean",
    "active_component_size",
    "active_degree_fraction",
    "remaining_fraction",
    "step_fraction",
    "tie_size_fraction",
]
RANK_FEATURES = BASE_FEATURES + CONTEXT_FEATURES


def write_tsv(path, rows, fieldnames=None):
    if not rows and not fieldnames:
        return
    names = fieldnames or list(rows[0])
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_groups(r2_path, feature_path, identity_model_path, n_barcodes, seed):
    import lightgbm as lgb

    booster = lgb.Booster(model_file=identity_model_path)
    available = []
    for barcode, lines, seqs, feature_rows in read_filter.iter_joined_groups(
            r2_path, feature_path):
        if not feature_rows:
            continue
        ids = [line.split("\t", 1)[0][read_filter.ID_START:] for line in lines]
        feature_by_id = {read_id: values for read_id, values in feature_rows}
        if any(read_id not in feature_by_id for read_id in ids):
            continue
        base = np.asarray([feature_by_id[read_id] for read_id in ids], dtype=float)
        available.append({
            "barcode": barcode,
            "read_ids": ids,
            "seqs": seqs,
            "base": base,
            "identity_prediction": np.asarray(booster.predict(base), dtype=float),
        })

    if n_barcodes > len(available):
        raise SystemExit("requested {} barcodes but only {} had complete feature rows".format(
            n_barcodes, len(available)))
    rng = random.Random(seed)
    selected = rng.sample(available, n_barcodes)
    selected.sort(key=lambda group: group["barcode"])
    return selected, len(available)


def active_degrees(conflicts, n, dropped):
    return {i: len(conflicts[i] - dropped)
            for i in range(n) if i not in dropped}


def tied_worst_candidates(conflicts, n, dropped):
    degrees = active_degrees(conflicts, n, dropped)
    worst_degree = max(degrees.values(), default=0)
    if worst_degree == 0:
        return [], degrees
    return [i for i, degree in degrees.items() if degree == worst_degree], degrees


def finish_index_policy(conflicts, n, dropped):
    dropped = set(dropped)
    while True:
        candidates, _ = tied_worst_candidates(conflicts, n, dropped)
        if not candidates:
            return dropped
        dropped.add(candidates[0])


def component_size(start, conflicts, n, dropped):
    seen = {start}
    pending = [start]
    while pending:
        node = pending.pop()
        for neighbor in conflicts[node] - dropped:
            if neighbor < n and neighbor not in seen:
                seen.add(neighbor)
                pending.append(neighbor)
    return len(seen)


def rank_features(group, conflicts, n, dropped, candidates, degrees):
    remaining = max(1, n - len(dropped))
    rows = []
    for candidate in candidates:
        neighbors = sorted(conflicts[candidate] - dropped)
        neighbor_scores = group["identity_prediction"][neighbors]
        neighbor_degrees = [degrees[index] for index in neighbors]
        context = [
            group["identity_prediction"][candidate],
            float(np.mean(neighbor_scores)),
            float(np.min(neighbor_scores)),
            float(np.max(neighbor_scores)),
            float(np.mean(neighbor_degrees)),
            float(component_size(candidate, conflicts, n, dropped)),
            degrees[candidate] / remaining,
            remaining / n,
            len(dropped) / n,
            len(candidates) / remaining,
        ]
        rows.append(np.concatenate([group["base"][candidate], context]))
    return np.asarray(rows, dtype=float)


def intervention_states(groups, same_molecule_id, min_overlap, max_reads):
    states = []
    tasks = []
    for group in groups:
        n, conflicts = read_filter.build_conflict_graph(
            group["seqs"], same_molecule_id, min_overlap, max_reads)
        group["n_checked"] = n
        group["conflicts"] = conflicts
        dropped = set()
        state_number = 0
        while True:
            candidates, degrees = tied_worst_candidates(conflicts, n, dropped)
            if not candidates:
                break
            if len(candidates) > 1:
                state_id = "{}__{:03d}".format(group["barcode"], state_number)
                features = rank_features(
                    group, conflicts, n, dropped, candidates, degrees)
                state_rows = []
                for row_number, candidate in enumerate(candidates):
                    final_drop = finish_index_policy(
                        conflicts, n, dropped | {candidate})
                    remaining = [seq for index, seq in enumerate(group["seqs"])
                                 if index not in final_drop]
                    query_id = "cf__{}__{:03d}__{:03d}".format(
                        group["barcode"], state_number, row_number)
                    tasks.append((query_id, remaining))
                    state_rows.append({
                        "query_id": query_id,
                        "barcode": group["barcode"],
                        "state_id": state_id,
                        "state_number": state_number,
                        "candidate": candidate,
                        "baseline_choice": int(candidate == candidates[0]),
                        "identity_prediction": group["identity_prediction"][candidate],
                        "features": features[row_number],
                        "n_remaining": len(remaining),
                    })
                states.append(state_rows)
                state_number += 1
            dropped.add(candidates[0])
    return states, tasks


def assemble_primary(task):
    query_id, seqs = task
    contigs = olc.assemble_umi(seqs)
    contigs = olc._polish_all(contigs, seqs, 20, 0.05, 10)
    primary = contigs[0] if contigs else ""
    return query_id, primary, len(seqs)


def assemble_tasks(tasks, processes):
    if processes == 1:
        return [assemble_primary(task) for task in tasks]
    with mp.Pool(processes) as pool:
        return list(pool.imap(assemble_primary, tasks, chunksize=4))


def write_fasta(path, assembled):
    with open(path, "w") as handle:
        for query_id, sequence, _n_reads in assembled:
            if sequence:
                handle.write(">{}\n{}\n".format(query_id, sequence))


def score_fasta(vsearch, fasta, reference, hits, threads):
    command = [
        vsearch, "--usearch_global", fasta, "--db", reference,
        "--id", "0.5", "--strand", "both", "--maxaccepts", "0",
        "--maxrejects", "0", "--top_hits_only", "--userout", hits,
        "--userfields", "query+target+id+ql+alnlen", "--threads", str(threads),
    ]
    subprocess.run(command, check=True)
    best = {}
    with open(hits) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            query = fields[0]
            identity = float(fields[2])
            if identity > best.get(query, (-1.0, 0, 0))[0]:
                best[query] = (identity, int(fields[3]), int(fields[4]))
    return best, command


def assign_rewards(states, scores, assembled):
    lengths = {query: len(sequence) for query, sequence, _n_reads in assembled}
    for state in states:
        unique_rewards = set()
        for row in state:
            identity, query_length, alignment_length = scores.get(
                row["query_id"], (0.0, lengths.get(row["query_id"], 0), 0))
            row["reward_identity"] = identity
            row["contig_length"] = query_length
            row["alignment_length"] = alignment_length
            unique_rewards.add(identity)
        ordered = {value: rank for rank, value in enumerate(sorted(unique_rewards))}
        for row in state:
            row["relevance"] = ordered[row["reward_identity"]]


def split_barcodes(groups, train_fraction, seed):
    barcodes = [group["barcode"] for group in groups]
    random.Random(seed + 1).shuffle(barcodes)
    n_train = max(1, min(len(barcodes) - 1,
                         int(round(len(barcodes) * train_fraction))))
    return set(barcodes[:n_train]), set(barcodes[n_train:])


def train_ranker(states, train_barcodes, out_model, seed):
    import lightgbm as lgb

    training_states = [state for state in states
                       if state[0]["barcode"] in train_barcodes]
    informative = [state for state in training_states
                   if len({row["reward_identity"] for row in state}) > 1]
    if not informative:
        raise SystemExit("no informative counterfactual training states")
    x = np.vstack([row["features"] for state in informative for row in state])
    y = np.asarray([row["relevance"] for state in informative for row in state])
    group_sizes = [len(state) for state in informative]
    ranker = lgb.LGBMRanker(
        objective="lambdarank", n_estimators=200, learning_rate=0.04,
        max_depth=4, num_leaves=15, min_child_samples=20,
        random_state=seed, verbose=-1,
    )
    ranker.fit(x, y, group=group_sizes)
    ranker.booster_.save_model(out_model)
    return ranker, len(training_states), len(informative)


def chosen_row(state, method, ranker=None):
    if method == "index":
        return min(state, key=lambda row: row["candidate"])
    if method == "identity":
        return min(state, key=lambda row: (row["identity_prediction"],
                                           row["candidate"]))
    if method == "oracle":
        return max(state, key=lambda row: (row["reward_identity"],
                                           -row["candidate"]))
    predictions = ranker.predict(np.vstack([row["features"] for row in state]))
    return max(zip(predictions, state), key=lambda pair: (pair[0],
                                                           -pair[1]["candidate"]))[1]


def decision_metrics(states, test_barcodes, ranker):
    from scipy.stats import binomtest, wilcoxon

    test_states = [state for state in states if state[0]["barcode"] in test_barcodes]
    rows = []
    chosen = {}
    for method in ("index", "identity", "counterfactual_ltr", "oracle"):
        rewards = []
        for state in test_states:
            row = chosen_row(state, method, ranker)
            rewards.append(row["reward_identity"])
            chosen[(state[0]["state_id"], method)] = row["reward_identity"]
        rows.append({
            "method": method,
            "n_states": len(rewards),
            "mean_counterfactual_identity": float(np.mean(rewards)),
            "below_90_pct": 100.0 * float(np.mean(np.asarray(rewards) < 90.0)),
            "oracle_choice_rate": float(np.mean([
                chosen[(state[0]["state_id"], method)] ==
                chosen[(state[0]["state_id"], "oracle")]
                for state in test_states
            ])) if method == "oracle" else float("nan"),
        })

    oracle = {state[0]["state_id"]: chosen_row(state, "oracle")["reward_identity"]
              for state in test_states}
    for result in rows:
        method = result["method"]
        values = [chosen[(state[0]["state_id"], method)] for state in test_states]
        result["mean_regret_to_oracle"] = float(np.mean([
            oracle[state[0]["state_id"]] - value
            for state, value in zip(test_states, values)
        ]))
        result["oracle_choice_rate"] = float(np.mean([
            value == oracle[state[0]["state_id"]]
            for state, value in zip(test_states, values)
        ]))

    pair_rows = []
    for state in test_states:
        identity = chosen_row(state, "identity")["reward_identity"]
        learned = chosen_row(state, "counterfactual_ltr", ranker)["reward_identity"]
        pair_rows.append(learned - identity)
    nonzero = [delta for delta in pair_rows if delta != 0]
    return rows, {
        "n_test_states": len(test_states),
        "ltr_vs_identity_wins": sum(delta > 0 for delta in pair_rows),
        "ltr_vs_identity_losses": sum(delta < 0 for delta in pair_rows),
        "ltr_vs_identity_ties": sum(delta == 0 for delta in pair_rows),
        "ltr_vs_identity_mean_delta": float(np.mean(pair_rows)),
        "ltr_vs_identity_sign_p": float(binomtest(
            sum(delta > 0 for delta in nonzero), len(nonzero)).pvalue),
        "ltr_vs_identity_wilcoxon_p": float(wilcoxon(pair_rows).pvalue),
    }


def policy_drop(group, method, ranker):
    conflicts = group["conflicts"]
    n = group["n_checked"]
    dropped = set()
    while True:
        candidates, degrees = tied_worst_candidates(conflicts, n, dropped)
        if not candidates:
            return dropped
        if len(candidates) == 1 or method == "index":
            choice = candidates[0]
        elif method == "identity":
            choice = min(candidates,
                         key=lambda index: (group["identity_prediction"][index], index))
        else:
            features = rank_features(
                group, conflicts, n, dropped, candidates, degrees)
            predictions = ranker.predict(features)
            choice = max(zip(predictions, candidates),
                         key=lambda pair: (pair[0], -pair[1]))[1]
        dropped.add(choice)


def evaluate_policies(groups, test_barcodes, ranker, out_dir, vsearch,
                      reference, processes):
    from scipy.stats import binomtest, wilcoxon

    test_groups = [group for group in groups if group["barcode"] in test_barcodes]
    tasks = []
    drop_counts = {}
    for method in ("index", "identity", "counterfactual_ltr"):
        for group in test_groups:
            dropped = policy_drop(group, method, ranker)
            drop_counts[(method, group["barcode"])] = len(dropped)
            remaining = [seq for index, seq in enumerate(group["seqs"])
                         if index not in dropped]
            tasks.append(("{}__{}".format(method, group["barcode"]), remaining))
    assembled = assemble_tasks(tasks, processes)
    fasta = os.path.join(out_dir, "heldout_policy_contigs.fa")
    hits = os.path.join(out_dir, "heldout_policy_hits.tsv")
    write_fasta(fasta, assembled)
    scores, command = score_fasta(vsearch, fasta, reference, hits, processes)

    lengths = {query: len(sequence) for query, sequence, _n_reads in assembled}
    per_barcode = []
    for group in test_groups:
        row = {"barcode": group["barcode"], "n_reads": len(group["seqs"])}
        for method in ("index", "identity", "counterfactual_ltr"):
            query = "{}__{}".format(method, group["barcode"])
            identity, _query_length, _alignment_length = scores.get(
                query, (float("nan"), lengths.get(query, 0), 0))
            row[method + "_identity"] = identity
            row[method + "_contig_length"] = lengths.get(query, 0)
            row[method + "_n_dropped"] = drop_counts[(method, group["barcode"])]
        per_barcode.append(row)

    metrics = []
    for method in ("index", "identity", "counterfactual_ltr"):
        identities = np.asarray([row[method + "_identity"] for row in per_barcode])
        valid = np.isfinite(identities)
        metrics.append({
            "method": method,
            "n_barcodes": len(per_barcode),
            "contig_hit_rate": float(np.mean(valid)),
            "mean_identity": float(np.mean(identities[valid])),
            "below_90_pct": 100.0 * float(np.mean(identities[valid] < 90.0)),
            "mean_contig_length": float(np.mean([
                row[method + "_contig_length"] for row in per_barcode])),
            "mean_reads_dropped": float(np.mean([
                row[method + "_n_dropped"] for row in per_barcode])),
        })
    paired = {}
    for baseline in ("index", "identity"):
        deltas = []
        for row in per_barcode:
            before = row[baseline + "_identity"]
            after = row["counterfactual_ltr_identity"]
            if math.isfinite(before) and math.isfinite(after):
                deltas.append(after - before)
        paired["ltr_vs_" + baseline] = {
            "n": len(deltas),
            "mean_delta": float(np.mean(deltas)),
            "wins": sum(delta > 0 for delta in deltas),
            "losses": sum(delta < 0 for delta in deltas),
            "ties": sum(delta == 0 for delta in deltas),
            "sign_p": float(binomtest(
                sum(delta > 0 for delta in deltas if delta != 0),
                sum(delta != 0 for delta in deltas)).pvalue),
            "wilcoxon_p": float(wilcoxon(deltas).pvalue),
        }
    return metrics, per_barcode, paired, command


def format_value(value):
    if isinstance(value, str):
        return value
    if not math.isfinite(float(value)):
        return "NA"
    return "{:.4f}".format(float(value))


def markdown_table(rows, fields):
    lines = ["| " + " | ".join(fields) + " |",
             "|" + "|".join(["---"] * len(fields)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(format_value(row[field]) for field in fields) + " |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--identity-model", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--vsearch", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-barcodes", type=int, default=600)
    parser.add_argument("--train-fraction", type=float, default=0.70)
    parser.add_argument("--same-molecule-id", type=float, default=0.90)
    parser.add_argument("--min-overlap", type=int, default=150)
    parser.add_argument("--max-reads", type=int, default=25)
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    groups, n_available = load_groups(
        args.r2, args.features, args.identity_model, args.n_barcodes, args.seed)
    train_barcodes, test_barcodes = split_barcodes(
        groups, args.train_fraction, args.seed)
    split_rows = [{"barcode": group["barcode"],
                   "split": "train" if group["barcode"] in train_barcodes else "test"}
                  for group in groups]
    write_tsv(os.path.join(args.out_dir, "barcode_split.tsv"), split_rows)

    states, intervention_tasks = intervention_states(
        groups, args.same_molecule_id, args.min_overlap, args.max_reads)
    assembled = assemble_tasks(intervention_tasks, args.processes)
    intervention_fasta = os.path.join(args.out_dir, "counterfactual_contigs.fa")
    intervention_hits = os.path.join(args.out_dir, "counterfactual_hits.tsv")
    write_fasta(intervention_fasta, assembled)
    counterfactual_scores, intervention_command = score_fasta(
        args.vsearch, intervention_fasta, args.reference,
        intervention_hits, args.processes)
    assign_rewards(states, counterfactual_scores, assembled)
    action_gaps = [
        max(row["reward_identity"] for row in state) -
        min(row["reward_identity"] for row in state)
        for state in states
    ]

    model_path = os.path.join(args.out_dir, "model_counterfactual_rank.lgb")
    ranker, n_train_states, n_informative_train_states = train_ranker(
        states, train_barcodes, model_path, args.seed)
    decision_rows, paired = decision_metrics(states, test_barcodes, ranker)
    policy_rows, policy_per_barcode, policy_paired, policy_command = evaluate_policies(
        groups, test_barcodes, ranker, args.out_dir, args.vsearch,
        args.reference, args.processes)

    dataset_rows = []
    for state in states:
        for row in state:
            out = {key: row[key] for key in (
                "query_id", "barcode", "state_id", "state_number", "candidate",
                "baseline_choice", "identity_prediction", "n_remaining",
                "reward_identity", "contig_length", "alignment_length", "relevance")}
            out.update({name: row["features"][i]
                        for i, name in enumerate(RANK_FEATURES)})
            out["split"] = "train" if row["barcode"] in train_barcodes else "test"
            dataset_rows.append(out)
    write_tsv(os.path.join(args.out_dir, "counterfactual_dataset.tsv"), dataset_rows)
    write_tsv(os.path.join(args.out_dir, "decision_metrics.tsv"), decision_rows)
    write_tsv(os.path.join(args.out_dir, "policy_metrics.tsv"), policy_rows)
    write_tsv(os.path.join(args.out_dir, "policy_per_barcode.tsv"), policy_per_barcode)
    importance_rows = sorted([
        {"feature": name, "importance_gain": value}
        for name, value in zip(RANK_FEATURES,
                               ranker.booster_.feature_importance(importance_type="gain"))
    ], key=lambda row: -row["importance_gain"])
    write_tsv(os.path.join(args.out_dir, "feature_importance.tsv"), importance_rows)

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "arguments": vars(args),
        "n_available_barcodes": n_available,
        "n_selected_barcodes": len(groups),
        "n_train_barcodes": len(train_barcodes),
        "n_test_barcodes": len(test_barcodes),
        "n_counterfactual_states": len(states),
        "n_counterfactual_actions": len(intervention_tasks),
        "n_train_states": n_train_states,
        "n_informative_train_states": n_informative_train_states,
        "n_informative_counterfactual_states": sum(gap > 0 for gap in action_gaps),
        "informative_counterfactual_state_fraction": float(np.mean(
            np.asarray(action_gaps) > 0)),
        "mean_action_gap_all_states": float(np.mean(action_gaps)),
        "mean_action_gap_informative_states": float(np.mean(
            [gap for gap in action_gaps if gap > 0])),
        "paired_decision_result": paired,
        "paired_policy_result": policy_paired,
        "intervention_vsearch_command": intervention_command,
        "policy_vsearch_command": policy_command,
    }
    with open(os.path.join(args.out_dir, "run.json"), "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    report = [
        "# A1 counterfactual learning-to-rank pilot", "",
        "Counterfactual labels are final primary-contig identities after a",
        "single forced tie action followed by the fixed index-order policy.",
        "Train and held-out barcodes are disjoint. This is an offline benchmark,",
        "not a production filter.", "", "## Decision-state evaluation", "",
        markdown_table(decision_rows, [
            "method", "n_states", "mean_counterfactual_identity",
            "below_90_pct", "mean_regret_to_oracle", "oracle_choice_rate"]),
        "", "- LTR vs identity wins/losses/ties: {}/{}/{}".format(
            paired["ltr_vs_identity_wins"], paired["ltr_vs_identity_losses"],
            paired["ltr_vs_identity_ties"]),
        "- LTR vs identity mean counterfactual delta: {:.4f}; sign p={:.4g}; "
        "Wilcoxon p={:.4g}".format(
            paired["ltr_vs_identity_mean_delta"],
            paired["ltr_vs_identity_sign_p"],
            paired["ltr_vs_identity_wilcoxon_p"]),
        "", "## Full held-out policy evaluation", "",
        markdown_table(policy_rows, [
            "method", "n_barcodes", "contig_hit_rate", "mean_identity",
            "below_90_pct", "mean_contig_length", "mean_reads_dropped"]),
        "", "- LTR vs identity wins/losses/ties: {}/{}/{}; mean delta {:.4f}; "
        "sign p={:.4g}; Wilcoxon p={:.4g}".format(
            policy_paired["ltr_vs_identity"]["wins"],
            policy_paired["ltr_vs_identity"]["losses"],
            policy_paired["ltr_vs_identity"]["ties"],
            policy_paired["ltr_vs_identity"]["mean_delta"],
            policy_paired["ltr_vs_identity"]["sign_p"],
            policy_paired["ltr_vs_identity"]["wilcoxon_p"]),
        "- LTR vs index wins/losses/ties: {}/{}/{}; mean delta {:.4f}; "
        "sign p={:.4g}; Wilcoxon p={:.4g}".format(
            policy_paired["ltr_vs_index"]["wins"],
            policy_paired["ltr_vs_index"]["losses"],
            policy_paired["ltr_vs_index"]["ties"],
            policy_paired["ltr_vs_index"]["mean_delta"],
            policy_paired["ltr_vs_index"]["sign_p"],
            policy_paired["ltr_vs_index"]["wilcoxon_p"]),
        "", "## Go/no-go", "",
        ("- **GO:** the counterfactual ranker significantly improves the current "
         "identity tie-break on held-out full-policy assemblies."
         if (policy_paired["ltr_vs_identity"]["mean_delta"] > 0 and
             policy_paired["ltr_vs_identity"]["wilcoxon_p"] < 0.05)
         else "- **NO-GO:** the counterfactual ranker does not significantly "
              "improve the current identity tie-break; keep the simpler incumbent."),
        "- This decision concerns replacing A1 only. Improvement over arbitrary",
        "  index-order tie resolution is not enough when the incumbent A1 is the",
        "  actual production comparator.",
        "", "## Top ranker features", "",
        markdown_table(importance_rows[:10], ["feature", "importance_gain"]), "",
    ]
    with open(os.path.join(args.out_dir, "report.md"), "w") as handle:
        handle.write("\n".join(report))

    print(json.dumps(metadata, indent=2, sort_keys=True))
    print("written={}".format(args.out_dir))


if __name__ == "__main__":
    main()
