from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
POSTALIGN = REPO_ROOT / "workflow" / "rules" / "postalign.smk"


def _rule_block(text: str, rule_name: str) -> str:
    marker = f"rule {rule_name}:"
    start = text.find(marker)
    assert start != -1, f"missing {marker}"
    next_rule = text.find("\nrule ", start + len(marker))
    if next_rule == -1:
        return text[start:]
    return text[start:next_rule]


def test_rule_all_requests_ucf_outputs_for_all_samples():
    snakefile = SNAKEFILE.read_text()

    assert (
        'f"{DATA_DIR}/Important_processed/ucf/{{sample}}.sorted.markd.ucf"'
        in snakefile
    )


def test_bam_to_ucf_rule_passes_configured_chrom_sizes():
    postalign = POSTALIGN.read_text()
    rule = _rule_block(postalign, "bam_to_ucf")

    assert 'bam = lambda wc: canonical_bam(wc.sample)' in rule
    assert 'chrom_sizes = config["CSIZES"]' in rule
    assert (
        'ucf = f"{DATA_DIR}/Important_processed/ucf/{{sample}}.sorted.markd.ucf"'
        in rule
    )
    assert 'ucf_cli = config.get("UCF_CLI", "ucf-cli")' in rule
    assert "{params.ucf_cli:q} build event-1d" in rule
    assert "--assay cut-and-tag" in rule
    assert "--input {input.bam:q}" in rule
    assert "--output {output.ucf:q}" in rule
    assert "--chrom-sizes {input.chrom_sizes:q}" in rule
    assert "--genome" not in rule
