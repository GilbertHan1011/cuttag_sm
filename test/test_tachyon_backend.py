from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text()


def test_configs_declare_opt_in_tachyon_backend():
    for name, backend in (("config/config.yml", "legacy"), ("config/config_test.yml", "legacy")):
        config = yaml.safe_load((ROOT / name).read_text())
        alignment = config["pipeline_steps"]["alignment"]
        assert alignment["backend"] == backend
        assert alignment["tachyon_upstream"]["executable"]
        assert alignment["tachyon_upstream"]["index"]
        assert alignment["tachyon_upstream"]["adapter"] == "illumina"
    overlay = yaml.safe_load((ROOT / "config/config_tachyon_test.yml").read_text())
    assert overlay["pipeline_steps"]["alignment"]["backend"] == "tachyon_upstream"
    fixture = yaml.safe_load((ROOT / "test/cuttag_tachyon_fixture.yml").read_text())
    assert fixture["pipeline_steps"]["alignment"]["backend"] == "tachyon_upstream"


def test_tachyon_producer_uses_raw_pairs_and_canonical_outputs():
    align = read("workflow/rules/align.smk")
    assert "rule align_tachyon_upstream:" in align
    assert "_paired_raw_fastqs" in align
    assert "--assay cut-and-tag" in align
    assert "--ucf" in align
    assert "--stats" in align
    assert "sorted.markd.bam.bai" in align
    tachyon_rule = align.split("rule align_tachyon_upstream:", 1)[1]
    assert "samtools index" not in tachyon_rule


def test_tachyon_backend_skips_legacy_producers_and_uses_native_ucf():
    assert "if ALIGNMENT_BACKEND == \"legacy\":" in read("workflow/rules/align.smk")
    assert "if config[\"alignment_backend\"] == \"legacy\":" in read("workflow/rules/preprocess.smk")
    postalign = read("workflow/rules/postalign.smk")
    assert "if ALIGNMENT_BACKEND == \"legacy\":" in postalign
    assert "canonical_bam" in read("workflow/rules/report.smk")


def test_paired_raw_fastqs_reject_incomplete_pairs():
    common = read("workflow/rules/common.smk")
    assert "def _paired_raw_fastqs" in common
    assert "missing paired FASTQ" in common


def test_peak_rules_use_local_macs2_environment():
    peaks = read("workflow/rules/peaks.smk")
    assert 'conda:\n        "../envs/macs2.yml"' in peaks
    assert "wrapper:" not in peaks
    assert (ROOT / "workflow/envs/macs2.yml").is_file()


def test_local_macs2_environment_pins_compatible_python():
    env = yaml.safe_load((ROOT / "workflow/envs/macs2.yml").read_text())
    assert "python=3.10" in env["dependencies"]
