import gzip
import os
import shutil
import stat
import subprocess
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
DEFAULT_UCF_CLI = Path("/home/gilberthan/disk1/projects/UCF/target/release/ucf-cli")


def _require_executable(path_or_name: str) -> str:
    path = Path(path_or_name)
    if path.is_absolute():
        assert path.exists() and os.access(path, os.X_OK), f"missing executable: {path}"
        return str(path)
    resolved = shutil.which(path_or_name)
    assert resolved, f"missing executable on PATH: {path_or_name}"
    return resolved


def _ucf_cli() -> str:
    return _require_executable(os.environ.get("CUTTAG_SMOKE_UCF_CLI", str(DEFAULT_UCF_CLI)))


def _run(command: list[str], *, cwd: Path = REPO_ROOT, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert completed.returncode == 0, (
        "command failed:\n"
        + " ".join(command)
        + "\n\noutput:\n"
        + completed.stdout
    )
    return completed


def _write_smoke_config(root: Path, sample_sheet: Path, ucf_cli: str) -> Path:
    db = root / "db"
    db.mkdir(parents=True)
    (db / "chrSmoke.chrom.sizes").write_text("chrSmoke\t1000\n")
    (db / "adapters.fa").write_text(">adapter\nAGATCGGAAGAGC\n")
    (db / "genes.bed").write_text("chrSmoke\t1\t100\tSMOKE\t0\t+\n")
    (db / "ref.fa").write_text(">chrSmoke\n" + "A" * 1000 + "\n")

    pep = root / "pep.yaml"
    pep.write_text(
        textwrap.dedent(
            f"""
            pep_version: "2.1.0"
            project:
              name: smoke
              raw_data: {root / "raw"}
            paths:
              process_dir: {root / "process"}
              database_dir: {db}
            """
        ).strip()
        + "\n"
    )

    config = root / "config.yml"
    config.write_text(
        textwrap.dedent(
            f"""
            project_info:
              name: smoke
              description: smoke test fixture
              genome: hg38
              samples:
                table: "{sample_sheet}"
                pep_config: "{pep}"
                pep_attributes: {{}}
              data:
                raw_data_root: "{root / "raw"}"
              outputs:
                root: "{root / "process"}/smoke"
                structure: {{}}
            workflow_resources:
              compute:
                threads_default: 1
                memory_gb_default: 1
              profiles:
                cluster_config: ""
            pipeline_steps:
              alignment:
                aligner: bowtie2
              peaks:
                use_igg: false
                igg_token: IgG
                n_intersects: 1
              peak_qval:
                H3K27ac: 0.001
              broad_marks: []
              broad_params: {{}}
              preprocessing:
                trim_adapters: false
              replicates:
                pval_thresh: 1e-3
                smooth_win: 150
                cap_num_peak: 300000
                pseudoreplication_random_seed: 0
                disable_tn5_shift: false
                idr_thresh: 0.05
                idr_rank: signal.value
                peak_combination_method: idr
                nonamecheck: false
            reference:
              genome_index: "{db / "smoke_index"}"
              genome_fasta: "{db / "ref.fa"}"
              genes_bed: "{db / "genes.bed"}"
              chrom_sizes: "{db / "chrSmoke.chrom.sizes"}"
              adapter_fasta: "{db / "adapters.fa"}"
              blacklist: ""
              genome_size_bp: hs
            UCF_CLI: "{ucf_cli}"
            """
        ).strip()
        + "\n"
    )
    return config


def _write_sample_sheet(path: Path, r1: Path, r2: Path) -> None:
    path.write_text(
        "sample,run,R1,R2,mark,condition,igg,gopeaks,replicate_sample_name\n"
        f"smoke,1,{r1},{r2},H3K27ac,smoke,smoke_IgG,-,smoke\n"
    )


def _write_tiny_bam(path: Path) -> None:
    samtools = _require_executable("samtools")
    path.parent.mkdir(parents=True, exist_ok=True)
    sam = path.with_suffix(".sam")
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chrSmoke\tLN:1000\n"
        f"read1\t99\tchrSmoke\t101\t60\t50M\t=\t151\t100\t{'A' * 50}\t{'I' * 50}\n"
        f"read1\t147\tchrSmoke\t151\t60\t50M\t=\t101\t-100\t{'T' * 50}\t{'I' * 50}\n"
    )
    _run([samtools, "view", "-bS", "-o", str(path), str(sam)])
    _run([samtools, "quickcheck", str(path)])


def _run_snakemake(config: Path, target: Path, *, env: dict[str, str] | None = None) -> str:
    snakemake = _require_executable("snakemake")
    completed = _run(
        [
            snakemake,
            "--workflow-profile",
            "none",
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(config),
            "-j",
            "1",
            str(target),
            "--printshellcmds",
            "--rerun-incomplete",
            "--nolock",
        ],
        env=env,
    )
    return completed.stdout


def _assert_valid_ucf(ucf_cli: str, path: Path) -> None:
    assert path.exists(), f"missing UCF output: {path}"
    assert path.stat().st_size > 0, f"empty UCF output: {path}"
    _run([ucf_cli, "validate", str(path)])
    inspect = _run([ucf_cli, "inspect", str(path)]).stdout
    assert "family=event_1d" in inspect
    assert "schema=atac_event_v1" in inspect
    assert "row_count=1" in inspect


def test_bam_to_ucf_smoke_runs_real_rule_on_minimal_bam(tmp_path: Path):
    root = tmp_path / "direct"
    raw = root / "raw"
    raw.mkdir(parents=True)
    r1 = raw / "smoke_R1.fastq.gz"
    r2 = raw / "smoke_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    sample_sheet = root / "samples.csv"
    _write_sample_sheet(sample_sheet, r1, r2)
    ucf_cli = _ucf_cli()
    config = _write_smoke_config(root, sample_sheet, ucf_cli)

    bam = root / "process" / "smoke" / "Important_processed" / "Bam" / "smoke.sorted.markd.bam"
    _write_tiny_bam(bam)
    ucf = root / "process" / "smoke" / "Important_processed" / "ucf" / "smoke.sorted.markd.ucf"

    output = _run_snakemake(config, ucf)

    assert "localrule bam_to_ucf" in output
    assert "--chrom-sizes" in output
    _assert_valid_ucf(ucf_cli, ucf)


def _write_fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("@read1/1\n" + "A" * 50 + "\n+\n" + "I" * 50 + "\n")


def _write_executable(path: Path, content: str) -> None:
    path.write_text(textwrap.dedent(content).strip() + "\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _write_tool_shims(bin_dir: Path) -> None:
    bin_dir.mkdir(parents=True)
    _write_executable(
        bin_dir / "fastp",
        r'''
        #!/usr/bin/env python3
        import pathlib
        import shutil
        import sys

        args = sys.argv[1:]
        def value(flag):
            return pathlib.Path(args[args.index(flag) + 1])

        for source_flag, output_flag in [('-i', '-o'), ('-I', '-O')]:
            source = value(source_flag)
            output = value(output_flag)
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, output)

        json_report = value('-j')
        html_report = value('-h')
        json_report.parent.mkdir(parents=True, exist_ok=True)
        html_report.parent.mkdir(parents=True, exist_ok=True)
        json_report.write_text('{"summary": "smoke"}\n')
        html_report.write_text('<html>smoke</html>\n')
        ''',
    )
    _write_executable(
        bin_dir / "bowtie2",
        r'''
        #!/usr/bin/env python3
        import sys

        sys.stderr.write('bowtie2 smoke shim\n')
        print('@HD\tVN:1.6\tSO:coordinate')
        print('@SQ\tSN:chrSmoke\tLN:1000')
        print('read1\t99\tchrSmoke\t101\t60\t50M\t=\t151\t100\t' + 'A' * 50 + '\t' + 'I' * 50)
        print('read1\t147\tchrSmoke\t151\t60\t50M\t=\t101\t-100\t' + 'T' * 50 + '\t' + 'I' * 50)
        ''',
    )
    _write_executable(
        bin_dir / "sambamba",
        r'''
        #!/usr/bin/env python3
        import pathlib
        import shutil
        import sys

        args = sys.argv[1:]
        assert args[0] == 'markdup', args
        positionals = []
        index = 1
        while index < len(args):
            arg = args[index]
            if arg == '-t':
                index += 2
            elif arg.startswith('--tmpdir='):
                index += 1
            elif arg.startswith('-'):
                index += 1
            else:
                positionals.append(pathlib.Path(arg))
                index += 1
        source, output = positionals[-2:]
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, output)
        ''',
    )


def test_end_to_end_smoke_generates_ucf_from_fastq_with_workflow_chain(tmp_path: Path):
    root = tmp_path / "e2e"
    raw = root / "raw"
    r1 = raw / "smoke_R1.fastq.gz"
    r2 = raw / "smoke_R2.fastq.gz"
    _write_fastq(r1)
    _write_fastq(r2)
    sample_sheet = root / "samples.csv"
    _write_sample_sheet(sample_sheet, r1, r2)
    ucf_cli = _ucf_cli()
    config = _write_smoke_config(root, sample_sheet, ucf_cli)
    shim_bin = root / "bin"
    _write_tool_shims(shim_bin)

    env = os.environ.copy()
    env["PATH"] = f"{shim_bin}:{env['PATH']}"
    ucf = root / "process" / "smoke" / "Important_processed" / "ucf" / "smoke.sorted.markd.ucf"

    output = _run_snakemake(config, ucf, env=env)

    assert "localrule fastp" in output
    assert "localrule align" in output
    assert "localrule markdup" in output
    assert "localrule bam_to_ucf" in output
    assert (root / "process" / "smoke" / "Report" / "fastp" / "smoke_run1_fastp.json").exists()
    assert (root / "process" / "smoke" / "Important_processed" / "Bam" / "smoke.sorted.markd.bam").exists()
    _assert_valid_ucf(ucf_cli, ucf)
