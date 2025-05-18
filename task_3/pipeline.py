from toil.common import Toil
from toil.job import Job
import subprocess
import os
import re
import sys

# Утилита для запуска shell-команд
def run_cmd(job, cmd, description):
    job.fileStore.logToMaster(f"==> {description}\n{cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        job.fileStore.logToMaster(f"[ERROR] {description} failed:\n{result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    job.fileStore.logToMaster(f"[OK] {description} succeeded.")
    return result.stdout

# FastQC
def fastqc(job, fastq_id, toil: Toil):
    fixed_fastq_path = os.path.join(job.tempDir, "sample.fastq.gz")
    job.fileStore.readGlobalFile(fastq_id, userPath=fixed_fastq_path)
    fastqc_report_path = os.path.join(job.tempDir, "sample_fastqc.html")
    run_cmd(job, f"fastqc {fixed_fastq_path}", f"Running FastQC on {fixed_fastq_path}")
    toil.export_file(
        job.fileStore.writeGlobalFile(fastqc_report_path), 
        "fastqc_report.html"
    )


# Индексация референса
def index_reference(job, ref_id):
    ref_path = job.fileStore.readGlobalFile(ref_id)
    index_path = os.path.join(job.tempDir, "ref.mmi")
    run_cmd(job, f"minimap2 -d {index_path} {ref_path}", f"Indexing reference")
    return job.fileStore.writeGlobalFile(index_path)

# Картирование
def align_reads(job, fastq_id, ref_index_id):
    fastq_path = job.fileStore.readGlobalFile(fastq_id)
    index_path = job.fileStore.readGlobalFile(ref_index_id)
    sam_path = os.path.join(job.tempDir, "aligned.sam")
    run_cmd(job, f"minimap2 -a {index_path} {fastq_path} > {sam_path}", f"Aligning reads")
    return job.fileStore.writeGlobalFile(sam_path)

# Конвертация в BAM
def sam_to_bam(job, sam_id):
    sam_path = job.fileStore.readGlobalFile(sam_id)
    bam_path = os.path.join(job.tempDir, "aligned.bam")
    run_cmd(job, f"samtools view -b {sam_path} > {bam_path}", f"Converting SAM to BAM")
    return job.fileStore.writeGlobalFile(bam_path)

# Статистика картирования
def flagstat(job, bam_id):
    bam_path = job.fileStore.readGlobalFile(bam_id)
    report_path = os.path.join(job.tempDir, "alignment_report.txt")
    run_cmd(job, f"samtools flagstat {bam_path} > {report_path}", "Generating alignment stats")
    return job.fileStore.writeGlobalFile(report_path)

# Оценка результата
def evaluate(job, report_id, toil: Toil):
    toil.export_file(report_id, "alignment_report.txt")

    report_path = job.fileStore.readGlobalFile(report_id)
    with open(report_path) as f:
        content = f.read()
    match = re.search(r"(\d+\.\d+)% mapped", content)
    if match:
        percent = float(match.group(1))
        job.fileStore.logToMaster(f"Mapped: {percent}%")
        if percent > 90:
            job.fileStore.logToMaster("Result: OK")
        else:
            job.fileStore.logToMaster("Result: Not OK")
    else:
        raise RuntimeError("Could not parse alignment percentage.")

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("fastq", help=".fastq.gz sequence file")
    parser.add_argument("reference", help=".fa reference file")

    options = parser.parse_args()
    options.clean = "always"
    options.logLevel = "INFO"
    options.writeGraph = "workflow_graph.dot"

    fastq_file = os.path.abspath(options.fastq)
    ref_file = os.path.abspath(options.reference)

    with Toil(options) as toil:
        root = Job()

        fastq_id = toil.importFile("file://" + fastq_file)
        ref_id = toil.importFile("file://" + ref_file)

        fastqc_job = root.addChild(Job.wrapJobFn(fastqc, fastq_id, toil))
        
        # index_job = root.addChild(Job.wrapJobFn(index_reference, ref_id))
        # align_job = index_job.addFollowOn(Job.wrapJobFn(align_reads, fastq_id, index_job.rv()))
        # bam_job = align_job.addFollowOn(Job.wrapJobFn(sam_to_bam, align_job.rv()))
        # flagstat_job = bam_job.addFollowOn(Job.wrapJobFn(flagstat, bam_job.rv(), toil))
        # flagstat_job.addFollowOn(Job.wrapJobFn(evaluate, flagstat_job.rv()))

        toil.start(root)

    try:
        subprocess.run("dot -Tpng workflow_graph.dot -o pipeline_graph.png", shell=True, check=True)
        print("Pipeline graph saved as pipeline_graph.png")
    except Exception as e:
        print("Failed to generate pipeline graph:", e)
