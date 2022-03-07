# Arriba Nexflow Workflow

Detect gene fusions from RNA sequencing data using Arriba.

## Quickstart

You need to have a working Conda installation. The first step is to download/generate the reference files. Note that this step requires 50 GB of RAM and 8 CPU threads.

```bash
data/download_references.sh
```

You can then run the workflow using the following command:

```bash
mkdir test_out/
nextflow run main.nf \
    -profile local,conda \
    -ansi-log \
    --fastq1=/path/to/reads_1.fq.gz \
    --fastq2=/path/to/reads_2.fq.gz \
    --outputDir=test_out
```

## Parameters

### Required parameters

- `--fastq1` / `--fastq2` / `--bam`: The sequencing data can be provided in gzipped FastQ format or BAM format. Paired-end FastQ files are supplied using the parameters `--fastq1` and `--fastq2`. A single-end FastQ file is supplied using just the parameter `--fastq1`; the parameter `--fastq2` is simply omitted. Multiple FastQ files can be specified as a comma-separated list. In this case, the FastQ files will be merged and analyzed as a single sample. Alternatively, the sequencing data can be supplied as a BAM file using the parameter `--bam`. The BAM file will be realigned. In contrast to FastQ files, only a single BAM file is supported.
- `--outputDir`: Output directory.

### Optional parameters

- `--preset`: A predefined combination of assembly + annotation to be used for alignment/annotation. The workflow automatically selects the appropriate reference files from the `data` directory. Available presets are: `hs37d5+GENCODE19`, `GRCh38+GENCODE28`, `GRCm38+GENCODEM25`. Alternatively, all reference files can be specified explicitly using the following parameters.
- `--starIndex`: Path to the STAR index directory.
- `--assembly`: Path to the reference genome assembly in FastA format. Arriba supports assemblies with coordinates that are compatible with `hg19/GRCh37/hs37d5`, `hg38/GRCh38`, and `mm10/GRCm38`.
- `--annotation`: Path to the gene model in GTF format. Any gene model can be used that is compatible with the assembly coordinates.
- `--blacklist`: Path to the Arriba blacklist.
- `--knownFusions`: Path to the Arriba known fusions file.
- `--proteinDomains`: Path to the Arriba protein domains file.
- `--cytobands`: Path to the Arriba cytobands file.
- `--sophiaSVs`: A file containing structural variants called by Sophia. The SVs are only used for the sake of annotating fusions with matched structural variants. They are not used for filtering.
- `--threads`: Number of threads to use for alignment and sorting.
- `--arribaMemory`: Upper limit of memory (in MB) consumed by Arriba.
- `--starMemory`: Memory (in MB) consumed by STAR index.
- `--sortMemory`: Memory (in MB) used for sorting of alignments. It should be the same as the STAR footprint, because sorting starts, when STAR finishes.
- `--arribaVersion`: The version of Arriba to use.

### Output

The following output files are generated by the workflow:

- `fusions.tsv`: List of fusions predicted by Arriba. The format is described in the [user manual](https://arriba.readthedocs.io/en/v2.1.0/output-files/#fusionstsv).
- `fusions.discarded.tsv.gz`: Compressed list of fusion candidates dismissed as artifacts by Arriba. The format is described in the [user manual](https://arriba.readthedocs.io/en/v2.1.0/output-files/#fusionsdiscardedtsv).
- `fusions.pdf`: Visualizations of fusion predictions.
- `virus_expression.tsv`: Quantification of viral genome expression. All RefSeq viral genomes are considered and any genome with at least one mapped read is reported. Only reads which map fully to a viral genome are counted. Partially mapping reads are ignored to reduce the impact of alignment artifacts. Other than that, the list is not filtered and may contain many irrelevant/false positive items.

## More Examples

### Run with Docker

A Dockerfile is provided. You will first have to build the container with

```bash
cd nf-arriba
docker build \
    --rm \
    --build-arg http_proxy=$HTTP_PROXY \
    --build-arg https_proxy=$HTTPS_PROXY \
    -t \
    nf-arriba \
    ./
```

Then to run the workflow locally with Docker you can do e.g.

```bash
nextflow run main.nf \
    -profile local,docker \
    -ansi-log \
    --fastq1=/path/to/reads_1.fq.gz \
    --fastq2=/path/to/reads_2.fq.gz \
    --outputDir=test_out
```

### Run with Singularity

To run the workflow with singularity, convert the previously build Docker container to Singularity (no native Singularity container, yet):

```bash
# Convert the Docker image to Singularity.
# Note that the image is stored in the current directory where it is then also expected by the singularity profile.
# Replace X.X.X with the respective workflow version.
singularity build nf-arriba_X.X.X.sif docker-daemon://nf-arriba:X.X.X

# Run with the singularity profile
nextflow run main.nf \
    -profile local,singularity \
    -ansi-log \
    --fastq1=/path/to/reads_1.fq.gz \
    --fastq2=/path/to/reads_2.fq.gz \
    --outputDir=test_out
```

## Environment and Execution

[Nextflow](https://www.nextflow.io/docs/latest/config.html#config-profiles)'s `-profile` parameter allows setting technical options for executing the workflow. You have already seen some of the profiles and that these can be combined. We conceptually separated the predefined profiles into two types, those concerning the "environment" and those for selecting the "executor".

The following "environment" profiles that define which environment will be used for executing the jobs are predefined in the `nextflow.config`:
  * `conda`
  * `docker`
  * `singularity`
  * `dkfzModules`: This environment uses the environment modules available in the DKFZ Cluster.

Currently, there are only two "executor" profiles that define the job execution method. These are
  * `local`: Just execute the jobs locally on the system that executes Nextflow.
  * `lsf`: Submit the jobs to an LSF cluster. Nextflow must be running on a cluster node on which `bsub` is available.
  * `dkfzLsf`: Submit jobs to the DKFZ Cluster, making use of locally attached storage as scratch space.

Here another example, if you want to run the workflow as Singularity containers in an LSF cluster:

```bash
nextflow run main.nf \
    -profile lsf,singularity \
    -ansi-log \
    --fastq1=/path/to/reads_1.fq.gz \
    --fastq2=/path/to/reads_2.fq.gz \
    --outputDir=test_out
```

Please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html) for defining other executors. Note that environments and executors cannot arbitrarily be combined. For instance, your LSF administrators may not allow Docker to be executed by normal users.

## Development

The integration tests can be run with

```bash
test/test1.sh test_out/
```

The integration tests are also run in Travis CI.

## Release Notes

- 0.1.0 (March 5, 2021)

  * Initial release based on Arriba 2.1.0 and STAR 2.7.8a.

## License & Contributors

See [LICENSE](LICENSE) and [CONTRIBUTORS](CONTRIBUTORS).