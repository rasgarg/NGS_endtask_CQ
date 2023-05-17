nextflow.enable.dsl = 2

params.accession = null
params.outdir = "SRA_data"
params.kmerlen = 71
params.cov_cutoff = 0

process prefetch{
    storeDir "${params.outdir}"

    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

    input:
    val accession
   
    output:
    path "sra/${accession}.sra" , emit: dotsra

    script:
    """
    prefetch $accession
    """
}

process convertfastq{
    storeDir "${params.outdir}/sra/${accession}/fastq"
    
    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

    input:
    val accession
    path srafile

    output:
    path "*.fastq" , emit: fastqfile_path

    script:
    """
    fastq-dump --split-files $srafile
    """
}

process qualityassess{
    storeDir "${params.outdir}/sra/${params.accession}/fastqc"

    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

    input:
    path fastqfile

    output:
    path "*.html", emit: fastQC_html
    path "*.zip", emit: fastQC_zip

    script:

    """
    fastqc $fastqfile 
    """
}

process trim{
    storeDir "${params.outdir}/sra/${params.accession}/fastp"

    container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"

    input:
    path fastqfile

    output:
    path "*.fastq", emit: fastp_file
    path "*.html", emit: report

    script:
    """
    fastp -i $fastqfile -o after_fastp.fastq -h after_fastp.html
    """
}

process qualityassessqc{
    storeDir "${params.outdir}/sra/${params.accession}/fastqc_after_fastp"

    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

    input:
    path fastq_trim

    output:
    path "*.html", emit: fastQC_html
    path "*.zip", emit: fastQC_zip

    script:
    """
    fastqc $fastq_trim
    """
}

process multiqc{
    storeDir "${params.outdir}/sra/${params.accession}/multiqc"

    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.13a--pyhdfd78af_1"

    input:
    path report_list
   
    output:
    path "*"

    script:
    """
    multiqc .
    """
}

process velveth{
    storeDir "${params.outdir}/sra/${params.accession}/velvet"

    container "https://depot.galaxyproject.org/singularity/velvet%3A1.2.10--h7132678_5"

    input:
    path fastq_trim

    output:
    path "hash"

    script:
    """
    velveth hash ${params.kmerlen} -fastq ${fastq_trim}
    """
}

process velvetg{
    storeDir "${params.outdir}/sra/${params.accession}/velvet"

    container "https://depot.galaxyproject.org/singularity/velvet%3A1.2.10--h7132678_5"

    input:
    path hash_data

    output:
    path "hash"

    script:
    """
    velvetg hash -cov_cutoff ${params.cov_cutoff}
    """
}

workflow{
    sra_file = prefetch(params.accession).dotsra
    fastq_file = convertfastq(params.accession, sra_file).fastqfile_path
    fastqc_file_flat = qualityassess(fastq_file)
    fastq_trim = trim(fastq_file).fastp_file
    fastqc_after_trim_flat = qualityassessqc(fastq_trim)
    reportQC = Channel.empty()
    reportQC_concat = reportQC.concat(fastqc_file_flat, fastqc_after_trim_flat)
    report_list = reportQC_concat.collect()
    multiqc(report_list)
    hash_data = velveth(fastq_trim)
    velvetg(hash_data)
}
