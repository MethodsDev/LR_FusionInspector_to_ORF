version 1.0

workflow LR_FusionInspector_to_ORF_wf {
    input {
        File igv_inputs_tar_gz
        File fusions_tsv

        Boolean HiFi = false
        Int? num_parallel_contigs
        Int blast_threads = 1
        String output_dir = "LR_FusionInspector_to_ORF_outdir"

        String docker = "lr_fusioninspector_to_orf:latest"
    }

    call LR_FusionInspector_to_ORF_task {
        input:
            igv_inputs_tar_gz = igv_inputs_tar_gz,
            fusions_tsv = fusions_tsv,
            HiFi = HiFi,
            num_parallel_contigs = num_parallel_contigs,
            blast_threads = blast_threads,
            output_dir = output_dir,
            docker = docker
    }

    output {
        File lraa_gtf = LR_FusionInspector_to_ORF_task.lraa_gtf
        File reference_peptides = LR_FusionInspector_to_ORF_task.reference_peptides
        File transdecoder_genome_gff3 = LR_FusionInspector_to_ORF_task.transdecoder_genome_gff3
        File transdecoder_pep = LR_FusionInspector_to_ORF_task.transdecoder_pep
        File transdecoder_cds = LR_FusionInspector_to_ORF_task.transdecoder_cds
        File transdecoder_bed = LR_FusionInspector_to_ORF_task.transdecoder_bed
        File fusion_orf_breakpoint_annotations = LR_FusionInspector_to_ORF_task.fusion_orf_breakpoint_annotations
        File output_dir_tar_gz = LR_FusionInspector_to_ORF_task.output_dir_tar_gz
    }
}

task LR_FusionInspector_to_ORF_task {
    input {
        File igv_inputs_tar_gz
        File fusions_tsv

        Boolean HiFi
        Int? num_parallel_contigs
        Int blast_threads
        String output_dir

        String docker
    }

    command <<<
        set -euo pipefail

        mkdir -p igv_inputs
        tar xzf ~{igv_inputs_tar_gz} -C igv_inputs

        genome=$(find igv_inputs -type f -name 'finspector.fa' | head -n 1)
        bam=$(find igv_inputs -type f -name 'finspector.junction_reads.bam' | head -n 1)
        gtf=$(find igv_inputs -type f -name 'finspector.gtf' | head -n 1)

        test -n "${genome}" || { echo "Missing finspector.fa in IGV inputs bundle" >&2; exit 1; }
        test -n "${bam}" || { echo "Missing finspector.junction_reads.bam in IGV inputs bundle" >&2; exit 1; }
        test -n "${gtf}" || { echo "Missing finspector.gtf in IGV inputs bundle" >&2; exit 1; }

        python3 /opt/LR_FusionInspector_to_ORF/LR_FusionInspector_to_ORF.py \
            --genome "${genome}" \
            --bam "${bam}" \
            --gtf "${gtf}" \
            --fusions_tsv ~{fusions_tsv} \
            --output_dir ~{output_dir} \
            --blast_threads ~{blast_threads} \
            ~{if HiFi then "--HiFi" else ""} \
            ~{if defined(num_parallel_contigs) then "--num_parallel_contigs " + num_parallel_contigs else ""}

        tar czf "~{output_dir}.tar.gz" "~{output_dir}"
    >>>

    output {
        File lraa_gtf = "~{output_dir}/LRAA.gtf"
        File reference_peptides = "~{output_dir}/finspector.gtf.pep"
        File transdecoder_genome_gff3 = "~{output_dir}/LRAA.cDNA.fasta.transdecoder.genome.gff3"
        File transdecoder_pep = "~{output_dir}/LRAA.cDNA.fasta.transdecoder.pep"
        File transdecoder_cds = "~{output_dir}/LRAA.cDNA.fasta.transdecoder.cds"
        File transdecoder_bed = "~{output_dir}/LRAA.cDNA.fasta.transdecoder.bed"
        File fusion_orf_breakpoint_annotations = "~{output_dir}/fusion_orf_breakpoint_annotations.tsv"
        File output_dir_tar_gz = "~{output_dir}.tar.gz"
    }

    runtime {
        docker: docker
        cpu: blast_threads + select_first([num_parallel_contigs, 1]) + 1
        memory: "32G"
        disks: "local-disk " + ceil(size(igv_inputs_tar_gz, "GB") * 8 + size(fusions_tsv, "GB") * 2 + 50) + " HDD"
    }
}
