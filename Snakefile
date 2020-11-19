
configfile: "config.yaml"

(IDS,PAIR) = glob_wildcards("data/fastq/{ip}_{pair}.fastq")

rule all:
     input:
        html= expand("result/fastqc/{ip}_{pair}.html", ip = IDS, pair = PAIR),
        zip = expand("result/fastqc/{ip}_{pair}_fastqc.zip", ip = IDS, pair = PAIR),
        idx1="result/index/reference.1.bt2",
        idx2="result/index/reference.2.bt2",
        idx3="result/index/reference.3.bt2",
        idx4="result/index/reference.4.bt2",
        idxrev1="result/index/reference.rev.1.bt2",
        idxrev2="result/index/reference.rev.2.bt2",
        quast_out = expand("result/quast/{ip}_quast/report.html",ip = IDS),
        metaphlan_out= expand("result/metaphlan2/{ip}_profiled.txt", ip=IDS),
        amr= expand("result/abricate/{ip}_amr_summary.txt",ip=IDS),
        plasmid= expand("result/abricate/{ip}_plasmid_summary.txt",ip=IDS),
        vf= expand("result/abricate/{ip}_vf_summary.txt",ip=IDS)

rule fastqc:
    input:
        "data/fastq/{ip}_{pair}.fastq"
    output:
        html = "result/fastqc/{ip}_{pair}.html",
        zip = "result/fastqc/{ip}_{pair}_fastqc.zip"# the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    conda:
        "env/fastqc.yaml"
    log:
        "logs/fastqc/{ip}_{pair}.log"
    priority: 50
    wrapper:
        "0.38.0/bio/fastqc"

rule trimmomatic_pe:
    input:
        r1="data/fastq/{ip}_1.fastq",
        r2="data/fastq/{ip}_2.fastq"
    output:
        r1_paired=temp("result/trimmed/{ip}_1_paired.fastq"),
        r2_paired=temp("result/trimmed/{ip}_2_paired.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired=temp("result/trimmed_unpaired/{ip}_1.unpaired.fastq"),
        r2_unpaired=temp("result/trimmed_unpaired/{ip}_2.unpaired.fastq")

    conda:
        "env/trimmomatic.yaml"

    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} SLIDINGWINDOW:5:20 MINLEN:55 LEADING:3 TRAILING:3"

rule bowtie2Build:
    input:
        config["fasta_path"]
    params:
        basename="result/index/reference"
    output:
        output1="result/index/reference.1.bt2",
        output2="result/index/reference.2.bt2",
        output3="result/index/reference.3.bt2",
        output4="result/index/reference.4.bt2",
        outputrev1="result/index/reference.rev.1.bt2",
        outputrev2="result/index/reference.rev.2.bt2"

    conda:
        "env/bowtie2.yaml"

    shell:
        "bowtie2-build {input} {params.basename}"

rule bowtie2Aln:
    input:
        sample = ["result/trimmed/{ip}_1_paired.fastq", "result/trimmed/{ip}_2_paired.fastq"],
        dummy_input = "result/index/reference.1.bt2"
    params:
        index = "result/index/reference",
        output_path = "result/mapping/{ip}_unmapped/unmapped.fastq"   # optional parameters
    output:
       temp("result/mapping/{ip}_unmapped/unmapped.1.fastq"),
       temp("result/mapping/{ip}_unmapped/unmapped.2.fastq")

    conda:
        "env/bowtie2.yaml"

    shell:
        "bowtie2 -x {params.index} -1 {input.sample[0]} -2 {input.sample[1]} --un-conc {params.output_path}"

rule assemble_spades_correction:
    input:
        forward = "result/mapping/{ip}_unmapped/unmapped.1.fastq",
        reverse = "result/mapping/{ip}_unmapped/unmapped.2.fastq"

    params:
        output_path = "result/assembly/spades/{ip}_spades_corrected"

    output: temp("result/assembly/spades/{ip}_spades_corrected/corrected/unmapped.1.00.0_0.cor.fastq.gz"),
            temp("result/assembly/spades/{ip}_spades_corrected/corrected/unmapped.2.00.0_0.cor.fastq.gz")

    conda:
        "env/spades.yaml"

    shell:
        "spades.py --meta --pe1-1 {input.forward} --pe1-2 {input.reverse} -o {params.output_path} -t 4 -m 10 --only-error-correction"

rule assemble_spades_assembly:
    input:
        forward = "result/assembly/spades/{ip}_spades_corrected/corrected/unmapped.1.00.0_0.cor.fastq.gz",
        reverse = "result/assembly/spades/{ip}_spades_corrected/corrected/unmapped.2.00.0_0.cor.fastq.gz"

    params:
        output_path = "result/assembly/spades/{ip}_spades_corrected/contigs"

    output: "result/assembly/spades/{ip}_spades_corrected/contigs/contigs.fasta"

    conda:
        "env/spades.yaml"

    shell:
        "spades.py --meta --pe1-1 {input.forward} --pe1-2 {input.reverse} -o {params.output_path} --only-assembler"

rule quast:
    input:
        quast1 = "result/assembly/spades/{ip}_spades_corrected/contigs/contigs.fasta"

    params:
        output_path = "result/quast/{ip}_quast"

    output: "result/quast/{ip}_quast/report.html"

    conda:
        "env/quast.yaml"

    shell:
        "metaquast.py {input.quast1} -o {params.output_path}"
#
rule metaphlan2_run:
   input:
      pe1= "result/trimmed/{ip}_1_paired.fastq",
      pe2= "result/trimmed/{ip}_2_paired.fastq"

   output:
      "result/metaphlan2/{ip}_profiled.txt"

   conda:
      "env/metaphlan2.yaml"

   shell:
      "metaphlan2.py {input.pe1} --input_type fastq > {output}"

rule merge_abundance:
    input: "result/metaphlan2/*_profiled.txt"
    output: "result/metaphlan2/merged_abundance_table.txt"
    shell: "scripts/python merge_metaphlan_tables.py {input} > {output}"

rule species_only:
    input: "result/metaphlan2/merged_abundance_table.txt"
    output: "result/metaphlan2/merged_abundance_table_species.txt"
    conda: "env/metaphlan2.yaml"
    shell: "./scripts/grep.sh {input} {output}"

rule heatmap:
    input: "result/metaphlan2/merged_abundance_table_species.txt"
    output: "result/metaphlan2/abundance_heatmap_species.png"
    conda: "env/metaphlan2.yaml"
    shell: "./scripts/hclust2.py -i {input} -o {output} --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300"

rule abricate_amr:
    input: "result/assembly/spades/{ip}_spades_corrected/contigs/contigs.fasta"

    output: "result/abricate/{ip}_abricate_amr.txt"
    log:
        "logs/abricate/{ip}.log"
    conda:
        "env/abricate.yaml"
    shell:
        "abricate --threads 2 --mincov 60 --db ncbi {input} > {output} 2> {log}"

rule abricate_plasmid:
    input: "result/assembly/spades/{ip}_spades_corrected/contigs/contigs.fasta"
    output:
        "result/abricate/{ip}_abricate_plasmid.txt"
    log:
        'logs/abricate/{ip}.log'
    conda:
        'env/abricate.yaml'
    shell:
        "abricate --threads 2 --mincov 60 --db plasmidfinder {input} > {output} 2> {log}"

rule abricate_vf:
    input: "result/assembly/spades/{ip}_spades_corrected/contigs/contigs.fasta"
    output:
        "result/abricate/{ip}_abricate_vf.txt"
    log:
        'logs/abricate_vf/{ip}.log'
    conda:
        'env/abricate.yaml'
    shell:
        "abricate --threads 2 --mincov 60 --db vfdb {input} > {output} 2> {log}"

rule abricate_amr_summary:
    input:
        "result/abricate/{ip}_abricate_amr.txt"
    output:
        "result/abricate/{ip}_amr_summary.txt"
    conda:
        "env/abricate.yaml"
    shell:
        "abricate --summary {input} > {output}"

rule abricate_vf_summary:
    input:
        "result/abricate/{ip}_abricate_vf.txt"
    output:
        "result/abricate/{ip}_vf_summary.txt"
    conda:
        "env/abricate.yaml"
    shell:
        "abricate --summary {input} > {output}"

rule abricate_plasmid_summary:
    input:
        "result/abricate/{ip}_abricate_plasmid.txt"
    output:
        "result/abricate/{ip}_plasmid_summary.txt"
    conda:
        "env/abricate.yaml"
    shell:
        "abricate --summary {input} > {output}"
