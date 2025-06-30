#strelka2 and annovar snakefile
#BLCA TAN
#Template made May 01, 2023
#Pushpa Itagi
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Java/11.0.2
ml picard/2.25.1-Java-11
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-11
ml tabix/0.2.6-GCCcore-8.3.0
ml annovar
### temp fix, since above version of python3 does not have pandas installed
ml fhPython/3.8.2-foss-2020a-Python-3.8.2

#ml Python/2.7.15-foss-2018b
#command to run snakemake (remove -np at end when done validating):
snakemake -s strelka2.snakefile  --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/config-strelka2.yaml"
configfile: "config/samples.yaml"


rule all:
	input:
		expand("results/{tumors}/runWorkflow.py", tumors = config["samples"]),
		expand("results/{tumors}/results/variants/somatic.snvs.vcf.gz", tumors = config["samples"]),
		expand("results/{tumors}/results/variants/somatic.snvs.hg38_multianno.vcf", tumors = config["samples"]),
		expand("results/{tumors}/results/variants/somatic.indels.hg38_multianno.vcf", tumors = config["samples"])

rule strelka2_configure:
	input:
		tumor_filepath = lambda wildcards: config["samples"][wildcards.tumors][0],
		normal_filepath = lambda wildcards: config["samples"][wildcards.tumors][2]
	output:
		pyout = protected("results/{tumors}/runWorkflow.py")
	params:
		reference_genome = config["reference_genome"],
		strelka2 = config["strelka2"],
		out_dir = "results/{tumors}/"
	log:
        	"logs/strelka2/{tumors}_process_strelka.txt"
	shell:
        	"({params.strelka2} \
		--normalBam={input.normal_filepath} \
		--tumorBam={input.tumor_filepath} \
		--referenceFasta={params.reference_genome} \
		--runDir={params.out_dir}) 2> {log}"

rule strelka2_submit:
	input:
		workflow = "results/{tumors}/runWorkflow.py"
	output:
		protected("results/{tumors}/results/variants/somatic.snvs.vcf.gz"),
		protected("results/{tumors}/results/variants/somatic.indels.vcf.gz")
	params:
		interpreter = config["interpreter"]
	log:
                "logs/strelka2/{tumors}_submit_strelka.txt"
	shell:
		"({params.interpreter} {input.workflow} -m local -j 5 ) 2> {log}"

rule runAnnovar_strelka2_snvs:
    input:
        input_vcfs = "results/{tumors}/results/variants/somatic.snvs.vcf.gz"
    output:
        output_annovar_vcf = protected("results/{tumors}/results/variants/somatic.snvs.hg38_multianno.vcf")
    params:
        annovar_python_script = ("config/{script}, script = config["annovar_python_script"]),
        interpreter_pd = config["interpreter_pd"]
    log:
        "logs/runAnnovar_strelka2/{tumors}_runAnnovar_snvs_strelka2.txt"
    shell:
        "({params.interpreter_pd} {params.annovar_python_script} --input_vcf_file_path {input.input_vcfs} ) 2> {log}"


rule runAnnovar_strelka2_indels:
    input:
        input_vcfs = "results/{tumors}/results/variants/somatic.indels.vcf.gz"
    output:
        output_annovar_vcf = protected("results/{tumors}/results/variants/somatic.indels.hg38_multianno.vcf")
    params:
        annovar_python_script = ("config/{script}, script = config["annovar_python_script"]),
        interpreter_pd = config["interpreter_pd"]
    log:
        "logs/runAnnovar_strelka2/{tumors}_runAnnovar_indels_strelka2.txt"
    shell:
        "({params.interpreter_pd} {params.annovar_python_script} --input_vcf_file_path {input.input_vcfs} ) 2> {log}"
