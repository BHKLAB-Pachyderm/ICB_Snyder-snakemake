from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Snyder-data/main/"

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript -e \
        '
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Snyder", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule format_snv:
    output:
        S3.remote(prefix + "processed/SNV.csv")
    input:
        S3.remote(prefix + "download/SNV.txt"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_expr:
    output:
        S3.remote(prefix + "processed/EXPR.csv")
    input:
        S3.remote(prefix + "download/EXPR.txt"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_EXPR.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_clin:
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    input:
        S3.remote(prefix + "download/CLIN.txt")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_CLIN.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_cased_sequenced:
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    input:
        S3.remote(prefix + "download/CLIN.txt")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_data:
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt"),
        S3.remote(prefix + "download/SNV.txt")
    resources:
        mem_mb=2000
    shell:
        """
        wget {data_source}CLIN.txt -O {prefix}download/CLIN.txt
        wget {data_source}EXPR.txt -O {prefix}download/EXPR.txt
        wget {data_source}SNV.txt -O {prefix}download/SNV.txt
        """ 