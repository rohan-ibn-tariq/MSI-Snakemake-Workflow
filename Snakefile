rule all:
    input:
        "raw_data/hello.txt"

rule hello:
    output:
        "raw_data/hello.txt"
    shell:
        "echo 'Hello, world' > {output}"
