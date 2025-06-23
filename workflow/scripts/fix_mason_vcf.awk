#!/usr/bin/awk -f

# Usage: awk -v format_tag="GT" -v description="Genotype" -f script.awk input.vcf

BEGIN { 
    OFS = "\t"
    
    format_tag = format_tag ? format_tag : "GT"
    description = description ? description : "Genotype"
    
    format_header = "##FORMAT=<ID=" format_tag ",Number=1,Type=String,Description=\"" description "\">"
    
    print "Using FORMAT tag: " format_tag > "/dev/stderr"
}

/^##/ { 
    print $0
    next
}

/^#CHROM/ {
    print format_header
    print $0
    next
}

{
    if ($9 == ".") {
        $9 = format_tag
    }
    print $0
}

END {
    close("/dev/stderr")
}
