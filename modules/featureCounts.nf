process featureCounts_mono{
        label 'big'
         
        publishDir "${params.outDir}/RUN/04_NUCS_READ_COUNTS/${sampleID}", mode: 'copy'


        input:
        file(saf)
        file(monoNucs)
        path genome

        output:
        tuple val(sampleID), file("${sampleID}_monoNucs_readCounts_wGC.csv")
        tuple val(sampleID), file("${sampleID}_monoNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o monoNucs_readCounts.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B --largestOverlap \
        $monoNucs



        #get GC content of called nucleosomes
        awk 'NR > 2 {print \$2"\\t"\$3"\\t"\$4}' ${sampleID}_monoNucs_readCounts.csv > tmp.bed
        bedtools nuc -fi ${genome} -bed tmp.bed > pre_${sampleID}_monoNucs_readCounts_wGC.csv
  

        #add nucID
        awk 'BEGIN {FS=OFS="\\t"} NR==FNR {if (FNR <= 2) next; key=\$2 FS \$3 FS \$4; geneid[key]=\$1; next} \\
        {key=\$1 FS \$2 FS \$3; if (key in geneid) print \$0, geneid[key]; else print \$0, "NA"}' \\
        ${sampleID}_monoNucs_readCounts.csv pre_${sampleID}_monoNucs_readCounts_wGC.csv > pre_${sampleID}_monoNucs_readCounts_wGC_ext.csv

        #prepare header
        numb=\$(cat ${sampleID}_monoNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        echo "GC_cont\\tnucID" | paste <(awk 'FNR == 2 {print}' ${sampleID}_monoNucs_readCounts.csv | cut -f 2-\$numb) - > header_all.csv
        cat header_all.csv pre_${sampleID}_monoNucs_readCounts_wGC_ext.csv > ${sampleID}_monoNucs_readCounts_wGC.csv

        #add header
        mv monoNucs_readCounts.csv.summary ${sampleID}_monoNucs_readCounts.csv.summary

        """
}


process featureCounts_sub{
        label 'big'
        
        publishDir "${params.outDir}/RUN/04_NUCS_READ_COUNTS/${sampleID}", mode: 'copy'


        input:
        file(saf)
        file(subNucs)
        path genome

        output:
        tuple val(sampleID), file("${sampleID}_subNucs_readCounts_wGC.csv")
        tuple val(sampleID), file("${sampleID}_subNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o subNucs_readCounts.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B --largestOverlap \
        $subNucs


        #get GC content of called nucleosomes
        awk 'NR > 2 {print \$2"\\t"\$3"\\t"\$4}' ${sampleID}_subNucs_readCounts.csv > tmp.bed
        bedtools nuc -fi ${genome} -bed tmp.bed > pre_${sampleID}_subNucs_readCounts_wGC.csv
  

        #add nucID
        awk 'BEGIN {FS=OFS="\\t"} NR==FNR {if (FNR <= 2) next; key=\$2 FS \$3 FS \$4; geneid[key]=\$1; next} \\
        {key=\$1 FS \$2 FS \$3; if (key in geneid) print \$0, geneid[key]; else print \$0, "NA"}' \\
        ${sampleID}_subNucs_readCounts.csv pre_${sampleID}_subNucs_readCounts_wGC.csv > pre_${sampleID}_subNucs_readCounts_wGC_ext.csv

        #prepare header
        numb=\$(cat ${sampleID}_subNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        echo "GC_cont\\tnucID" | paste <(awk 'FNR == 2 {print}' ${sampleID}_subNucs_readCounts.csv | cut -f 2-\$numb) - > header_all.csv
        cat header_all.csv pre_${sampleID}_subNucs_readCounts_wGC_ext.csv > ${sampleID}_subNucs_readCounts_wGC.csv

        #add header
        mv subNucs_readCounts.csv.summary ${sampleID}_subNucs_readCounts.csv.summary

        """
}
