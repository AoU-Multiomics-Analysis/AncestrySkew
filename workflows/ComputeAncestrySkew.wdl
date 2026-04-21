version 1.0

workflow ComputeAncestrySkew {
    input {
        File AnnotationData
        String OutputFile
    }
    
    call GetAncestrySkew {
        input:
            AnnotationData = AnnotationData,
            OutputFile = OutputFile
    }

    output {
        File Output = GetAncestrySkew.AncestrySkewOutput
    }
}



task GetAncestrySkew {
    input {
        File AnnotationData
        String OutputFile
    }
    command <<< 
    Rscript /tmp/ComputeAncestrySkew.R \
       --AnnotationData ~{AnnotationData} \
       --OuputFile ~{OutputFile}
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/ComputeAncestrySkew:main"
        memory: "96G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    
    output {
        File AncestrySkewOutput = "~{OutputFile}"
    }
}
