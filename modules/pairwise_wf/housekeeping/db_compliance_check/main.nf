process DB_COMPLIANCE_CHECK {
    tag "$params.runID"
    
    publishDir "${params.outDir}/reports", mode: 'copy', pattern: "*.{txt,log}"
    
    input:
        val(sampleID_list)
    
    output:
        val("PASS")
        path("db_integrity_report.txt"), emit: report
        path("db_compliance_check.log"), emit: log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Setup logging
    exec > >(tee -a db_compliance_check.log)
    exec 2>&1
    
    echo "=== DATABASE COMPLIANCE CHECK ===" | tee db_integrity_report.txt
    echo "Run ID: $params.runID" | tee -a db_integrity_report.txt
    echo "Timestamp: \$(date)" | tee -a db_integrity_report.txt
    echo "Database path: $params.outDir/db/samples/" | tee -a db_integrity_report.txt
    echo "" | tee -a db_integrity_report.txt
    
    # Check if the database samples directory exists
    if [[ ! -d "$params.outDir/db/samples/" ]]; then
        echo "ERROR: Database samples directory not found: $params.outDir/db/samples/" | tee -a db_integrity_report.txt
        exit 1
    fi
    
    # Generate a list of sample IDs from the db samples directory
    ls "$params.outDir/db/samples/" > sample_ids.txt
    
    # Check if any samples exist
    if [[ ! -s sample_ids.txt ]]; then
        echo "ERROR: No samples found in database directory" | tee -a db_integrity_report.txt
        exit 1
    fi
    
    total_samples=\$(wc -l < sample_ids.txt)
    echo "Found \$total_samples samples in database" | tee -a db_integrity_report.txt
    echo "" | tee -a db_integrity_report.txt
    
    # Initialize counters and arrays for detailed reporting
    valid_samples=0
    invalid_samples=0
    missing_subdirs=0
    declare -a invalid_naming=()
    declare -a missing_subdirs_list=()
    
    # Required subdirectories
    required_subdirs=("mtbseq" "tbprofiler" "snippy")
    
    echo "DETAILED SAMPLE ANALYSIS:" | tee -a db_integrity_report.txt
    echo "=========================" | tee -a db_integrity_report.txt
    
    # Check naming convention and directory structure
    while read -r sample; do
        # Remove any trailing slashes just in case
        sample="\${sample%/}"
        
        echo "Checking sample: \$sample"
        
        # Check naming convention (expected format: prefix_suffix)
        if [[ "\$sample" =~ ^[^_]+_[^_]+\$ ]]; then
            echo "  ✓ Naming convention: VALID"
            ((valid_samples++))
            
            # Check if sample directory contains required subdirectories
            sample_path="$params.outDir/db/samples/\$sample"
            
            if [[ ! -d "\$sample_path" ]]; then
                echo "  ✗ Sample directory not found: \$sample_path" | tee -a db_integrity_report.txt
                invalid_naming+=("\$sample")
                ((invalid_samples++))
                continue
            fi
            
            # Check for required subdirectories
            missing_count=0
            missing_list=()
            
            for subdir in "\${required_subdirs[@]}"; do
                if [[ -d "\$sample_path/\$subdir" ]]; then
                    echo "  ✓ Subdirectory found: \$subdir"
                else
                    echo "  ✗ Missing subdirectory: \$subdir"
                    missing_list+=("\$subdir")
                    ((missing_count++))
                fi
            done
            
            if [[ \$missing_count -gt 0 ]]; then
                echo "  ✗ Sample \$sample is missing \$missing_count required subdirectories: \${missing_list[*]}" | tee -a db_integrity_report.txt
                missing_subdirs_list+=("\$sample: \${missing_list[*]}")
                ((missing_subdirs++))
            else
                echo "  ✓ All required subdirectories present"
            fi
            
        else
            echo "  ✗ Naming convention: INVALID (expected format: prefix_suffix)" | tee -a db_integrity_report.txt
            invalid_naming+=("\$sample")
            ((invalid_samples++))
        fi
        
        echo ""
        
    done < sample_ids.txt
    
    # Generate detailed report
    echo "" | tee -a db_integrity_report.txt
    echo "SUMMARY REPORT:" | tee -a db_integrity_report.txt
    echo "===============" | tee -a db_integrity_report.txt
    echo "Total samples checked: \$total_samples" | tee -a db_integrity_report.txt
    echo "Valid samples: \$valid_samples" | tee -a db_integrity_report.txt
    echo "Invalid naming convention: \$invalid_samples" | tee -a db_integrity_report.txt
    echo "Samples with missing subdirectories: \$missing_subdirs" | tee -a db_integrity_report.txt
    echo "" | tee -a db_integrity_report.txt
    
    # List problematic samples
    if [[ \${#invalid_naming[@]} -gt 0 ]]; then
        echo "SAMPLES WITH INVALID NAMING:" | tee -a db_integrity_report.txt
        printf '%s\\n' "\${invalid_naming[@]}" | tee -a db_integrity_report.txt
        echo "" | tee -a db_integrity_report.txt
    fi
    
    if [[ \${#missing_subdirs_list[@]} -gt 0 ]]; then
        echo "SAMPLES WITH MISSING SUBDIRECTORIES:" | tee -a db_integrity_report.txt
        printf '%s\\n' "\${missing_subdirs_list[@]}" | tee -a db_integrity_report.txt
        echo "" | tee -a db_integrity_report.txt
    fi
    
    # Determine overall result
    if [[ \$invalid_samples -eq 0 && \$missing_subdirs -eq 0 ]]; then
        echo "✓ DATABASE INTEGRITY CHECK: PASSED" | tee -a db_integrity_report.txt
        echo "All samples follow the correct naming convention and contain required subdirectories (mtbseq, tbprofiler, snippy)" | tee -a db_integrity_report.txt
    else
        echo "✗ DATABASE INTEGRITY CHECK: FAILED" | tee -a db_integrity_report.txt
        if [[ \$invalid_samples -gt 0 ]]; then
            echo "  - \$invalid_samples samples have invalid naming convention" | tee -a db_integrity_report.txt
        fi
        if [[ \$missing_subdirs -gt 0 ]]; then
            echo "  - \$missing_subdirs samples are missing required subdirectories" | tee -a db_integrity_report.txt
        fi
        echo "" | tee -a db_integrity_report.txt
        echo "Please review the detailed analysis above and fix the identified issues." | tee -a db_integrity_report.txt
        exit 1
    fi
    
    echo "" | tee -a db_integrity_report.txt
    echo "Check completed at: \$(date)" | tee -a db_integrity_report.txt
    """
}