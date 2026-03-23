
MINIWDL_RUN_DIR := miniwdl-run

test:
	../LR_FusionInspector_to_ORF.py --genome FI/IGV_inputs/finspector.fa --bam FI/IGV_inputs/finspector.junction_reads.bam --gtf FI/IGV_inputs/finspector.gtf --fusions_tsv FI/finspector.FusionInspector.fusions.abridged.tsv

FI/IGV_inputs.tar.gz:
	tar czf $@ -C FI IGV_inputs

test-miniwdl: FI/IGV_inputs.tar.gz
	mkdir -p $(MINIWDL_RUN_DIR)
	miniwdl run ../WDL/LR_FusionInspector_to_ORF.wdl -i miniwdl.inputs.json -d $(MINIWDL_RUN_DIR) -o $(MINIWDL_RUN_DIR)/outputs.json --error-json --verbose


clean:
	rm -rf ./LR_FusionInspector_to_ORF_outdir
	rm -rf ./$(MINIWDL_RUN_DIR)
