# LR_FusionInspector_to_ORF
model ORFs for fusion transcripts

Run `LR_FusionInspector_to_ORF.py` with `--fusions_tsv` pointing to
`finspector.FusionInspector.fusions.abridged.tsv` to generate
`fusion_orf_breakpoint_annotations.tsv`, which reports whether each predicted
ORF contains coding sequence on both sides of the matched fusion breakpoint.
