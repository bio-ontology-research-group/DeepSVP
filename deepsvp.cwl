cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["sh","deepsvp.sh"]
inputs:
  path_data:
    type: Directory
    doc: |
      string(s): list files in a directory.
    inputBinding:
      position: 1
  input_file:
    type: File
  input_pheno:
    type: string
  onto_type:
    type: string
  agg_method:
    type: string 
  maf_value:
    type: float?
  output_file:
    type: string
outputs:
  console_out: stdout
  anno_file_out: 
    type: File
    outputBinding:
      glob: $(inputs.output_file)
arguments:
  - valueFrom: $(inputs.path_data)
  - valueFrom: $(inputs.input_file.basename)
  - valueFrom: $(inputs.input_pheno)
  - valueFrom: $(inputs.maf_value)
  - valueFrom: $(inputs.onto_type)
  - valueFrom: $(inputs.agg_method)
  - valueFrom: $(inputs.output_file)
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: deepsvp.sh
        entry: |-
          deepsvp -d $(inputs.path_data.path)'/' -i $(inputs.input_file.basename) -p input_pheno -maf $(inputs.maf_value)  -m $(inputs.onto_type) -ag $(inputs.agg_method) -o $(inputs.output_file)

