cwlVersion: v1.0
class: CommandLineTool

baseCommand: deepsvp
requirements:
  - class: DockerRequirement
    dockerPull: 'coolmaksat/deepsvp'
  - class: InlineJavascriptRequirement

inputs:
  path_data:
    type: Directory
    doc: |
      string(s): list files in a directory.
    inputBinding:
      position: 1
      prefix: '-d'
  input_file:
    type: File
    inputBinding:
      prefix: '-i'
  input_pheno:
    type: string
    inputBinding:
      prefix: '-p'
  onto_type:
    type: string
    inputBinding:
      prefix: '-m'
  agg_method:
    type: string 
    inputBinding:
      prefix: '-ag'
  maf_value:
    type: float?
    inputBinding:
      prefix: '-maf'
  output_file:
    type: string
    inputBinding:
      prefix: '-o'
outputs:
  console_out: stdout
  anno_file_out: 
    type: File
    outputBinding:
      glob: $(inputs.output_file)
