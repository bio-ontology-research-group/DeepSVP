cwlVersion: v1.0
class: CommandLineTool

baseCommand: /AnnotSV/bin/AnnotSV
requirements:
  - class: DockerRequirement
    dockerPull: 'coolmaksat/deepsvp'
  - class: InlineJavascriptRequirement

inputs:
  annotations_dir:
    type: Directory
    doc: |
      string(s): list files in a directory.
    inputBinding:
      position: 1
      prefix: '-annotationsDir'
  input_file:
    type: File
    inputBinding:
      prefix: '-SVinputFile'
  genome_build:
    type: string 
    inputBinding:
      prefix: '-genomeBuild'
  output_file:
    type: string
    inputBinding:
      prefix: '-outputFile'
outputs:
  anno_file_out: 
    type: File
    outputBinding:
      glob: $(inputs.output_file)
