# miRDeepResultsToGFF3

this scripts using the prediction output of miRDeep2 (result_08_10_2021_t_09_57_05.csv) to create GFF3 file.

* Dependency:
  * Python 3.x
  * pandas

* Manual:

  `-i <path> : miRDeep2 prediction output path, like 'result_08_10_2021_t_09_57_05.csv'.`
  
  `-o <path> : GFF3 output path.`
  
  `-seed <path> : classify the reads by seed file, should be separated by tab with columns [miRBase_name, seed], default: None.`
  
  `--create-fasta <path>: create fasta file from the gff3 table. default: None`
  
  `--filter-tp <float> : threshold for the true positive estimate, any value between 0 - 100, default: None.`
  
  `--filter-s <float> : threshold for score, default: None.`
  
  `--exclude-c <int> : term to ignore the score filter threshold if total counts are higher, default: None.`
  
  `--csv-save : will save the inner tables of miRDeep2 output results as csv.`

* Example Run:

  `python miRDeepResultsToGFF3.py --csv-save -i result_08_10_2021_t_09_57_05.csv -o mirdeep_output.gff3`
