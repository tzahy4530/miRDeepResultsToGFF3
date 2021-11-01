# miRDeepResultsToGFF3

this scripts using the prediction output of miRDeep2 (result_08_10_2021_t_09_57_05.csv) to create GFF3 file.

* Dependency:
  * Python 3.x
  * pandas

* Manual:

  `-i <path> : miRDeep2 prediction output path, like 'result_08_10_2021_t_09_57_05.csv'.`
  
  `-o <path> : GFF3 output path.`
  
  `-t <float> : threshold for the true positive estimate, any value between 0 - 100, default: None.`
  
  `--csv-save : will save the inner tables of miRDeep2 output results as csv.`

* Example Run:

  `python miRDeepResultsToGFF3.py --csv-save -i result_08_10_2021_t_09_57_05.csv -o mirdeep_output.gff3`
