#!/usr/bin/python

import sys
import pandas as pd
import io
pd.options.mode.chained_assignment = None

def filterInputs(inputs_arr, score_threshold, true_positive_threshold, exclude_counts):
    """
    This Function filtering the inputs Dataframe by threshold
    :param inputs_arr: inputs array
    :param score_threshold: Float - threshold for score.
    :param true_positive_threshold: Float - threshold for the true positive estimate.
    :return: filtered dataframes array
    """
    if exclude_counts is None:
        exclude_counts = float('inf')
    if score_threshold is None and true_positive_threshold is None:
        return inputs_arr

    filtered_inputs_arr = []


    for input in inputs_arr:
        if score_threshold is not None:
            input['total read count'] = pd.to_numeric(input['total read count'])
            input['miRDeep2 score'] = pd.to_numeric(input['miRDeep2 score'])
            input = input[(input['total read count'] >= exclude_counts) | (input['miRDeep2 score'] >= score_threshold)]
        if true_positive_threshold is not None:
            try:
                # Novel MicroRNA
                input['true positive probability'] = pd.to_numeric(input.apply(
                    lambda row: row['estimated probability that the miRNA candidate is a true positive'].split(' ')[0],
                    axis=1))
            except:
                # Known MicroRNA
                input['true positive probability'] = pd.to_numeric(input.apply(
                    lambda row: row['estimated probability that the miRNA is a true positive'].split(' ')[0],
                  axis=1))

            input = input[input['true positive probability'] >= true_positive_threshold]
        filtered_inputs_arr.append(input)

    return filtered_inputs_arr


def readMirbaseResults(input_path):
    """
    This function read the miRDeep2 results file (results_<date>.csv)
    and extract the inner tables of the predictions to DataFrames.
    :param input_path: String - the path of the results_<date>.csv.
    :return: Array - of Dataframes - which contain the prediction tables.
    """
    novel_string = ''
    mirbase_string = ''

    read_novel = False
    read_mirbase = False

    with open(input_path) as mirdeep_results:
        lines = mirdeep_results.readlines()

        for line in lines:
            if not read_mirbase:
                if read_novel and len(line) < 100:
                    read_novel = False
                if read_novel:
                    novel_string += line
                if 'novel miRNAs predicted by miRDeep2' in line:
                    read_novel = True

            if not read_novel:
                if read_mirbase and len(line) < 100:
                    read_mirbase = False
                if read_mirbase:
                    mirbase_string += line
                if 'mature miRBase miRNAs detected by miRDeep2' in line:
                    read_mirbase = True

    inputs = []
    try:
        novel_data = io.StringIO(novel_string)
        novel_df = pd.read_csv(novel_data, sep="\t")
        inputs.append(novel_df)
    except:
        pass

    try:
        mirbase_data = io.StringIO(mirbase_string)
        mirbase_df = pd.read_csv(mirbase_data, sep="\t")
        inputs.append(mirbase_df)
    except:
        pass

    return inputs


def run(inputs, output, threshold_tp, threshold_s, exclude_c):
    f"""
    This function create gff3 file from the inputs, This function will add to the miRNA ID -m / -s if its mature/star,
    -number which the number its the frequency of this seq among its type (mature/star)
    :param exclude_c: Int - exclude rows from the filtering when total counts higher than exclude_c .
    :param threshold_s: Float - threshold for score.
    :param threshold_tp: Float - threshold for the true positive estimate.
    :param inputs: Array - with DataFrames of miRDeep results.csv predict tables.
    :param output: String - output path of the GFF formatted file.
    :return: None, at the end of the function gff3 file will created.
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)

    filtered_input = filterInputs(inputs, threshold_s, threshold_tp, exclude_c)

    for input in filtered_input:
        for index, row in input.iterrows():
            details = row['precursor coordinate']
            name = details.split(':')[0]  # *
            positions = details.split(':')[1]
            strand = details.split(':')[2]  # *

            try:
                seq_id = row['provisional id']  # *

            except:
                # That's mean that belong to known microRNA
                seq_id = row["mature miRBase miRNA"]
                seq_id = seq_id.replace('-3p', '')
                seq_id = seq_id.replace('-5p', '')

            star_seq = row['consensus star sequence']
            mature_seq = row['consensus mature sequence']
            hairpin = row['consensus precursor sequence']  # *

            star_position = hairpin.index(star_seq)
            mature_position = hairpin.index(mature_seq)

            seq5p_id = seq_id + '-5p'  # *
            seq3p_id = seq_id + '-3p'  # *

            if star_position > mature_position:
                seq5p = row['consensus mature sequence']  # *
                seq3p = row['consensus star sequence']  # *
                seq5p_freq = len(input[input['consensus mature sequence'] == seq5p])
                seq3p_freq = len(input[input['consensus star sequence'] == seq3p])
                seq5p_id += f'-m-{seq5p_freq}'
                seq3p_id += f'-s-{seq3p_freq}'

            else:
                seq5p = row['consensus star sequence']  # *
                seq3p = row['consensus mature sequence']  # *
                seq5p_freq = len(input[input['consensus star sequence'] == seq5p])
                seq3p_freq = len(input[input['consensus mature sequence'] == seq3p])
                seq5p_id += f'-m-{seq5p_freq}'
                seq3p_id += f'-s-{seq3p_freq}'

            start = int(positions.split('..')[0]) + 1
            end = int(positions.split('..')[1])

            gff_row = [[name, '.', 'pre_miRNA', start, end, '.', strand, '.', f'ID={seq_id}']]

            if strand == '+':
                try:
                    if seq5p == '-':
                        raise ('5p Not Exists')
                    offset5p = len(hairpin.split(seq5p)[0])
                    start5p = start + offset5p
                    end5p = start + offset5p + len(seq5p) - 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seq5p_id}'])
                except:
                    pass

                try:
                    if seq3p == '-':
                        raise ('3p Not Exists')
                    offset3p = len(hairpin.split(seq3p)[0])
                    start3p = start + offset3p
                    end3p = start + offset3p + len(seq3p) - 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seq3p_id}'])
                except:
                    pass

            else:
                try:
                    offset5p = len(hairpin.split(seq5p)[0])
                    end5p = end - offset5p
                    start5p = end - offset5p - len(seq5p) + 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seq5p_id}'])
                except:
                    pass

                try:
                    offset3p = len(hairpin.split(seq3p)[0])
                    end3p = end - offset3p
                    start3p = end - offset3p - len(seq3p) + 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seq3p_id}'])
                except:
                    pass

            miRNAs = pd.DataFrame(gff_row, columns=gff3_columns)

            gff3 = gff3.append(miRNAs)

    with open(output, 'w') as file:
        file.write(version)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    input = None
    output = None
    csv_save = False
    threshold_tp = None
    threshold_s = None
    exclude_c = None
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '--filter-tp':
            threshold_tp = sys.argv[i + 1]
        elif arg == '--filter-s':
            threshold_s = sys.argv[i + 1]
        elif arg == '--exclude-c':
            exclude_c = sys.argv[i + 1]
        elif arg == '--csv-save':
            csv_save = True
            i += 1
            continue

        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path> : miRDeep2 prediction output path, like result_08_10_2021_t_09_57_05\n'
                  f' -o <path> : output path.\n'
                  f' --filter-tp <float> : threshold for the true positive estimate, any value between 0 - 100, '
                  f'default: None.\n '
                  f' --filter-s <float> : threshold for score, default: None.\n'
                  f' --exclude-c <int> : term to ignore the score filter threshold if total counts are higher, default: '
                  f'None.\n '
                  f' --csv-save : will save the inner tables of miRDeep2 output results as csv.\n')
            sys.exit()
        i += 2

    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')

    inputs = readMirbaseResults(input)
    if csv_save is not None:
        count = 1
        for input in inputs:
            input.to_csv(f'table{count}.csv', sep='\t')
            count += 1

    if threshold_tp is not None:
        threshold_tp = float(threshold_tp)

    if threshold_s is not None:
        threshold_s = float(threshold_s)

    if exclude_c is not None:
        exclude_c = int(exclude_c)

    run(inputs, output, threshold_tp, threshold_s, exclude_c)
