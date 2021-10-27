#!/usr/bin/python

import sys
import pandas as pd
import io


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


def run(inputs, output):
    """
    This function create gff3 file from the inputs.
    :param inputs: Array - with DataFrames of miRDeep results.csv predict tables.
    :param output: String - output path of the GFF formatted file.
    :return: None, at the end of the function gff3 file will created.
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    for input in inputs:
        for index, row in input.iterrows():
            details = row['precursor coordinate']
            name = details.split(':')[0]  # *
            positions = details.split(':')[1]
            strand = details.split(':')[2]  # *

            star_seq = row['consensus star sequence']
            mature_seq = row['consensus mature sequence']
            hairpin = row['consensus precursor sequence']  # *

            star_position = hairpin.index(star_seq)
            mature_position = hairpin.index(mature_seq)

            if star_position > mature_position:
                seq5p = row['consensus mature sequence']  # *
                seq3p = row['consensus star sequence']  # *

            else:
                seq5p = row['consensus star sequence']  # *
                seq3p = row['consensus mature sequence']  # *

            try:
                seqId = row['provisional id']  # *

            except:
                seqId = row['tag id']
                name = row["mature miRBase miRNA"]
                name = name.replace('-3p', '')
                name = name.replace('-5p', '')

            seqId5p = name + '-5p'  # *
            seqId3p = name + '-3p'  # *

            start = int(positions.split('..')[0])+1  # *
            end = int(positions.split('..')[1])  # *

            gff_row = [[name, '.', 'pre_miRNA', start, end, '.', strand, '.', f'ID={seqId}']]

            if strand == '+':
                try:
                    if seq5p == '-':
                        raise ('5p Not Exists')
                    offset5p = len(hairpin.split(seq5p)[0])
                    start5p = start + offset5p
                    end5p = start + offset5p + len(seq5p) - 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seqId5p}'])
                except:
                    pass

                try:
                    if seq3p == '-':
                        raise ('5p Not Exists')
                    offset3p = len(hairpin.split(seq3p)[0])
                    start3p = start + offset3p
                    end3p = start + offset3p + len(seq3p) - 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seqId3p}'])
                except:
                    pass

            else:
                try:
                    offset5p = len(hairpin.split(seq5p)[0])
                    end5p = end - offset5p
                    start5p = end - offset5p - len(seq5p) + 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seqId5p}'])
                except:
                    pass

                try:
                    offset3p = len(hairpin.split(seq3p)[0])
                    end3p = end - offset3p
                    start3p = end - offset3p - len(seq3p) + 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seqId3p}'])
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
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '--csv-save':
            csv_save = True
            i += 1
            continue

        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path> : miRDeep2 prediction output path, like result_08_10_2021_t_09_57_05\n'
                  f' -o <path> : output path.\n'
                  f' --csv-save : will save the inner tables of miRDeep2 output results as csv.\n')
            sys.exit()
        i += 2


    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')

    inputs = readMirbaseResults(input)
    if csv_save:
        count = 1
        for input in inputs:
            input.to_csv(f'table{count}.csv',sep='\t')
            count += 1

    run(inputs, output)

