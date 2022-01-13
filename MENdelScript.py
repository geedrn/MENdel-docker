import subprocess
import os, sys, getopt, time, shutil, glob
import argparse
import pandas as pd
from Bio.Seq import Seq


def runMenthu(args):
    # subprocess.run(["Rscript", "ExampleRScript.R"])
    outdir = "MENdel_Output"
    menthudir = ''
    lindeldir = ''
    for i in os.listdir(".." + os.sep):
        if "MENTHU" in i:
            menthudir = i
        elif "Lindel" in i:
            lindeldir = i

    if not os.path.exists(".." + os.sep + outdir):
        os.mkdir(".." + os.sep + outdir)
    # os.chdir(".." + os.sep + outdir)
    os.chdir(".." + os.sep + menthudir)
    # print(menthudir)
    # print(os.getcwd())
    subprocess.run(["Rscript", "menthu.R", ".."+os.sep+outdir+os.sep+"MENTHUSummary_"+args.Outfile, args.crispr,
                    args.pam, args.dist2dsb, args.overhang,
                    args.talenOption, args.talenScheme, args.genInputType, args.genInput, args.scoreThreshold,
                    args.T7opt, args.verbose, args.validate])
    time_to_wait = 20
    time_counter = 0
    while not os.path.exists(args.Outfile):
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            break
    # shutil.move(args.Outfile, ".." + os.sep + args.Outfile)
    os.chdir("..")
    runLindel("MENTHUSummary_"+args.Outfile, args.Outfile, lindeldir)


def runLindel(menthuoutfile, outfile, lindeldir):
    # shutil.move(outfile, lindeldir + os.sep + outfile)
    outdir = "MENdel_Output"
    # os.chdir(lindeldir)
    os.chdir(outdir)
    df = pd.read_csv(menthuoutfile)
    i = 0
    lindelScore = []
    lindelPrediction = []
    dna_string = ""
    for j in df['Tool_Type']:
        # print(df.iloc[i])
        # print(j)
        if j == 'NGG' and df['Strand'][i] == 'forward':
            dna_string = df['Context'][i][10::]
            dna_string = dna_string[:-10]
            # print(dna_string)
            # print(df['MENTHU_Score'][i])
        elif j == 'NGG' and df['Strand'][i] == 'complement':
            dna_string = Seq(df['Context'][i]).reverse_complement()[10::]
            dna_string = dna_string[:-10]
            # print(dna_string)
            # print(df['MENTHU_Score'][i])
        subprocess.run(["python", ".."+os.sep+lindeldir+os.sep+"Lindel_prediction.py", str(dna_string), "test_out"])
        txt = glob.glob(os.getcwd() + os.sep + "test_out*.txt")
        for txtfile in txt:
            shutil.move(txtfile, str(outfile[:-4]) + "_" + str(i + 1) + ".tsv")
            fho = open(str(outfile[:-4]) + "_" + str(i + 1) + ".tsv", "r")
            fho.readline()
            for lindelLine in fho:
                # print(lindelLine)
                if "I" in lindelLine.rstrip().split("\t")[2]:
                    lindelScore.append(lindelLine.rstrip().split("\t")[1])
                    lindelPrediction.append(lindelLine.rstrip().split("\t")[2])
                    # print(lindelLine.rstrip().split("\t")[2])
                    break
            fho.close()
            fh = open(str(outfile[:-4]) + "_" + str(i + 1) + ".tsv", "r+")
            content = fh.read()
            fh.seek(0, 0)
            pd.set_option('display.max_rows', None)
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', None)
            pd.set_option('display.max_colwidth', None)
            fh.write(str(df.iloc[i]) + "\n" + str(content))
            fh.close()
        # shutil.move(str(outfile[:-4]) + "_" + str(i + 1) + ".tsv",
        #             ".." + os.sep + str(outfile[:-4]) + "_" + str(i + 1) + ".tsv")
        i += 1
    # shutil.move(outfile, ".." + os.sep + outfile)
    outputProcessing(menthuoutfile, outfile, lindelScore, lindelPrediction)


def outputProcessing(menthuoutfile, outfile, lindelScore, lindelPrediction):
    # os.chdir("..")
    fhs = open("MENdelSummary_" + str(outfile[:-4]) + ".tsv", "w")
    index = 1
    fhs.write(
        "Index\tTarget Sequence\tMENdel SMO\tSMO via MENTHU\tMENTHU Score\tMENTHU Microhomology\tMENTHU Frameshift\t"
        "SMO via Lindel\tLindelScore\tLindel Inserted base\tMENdel Frameshift\n")
    fhr = open(os.getcwd() + os.sep + str(menthuoutfile[:-4]) + ".csv")
    i = 0
    fhr.readline()
    for line in fhr:
        menthuOut = line.rstrip().split(",")
        fhs.write(str(index) + "\t" + str(menthuOut[0]) + "\t")
        index += 1
        # print(menthuOut[1], lindelScore[i])
        if float(menthuOut[1]) > 1.5 or float(lindelScore[i]) > 50.00:
            #fhs.write("1\t")
            fhs.write("Yes\t")
            if float(menthuOut[1]) > 1.5:
                fhs.write("Yes\t" + str(menthuOut[1]) + "\t" + str(menthuOut[7]) + "\t" + menthuOut[2] + "\t")
            else:
                fhs.write("No\t" + str(menthuOut[1]) + "\tNA\t" + menthuOut[2] + "\t")

            if float(lindelScore[i]) > 50.00:
                fhs.write("Yes\t" + str(lindelScore[i]) + "\t" + lindelPrediction[i] + "\t" + "Yes" + "\n")
            else:
                fhs.write("No\t" + str(lindelScore[i]) + "\t" + lindelPrediction[i] + "\t" + menthuOut[2] + "\n")
        else:
            fhs.write("No\tNo\t" + str(menthuOut[1]) + "\tNA\t" + menthuOut[2] + "\tNo\t" + str(lindelScore[i]) + "\t" +
                      lindelPrediction[i] + "\t" + menthuOut[2] + "\n")
        # fhs.write(menthuOut[2] + "\n")
        i += 1
    fhr.close()
    fhs.close()
    txt = glob.glob(os.getcwd() + os.sep + str(outfile[:-4]) + "*.tsv")
    # print(str(outfile[:-4]))


def main():
    parser = argparse.ArgumentParser(description='Inputs arguments for input to Menthu')
    parser.add_argument('-o', '--Outfile', type=str, default='MENdel_outfile.csv',
                        help='character string file name The name of the file to output your results to. If using a '
                             'fasta file with multiple sequences, multiple files will be created, using this as a '
                             'prefix')
    parser.add_argument('-c', '--crispr', type=str, default='T',
                        help='T or F Flags the system to use CRISPR nuclease processing. If this option is T (true), '
                             '"TALEN Option" must be F (false)')
    parser.add_argument('-p', '--pam', type=str, default='NGG',
                        help='A PAM sequence The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are '
                             'allowed. Using N will scan every possible cut site in the target sequence. This '
                             'parameter must be present, but is not used, if "CRISPR Option" is false (i.e., '
                             'you can put a 0 or NA in this spot.)')
    parser.add_argument('-d', '--dist2dsb', type=str, default='-3',
                        help='Integer The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, '
                             'use a negative value (e.g., -3 for SpCa9); for downstream, use a positive value (e.g., '
                             '18 for Cas12a.) This parameter must be present, but is not used, if "CRISPR Option" is '
                             'false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-oh', '--overhang', type=str, default='0',
                        help='Integer >= 0 The length of 5\' overhang produced by the nuclease (e.g., 5 for Cas12a). '
                             'Use 0 for blunt-cutting nucleases. This parameter must be present, but is not used, '
                             'if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-to', '--talenOption', type=str, default='F',
                        help='T or F Flags the system to use TALEN processing. If this option is T (true), "CRISPR '
                             'Option" must be F (false)')
    parser.add_argument('-ts', '--talenScheme', type=str, default='0',
                        help='15-18/14 or 16/15-18 The left arm length, spacer length, and right arm length to use '
                             'when searching for TALEN locations. E.g., for a TALEN with arms 15 nt long, with spacer '
                             '14 nt, use 15/14/15. TALEN arms can be 15-18 nt in length; the spacer should be 14 OR '
                             '16 nt in length (15 is not allowed for the spacer) This parameter must be present, '
                             'but is not used, if "TALEN Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-g', '--genInputType', type=str, default='gb',
                        help='gb ens seq file Flags the system to get a GenBank/RefSeq ID (gb), Ensembl ID (ens), '
                             'DNA sequence (seq), or to expect a FASTA file (file)')
    parser.add_argument('-i', '--genInput', type=str, default='AY214391.1',
                        help='See explanation Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for '
                             '"file" option, or DNA sequence for "seq". If the file name has spaces in it, '
                             'put this parameter in quotes.')
    parser.add_argument('-st', '--scoreThreshold', type=str, default='1',
                        help='Positive number Only output results with MENTHU score above this threshold. Default is '
                             '1.0. We recommend to only use sites with score >= 1.5')
    parser.add_argument('-t7', '--T7opt', type=str, default='F',
                        help='T or F If T (true), only displays results where the gRNA is compatible with T7-cloning.')
    parser.add_argument('-v', '--verbose', type=str, default='F',
                        help='T or F If T (true), outputs progress messages to the console.')
    parser.add_argument('-va', '--validate', type=str, default='F',
                        help='T or F If T (true), checks the command line arguments to make sure they are all valid ('
                             'this may take some time); if F, skip validation checks')
    args = parser.parse_args()
    runMenthu(args)


if __name__ == '__main__':
    main()
