import subprocess as sp
import os
import shlex
import pandas as pd
import numpy as np
import skbio
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.Seq import Seq
from tqdm import tqdm
import re
from optparse import OptionParser
import pickle

FORWARD = True

def Init_Folder(path, user, project):

    if not os.path.exists(path+"/Input/Fastq/split"): 
        os.makedirs(path+"/Input/Fastq/split")
    if not os.path.exists(path+"/Output"):
        os.makedirs(path+"/Output")
    print("makedir")

def FileSplit(path, input, size): #split NGS fastq file as n sequence for each file
    
    print(f"split {path}/Fastq/{input} -l {4 * size} -d -a 6 --additional-suffix=.fastq 'split/split_'")
    sp.run(
        shlex.split(
                f"split {path}/Input/Fastq/{input} -l {4 * size} -d -a 6 --additional-suffix=.fastq '{path}/Input/Fastq/split/split_'"
            )
        )

def ImportBarcode(file, match_type):
    f_barcode = pd.read_csv(file, names = ['Barcode_F', 'Barcode_R', 'Unedit', 'Edit'], skiprows=[0])
    if match_type:
        return f_barcode
    else:
        fwd_dict = dict()
        for fwd in f_barcode['Barcode_F'].dropna():
            rev_dict = dict()
            for rev in f_barcode['Barcode_R'].dropna():
                edit_dict = {'edit':[f_barcode['Unedit'][0],0], 'unedit':[f_barcode['Edit'][0],0],'others':0}
                rev_dict[rev]=edit_dict
            fwd_dict[fwd]=rev_dict
        
        return fwd_dict
    

def CountBarcode(file, barcode):
    seqs = skbio.io.read(file, format = 'fastq', verify = True, variant = 'illumina1.8')
    n=0
    count_edit = 0
    count_unedit = 0
    count_others = 0
    others = []
    print(len(barcode))
    #print(result)
    for seq in tqdm(seqs):
        n+=1
        
        if not FORWARD :
            seq = str(seq)
            seq = Seq.reverse_complement(Seq(seq))
        for fwd in barcode.keys():
            if fwd in seq:
                for rev in barcode[fwd].keys():
                    if rev in seq:
                        
                        if barcode[fwd][rev]['edit'][0] in seq:
                            barcode[fwd][rev]['edit'][1] = barcode[fwd][rev]['edit'][1] +1
                            count_edit+=1
                        elif barcode[fwd][rev]['unedit'][0] in seq:
                            barcode[fwd][rev]['unedit'][1] = barcode[fwd][rev]['unedit'][1] +1
                            count_unedit+=1
                        else:
                            barcode[fwd][rev]['others'] = barcode[fwd][rev]['others']+1
                            count_others+=1
            #if re.search(barcode[i][0], str(seq)) != None:
            #if barcode[i][0] in seq:
            #    result[i][0]+=1
            #    if barcode[i][1] in seq:
            #        result[i][1]+=1
                #elif barcode[i][2] in seq:
                #    result[i][2]+=1
                #else:
                #    result[i][3]+=1
            #    break
    #np.savetxt(file+'_result.csv', result, delimiter=',')
    print(n, count_edit, count_unedit, count_others)
    return barcode


def dict_sum(total, result):
    count=0
    for fwd in total.keys():
        for rev in total[fwd].keys():
            total[fwd][rev]['edit'][1]+=result[fwd][rev]['edit'][1]
            total[fwd][rev]['unedit'][1]+=result[fwd][rev]['unedit'][1]
            total[fwd][rev]['others']+=result[fwd][rev]['others']
            count = count+total[fwd][rev]['edit'][1]+total[fwd][rev]['unedit'][1]+total[fwd][rev]['others']
    print(count)
    return total

def dict_output(result):
    output = []
    for fwd in result.keys():
        for rev in result[fwd].keys():
            output.append([fwd,rev,result[fwd][rev]['edit'][0],result[fwd][rev]['unedit'][0],result[fwd][rev]['edit'][1],result[fwd][rev]['unedit'][1],result[fwd][rev]['others']])
    return output

def SplitedFileLoad(barcode, path, user, project):
    file_list = os.listdir(path + "/Input/Fastq/split")
    file_num = len(file_list)
    for i in range(len(file_list)):
        file_list[i] = path + "/Input/Fastq/split/"+file_list[i]
    result = {}
    print(file_num)
    with ProcessPoolExecutor(max_workers=32) as executor:
        future_to_output = {executor.submit(CountBarcode, file_list[i], barcode): i for i in range(len(file_list))}

        for future in as_completed(future_to_output):
            output = future_to_output[future]
            try:
                if len(result)==0:
                    result = future.result()

                else:
                    result = dict_sum(result, future.result())

                #result = result + future.result()
            except Exception as exc:
                print("%r generated an exception : %s"%(output, exc))
            else:
                with open(path+"/Output/temp.pkl", 'wb') as f: pickle.dump(result,f)
    output = dict_output(result)
    return pd.DataFrame(output, columns=['5-barcode', '3-barcode', 'edit', 'unedit', 'edit_count', 'unedit_count', 'others_count'])
    
def main():
    parser = OptionParser(
        "Extractor_V2 SKKU-GE"
    )
    parser.add_option(
        "-t",
        "--thread",
        default="1",
        type="int",
        dest="multicore",
        help="multiprocessing number, recommendation:t<16",
    )
    parser.add_option(
        "-c",
        "--chunk_number",
        default="100000",
        type="int",
        dest="chunk_number",
        help="split FASTQ, must be multiples of 4. file size < 1G recommendation:40000, size > 1G recommendation:400000",
    )
    parser.add_option(
        "-q",
        "--base_quality",
        default="20",
        dest="base_quality",
        help="NGS read base quality",
    )
    parser.add_option("--user", dest="user_name", help="The user name with no space", default="LBJ")
    parser.add_option(
        "--project", dest="project_name", help="The project name with no space", default="test"
    )
    options, args = parser.parse_args()

    user = options.user_name
    project = options.project_name
    path = os.path.join(os.getcwd(),user,project)

    Init_Folder(path, user, project)
    
    try :
        print(os.listdir(path + "/Input/Fastq/"))
        input_file = [file for file in os.listdir(path + "/Input/Fastq/") if file.endswith('fastq')][0]
        print(input_file)
        barcode = ImportBarcode(path + "/Input/" + project + "_barcode.csv", False)
        print(len(barcode))
        #output_init = np.zeros((barcode.shape[0], barcode.shape[1]+1))
        FileSplit(path, input_file, options.chunk_number)
        output = SplitedFileLoad(barcode, path, user, project)
        output.to_csv(path+"/Output/"+project+".csv")
        #np.savetxt(path + "/" + user + "/" + project + "/Output/PJS_AC_5.csv", output, delimiter=',') #
        
    except:
        print('path initialization')
        
if __name__ == "__main__":
    main()